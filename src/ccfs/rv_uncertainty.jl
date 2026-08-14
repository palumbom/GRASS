# variance-weighted Gaussian fits to a CCF, yielding an RV and its uncertainty
import LsqFit: LsqFitResult
# bindings of the RVFromCCF submodule, not re-exported by EchelleCCFs: importing
# them from the parent module yields an unassigned binding rather than an error
import EchelleCCFs.RVFromCCF: default_frac_of_width_to_fit, default_measure_width_at_frac_depth

@. gaussian_line_helper(x, p) = p[4] + p[3] * exp(-0.5*((x-p[1])/p[2])^2)

"""
    est_full_width(vels, ccf; measure_width_at_frac_depth=0.5)

Estimate the full width of a CCF at a given fraction of its depth.

Returns `NaN` when the CCF contains `NaN`, which makes every depth comparison
false. A finite CCF always yields a width, a flat one the full span of `vels`.
"""
function est_full_width(vels::AA{T1,1}, ccf::AA{T2,1};
                        measure_width_at_frac_depth::Real=default_measure_width_at_frac_depth
                        ) where {T1<:AF, T2<:Real}
    minccf, maxccf = extrema(ccf)
    depth = maxccf - minccf
    target_val = minccf + (1-measure_width_at_frac_depth) * depth
    ind1 = findfirst(ccf .<= target_val)
    ind2 = findlast(ccf .<= target_val)
    if isnothing(ind1) || isnothing(ind2)
        return T1(NaN)
    end
    return vels[ind2] - vels[ind1]
end

"""
    find_idx_at_and_around_minimum(vels, ccf; frac_of_width_to_fit=0.5,
                                   measure_width_at_frac_depth=0.5)

Locate the CCF minimum and the index range spanning `frac_of_width_to_fit` of the
CCF width on either side of it.

Returns `(NaN, 1:length(vels))` when the width cannot be estimated; callers detect
this with `isnan` on the first element.
"""
function find_idx_at_and_around_minimum(vels::AA{T1,1}, ccf::AA{T2,1};
                                        frac_of_width_to_fit::Real=default_frac_of_width_to_fit,
                                        measure_width_at_frac_depth::Real=default_measure_width_at_frac_depth
                                        ) where {T1<:AF, T2<:Real}
    # do a prelim fit to get the width
    full_width = est_full_width(vels, ccf, measure_width_at_frac_depth=measure_width_at_frac_depth)
    if isnan(full_width)
       return (NaN, 1:length(vels))
    end

    # find the min and fit only that
    amin = argmin(ccf)
    if amin == 1 || amin==length(vels)
        offset = max(1,floor(Int64,length(vels)//4))
        amin = argmin(view(ccf,offset:(length(vels)-offset)))
        amin += offset-1
    end
    # convert to T1 so the bounds share the element type of vels: searchsortednearest
    # in utils.jl dispatches on the needle and the haystack eltype being identical
    lend = vels[amin] - T1(frac_of_width_to_fit) * full_width
    rend = vels[amin] + T1(frac_of_width_to_fit) * full_width

    # get the indices
    lind = searchsortednearest(view(vels,1:amin), lend)
    rind = amin + searchsortednearest(view(vels,(amin+1):length(vels)), rend)
    inds = lind:rind

    return (amin, inds)
end

function ccf_vcov(fit::LsqFitResult, ccf_var::AA{T,1}) where T<:Real
    # covariance matrix of the fit parameters
    J = fit.jacobian
    covar = pinv(J' * J) * mean(ccf_var)
    return covar
end

function ccf_stderror(fit::LsqFitResult, ccf_var::AA{T,1}) where T<:Real
    covar = ccf_vcov(fit, ccf_var)
    # standard errors are the sqrt of the diagonal
    vars = diag(covar)
    return sqrt.(abs.(vars))
end

"""
    measure_rv_from_ccf_gaussian(vels, ccf, ccf_var, mrv)

Fit a Gaussian to the core of a CCF, weighting by the inverse CCF variance, and
return `(rv=..., σ_rv=...)`. Falls back to a variance-weighted quadratic fit when
the Gaussian fit does not converge.

`ccf` and `ccf_var` must be on a common absolute scale, as returned together by
`calc_ccf(λs, flux, var, ...)`; normalizing one without the other changes σ_rv.
`mrv` supplies `frac_of_width_to_fit`, `measure_width_at_frac_depth`, and
`init_guess_ccf_σ`, so it must be a Gaussian measurement type.

Weighting by `1 ./ ccf_var` treats the velocity lags as independent, which they are
not, so σ_rv scales with the CCF grid spacing and understates the true velocity
scatter. See the `calc_rvs_from_ccf(v_grid, ccf, ccf_var)` docstring before relying
on it as an error bar.
"""
function measure_rv_from_ccf_gaussian(vels::AA{T1,1}, ccf::AA{T2,1},
                                      ccf_var::AA{T3,1}, mrv
                                      ) where {T1<:AF, T2<:Real, T3<:Real}
    if all(ccf .== zero(eltype(ccf))) return (rv=NaN, σ_rv=NaN) end

    # find the min and fit only the part near the minimum of the CCF
    amin, inds = find_idx_at_and_around_minimum(vels, ccf,
                                                frac_of_width_to_fit=mrv.frac_of_width_to_fit,
                                                measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)
    if isnan(amin) return (rv=NaN, σ_rv=NaN) end

    # make initial guess parameters
    μ = vels[amin]
    σ = mrv.init_guess_ccf_σ
    minccf, maxccf = extrema(ccf)
    amp = minccf - maxccf
    y0 = maxccf
    p0 = [μ, σ, amp, y0]

    # a zero variance would give infinite weight; substitute on a copy so the
    # caller's array is left untouched
    var_safe = copy(ccf_var)
    var_safe[var_safe .== 0.0] .= eps()

    local rvfit
    result = curve_fit(gaussian_line_helper, view(vels,inds), view(ccf,inds),
                       (1.0 ./ view(var_safe,inds)), p0)
    if result.converged
        rv = coef(result)[1]
        sigma_rv = ccf_stderror(result, view(var_safe ./ maximum(var_safe), inds))[1]
        rvfit = (rv=rv, σ_rv=sigma_rv)
    else
        @warn "Fit of Gaussian to CCF did not converge. Reverting to fit quadratic to CCF."
        quad_fit_to_ccf = QuadraticFit(frac_of_width_to_fit=mrv.frac_of_width_to_fit,
                                       measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)
        rvfit = quad_fit_to_ccf(vels, ccf, var_safe)
    end
    return rvfit
end
