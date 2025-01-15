# wrapper functions for measuring velocities with EchelleCCFs
import EchelleCCFs
import EchelleCCFs: BasicCCFPlan, calc_doppler_factor, project_mask!
import EchelleCCFs: TopHatCCFMask as TopHatMask
import EchelleCCFs: MeasureRvFromCCFGaussian as GaussianFit
import EchelleCCFs: MeasureRvFromCCFQuadratic as QuadraticFit
import EchelleCCFs: AbstractCCFMaskShape as MaskShape
import EchelleCCFs: AbstractMeasureRvFromCCF as FitType
import RvSpectMLBase: searchsortednearest

"""
    calc_ccf(lambdas, flux, lines, depths, resolution; normalize=true)

Compute the cross correlation function from a spectrum (λs and flux) with a mask computed from a line list (lines).

# Arguments
- `lambdas::AbstractArray{Float64,1}`: List of wavelengths in spectrum.
- `flux::AbstractArray{Float64,1}`: List of flux in spectrum.
- `lines::AbstractArray{Float64,1}`: List of line centers for use in CCF mask.
- `depths::AbstractArray{Float64,1}`: List of line depths for use as weights in CCF mask.
- `resolution::Float64`: Spectral resolution of spcectrum.
"""
function calc_ccf(λs::AA{T1,1}, flux::AA{T1,1}, var,
                  lines::AA{T1,1}, depths::AA{T1,1},
                  resolution::T1; normalize::Bool=true,
                  mask_width::T1=c_ms/resolution,
                  Δv_max::T1=15e3, Δv_step::Float64=100.0,
                  mask_type::Type{T2}=TopHatMask) where {T1<:AF, T2<:MaskShape}
    # make sure lines are sorted
    if !issorted(lines)
        idx = sortperm(lines)
        lines = view(lines, idx)
        depths = view(depths, idx)
    end

    # make a line list
    line_list = EchelleCCFs.BasicLineList(lines, depths)

    # set mask width
    mask_shape = T2(mask_width)

    # make ccf_plan
    ccf_plan = BasicCCFPlan(line_list=line_list, mask_shape=mask_shape, step=Δv_step, max=Δv_max)

    # get the velocity grid for the ccfs
    v_grid = EchelleCCFs.calc_ccf_v_grid(ccf_plan)

    # compute the mask projection
    ccf = zeros(T1, length(v_grid))
    ccf_var_out = zeros(T1, length(v_grid))
    projection = zeros(T1, length(λs), 1)
    for i in eachindex(v_grid)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask!(projection, λs, ccf_plan, shift_factor=doppler_factor)

        # compute the ccf value at the current velocity shift
        ccf[i] = ccf_plan.allow_nans ? nansum(projection .* flux) : sum(projection .* flux)
        ccf_var_out[i] = ccf_plan.allow_nans ? nansum(var .* projection.^2) : sum(var .* projection.^2)
    end
    return v_grid, ccf, ccf_var_out
end

function calc_ccf(λs::AA{T1,1}, flux::AA{T1,1},
                  lines::AA{T1,1}, depths::AA{T1,1},
                  resolution::T1; normalize::Bool=true,
                  mask_width::T1=c_ms/resolution,
                  Δv_max::T1=15e3, Δv_step::Float64=100.0,
                  mask_type::Type{T2}=TopHatMask) where {T1<:AF, T2<:MaskShape}
    # make sure lines are sorted
    if !issorted(lines)
        idx = sortperm(lines)
        lines = view(lines, idx)
        depths = view(depths, idx)
    end

    # make a line list
    line_list = EchelleCCFs.BasicLineList(lines, depths)

    # set mask width
    mask_shape = T2(mask_width)

    # make ccf_plan
    ccf_plan = BasicCCFPlan(line_list=line_list, mask_shape=mask_shape, step=Δv_step, max=Δv_max)

    # get the velocity grid for the ccfs
    v_grid = EchelleCCFs.calc_ccf_v_grid(ccf_plan)

    # compute the mask projection
    ccf = zeros(T1, length(v_grid))
    projection = zeros(T1, length(λs), 1)
    for i in eachindex(v_grid)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask!(projection, λs, ccf_plan, shift_factor=doppler_factor)

        # compute the ccf value at the current velocity shift
        ccf[i] = ccf_plan.allow_nans ? nansum(projection .* flux) : sum(projection .* flux)
    end

    # normalize if normalize==true
    if normalize
        ccf ./= maximum(ccf)
    end
    return v_grid, ccf
end

function calc_ccf(λs::AA{T1,1}, flux::AA{T1,2},
                  lines::AA{T1,1}, depths::AA{T1,1},
                  resolution::T1; normalize::Bool=true,
                  mask_width::T1=c_ms/resolution,
                  Δv_max::T1=15e3, Δv_step::T1=100.0,
                  mask_type::Type{T2}=TopHatMask) where {T1<:AF, T2<:MaskShape}
    # make sure lines are sorted
    if !issorted(lines)
        idx = sortperm(lines)
        lines = view(lines, idx)
        depths = view(depths, idx)
    end

    # make a line list
    line_list = EchelleCCFs.BasicLineList(lines, depths)

    # set mask width
    mask_shape = T2(mask_width)

    # make ccf_plan
    ccf_plan = BasicCCFPlan(line_list=line_list, mask_shape=mask_shape, step=Δv_step, max=Δv_max)

    # get the velocity grid for the ccfs
    v_grid = EchelleCCFs.calc_ccf_v_grid(ccf_plan)

    # loop over velocity shift
    ccf = zeros(T1, length(v_grid), size(flux, 2))
    projection = zeros(T1, length(λs), 1)
    proj_flux = zeros(T1, length(λs))
    for i in 1:size(ccf, 1)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask!(projection, λs, ccf_plan, shift_factor=doppler_factor)

        # loop over epochs of lfux
        for j in 1:size(ccf, 2)
            # compute the ccf value at the current velocity shift
            proj_flux .= projection .* view(flux, :, j)
            ccf[i,j] = ccf_plan.allow_nans ? nansum(proj_flux) : sum(proj_flux)
        end
    end

    # normalize if normalize==true
    if normalize
        ccf ./= maximum(ccf)
    end
    return v_grid, ccf
end

function calc_ccf!(v_grid::AA{T1,1}, projection::AA{T1,2}, proj_flux::AA{T1,1},
                   ccf::AA{T1,2}, λs::AA{T1,1}, flux::AA{T1,2}, lines::AA{T1,1},
                   depths::AA{T1,1}, resolution::T1; normalize::Bool=true,
                   mask_width::T1=c_ms/resolution, Δv_max::T1=15e3,
                   Δv_step::T1=100.0, mask_type::Type{T2}=TopHatMask) where {T1<:AF, T2<:MaskShape}
    # make sure lines are sorted
    if !issorted(lines)
        idx = sortperm(lines)
        lines = view(lines, idx)
        depths = view(depths, idx)
    end

    # make a line list
    line_list = EchelleCCFs.BasicLineList(lines, depths)

    # set mask width
    mask_shape = T2(mask_width)

    # make ccf_plan
    ccf_plan = BasicCCFPlan(line_list=line_list, mask_shape=mask_shape, step=Δv_step, max=Δv_max)

    # check the size of the ccf
    @assert size(ccf) == (length(v_grid), size(flux, 2))
    @assert size(projection) == (length(λs), 1)
    @assert length(proj_flux) == length(λs)

    ccf .= 0.0
    projection .= 0.0
    proj_flux .= 0.0

    # loop over velocity shift
    for i in 1:size(ccf, 1)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask!(projection, λs, ccf_plan, shift_factor=doppler_factor)

        # loop over epochs of lfux
        for j in 1:size(ccf, 2)
            # compute the ccf value at the current velocity shift
            proj_flux .= projection .* view(flux, :, j)
            ccf[i,j] = ccf_plan.allow_nans ? nansum(proj_flux) : sum(proj_flux)
        end
    end

    # normalize if normalize==true
    if normalize
        ccf ./= maximum(ccf)
    end
    return nothing
end

function calc_ccf(λs::AA{T,1}, flux::AA{T,1}, var::Vector{Float32},
                  spec::SpecParams{T}; kwargs...) where {T<:Float64}
    return calc_ccf(λs, flux, var, spec.lines, spec.depths, spec.resolution; kwargs...)
end

function calc_ccf(λs::AA{T,1}, flux::AA{T,1},
                  spec::SpecParams{T}; kwargs...) where {T<:Float64}
    return calc_ccf(λs, flux, spec.lines, spec.depths, spec.resolution; kwargs...)
end

function calc_ccf(λs::AA{T,1}, flux::AA{T,2},
                  spec::SpecParams{T}; kwargs...) where {T<:Float64}
    # func = x -> calc_ccf(λs, x, spec; kwargs...)
    # outs = mapslices(func, flux, dims=1)
    # return outs[1][1], cat([x[2] for x in outs]..., dims=2)
    return calc_ccf(λs, flux, spec.lines, spec.depths, spec.resolution; kwargs...)
end


"""
    calc_rvs_from_ccf(v_grid, ccf)

Calculate apparent radial velocity from a CCF and velocity grid.

# Arguments
- `v_grid::AbstractArray{Float64,1}`: List of velocities returned by calc_ccf.
- `ccf::AbstractArray{Float64,1}`: CCF values returned by calc_ccf.
"""
function calc_rvs_from_ccf(v_grid::AA{Float64,1}, ccf::AA{Float64,1};
                            frac_of_width_to_fit::Float64=0.75,
                            fit_type::Type{T}=GaussianFit) where {T<:FitType}
    mrv = T(frac_of_width_to_fit=frac_of_width_to_fit)
    return mrv(v_grid, ccf)
end

function calc_rvs_from_ccf(v_grid::AA{T,1}, ccf::AA{T,2}; kwargs...) where {T<:Float64}
    func = x -> calc_rvs_from_ccf(v_grid, x; kwargs...)
    out = mapslices(func, ccf, dims=1)
    return vec.(unzip(out))
end   

function calc_rvs_from_ccf(v_grid::AbstractArray{Float64,1}, ccf::AbstractArray{Float64,1}, var::AbstractArray{Float64,1};
                           frac_of_width_to_fit::Float64=0.75,
                           fit_type::Type{T}=GaussianFit) where {T<:FitType}
    mrv = T(frac_of_width_to_fit=frac_of_width_to_fit)
    return MeasureRvFromCCFGaussian_New(v_grid, ccf, var, mrv)
end

@. gaussian_line_helper(x, p) = p[4] + p[3] * exp(-0.5*((x-p[1])/p[2])^2)

function MeasureRvFromCCFGaussian_New(vels::A1, ccf::A2, ccf_var::A3, mrv) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1} }
    if all(ccf .== zero(eltype(ccf))) return (rv=NaN, σ_rv=NaN)     end
    # find the min and fit only the part near the minimum of the CCF
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)
    if isnan(amin)  return (rv=NaN, σ_rv=NaN)     end

    # make initial guess parameters
    μ = vels[amin]
    σ = mrv.init_guess_ccf_σ                   
    minccf, maxccf = extrema(ccf)
    amp = minccf - maxccf
    y0 = maxccf
    p0 = [μ, σ, amp, y0]

    local rvfit
   # fit and return the mean of the distribution
    ccf_var[ccf_var .== 0.0] .= eps()
    result = curve_fit(gaussian_line_helper, view(vels,inds), view(ccf,inds), (1.0 ./ view(ccf_var,inds)),  p0)
    if result.converged
          rv = coef(result)[1]
          sigma_rv = stderror_new(result, view(ccf_var ./ maximum(ccf_var),inds))[1]
          rvfit = (rv=rv, σ_rv=sigma_rv)
    else
        @warn "Fit of Gaussian to CCF did not converge.  Reverting to fit quadratic to CCF."
        quad_fit_to_ccf = QuadraticFit(frac_of_width_to_fit=mrv.frac_of_width_to_fit,measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)
        rvfit = quad_fit_to_ccf(vels,ccf,ccf_var)
    end
    return rvfit
end

function vcov_new(fit, ccf_var)
    # computes covariance matrix of fit parameters
    J = fit.jacobian
    covar = pinv(J' * J) * mean(ccf_var)
    return covar
end

function stderror_new(fit, ccf_var; rtol::Real=NaN, atol::Real=0)
    covar = vcov_new(fit, ccf_var)
    # then the standard errors are given by the sqrt of the diagonal
    vars = diag(covar)
    return sqrt.(abs.(vars))
end

function est_full_width(vels::A1, ccf::A2; measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth ) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    minccf, maxccf = extrema(ccf)
    depth = maxccf - minccf
    target_val = minccf + (1-measure_width_at_frac_depth) * depth
    ind1 = findfirst(ccf .<= target_val)
    ind2 = findlast(ccf .<= target_val)
    if isnothing(ind1) || isnothing(ind2)
        return NaN
        println("ccf = ",ccf)
        println("minccf= ", minccf, " maxccf= ", maxccf, " depth= ", depth, " measure_width_at_frac_depth= ", measure_width_at_frac_depth, " targetval= ",target_val, " ind1= ", ind1, " ind2= ", ind2)
        @error "est_full_width failed."
    end
    return vels[ind2] - vels[ind1]
end

function find_idx_at_and_around_minimum(vels::A1, ccf::A2; frac_of_width_to_fit::Real = default_frac_of_width_to_fit, measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
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
    lend = vels[amin] - frac_of_width_to_fit * full_width
    rend = vels[amin] + frac_of_width_to_fit * full_width
    # get the indices
    lind = searchsortednearest(view(vels,1:amin), lend)
    rind = amin + searchsortednearest(view(vels,(amin+1):length(vels)), rend)
    inds = lind:rind

    return (amin, inds)
end
