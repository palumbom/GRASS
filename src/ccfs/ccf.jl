# wrapper functions for measuring velocities with EchelleCCFs
import EchelleCCFs
import EchelleCCFs: BasicCCFPlan, calc_doppler_factor, project_mask!
import EchelleCCFs: TopHatCCFMask as TopHatMask
import EchelleCCFs: MeasureRvFromCCFGaussian as GaussianFit
import EchelleCCFs: MeasureRvFromCCFQuadratic as QuadraticFit
import EchelleCCFs: AbstractCCFMaskShape as MaskShape
import EchelleCCFs: AbstractMeasureRvFromCCF as FitType

"""
    calc_ccf(λs, flux, var, lines, depths, resolution; mask_width=c_ms/resolution,
             Δv_max=15e3, Δv_step=100.0, mask_type=TopHatMask)

Compute the cross-correlation function (CCF) and its variance from a spectrum
(λs, flux) and a per-pixel flux variance, using a line-list mask.

The CCF and its variance are returned unnormalized and on a common absolute scale.
`calc_rvs_from_ccf(v_grid, ccf, ccf_var)` requires that; rescaling one without the
other changes the reported σ_rv.

# Arguments
- `λs::AbstractArray{Float64,1}`: wavelength grid for the spectrum.
- `flux::AbstractArray{Float64,1}`: flux values on the wavelength grid.
- `var::AbstractArray{<:Real,1}`: per-pixel variance of `flux`.
- `lines::AbstractArray{Float64,1}`: line centers for the CCF mask.
- `depths::AbstractArray{Float64,1}`: line depths used as mask weights.
- `resolution::Float64`: spectral resolution of the spectrum.

# Keyword Arguments
- `mask_width::Float64=c_ms/resolution`: mask width in velocity units.
- `Δv_max::Float64=15e3`: maximum velocity shift (m/s) for the CCF grid.
- `Δv_step::Float64=100.0`: velocity step size (m/s) for the CCF grid.
- `mask_type::Type=TopHatMask`: mask shape type from `EchelleCCFs`.
"""
function calc_ccf(λs::AA{T1,1}, flux::AA{T1,1}, var::AA{T3,1},
                  lines::AA{T1,1}, depths::AA{T1,1},
                  resolution::T1;
                  mask_width::T1=c_ms/resolution,
                  Δv_max::T1=15e3, Δv_step::Float64=100.0,
                  mask_type::Type{T2}=TopHatMask) where {T1<:AF, T2<:MaskShape, T3<:Real}
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

"""
    calc_ccf(λs, flux, lines, depths, resolution; normalize=true, mask_width=c_ms/resolution,
             Δv_max=15e3, Δv_step=100.0, mask_type=TopHatMask)

Compute the cross-correlation function (CCF) from a spectrum (λs, flux) using a line-list mask.

# Arguments
- `λs::AbstractArray{Float64,1}`: wavelength grid for the spectrum.
- `flux::AbstractArray{Float64,1}`: flux values on the wavelength grid.
- `lines::AbstractArray{Float64,1}`: line centers for the CCF mask.
- `depths::AbstractArray{Float64,1}`: line depths used as mask weights.
- `resolution::Float64`: spectral resolution of the spectrum.

# Keyword Arguments
- `normalize::Bool=true`: normalize the CCF to its maximum value.
- `mask_width::Float64=c_ms/resolution`: mask width in velocity units.
- `Δv_max::Float64=15e3`: maximum velocity shift (m/s) for the CCF grid.
- `Δv_step::Float64=100.0`: velocity step size (m/s) for the CCF grid.
- `mask_type::Type=TopHatMask`: mask shape type from `EchelleCCFs`.
"""
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

function calc_ccf(λs::AA{T,1}, flux::AA{T,1}, var::AA{<:Real,1},
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
    calc_rvs_from_ccf(v_grid, ccf; frac_of_width_to_fit=0.75, fit_type=GaussianFit)

Calculate apparent radial velocity from a CCF and velocity grid.

# Arguments
- `v_grid::AbstractArray{Float64,1}`: velocity grid returned by `calc_ccf`.
- `ccf::AbstractArray{Float64,1}`: CCF values returned by `calc_ccf`.

# Keyword Arguments
- `frac_of_width_to_fit::Float64=0.75`: fraction of the CCF width used in the fit.
- `fit_type::Type=GaussianFit`: fit model from `EchelleCCFs`.
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

"""
    calc_rvs_from_ccf(v_grid, ccf, ccf_var; frac_of_width_to_fit=0.75, fit_type=GaussianFit)

Calculate apparent radial velocity and its uncertainty from a CCF, the CCF variance,
and a velocity grid.

`ccf` and `ccf_var` must be on a common absolute scale, as returned together by
`calc_ccf(λs, flux, var, ...)`.

!!! warning "σ_rv is not a calibrated uncertainty"
    `calc_ccf` returns only the diagonal of the CCF covariance, and the fit weights
    by `1 ./ ccf_var`, so neighbouring velocity lags are treated as independent when
    they are not: the tophat mask spans `c_ms/resolution` in velocity, several times
    the default `Δv_step`. The returned `σ_rv` therefore shrinks as `sqrt(Δv_step)`
    without bound and understates the true scatter of recovered velocities. Use it
    as a relative weight between epochs sharing one `Δv_step`, not as an absolute
    error bar. Handling this properly needs the off-diagonal covariance
    (`EchelleCCFs.RVFromCCF` provides `ccf_sample_covar.jl` for that).

# Arguments
- `v_grid::AbstractArray{Float64,1}`: velocity grid returned by `calc_ccf`.
- `ccf::AbstractArray{Float64,1}`: CCF values returned by `calc_ccf`.
- `ccf_var::AbstractArray{Float64,1}`: CCF variance returned by `calc_ccf`.

# Keyword Arguments
- `frac_of_width_to_fit::Float64=0.75`: fraction of the CCF width used in the fit.
- `fit_type::Type=GaussianFit`: supplies the fit hyperparameters; the Gaussian path
  is always used.
"""
function calc_rvs_from_ccf(v_grid::AA{Float64,1}, ccf::AA{Float64,1}, ccf_var::AA{Float64,1};
                           frac_of_width_to_fit::Float64=0.75,
                           fit_type::Type{T}=GaussianFit) where {T<:FitType}
    mrv = T(frac_of_width_to_fit=frac_of_width_to_fit)
    return measure_rv_from_ccf_gaussian(v_grid, ccf, ccf_var, mrv)
end
