# wrapper functions for measuring velocities with EchelleCCFs
import EchelleCCFs
import EchelleCCFs: BasicCCFPlan
import EchelleCCFs: TopHatCCFMask as TopHatMask
import EchelleCCFs: MeasureRvFromCCFGaussian as GaussianFit
import EchelleCCFs: MeasureRvFromCCFQuadratic as QuadraticFit
import EchelleCCFs: AbstractCCFMaskShape as MaskShape
import EchelleCCFs: AbstractMeasureRvFromCCF as FitType


"""
    calc_ccf(lambdas, intensities, lines, depths, resolution; normalize=true)

Compute the cross correlation function from a spectrum (lambdas and intensities) with a mask computed from a line list (lines).

# Arguments
- `lambdas::AbstractArray{Float64,1}`: List of wavelengths in spectrum.
- `intensities::AbstractArray{Float64,1}`: List of intensities in spectrum.
- `lines::AbstractArray{Float64,1}`: List of line centers for use in CCF mask.
- `depths::AbstractArray{Float64,1}`: List of line depths for use as weights in CCF mask.
- `resolution::Float64`: Spectral resolution of spcectrum.
"""
function calc_ccf(lambdas::AA{Float64,1}, intensities::AA{Float64,1},
                  lines::AA{Float64,1}, depths::AA{Float64,1},
                  resolution::Float64; normalize::Bool=true, Δv_max=15e3,
                  mask_type::Type{T}=TopHatMask) where {T<:MaskShape}
    # make a line list
    line_list = EchelleCCFs.BasicLineList(lines, depths)

    # set mask properties
    speed_of_light = c_ms
    mask_width = speed_of_light/resolution
    mask_shape = T(mask_width)

    Δv_step = 200.0

    # make ccf_plan
    ccf_plan = BasicCCFPlan(line_list=line_list, mask_shape=mask_shape, step=Δv_step, max=Δv_max)

    # calculate the ccf
    v_grid = EchelleCCFs.calc_ccf_v_grid(ccf_plan)
    ccf = EchelleCCFs.ccf_1D(lambdas, intensities, ccf_plan)

    # normalize if normalize==true
    if normalize
        ccf ./= maximum(ccf)
    end
    return v_grid, ccf
end

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,2},
                  lines::AA{T,1}, depths::AA{T,1},
                  resolution::T; kwargs...) where {T<:Float64}
    func = x -> calc_ccf(lambdas, x, lines, depths, resolution; kwargs...)
    out = mapslices(func, intensities, dims=1)
    return out[1][1], cat([x[2] for x in out]..., dims=2)
end

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,1},
                  spec::SpecParams{T}; kwargs...) where {T<:Float64}
    return calc_ccf(lambdas, intensities, spec.lines, spec.depths, spec.resolution; kwargs...)
end

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,2},
                  spec::SpecParams{T}; kwargs...) where {T<:Float64}
    func = x -> calc_ccf(lambdas, x, spec; kwargs...)
    outs = mapslices(func, intensities, dims=1)
    return outs[1][1], cat([x[2] for x in outs]..., dims=2)
end


"""
    calc_rvs_from_ccf(v_grid, ccf)

Calculate apparent radial velocity from a CCF and velocity grid.

# Arguments
- `v_grid::AbstractArray{Float64,1}`: List of velocities returned by calc_ccf.
- `ccf::AbstractArray{Float64,1}`: CCF values returned by calc_ccf.
"""
function calc_rvs_from_ccf(v_grid::AA{Float64,1}, ccf::AA{Float64,1};
                           frac_of_width_to_fit::Float64=0.5,
                           fit_type::Type{T}=GaussianFit) where {T<:FitType}
    mrv = T(frac_of_width_to_fit=frac_of_width_to_fit)
    return mrv(v_grid, ccf)
end

function calc_rvs_from_ccf(v_grid::AA{T,1}, ccf::AA{T,2}; kwargs...) where {T<:Float64}
    func = x -> calc_rvs_from_ccf(v_grid, x; kwargs...)
    out = mapslices(func, ccf, dims=1)
    return vec.(unzip(out))
end
