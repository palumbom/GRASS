# wrapper functions for measuring velocities with EchelleCCFs
import EchelleCCFs

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
                  resolution::Float64; normalize::Bool=true,
                  mask_type::Type{T}=EchelleCCFs.TopHatCCFMask) where {T<:EchelleCCFs.AbstractCCFMaskShape}
    # make a line list
    line_list = EchelleCCFs.BasicLineList(lines, depths)

    # set mask properties
    speed_of_light = c_ms
    mask_width = speed_of_light/resolution
    mask_shape = T(mask_width)
    midpoint = 0.0  # just by eye so things aren't wildly shifted
    max_bc = 0.0    #  ~ 2pi AU/year  in m/s
    Δv_step = 400   # m/s arbitrary probably smaller/slower than you need
    Δv_max = 30e3   # m/s arbitrary

    # make ccf_plan
    ccf_plan = EchelleCCFs.BasicCCFPlan(line_list=line_list, mask_shape=mask_shape,
                                        step=Δv_step, max=Δv_max, midpoint=midpoint,
                                        range_no_mask_change=max_bc)

    # calculate the ccf
    v_grid = EchelleCCFs.calc_ccf_v_grid(ccf_plan)
    ccf = EchelleCCFs.ccf_1D(lambdas, intensities, ccf_plan)

    # normalize if normalize==true
    if normalize
        ccf ./= maximum(ccf)
    end
    return v_grid, ccf
end

function calc_ccf(lambdas::AA{Float64,1}, intensities::AA{Float64,2},
                  lines::AA{Float64,1}, depths::AA{Float64,1},
                  resolution::Float64; normalize::Bool=true,
                  mask_type::Type{T}=EchelleCCFs.TopHatCCFMask) where {T<:EchelleCCFs.AbstractCCFMaskShape}
    func = x -> calc_ccf(lambdas, x, lines, depths, resolution, normalize=normalize, mask_type=mask_type)
    out = mapslices(func, intensities, dims=1)
    return out[1][1], cat([x[2] for x in out]..., dims=2)
end

function calc_ccf(lambdas::AA{Float64,1}, intensities::AA{Float64,1},
                  spec::SpecParams{Float64}; normalize::Bool=true,
                  mask_type::Type{T}=EchelleCCFs.TopHatCCFMask) where {T<:EchelleCCFs.AbstractCCFMaskShape}
    return calc_ccf(lambdas, intensities, spec.lines, spec.depths, spec.resolution, normalize=normalize, mask_type=mask_type)
end

function calc_ccf(lambdas::AA{Float64,1}, intensities::AA{Float64,2},
                  spec::SpecParams{Float64}; normalize::Bool=true,
                  mask_type::Type{T}=EchelleCCFs.TopHatCCFMask) where {T<:EchelleCCFs.AbstractCCFMaskShape}
    func = x -> calc_ccf(lambdas, x, spec, normalize=normalize, mask_type=mask_type)
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
                           frac_of_width_to_fit::Float64=0.65,
                           fit_type::Type{T}=EchelleCCFs.MeasureRvFromCCFGaussian) where {T<:EchelleCCFs.AbstractMeasureRvFromCCF}
    mrv = T(frac_of_width_to_fit=frac_of_width_to_fit)
    return mrv(v_grid, ccf)
end

function calc_rvs_from_ccf(v_grid::AA{Float64,1}, ccf::AA{Float64,2};
                           frac_of_width_to_fit::Float64=0.65,
                           fit_type::Type{T}=EchelleCCFs.MeasureRvFromCCFGaussian) where {T<:EchelleCCFs.AbstractMeasureRvFromCCF}
    func = x -> calc_rvs_from_ccf(v_grid, x, frac_of_width_to_fit=frac_of_width_to_fit, fit_type=fit_type)
    out = mapslices(func, ccf, dims=1)
    return vec.(unzip(out))
end
