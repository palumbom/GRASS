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
function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,1}, lines::AA{T,1},
                  depths::AA{T,1}, resolution::T; normalize::Bool=true) where T<:AF
    # make a line list
    line_list = EchelleCCFs.BasicLineList(lines, depths)

    # set mask properties
    speed_of_light = c_ms::T
    mask_width = speed_of_light/resolution
    mask_shape = EchelleCCFs.TopHatCCFMask(mask_width)
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

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,2}, lines::AA{T,1},
                  depths::AA{T,1}, resolution::T; normalize::Bool=true) where T<:AF
    func = x -> calc_ccf(lambdas, x, lines, depths, resolution, normalize=normalize)
    out = mapslices(func, intensities, dims=1)
    return out[1][1], cat([x[2] for x in out]..., dims=2)
end

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,1}, spec::SpecParams{T}; normalize::Bool=true) where T<:AF
    return calc_ccf(lambdas, intensities, spec.lines, spec.depths, spec.resolution, normalize=normalize)
end

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,2}, spec::SpecParams{T}; normalize::Bool=true) where T<:AF
    func = x -> calc_ccf(lambdas, x, spec, normalize=normalize)
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
function calc_rvs_from_ccf(v_grid::AA{T,1}, ccf::AA{T,1}; frac_of_width_to_fit::T=0.65) where T<:AF
    mrv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=frac_of_width_to_fit)
    return mrv(v_grid, ccf)
end

function calc_rvs_from_ccf(v_grid::AA{T,1}, ccf::AA{T,2}; frac_of_width_to_fit::T=0.65) where T<:AF
    func = x -> calc_rvs_from_ccf(v_grid, x, frac_of_width_to_fit=frac_of_width_to_fit)
    out = mapslices(func, ccf, dims=1)
    return vec.(unzip(out))#(getfield.(out, x) for x in fieldnames(eltype(out)))
end
