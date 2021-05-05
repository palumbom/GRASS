# wrapper functions for measuring velocities with EchelleCCFs
import EchelleCCFs

"""

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

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,1}, spec::SpecParams{T}; normalize::Bool=true) where T<:AF
    return calc_ccf(lambdas, intensities, spec.lines, spec.depths, spec.resolution, normalize=normalize)
end

function calc_ccf(lambdas::AA{T,1}, intensities::AA{T,2}, spec::SpecParams{T}; normalize::Bool=true) where T<:AF
    func = x -> calc_ccf(lambdas, x, spec, normalize=normalize)
    out = map(func, intensities[:,i] for i in 1:size(intensities, 2))
    return out[1][1], cat([x[2] for x in out]..., dims=2)
end


"""

"""
function calc_rvs_from_ccf(v_grid::AA{T,1}, ccf1::AA{T,1}; frac_of_width_to_fit::T=0.65) where T<:AF
    mrv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=frac_of_width_to_fit)
    return mrv(v_grid, ccf1)
end

function calc_rvs_from_ccf(v_grid::AA{T,1}, ccf1::AA{T,2}; frac_of_width_to_fit::T=0.65) where T<:AF
    func = x -> calc_rvs_from_ccf(v_grid, x, frac_of_width_to_fit=frac_of_width_to_fit)
    out = map(func, ccf1[:,i] for i in 1:size(ccf1, 2))
    return cat([x[1] for x in out]..., dims=1), cat([x[2] for x in out]..., dims=1)
end
