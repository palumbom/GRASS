struct GPUSolarData{T1<:AF}
    bis::CuArray{T1,3}
    int::CuArray{T1,3}
    wid::CuArray{T1,3}
    dep_contrast::CuArray{T1,1}
    cbs::CuArray{T1,1}
    len::CuArray{Int32,1}
    ax::CuArray{Int32,1}
    mu::CuArray{T1,1}
end

function GPUSolarData(soldata::SolarData{T}; precision::DataType=Float64) where T<:AF
    # collect attributes that are 1D arrays
    len = collect(values(soldata.len))
    cbs = collect(values(soldata.cbs))
    dep_contrast = collect(values(soldata.dep_contrast))

    # allocate memory for bisector + width data
    npositions = length(len)
    bis = zeros(100, maximum(len), npositions)
    int = zeros(100, maximum(len), npositions)
    wid = zeros(100, maximum(len), npositions)
    for (ind,(key,val)) in enumerate(soldata.len)
        bis[:, 1:val, ind] .= soldata.bis[key]
        int[:, 1:val, ind] .= soldata.int[key]
        wid[:, 1:val, ind] .= soldata.wid[key]
    end

    # get the value of mu and ax codes
    disc_ax = GRASS.parse_ax_string.(getindex.(keys(soldata.len),1))
    disc_mu = GRASS.parse_mu_string.(getindex.(keys(soldata.len),2))

    # get indices to sort by mus
    inds_mu = sortperm(disc_mu)
    disc_mu .= disc_mu[inds_mu]
    disc_ax .= disc_ax[inds_mu]

    # get the arrays in mu sorted order
    len .= len[inds_mu]
    cbs .= cbs[inds_mu]
    dep_contrast .= dep_contrast[inds_mu]
    bis .= view(bis, :, :, inds_mu)
    int .= view(int, :, :, inds_mu)
    wid .= view(wid, :, :, inds_mu)

    # get indices to sort by axis within mu sort
    for mu_val in unique(disc_mu)
        inds1 = (disc_mu .== mu_val)
        inds2 = sortperm(disc_ax[inds1])
        disc_mu[inds1] .= disc_mu[inds1][inds2]
        disc_ax[inds1] .= disc_ax[inds1][inds2]

        len[inds1] .= len[inds1][inds2]
        cbs[inds1] .= cbs[inds1][inds2]
        dep_contrast[inds1] .= dep_contrast[inds1][inds2]
        bis[:, :, inds1] .= bis[:, :, inds1][:, :, inds2]
        int[:, :, inds1] .= int[:, :, inds1][:, :, inds2]
        wid[:, :, inds1] .= wid[:, :, inds1][:, :, inds2]
    end

    # copy to GPU and return composite type
    bis_gpu = CuArray{precision}(bis)
    int_gpu = CuArray{precision}(int)
    wid_gpu = CuArray{precision}(wid)
    dep_contrast_gpu = CuArray{precision}(dep_contrast)
    cbs_gpu = CuArray{precision}(cbs)
    len_gpu = CuArray{Int32}(len)
    disc_ax_gpu = CuArray{Int32}(disc_ax)
    disc_mu_gpu = CuArray{precision}(disc_mu)

    return GPUSolarData(bis_gpu, int_gpu, wid_gpu, dep_contrast_gpu,
                        cbs_gpu, len_gpu, disc_ax_gpu, disc_mu_gpu)
end
