function measure_convective_blueshifts(fname)
    h5open(fname, "r+") do f
        # get rest wavelength
        λrest = read(HDF5.attributes(f)["air_wavelength"])

        # loop over disk positions
        for k in keys(f)
            # get the axis and mu values
            attr = HDF5.attributes(f[k])
            ax = read(attr["axis"])
            mu = read(attr["mu"])

            # loop over times at disk position
            vlos = zeros(length(keys(f[k])))
            for (i, t) in enumerate(keys(f[k]))
                # get bisectors
                top = read(f[k][t]["top_ints"])
                bis = read(f[k][t]["bisectors"])
                int = read(f[k][t]["bis_intensities"])
                wid = read(f[k][t]["widths"])
                int2 = read(f[k][t]["wid_intensities"])

                # identify bad columns and strip them out
                badcols = identify_bad_cols(bis, int, wid, int2)

                if sum(badcols) > 0
                    bis = strip_columns(bis, badcols)
                    int = strip_columns(int, badcols)
                    wid = strip_columns(wid, badcols)
                    top = strip_columns(top, badcols)
                end

                # loop over time
                for j in 1:size(bis, 2)
                    # get views for time slice
                    bist = view(bis, :, j)
                    intt = view(int, :, j)
                    topt = top[j]

                    # find the depth
                    bot = minimum(intt)
                    dep = 1.0 - bot

                    # average of bisector velocity above bottom and below top
                    idx1 = 5 # the lowest-most measurements are usually janky
                    idx2 = findfirst(x -> x .>= topt, intt)

                    # take the mean λ of bisectors
                    vlos[i] += mean(view(bist, idx1:idx2))
                end

                # take the mean over time
                vlos[i] /= size(bis, 2)

                # now turn that λ into a velocity
                # gravitational redshift from Stief et al. (2019)
                vlos[i] = c_ms * (vlos[i] - λrest) / λrest - 633.5
            end

            # write it as an attribute for this disk position
            attr["vconv"] = mean(vlos)
        end
    end
    return nothing
end

function retrieve_vconvs(fname)
    mus = []
    vconvs = []


    h5open(fname, "r+") do f
        for k in keys(f)
            # get the axis and mu values
            attr = HDF5.attributes(f[k])
            ax = read(attr["axis"])
            push!(mus, parse_mu_string(read(attr["mu"])))
            push!(vconvs, read(attr["vconv"]))
        end
    end
    return mus, vconvs
end
