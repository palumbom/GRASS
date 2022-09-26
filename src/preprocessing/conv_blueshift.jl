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
            vlos = zeros(len(keys(f[k])))
            for (i, t) in enumerate(keys(f[k]))
                # get bisectors
                top = read(f[k][t]["top_ints"])
                bis = read(f[k][t]["bisectors"])
                int = read(f[k][t]["intensities"])
                wid = read(f[k][t]["widths"])

                # identify bad columns and strip them out
                badcols = identify_bad_cols(bis, int, wid)

                if sum(badcols) > 0
                    bis = strip_columns(bis, badcols)
                    int = strip_columns(int, badcols)
                    wid = strip_columns(wid, badcols)
                    top = strip_columns(top, badcols)
                end

                plt.plot(bis, int)
                plt.show()

                # loop over time
                for j in 1:size(bis,2)
                    # get views for time slice
                    bist = view(bis, : ,j)
                    intt = view(int, : ,j)

                    # find the depth
                    bot = minimum(intt)
                    dep = 1.0 - bot

                    # find the lower 5% of the bisector curve
                    idx1 = 2 # the lowest measurement is usually janky
                    idx2 = findfirst(x -> x .>= bot + 0.1 * dep, intt)

                    # take the mean Δλ
                    vlos[i] += mean(view(idx1:idx2))
                end

                # take the average
                vlos[i] /= size(bis,2)

                # now turn that Δλ into a velocity
                # gravitational redshift from Stief et al. (2019)
                vlos[i] = c_ms * vlos[i]/λrest - 633.5
            end

            # write it as an attribute for this disk position
            attr["vconv"] = mean(vlos)
        end
    end
    return nothing
end
