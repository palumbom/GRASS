# Earth body-fixed frame for the observatory position.
#
# ITRF93 needs a binary Earth PCK that covers the epoch. IAU_EARTH (text PCK) never
# runs out of coverage but is off by about 0.08 deg, i.e. ~7 km at the surface,
# ~0.5 m/s of observatory velocity and ~4 arcsec of lunar position, so it is used
# only as an explicit last resort.

const earth_pck_refresh_attempted = Ref(false)

function itrf93_covers(et::Real)
    try
        pxform("ITRF93", "J2000", et)
        return true
    catch err
        err isa SpiceError || rethrow()
        return false
    end
end

# re-download the stable-name high-precision Earth PCK and re-furnish it; at most one
# attempt per session. Returns true if a fresh copy is now loaded.
function refresh_earth_pck!()
    earth_pck_refresh_attempted[] && return false
    earth_pck_refresh_attempted[] = true

    # a copy less than a day old cannot gain coverage; NAIF updates twice a week
    earth_pck_stale(max_age_days=1) || return false
    path = earth_pck_path()
    isfile(path) && unload(path)
    try
        download_earth_pck!()
        furnsh(path)
        return true
    catch err
        @warn "Could not refresh $(basename(path)) from NAIF; keeping the local copy" exception=err
        isfile(path) && furnsh(path)
        return false
    end
end

function earth_frame(et::Real)
    itrf93_covers(et) && return "ITRF93"
    refresh_earth_pck!() && itrf93_covers(et) && return "ITRF93"
    @warn string("No Earth PCK covers ", et2utc(et, "ISOC", 0),
                 "; using IAU_EARTH for the observatory (about 7 km / 0.5 m/s error). ",
                 "For epochs beyond the high-precision file, place a NAIF earth_*_predict.bpc ",
                 "or earth_*_combined.bpc in ", datdir) maxlog=1
    return "IAU_EARTH"
end
