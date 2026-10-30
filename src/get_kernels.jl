# URLs for ephemerides 
const KERNELS = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/"
const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
const SPK_JUP = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp"
const BPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc"
const EARTH_BPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc"
const EARTH_default = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf"
const TPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"

import Downloads.download

# NAIF regenerates earth_latest_high_prec.bpc twice a week and its predicted span ends
# about three months after generation, so a local copy older than this is re-downloaded
const earth_pck_max_age_days = 30

earth_pck_path() = joinpath(datdir, "earth_latest_high_prec.bpc")

function earth_pck_stale(; max_age_days=earth_pck_max_age_days)
    path = earth_pck_path()
    return !isfile(path) || (time() - mtime(path)) > max_age_days * 86400
end

# download the stable-name Earth PCK to a temporary file and move it into place, so a
# failed download never clobbers the existing copy (no SPICE dependency)
function download_earth_pck!()
    path = earth_pck_path()
    tmp = path * ".part"
    try
        download(EARTH_BPC, tmp; timeout=120)
        mv(tmp, path; force=true)
    finally
        isfile(tmp) && rm(tmp)
    end
    return path
end

# download the SPICE kernels (no SPICE dependency, so deps/build.jl can include this file)
function download_kernels()
    if !isfile(joinpath(datdir, "de440.bsp")); download(SPK, joinpath(datdir, "de440.bsp")); end
    if !isfile(joinpath(datdir, "naif0012.tls")); download(LSK, joinpath(datdir, "naif0012.tls")); end
    if !isfile(joinpath(datdir, "pck00010.tpc")); download(TPC, joinpath(datdir, "pck00010.tpc")); end
    if !isfile(joinpath(datdir, "jup365.bsp")); download(SPK_JUP, joinpath(datdir, "jup365.bsp")); end
    if !isfile(joinpath(datdir, "moon_pa_de440_200625.bpc")); download(BPC, joinpath(datdir, "moon_pa_de440_200625.bpc")); end
    if !isfile(joinpath(datdir, "earth_latest_high_prec.bpc")); download(EARTH_BPC, joinpath(datdir, "earth_latest_high_prec.bpc")); end
    if !isfile(joinpath(datdir, "earth_assoc_itrf93.tf")); download(EARTH_default, joinpath(datdir, "earth_assoc_itrf93.tf")); end

    return nothing
end

# furnish the (already-downloaded) kernels into the SPICE kernel pool; must run at runtime
function furnsh_kernels()
    furnsh(joinpath(datdir, "naif0012.tls"))
    furnsh(joinpath(datdir, "de440.bsp"))
    furnsh(joinpath(datdir, "moon_pa_de440_200625.bpc"))
    furnsh(joinpath(datdir, "pck00010.tpc"))

    # optional long-range Earth PCKs (NAIF earth_*_predict.bpc or earth_*_combined.bpc,
    # placed in data/ by hand; their names are not stable enough to download here).
    # Furnished first so the high-precision file below takes precedence where both cover.
    for f in sort(filter(f -> occursin(r"^earth_.*_(predict|combined)\.bpc$", f), readdir(datdir)))
        furnsh(joinpath(datdir, f))
    end

    furnsh(earth_pck_path())
    furnsh(joinpath(datdir, "earth_assoc_itrf93.tf"))
    furnsh(joinpath(datdir, "jup365.bsp"))

    return nothing
end
