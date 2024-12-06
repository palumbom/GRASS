# URLs for ephemerides 
const KERNELS = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/"
const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
const SPK_JUP = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp"
const BPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc"
const EARTH_BPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc"
const EARTH_default = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf"
const TPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"

# function to download and load kernels
function get_kernels()
    # Download kernels
    if !isfile(joinpath(datdir, "de440.bsp")); download(SPK, joinpath(datdir, "de440.bsp")); end
    if !isfile(joinpath(datdir, "naif0012.tls")); download(LSK, joinpath(datdir, "naif0012.tls")); end
    if !isfile(joinpath(datdir, "pck00010.tpc")); download(TPC, joinpath(datdir, "pck00010.tpc")); end
    if !isfile(joinpath(datdir, "jup365.bsp")); download(SPK_JUP, joinpath(datdir, "jup365.bsp")); end
    if !isfile(joinpath(datdir, "moon_pa_de440_200625.bpc")); download(BPC, joinpath(datdir, "moon_pa_de440_200625.bpc")); end
    if !isfile(joinpath(datdir, "earth_latest_high_prec.bpc")); download(EARTH_BPC, joinpath(datdir, "earth_latest_high_prec.bpc")); end
    if !isfile(joinpath(datdir, "earth_assoc_itrf93.tf")); download(EARTH_default, joinpath(datdir, "earth_assoc_itrf93.tf")); end

    # load into memory
    furnsh(joinpath(datdir, "naif0012.tls"))
    furnsh(joinpath(datdir, "de440.bsp"))
    furnsh(joinpath(datdir, "moon_pa_de440_200625.bpc"))
    furnsh(joinpath(datdir, "pck00010.tpc"))
    furnsh(joinpath(datdir, "earth_latest_high_prec.bpc"))
    furnsh(joinpath(datdir, "earth_assoc_itrf93.tf"))
    furnsh(joinpath(datdir, "jup365.bsp"))
    furnsh(joinpath(datdir, "eclipse_SPK.bsp"))
    furnsh(joinpath(datdir, "eclipse_ID_mapping_SPK.txt"))

    return nothing
end

get_kernels()
