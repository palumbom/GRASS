# URLs for ephemerides 
const KERNELS = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/"
const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
const BPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc"
const TPC = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"
const moddir = abspath(joinpath(@__DIR__, ".."))
const datdir = joinpath(moddir, "data/")
 
if !isdir(datdir)
    mkdir(datdir)
end

# function to download and load kernels
function get_kernels()
    # Download kernels
    if !isfile(joinpath(datdir, "de440.bsp")); download(SPK, datdir * "de440.bsp"); end
    if !isfile(datdir * "naif0012.tls"); download(LSK, datdir * "naif0012.tls"); end
    if !isfile(datdir * "pck00010.tpc"); download(TPC, datdir * "pck00010.tpc"); end
    if !isfile(datdir * "moon_pa_de440_200625.bpc"); download(BPC, datdir * "moon_pa_de440_200625.bpc"); end

    # load into memory
    furnsh(datdir * "naif0012.tls")
    furnsh(datdir * "de440.bsp")
    furnsh(datdir * "moon_pa_de440_200625.bpc")
    furnsh(datdir * "pck00010.tpc")
    return nothing
end

get_kernels()
