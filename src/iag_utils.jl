using HTTP
using GZip
using CSV
using DataFrames
using EchelleCCFs: λ_air_to_vac, λ_vac_to_air

function download_iag()
    println(">>> Downloading IAG atlas...")
    url = "http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/txt.gz?J/A+A/587/A65/spvis.dat.gz"
    file = HTTP.download(url, moddir * "input_data/spvis.dat.gz", update_period=Inf)
    println(">>> IAG atlas downloaded!")
    return nothing
end

function read_iag(; isolate::Bool=false, airwav::Float64=5434.5)
    # download the IAG atlas
    file = GRASS.moddir * "input_data/spvis.dat.gz"
    if !isfile(file) || mtime(file) < 1.66e9
        download_iag()
    end

    # read in the IAG atlas
    iag = GZip.open(file, "r") do io
        CSV.read(io, DataFrame, ignorerepeated=true, delim="|", skipto=5,
                 footerskip=1, header=["wavenum", "nflux", "flux"])
    end

    # convert wavenumber to wavelength in angstroms
    wavs = (1.0 ./ iag.wavenum) * 1e8

    # reverse to deal with conversion of units
    reverse!(wavs)
    reverse!(iag.nflux)

    # isolate region around line
    if isolate
        vacwav = λ_air_to_vac(airwav)
        ind1 = findfirst(x -> x .> vacwav-1.0, wavs)
        ind2 = findfirst(x -> x .> vacwav+1.0, wavs[ind1:end]) + ind1
        return λ_vac_to_air.(view(wavs,ind1:ind2)), view(iag.nflux, ind1:ind2)
    else
        return λ_vac_to_air.(wavs), iag.nflux
    end
end
