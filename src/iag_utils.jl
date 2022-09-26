using CSV
using HTTP
using DataFrames
using EchelleCCFs: λ_air_to_vac, λ_vac_to_air

function download_iag_atlas()
    # download the file
    println(">>> Downloading IAG atlas...")
    url = "http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/txt.gz?J/A+A/587/A65/spvis.dat.gz"
    file = HTTP.download(url, datdir * "spvis.dat.gz", update_period=Inf)

    # decompress it
    @assert isfile(file)
    try
        run(`gunzip -q $file`)
    catch e
        nothing
    end

    println(">>> IAG atlas downloaded to " * file)
    return nothing
end

function read_iag_atlas(;isolate::Bool=true, airwav::Float64=5434.5232, buffer::Float64=1.0)
    # download the IAG atlas
    file = datdir * "spvis.dat"
    if !isfile(file) || mtime(file) < 1.66e9
        download_iag_atlas()
    end

    # read in the IAG atlas
    iag = CSV.read(file, DataFrame, ignorerepeated=true, delim="|", skipto=5,
                   footerskip=1, header=["wavenum", "nflux", "flux"])

    # convert wavenumber to wavelength in angstroms
    wavs = (1.0 ./ iag.wavenum) * 1e8

    # reverse to deal with conversion of units
    reverse!(wavs)
    reverse!(iag.nflux)

    # isolate region around line
    if isolate
        vacwav = λ_air_to_vac(airwav)
        ind1 = findfirst(x -> x .> vacwav-buffer, wavs)
        ind2 = findfirst(x -> x .> vacwav+buffer, wavs[ind1:end]) + ind1
        return λ_vac_to_air.(view(wavs,ind1:ind2)), view(iag.nflux, ind1:ind2)
    else
        return λ_vac_to_air.(wavs), iag.nflux
    end
end

function download_iag_blueshifts()
    # download the file
    println(">>> Downloading IAG convective blueshift data...")
    url = "http://cdsarc.u-strasbg.fr/ftp/cats/A&A/587/A65/tablea1.dat"
    file = HTTP.download(url, datdir * "iag_blueshifts.dat", update_period=Inf)

    println(">>> IAG blueshifts downloaded to " * file)
    return nothing
end

function read_iag_blueshifts()
     # download the IAG blueshifts
    file = datdir * "iag_blueshifts.dat"
    if !isfile(file)
        download_iag_blueshifts()
    end

    # read in the IAG atlas
    iag = CSV.read(file, DataFrame, ignorerepeated=true, delim=" ", types=Float64,
                   missingstring="---", header=["wavelength", "EW", "depth", "blueshift"])

    # convert missing to nan
    iag.EW[ismissing.(iag.EW)] .= NaN
    iag[!,:EW] = convert.(Float64, iag[!,:EW])

    iag.blueshift[ismissing.(iag.blueshift)] .= NaN
    iag[!,:blueshift] = convert.(Float64, iag[!,:blueshift])

    # convert to angstroms
    iag.wavelength .*= 10.0
    return iag
end
