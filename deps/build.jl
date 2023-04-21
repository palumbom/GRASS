using Tar, CodecZlib

download_url = "https://zenodo.org/record/7853524/files/input_data.tar.gz?download=1"
download_loc = joinpath("..", "input_data.tar.gz")

download(download_url, download_loc)

open(GzipDecompressorStream, download_loc) do io
    Tar.extract(io, joinpath("..", "tmp"))
end

mv(joinpath("..", "tmp", "input_data"), joinpath("..", "input_data"), force=true)
