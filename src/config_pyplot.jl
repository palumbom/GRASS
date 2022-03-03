# plots from python
using PyCall
import PyPlot; const plt = PyPlot; plt.ioff()

const pycopy = PyNULL()
const mpl = PyNULL()
function __init__()
    copy!(pycopy, pyimport("copy"))
    copy!(mpl, pyimport("matplotlib"))
    return nothing
end
