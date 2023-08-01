function rotation_period_gpu(sin_lat::T, A::T, B::T, C::T) where T<:AF
    return convert(T, 360.0)/(A + B * sin_lat^2 + C * sin_lat^4)
end

function patch_velocity_los_gpu(x::T, y::T, rstar::T, polex::T, poley::T, polez::T) where T<:AF
    v0 = convert(T, 0.000168710673) # in (Rsol/day)/speed_of_light
    z = CUDA.sqrt(rstar - calc_r2(x,y))
    sin_lat = (x*polex) + (y*poley) + (z*polez)
    vmax = v0*rstar/rotation_period_gpu(sin_lat, convert(T, 14.713), convert(T, -2.396), convert(T, -1.787))
    return vmax * (polex*y - poley*x)
end
