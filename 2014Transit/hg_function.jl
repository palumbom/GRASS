import PyPlot
plt = PyPlot

function heney_greenstein_CB(B, h, theta)
    x = tan.(theta ./ 2) ./ h
    delta = (B/2) .* ((1 .+ ((1 .- exp.(-x)) ./ (x))) ./ ((1 .+ x).^2)) 
    return (1 .+ delta)
end

function heney_greenstein_SH(B, h, theta)
    delta = B ./ (1 .+ (1/h) .* tan.(theta ./ 2)) 
    return (1 .+ delta)
end

B_CB = 0.45
h_CB = 0.00045
B_SH = 0.65
B_SH_combined = 0.85
h_SH = 0.00033

# Generate a range of theta values from 0 to π (scattering angles in radians)
theta_values = (LinRange(-0.0005, 1, 100))

# Calculate the phase function for each theta value
phase_values_CB = heney_greenstein_CB(B_CB, h_CB, deg2rad.(theta_values))
phase_values_SH = heney_greenstein_SH(B_SH, h_SH, deg2rad.(theta_values))
phase_values_SH_combined = heney_greenstein_SH(B_SH_combined, h_SH, deg2rad.(theta_values))

# Plot the phase function
plt.plot(theta_values, phase_values_CB, label="CB")
plt.plot(theta_values, phase_values_SH, label="SH")
plt.plot(theta_values, 0.15 .* phase_values_SH_combined .+ phase_values_CB, label="Combined")
plt.xlabel("Scattering Angle (θ, degrees)")
plt.ylabel("Phase Function P(θ)")
plt.legend()
plt.savefig("test.png")
plt.show()
plt.clf()