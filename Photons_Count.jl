using Statistics, Plots, NumericalIntegration
using Planck
using Unitful
using Unitful: h, c0, k

# test?


wavelength_arr = collect(1:20000).*1e-9u"m"

star_temp = 5800.0
star_temp *= u"K"

#plot(1:length(wavelength_arr), Planck.blackbody.(u"W/m^3/sr", wavelength_arr, star_temp), 
#    title="Temprature is $(star_temp)", legend="$(star_temp)",
#    lw=4.0, color="black")




function blackbody_radiation(T_S, R_S)
    function to_integrate(λ)
        return (4.0 * pi^2 * R_S^2.0 * 2.0 * h * c0^2) / (λ^5 * expm1(h * c0 / (λ * k * T_S)))
    end
end

f = blackbody_radiation(star_temp, R_S)

y = f.(wavelength_arr)

all_flux = integrate(wavelength_arr, y)

(all_flux) / (4.0 * pi * distance^2)


plot(wavelength_arr, y)



function stefan_boltzmann_power(R_S, T_S)
    return 4.0 * pi * R_S^2.0 * (5.670e-8u"W/m^2/K^4") * T_S^4.0
end

function measured_power(measurement_area, distance, R_S, T_S)
    return stefan_boltzmann_power(R_S, T_S) * (measurement_area / (4.0 * pi * distance^2.0))
end

radius_measurement = 1.0
measurement_area = 1.0u"m^2" # pi * r_measurement^2.0
distance = 150.0e9u"m"
R_S = 696340e3u"m"
T_S = 5772.0u"K"

power1 = measured_power(measurement_area, distance, R_S, T_S)


Unitful.σ

# constants in SI units
const _h = ustrip(u"J*s", h)
const _c0 = ustrip(u"m/s", c0)
const _k = ustrip(u"J/K", k)

using Unitful: 𝐋, 𝐓



