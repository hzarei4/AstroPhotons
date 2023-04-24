using Statistics, Plots, NumericalIntegration
using Unitful
using Unitful: h, c0, k

# test?


wavelength_arr = collect(400:700).*1e-9u"m"

radius_measurement = 1.0
measurement_area = 1.0u"m^2" # pi * r_measurement^2.0
distance = 150.0e9u"m"
R_S = 696340e3u"m"
T_S = 5772.0u"K"


Base.@kwdef mutable struct Star
    distance # in ly
    R_S # in R_sun
    T_S # in K
end

Vega = Star(2500.0, 2.5, 9600.0u"K")

Vega.T_S
#plot(1:length(wavelength_arr), Planck.blackbody.(u"W/m^3/sr", wavelength_arr, star_temp), 
#    title="Temprature is $(star_temp)", legend="$(star_temp)",
#    lw=4.0, color="black")


function blackbody_radiation(T_S, R_S)
    R_S *= 6.957e8u"m"
    
    function to_integrate(λ)
        return (4.0 * pi^2 * R_S^2.0 * 2.0 * h * c0^2) / (λ^5 * expm1(h * c0 / (λ * k * T_S)))
    end
end

f = blackbody_radiation(Vega.T_S, Vega.R_S)

y = f(mean(wavelength_arr))
y1 = y / (4.0 * pi * (9.46e12u"m" * Vega.distance)^2) .* ((pi * 0.2^2)u"m^2")

all_flux = y1 * mean(wavelength_arr)

all_flux = integrate(wavelength_arr, y1);

(all_flux) 

# all flux in the v band of the spectrum of the Vega star is 3.636e-23 J s^-1 m^-2 Hz^-1


n_photons = all_flux / (h * c0 / mean(wavelength_arr))

plot(wavelength_arr, y)



function stefan_boltzmann_power(R_S, T_S)
    return 4.0 * pi * R_S^2.0 * (5.670e-8u"W/m^2/K^4") * T_S^4.0
end

function measured_power(measurement_area, distance, R_S, T_S)
    return stefan_boltzmann_power(R_S, T_S) * (measurement_area / (4.0 * pi * distance^2.0))
end



power1 = measured_power(measurement_area, distance, R_S, T_S)


Unitful.σ

# constants in SI units
const _h = ustrip(u"J*s", h)
const _c0 = ustrip(u"m/s", c0)
const _k = ustrip(u"J/K", k)

using Unitful: 𝐋, 𝐓



