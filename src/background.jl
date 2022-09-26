"""
Code taken from Cosmology.jl (https://github.com/marius311/Cosmology.jl/tree/master)
"""

export Cosmology
export σT, mH, mHe, mp



const eV = 1.782661907e-36  # [kg] - eV/c^2
const mpc = 3.085677581491367e22  # [m]
const speed_of_light_km_s = 299792.458
const G = 6.67408e-11  # [m^3/kg/s^2]
const aRad = 4.0*5.670400e-8 # radiation constant
const σT = 6.6524587158e-29  # [m**2]
const mp = 938.2720813  # [MeV/c**2]
const msun = 1.98848e30  # [kg]
const ħ = 6.582119514e-16  # eV s
const evc² = 1.782661907e-36  # [kg] - eV/c^2
const kb = 8.6173303e-5  # [eV/K]


const H0units = 1. #km/second/Mpc
const ρx_over_ωx = 3(100H0units)^2/(8π)
const π² = π^2


@with_kw mutable struct Cosmology{T<:Real}

    Yp::T = 0.24
    ω_b::T = 0.0225
    ω_c::T = 0.12
    h::T = 0.67
    Neff::T = 3.044
    Ω_k₀::T = 0.
    T_cmb::T = 2.725
    w0::T = -1
    wa::T = 0
    
    # Derived params
    ρ_b₀::T = NaN
    ρ_c₀::T = NaN
    ρ_Λ₀::T = NaN
    ρ_ν₀::T = NaN
    ρ_γ₀::T = NaN
    ρ_k₀::T = NaN
    ρ_crit::T = NaN

    ω_ν::T = NaN
    ω_γ::T = NaN
    ω_Λ::T = NaN
    ω_k::T = NaN

    Ω_b₀::T = NaN
    Ω_c₀::T = NaN
    Ω_Λ₀::T = NaN
    Ω_ν₀::T = NaN
    Ω_γ₀::T = NaN
    Ω_m₀::T = NaN
    Ω_r₀::T = NaN

    T_γ₀::T = NaN
    T_ν₀::T = NaN

    h²::T = NaN

    rh::T = NaN



end #struct

function derive_params!(c)
    c.h² = c.h^2
    H₀ = c.h * 100
    c.ρ_crit = 3 * H₀^2 / (8π * G) * 1e6 * mpc / msun / c.h² #[h^2 M_sun Mpc^-3]
    c.rh = speed_of_light_km_s * c.h / H₀ # Mpc / h
    
    # Photons
    c.T_γ₀ = c.T_cmb #* Kelvin
    c.ρ_γ₀ = 3 * 100^2 / (8π * G) * ħ^3 * (1. / (1e-3 * mpc))^2 * (1. / evc²) * (speed_of_light_km_s * 1e3)^3 # eV
    c.Ω_γ₀ = (π² / 15) * (2.725 * kb)^4 / c.ρ_γ₀ * (c.T_cmb / 2.725)^4 / c.h²
    c.ω_γ = c.Ω_γ₀ * c.h²
    # Neutrinos
    c.Ω_ν₀ = c.Neff * (7 / 8 * (4 / 11)^(4 / 3)) * c.Ω_γ₀
    c.ω_ν = c.Ω_ν₀ * c.h²
    c.ρ_ν₀ = c.Ω_ν₀ * c.ρ_crit #[h^2 M_sun Mpc^-3]
    c.T_ν₀ = c.T_γ₀ * (4 / 11)^(1 / 3)
    # baryons
    c.ρ_b₀ = c.ω_b * ρx_over_ωx
    c.Ω_b₀ = c.ω_b / c.h²
    # CDM
    c.ρ_c₀ = c.ω_c * ρx_over_ωx
    c.Ω_c₀ = c.ω_c / c.h²
    # Curvature
    c.ω_k = c.Ω_k₀ * c.h²
    c.ρ_k₀ = c.ω_k * ρx_over_ωx
    # Dark energy
    c.Ω_Λ₀ = 1 - c.Ω_k₀ - c.Ω_b₀ - c.Ω_c₀ - c.Ω_ν₀ - c.Ω_γ₀
    c.ω_Λ = c.Ω_Λ₀ * c.h²
    c.ρ_Λ₀ = c.ω_Λ * ρx_over_ωx

    c.Ω_m₀ = c.Ω_b₀ + c.Ω_c₀
    c.Ω_r₀ = c.Ω_γ₀ + c.Ω_ν₀


end #func

function Cosmology(args...)
    c = Cosmology(promote(map(float, args)...)...)
    derive_params!(c)
    c
end
function Cosmology(;kwargs...)
    c = Cosmology(;promote(map(float, kwargs)...)...)
    derive_params!(c)
    c
end

# Density evolution

# Radiation
ρ_ν(c::Cosmology, z) = c.ρ_ν₀ * (1 + z)^4
ρ_γ(c::Cosmology, z) = c.ρ_γ₀ * (1 + z)^4
Ω_ν(c::Cosmology, z) = c.Ω_ν₀ * (1 + z)^4
Ω_γ(c::Cosmology, z) = c.Ω_γ₀ * (1 + z)^4
# Matter
ρ_b(c::Cosmology, z) = c.ρ_b₀ * (1 + z)^3
ρ_c(c::Cosmology, z) = c.ρ_c₀ * (1 + z)^3
Ω_b(c::Cosmology, z) = c.Ω_b₀ * (1 + z)^3
Ω_c(c::Cosmology, z) = c.Ω_c₀ * (1 + z)^3
# Curvature
ρ_k(c::Cosmology, z) = c.ρ_k₀ * (1 + z)^2
Ω_k(c::Cosmology, z) = c.Ω_k₀ * (1 + z)^2
# Dark energy
w(c::Cosmology, z) = c.w0 + (1 - (1+z)^-1) * c.wa
f_DE(c::Cosmology, a) = (-3 * (1 + c.w0) + 3 * c.wa * ((a - 1) / log(a) - 1))
ρ_Λ(c::Cosmology, z) = c.ρ_Λ₀ * (1 + z)^(-f_DE(c, (1+z)^-1))
Ω_Λ(c::Cosmology, z) = c.Ω_Λ₀ * (1 + z)^(-f_DE(c, (1+z)^-1))


# Cosmological quanities
E(c::Cosmology, z) = sqrt(Ω_ν(c, z) + Ω_γ(c, z) + Ω_b(c,z) + Ω_c(c,z) + Ω_k(c,z) +  Ω_Λ(c, z))
Eρ(c::Cosmology, z) = sqrt(ρ_ν(c, z) + ρ_γ(c, z) + ρ_b(c,z) + ρ_c(c,z) + ρ_k(c,z) +  ρ_Λ(c, z))
H(c::Cosmology, z) = c.h * 100 * E(c, z)
χ(c::Cosmology, z1, z2) = quadgk(z -> 1 / H(c, z), z1, z2, rtol=1e-8)[1] # Mpc s / km
χ_a(c::Cosmology, a1, a2) = quadgk(a -> 1 / H(c, 1 / a - 1) / a^2, a1, a2, rtol=1e-8)[1] # Mpc s / km
χ(c::Cosmology, z) = χ(c, 0, z) # Mpc s / km
d_hubble(c::Cosmology) = 1 / (c.h * 100)
hubble_distance(c::Cosmology) = speed_of_light_km_s * d_hubble(c)
comoving_distance(c::Cosmology, z) = speed_of_light_km_s * χ(c, z)  # Mpc
η(c::Cosmology, z1, z2) = quadgk(z -> 1 / H(c, z), z1, z2, rtol=1e-8)[1] * c.h * 100 * c.rh
η(c::Cosmology, z) = η(c::Cosmology, z, Inf)
η_a(c::Cosmology, a1, a2) = quadgk(a -> 1 / H(c, 1 / a - 1) / a^2, a1, a2, rtol=1e-8)[1] * c.h * 100 * c.rh
η_a(c::Cosmology, a) = a < 1e-20 ? 0 : η_a(c::Cosmology, 0, a)









