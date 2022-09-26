using Plots
using ModelingToolkit
using Pkg
using Boltzmann
using Latexify
using QuadGK


c = Boltzmann.Cosmology(h = 0.7, ω_c = 0.25 * (0.7)^2, ω_b = 0.05 * (0.7)^2, Neff = 0, T_cmb = 2.7255)
Boltzmann.derive_params!(c)

const ℓ_max = 3
const speed_of_light_km_s = 299792.458

@variables a z t δ(t) Φ(t) u(t) Ψ(t) δ_b(t) u_b(t) Θ_ℓ(t)[1:ℓ_max + 1] Θ_Pℓ(t)[1:ℓ_max + 1] N_ℓ(t)[1:ℓ_max + 1] τ(a) H(a) cₛ(a) R(a) η(a)
@parameters h Ω_γ Ω_ν Ω_b Ω_c k rh Ω_r Ω_m
function symbolic_boltzmann_ics()

    H₀ = 100 * h
    η₀_approx = minimum((1e-3 / k, 1e-1 * h))
    adotrad = √Ω_r / rh
    a₀ = η₀_approx * adotrad
    η₀ = substitute(η, Dict(a => a₀))
    H_a₀ = substitute(H, Dict(a => a₀))
    ha = H_a₀ / H₀ / rh
    a′_over_a = a₀ * ha
    ω = Ω_m / Ω_r^0.5 / rh

    # ω=H0*Ω_m/sqrt(Ω_r) [h Mpc^-1]
    ω = Ω_m / Ω_r^0.5 / rh
    # neutrino/radiation (constant) ratio [1]
    F_ν = Ω_ν / Ω_r
    # dark matter/matter (constant) ratio[1]
    F_c = Ω_c / Ω_m
    # baryon/matter (constant) ratio [1]
    F_b = 1.0 - F_c
    # photon/radiation (constant) ratio [1]
    F_γ = 1.0 - F_ν
    # matter to radiation energy density (time-dependent)ratio [1]
    ρ_m_over_ρ_r = Ω_m / Ω_r * a₀
    kτ = k * η₀
    kτ_two = kτ^2
    kτ_three = kτ^3
    
    # initial perturbations in synchronous gauge and in Ma&Bertschinger95 notation
    δ_γ = -kτ_two / 3.0 * (1.0 - ω * η₀ / 5.0)  # photon density
    θ_γ = (
        -k
        * kτ_three
        / 36.0
        * (
            1.0
            - 3.0 * (1.0 + 5.0 * F_b - F_ν) / 20.0 / (1.0 - F_ν) * ω * η₀
        )
    )  # photon velocity
    δ_b = 3.0 / 4.0 * δ_γ  # baryon density
    θ_b = θ_γ  # baryon velocit

    δ_c = (
        3.0 / 4.0 * δ_γ
    )  # dm density - note: dm velocity=0 in synchronous gauge
    δ_ν = δ_γ  # neutrino (massless) density
    # neutrino velocity
    θ_ν = (
        -k
        * kτ_three
        / 36.0
        / (4.0 * F_ν + 15.0)
        * (
            4.0 * F_ν
            + 11.0
            + 12.0
            - 3.0
            * (8.0 * F_ν * F_ν + 50.0 * F_ν + 275.0)
            / 20.0
            / (2.0 * F_ν + 15.0)
            * η₀
            * ω
        )
    )
    # neutrino shear
    σ_ν = (
        kτ_two
        / (45.0 + 12.0 * F_ν)
        * (3.0 - 1.0)
        * (1.0 + (4.0 * F_ν - 5.0) / 4.0 / (2.0 * F_ν + 15.0) * η₀ * ω)
    )
    l3_ν = kτ_three * 2.0 / 7.0 / (12.0 * F_ν + 45.0)  # l=3 neutrino mωent - TBC
    # metric perturbation in synchronous gauge
    η_sync = 1.0 - kτ_two / 12.0 / (15.0 + 4.0 * F_ν) * (
        5.0
        + 4.0 * F_ν
        - (16.0 * F_ν^2 + 280.0 * F_ν + 325.0)
        / 10.0
        / (2.0 * F_ν + 15.0)
        * η₀
        * ω
    )

    # cωpute factor α to convert frω synchronous to newtonian gauge
    δ_tot = (
        F_γ * δ_γ
        + F_ν * δ_ν
        + ρ_m_over_ρ_r * (F_b * δ_b + F_c * δ_c)
    ) / (1.0 + ρ_m_over_ρ_r)
    v_tot = (
        (4.0 / 3.0) * (F_γ * θ_γ + F_ν * θ_ν)
        + ρ_m_over_ρ_r * F_b * θ_b
    ) / (1.0 + ρ_m_over_ρ_r)
    α = (
        η_sync
        + 3.0
        / 2.0
        * a′_over_a
        * a′_over_a
        / k
        / k
        * (δ_tot + 3.0 * a′_over_a / k / k * v_tot)
    ) / a′_over_a

    # convert to newtonian gauge
    # newtonian potential perturbation (differs by minus sign in MB95)
    Φ_mb95 = η_sync - a′_over_a * α
    δ_γ = δ_γ - 4.0 * a′_over_a * α
    θ_γ = θ_γ + k * k * α
    δ_b = δ_b - 3.0 * a′_over_a * α
    θ_b = θ_b + k * k * α
    δ_c = δ_c - 3.0 * a′_over_a * α
    θ_c = k * k * α
    δ_ν = δ_ν - 4.0 * a′_over_a * α
    θ_ν = θ_ν + k * k * α
    # σ_ν and l3_ν are gauge invariant

    # convert to initial condition vector
    y = Vector(undef, 5 + 3 * (ℓ_max + 1))

    Θ_ℓ = view(y, 6:6+ℓ_max)
    Θ_Pℓ = view(y, 6 + ℓ_max + 1:6 + 2ℓ_max + 1)
    N_ℓ = view(y, 6 + 2ℓ_max + 2:6 + 3ℓ_max + 2)

    y[1] = -Φ_mb95  # Φ
    y[2] = δ_c  # delta
    y[3] = θ_c / k  # u
    y[4] = δ_b  # deltab
    y[5] = θ_b / k  # ub
    Θ_ℓ[1] = δ_γ / 4.0  # thη₀
    N_ℓ[1] = δ_ν / 4.0  # N_0
    Θ_ℓ[2] = θ_γ / k / 3.0  # theta_1

    N_ℓ[2] = θ_ν / k / 3.0  # N_1
    N_ℓ[3] = σ_ν / 2.0  # N_2
    N_ℓ[4] = l3_ν / 4.0 # N_3    l3_ν=F_nu,3 in MB95 notation = N_3*4 in Dodelson notation
    
    y
end #func
substitute(symbolic_boltzmann_ics()[1], Dict(h=>0.67, k=>1e-2, Ω_m=>0.3, rh=>speed_of_light_km_s * h / (100h)))