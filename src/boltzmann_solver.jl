

export boltzmann_derivatives!, solve_boltzmann

function boltzmann_derivatives!(dy, y, p, t)
    k, c, ℓ_max, cₛ_fun, τ_dot_fun = p
    lna = t #Our time coordinate is lna
    a = exp(lna)
    z = 1 / a - 1
    
    # Equations from Appendix A in https://arxiv.org/pdf/2112.08395.pdf
    # Derivatives wrt lna
    Hz = H(c, z)
    R = 3Ω_b(c,z)/4Ω_γ(c,z) * a
    cₛ² = cₛ_fun(c, a)^2
    τ_dot = τ_dot_fun(c,a)
    ηz = η(c, z) 


    Φ, δ, u, δ_b, u_b = y[1:5]
    Θ_ℓ = view(y, 6:6+ℓ_max)
    Θ_Pℓ = view(y, 6 + ℓ_max + 1:6 + 2ℓ_max + 1)
    N_ℓ = view(y, 6 + 2ℓ_max + 2:6 + 3ℓ_max + 2)

    Θ_ℓ′ = view(dy, 6:6+ℓ_max)
    Θ_Pℓ′ = view(dy, 6 + ℓ_max + 1:6 + 2ℓ_max + 1)
    N_ℓ′ = view(dy, 6 + 2ℓ_max + 2:6 + 3ℓ_max + 2)
    
    #@show a
    #@show z
    #@show τ_dot
    #@show cₛ²

    # Einstein  Equations

    #Ψ = -Φ - 12 * (100 * c.h / (k * a))^2 * (Ω_γ(c, z) * Θ_ℓ[3] + Ω_ν(c, z) * N_ℓ[3])
    Ψ = -Φ - 12 * (100 * c.h / (k * a))^2 * (c.Ω_γ₀ * Θ_ℓ[3] + c.Ω_ν₀ * N_ℓ[3])
    Π = Θ_ℓ[3] + Θ_Pℓ[1] + Θ_Pℓ[3]
    Φ′ = Ψ - (k / (a*Hz))^2 * Φ / 3 + 0.5 * (c.h * 100 / Hz)^2 * ((c.Ω_c₀ * δ + c.Ω_b₀ * δ_b) * a^(-3) + 4 * a^(-4) * (c.Ω_γ₀ * Θ_ℓ[1] + c.Ω_ν₀ * N_ℓ[1]))

    # Dark Matter
    δ′ = - k * u / (a*Hz) - 3 * Φ′
    u′ = -u + k * Ψ / (a*Hz)

    # Baryonic Matter
    δ_b′ = - k * u_b / (a*Hz) - 3 * Φ′
    u_b′ = -u_b + k * Ψ / (a*Hz) + τ_dot / (R * a * Hz) * (u_b - 3*Θ_ℓ[2]) + k * cₛ² * δ_b / (a*Hz)
    
    # Photon temperature
    Θ_ℓ′[1] = -k * Θ_ℓ[2] / (a*Hz) - Φ′
    Θ_ℓ′[2] = k * (Θ_ℓ[1] - 2Θ_ℓ[3] + Ψ) / (3a*Hz) + τ_dot * (Θ_ℓ[2] - u_b / 3) / (a*Hz)
    Θ_ℓ′[3] = k * (2Θ_ℓ[2] - 3Θ_ℓ[4]) / (5a*Hz) + τ_dot * (Θ_Pℓ[1] - Π / 2) / (a*Hz) #Eq in paper, equivalent to the one in PyCosmo source
    for ℓ_ind in 4:ℓ_max
        ℓ = ℓ_ind -1 # ℓ = 3 is at index 4 and so on
        Θ_ℓ′[ℓ_ind] = k * (ℓ * Θ_ℓ[ℓ_ind - 1] - (ℓ + 1) * Θ_ℓ[ℓ_ind+1]) / (a * Hz * (2ℓ + 1)) + τ_dot * Θ_ℓ[ℓ_ind] / (a*Hz)
    end #for

    # Photon polarization

    Θ_Pℓ′[1] = -k * Θ_Pℓ[2] / (a*Hz) + τ_dot * (Θ_Pℓ[1] - Π / 2) / (a*Hz)
    Θ_Pℓ′[2] = k * (Θ_Pℓ[1] - 2Θ_Pℓ[3]) / (3a*Hz) + τ_dot * Θ_Pℓ[2] / (a*Hz)
    Θ_Pℓ′[3] = k * (2Θ_Pℓ[2] - 3Θ_Pℓ[4]) / (5a*Hz) + τ_dot * (Θ_Pℓ[3] - Π / 10) / (a*Hz)

    for ℓ_ind in 4:ℓ_max
        ℓ = ℓ_ind -1
        Θ_Pℓ′[ℓ_ind] = k * (ℓ * Θ_Pℓ[ℓ_ind - 1] - (ℓ + 1) * Θ_Pℓ[ℓ_ind+1]) / (a * Hz * (2ℓ + 1)) + τ_dot * Θ_Pℓ[ℓ_ind] / (a*Hz)
    end #for

    # Massless neutrinos

    N_ℓ′[1] = -k * N_ℓ[2] / (a*Hz) - Φ′
    N_ℓ′[2] = k * (N_ℓ[1] - 2N_ℓ[3] + Ψ) / (3a*Hz)

    for ℓ_ind in 3:ℓ_max
        ℓ = ℓ_ind -1
        N_ℓ′[ℓ_ind] = k * (ℓ * N_ℓ[ℓ_ind - 1] - (ℓ + 1) * N_ℓ[ℓ_ind+1]) / (a * Hz * (2ℓ + 1))
    end #for
    

    # truncations for ℓ_max
    Θ_ℓ′[ℓ_max + 1] = (k * Θ_ℓ[ℓ_max] - ((ℓ_max + 1) / ηz - τ_dot) * Θ_ℓ[ℓ_max+1]) / (a*Hz)
    Θ_Pℓ′[ℓ_max + 1] = (k * Θ_Pℓ[ℓ_max] - ((ℓ_max + 1) / ηz - τ_dot) * Θ_Pℓ[ℓ_max+1]) / (a*Hz)
    N_ℓ′[ℓ_max + 1] = (k * N_ℓ[ℓ_max] - (ℓ_max + 1) * N_ℓ[ℓ_max+1] / ηz) / (a*Hz)

    dy[1] = Φ′
    dy[2] = δ′
    dy[3] = u′
    dy[4] = δ_b′
    dy[5] = u_b′   


    #@show @view dy[1:5]
    #@show Θ_ℓ′
    #@show Θ_Pℓ′
    #@show N_ℓ′

end #func


function boltzmann_ics(c::Cosmology, k, ℓ_max)
    @assert ℓ_max > 2
    
    H₀ = c.h * 100
    Ω_m = c.Ω_c₀ + c.Ω_b₀
    Ω_r = c.Ω_ν₀ + c.Ω_γ₀
    rh = hubble_distance(c) * c.h #h^-1 Mpc
    η₀_approx = minimum((1e-3 / k, 1e-1 * c.h))
    a₀ = η₀_approx * sqrt(Ω_r) / rh
    z₀ = 1 / a₀ - 1
    H_a₀ = H(c, z₀)
    ha = H_a₀ / H₀ / rh
    

    # initial conformal time [h^-1 Mpc] - warning: this is in MB95 notation and thus not
    # the optical depth
    η₀ = η(c, z₀)
    # da/d(eta)/a [h Mpc^-1]
    a′_over_a = a₀ * ha
    # ω=H0*Ω_m/sqrt(Ω_r) [h Mpc^-1]
    ω = Ω_m / Ω_r^0.5 / rh
    # neutrino/radiation (constant) ratio [1]
    F_ν = c.Ω_ν₀ / Ω_r
    # dark matter/matter (constant) ratio[1]
    F_c = c.Ω_c₀ / Ω_m
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
    y = zeros(Float64, 5 + 3 * (ℓ_max + 1))

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
    
    #@show y
    #@show Θ_ℓ
    #@show Θ_Pℓ
    #@show N_ℓ
    #@show a₀
    #@show log(a₀)

    return y, log(a₀)

end #func

function solve_boltzmann(c::Cosmology, k, ℓ_max)

    z, x_He, x_H, T_m, x_e = precompute_recfast(c::Cosmology)
    cₛ_fun, τ_dot_fun =  recombination_functions(c, z, x_He, x_H, T_m, x_e; autodiff = false)
    p = k, c, ℓ_max, cₛ_fun, τ_dot_fun
    y0, loga0 = boltzmann_ics(c::Cosmology, k, ℓ_max)
    tspan = (loga0, log(1.))
    prob = ODEProblem(boltzmann_derivatives!, y0, tspan, p)
    alg = TRBDF2()
    #alg = CVODE_BDF()
    #alg = Tsit5()
    sol = solve(prob, alg, reltol=1e-5, abstol=1e-5)
    sol
    

end #func

