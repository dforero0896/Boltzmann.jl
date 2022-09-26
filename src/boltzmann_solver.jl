

export boltzmann_derivatives!, solve_boltzmann

function einstein_zero(c, y, a, k, ℓ_max)
    z = 1 / a -1
    rh = hubble_distance(c) * c.h #h^-1 Mpc
    H₀ = 100c.h
    H₀_Mpc = 1 / rh
    H_Mpc = H(c, z) / (100c.h) * H₀_Mpc
    Φ, δ, u, δ_b, u_b = y[1:5]
    Θ_ℓ = view(y, 6:6+ℓ_max)
    Θ_Pℓ = view(y, 6 + ℓ_max + 1:6 + 2ℓ_max + 1)
    N_ℓ = view(y, 6 + 2ℓ_max + 2:6 + 3ℓ_max + 2)
    (
        (-2 / 3 * (k / (a * H₀))^2 * Φ) + 
        (
            (c.Ω_c₀ * δ + c.Ω_b₀ * δ_b) / a^3
             + 4 * (c.Ω_γ₀ * Θ_ℓ[1] + c.Ω_ν₀ * N_ℓ[1]) / a^4
        )
        + (3 * a * H_Mpc / k)
        * (
            (c.Ω_c₀ * u + c.Ω_b₀ * u_b) / a^3
            + 4 * (c.Ω_γ₀ * Θ_ℓ[2] + c.Ω_ν₀ * N_ℓ[2]) / a^4
        )
    ) / ((c.Ω_c₀ + c.Ω_b₀) / a^3 + (c.Ω_γ₀ + c.Ω_ν₀) / a^4)
end #func

function boltzmann_derivatives_dae!(out, dy, y, p, t)
    k, c, ℓ_max, cₛ_fun, τ_dot_fun = p
    lna = t #Our time coordinate is lna
    a = exp(lna)
    z = 1 / a - 1

    boltzmann_derivatives!(dy, y, p, t)
    zero = einstein_zero(c, y, a, k, ℓ_max)
    out .= y .- dy
    out[end] = zero

end #func

function boltzmann_derivatives!(dy, y, p, lna)
    k, c, ℓ_max, cₛ_fun, τ_dot_fun = p
    a = exp(lna)
    z = 1 / a - 1
    
    
   
    # Equations from Appendix A in https://arxiv.org/pdf/2112.08395.pdf
    # Derivatives wrt lna
    rh = hubble_distance(c) * c.h #h^-1 Mpc
    H₀_Mpc = 1 / rh
    H_Mpc = H(c, z) / (100c.h) * H₀_Mpc
    k_aH = k / (a * H_Mpc)
    
    #R = 3Ω_b(c,z)/4Ω_γ(c,z) * a
    R = 3 * c.Ω_b₀ / (4 * c.Ω_γ₀) * a
    cₛ² = cₛ_fun(c, a)^2
    τ_dot = τ_dot_fun(c,a)
    ηz = η(c, z) 


    Φ, δ, u, δ_b, u_b = y[1:5]
    Θ_ℓ = view(y, 6:6+ℓ_max)
    Θ_Pℓ = view(y, 6 + ℓ_max + 1:6 + 2ℓ_max + 1)
    N_ℓ = view(y, 6 + 2ℓ_max + 2:6 + 3ℓ_max + 2)

    #Θ_ℓ = @view y[6:ℓ_max:end]
    #Θ_Pℓ = @view y[7:ℓ_max:end]
    #N_ℓ = @view y[8:ℓ_max:end]

    Θ_ℓ′ = view(dy, 6:6+ℓ_max)
    Θ_Pℓ′ = view(dy, 6 + ℓ_max + 1:6 + 2ℓ_max + 1)
    N_ℓ′ = view(dy, 6 + 2ℓ_max + 2:6 + 3ℓ_max + 2)
    #Θ_ℓ′ = @view y[6:ℓ_max:end]
    #Θ_Pℓ′ = @view y[7:ℓ_max:end]
    #N_ℓ′ = @view y[8:ℓ_max:end]
    

    # Einstein  Equations
    
    #Ψ = -Φ - 12 * (100 * c.h / (k * a))^2 * (Ω_γ(c, z) * Θ_ℓ[3] + Ω_ν(c, z) * N_ℓ[3])
    Ψ = -Φ - 12 * (H₀_Mpc / (k * a))^2 * (c.Ω_γ₀ * Θ_ℓ[3] + c.Ω_ν₀ * N_ℓ[3])
    Π = Θ_ℓ[3] + Θ_Pℓ[1] + Θ_Pℓ[3]
    Φ′ = (
        Ψ 
        - (k_aH)^2 / 3 * Φ 
        + (H₀_Mpc / H_Mpc)^2 
        /2 
        * (
            (c.Ω_c₀ * δ + c.Ω_b₀ * δ_b) * a^(-3) 
            + 4 * (c.Ω_γ₀ * Θ_ℓ[1] + c.Ω_ν₀ * N_ℓ[1]) * a^(-4)
            )
            )
    
    

    # Dark Matter
    δ′ = - k_aH * u - 3 * Φ′
    u′ = -u + k_aH * Ψ 

    # Baryonic Matter
    δ_b′ = - k_aH * u_b - 3 * Φ′
    u_b′ = (
        -u_b 
        + k_aH * Ψ 
        + τ_dot / (R * a * H_Mpc) * (u_b - 3 * Θ_ℓ[2]) 
        + k_aH * cₛ² * δ_b
        )
    
    # Photon temperature
    Θ_ℓ′[1] = -k_aH * Θ_ℓ[2] - Φ′
    
    Θ_ℓ′[2] = k_aH / 3 * (Θ_ℓ[1] - 2Θ_ℓ[3] + Ψ) + τ_dot / (
        a * H_Mpc
        ) * (
            Θ_ℓ[2] - u_b / 3
            )
    Θ_ℓ′[3] = k_aH / 5 * (2Θ_ℓ[2] - 3Θ_ℓ[4]) + τ_dot / (
        a * H_Mpc
        ) * (
            Θ_ℓ[3] - Π / 10
            )
    for ℓ_ind in 4:ℓ_max
        ℓ = ℓ_ind - 1 # ℓ = 3 is at index 4 and so on
        Θ_ℓ′[ℓ_ind] = (
            k_aH  / (2ℓ + 1) * (ℓ * Θ_ℓ[ℓ_ind - 1] - (ℓ + 1) * Θ_ℓ[ℓ_ind+1]) 
            + τ_dot / (a * H_Mpc) * Θ_ℓ[ℓ_ind])
    end #for

    # Photon polarization

    Θ_Pℓ′[1] = k_aH / 1 * -Θ_Pℓ[2] + τ_dot / (a * H_Mpc) * (
            Θ_Pℓ[1] - Π / 2
            )
    Θ_Pℓ′[2] = (
        k_aH / 3 * (Θ_Pℓ[1] - 2Θ_Pℓ[3]) 
        + τ_dot / (a * H_Mpc) * Θ_Pℓ[2]
        )
    Θ_Pℓ′[3] = k_aH / 5 * (2Θ_Pℓ[2] - 3Θ_Pℓ[4]) + τ_dot / (a * H_Mpc) * (
        Θ_Pℓ[3] - Π / 10
        )

    for ℓ_ind in 4:ℓ_max
        ℓ = ℓ_ind -1
        Θ_Pℓ′[ℓ_ind] = (
            k_aH / (2ℓ + 1) * (ℓ * Θ_Pℓ[ℓ_ind - 1] - (ℓ + 1) * Θ_Pℓ[ℓ_ind+1]) 
            + τ_dot / (a * H_Mpc) * Θ_Pℓ[ℓ_ind]
            )
    end #for

    # Massless neutrinos

    N_ℓ′[1] = -k_aH * N_ℓ[2] - Φ′
    N_ℓ′[2] = k_aH / 3 * (N_ℓ[1] - 2N_ℓ[3] + Ψ)

    for ℓ_ind in 3:ℓ_max
        ℓ = ℓ_ind -1
        N_ℓ′[ℓ_ind] = k_aH  / (2ℓ + 1) * (ℓ * N_ℓ[ℓ_ind - 1] - (ℓ + 1) * N_ℓ[ℓ_ind+1])
    end #for
    

    # truncations for ℓ_max
    ℓ_max_ind = ℓ_max + 1
    Θ_ℓ′[ℓ_max_ind] = (
        1
        / (a * H_Mpc)
        * (
            k * Θ_ℓ[ℓ_max_ind - 1] 
            - ((ℓ_max + 1) / ηz - τ_dot) * Θ_ℓ[ℓ_max_ind]
            )
            )
    Θ_Pℓ′[ℓ_max_ind] = (
        1
        / (a * H_Mpc)
        * (
            k * Θ_Pℓ[ℓ_max_ind - 1] 
            - ((ℓ_max + 1) / ηz - τ_dot) * Θ_Pℓ[ℓ_max_ind]
            )
            )
    N_ℓ′[ℓ_max_ind] = (
        1 / (a * H_Mpc) * (k * N_ℓ[ℓ_max_ind - 1] - (ℓ_max + 1) / ηz * N_ℓ[ℓ_max_ind])
        )

    dy[1] = Φ′
    dy[2] = δ′
    dy[3] = u′
    dy[4] = δ_b′
    dy[5] = u_b′   

    

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

    #Θ_ℓ = @view y[6:ℓ_max:end]
    #Θ_Pℓ = @view y[7:ℓ_max:end]
    #N_ℓ = @view y[8:ℓ_max:end]
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

    return y, log(a₀)

end #func

function solve_boltzmann(c::Cosmology, k::Real, ℓ_max::Int, a_grid::Vector{<:AbstractFloat})

    z, x_He, x_H, T_m, x_e = precompute_recfast(c::Cosmology)
    cₛ_fun, τ_dot_fun =  recombination_functions(c, z, x_He, x_H, T_m, x_e; autodiff = false)
    p = k, c, ℓ_max, cₛ_fun, τ_dot_fun
    y0, loga0 = boltzmann_ics(c::Cosmology, k, ℓ_max)
    tspan = (loga0, 0)
    prob = ODEProblem(boltzmann_derivatives!, y0, tspan, p)
    alg = KenCarp4()
    sol = solve(prob, alg, reltol=1e-5, abstol=1e-5, saveat = log.(a_grid))
    sol
    

end #func

function unpack_fields(sol, ℓ_max)
    fields = Dict()
    fields[:Φ] = [s[1] for s in sol.u]
    fields[:δ] = [s[2] for s in sol.u]
    fields[:u] = [s[3] for s in sol.u]
    fields[:δ_b] = [s[4] for s in sol.u]
    fields[:u_b] = [s[5] for s in sol.u]
    fields[:Θ_ℓ] = [s[6+i] for s in sol.u, i in 0:ℓ_max]
    fields[:Θ_Pℓ] = [s[6 + ℓ_max + 1 + i] for s in sol.u, i in 0:ℓ_max]
    fields[:N_ℓ] = [s[6 + 2ℓ_max + 2 + i] for s in sol.u, i in 0:ℓ_max]
    fields
end #func