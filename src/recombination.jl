


function precompute_recfast(c::Cosmology)

    params = Recfast.Params(Yp = c.Yp,
                    T0 = c.T_cmb,
                    Omega_M = c.Ω_m₀,
                    Omega_B = c.Ω_b₀,
                    Omega_K = c.Ω_k₀,
                    h100 = c.h,
                    n_eff = c.Neff,
                    )

    sol = Evaluate_recombination(params, logzstart = 6., logzend = -5,dt = 1e-2, dtmin = 1e-30)
    x_He = [s[1] for s in sol.u]
    x_H = [s[2] for s in sol.u]
    T_m = [s[3] for s in sol.u]
    x_e = x_H .+ x_He
    
    return sol.t, x_He, x_H, T_m, x_e
end #func

function interpolate_recfast(z::Vector{T}, x_He::Vector{T}, x_H::Vector{T}, T_m::Vector{T}, x_e::Vector{T}) where T <: AbstractFloat

    a = 1 ./ (1 .+ z)
    loga = log.(a)
    x_He_i = extrapolate(interpolate((loga,), x_He, Gridded(Linear())), Flat())
    x_H_i = extrapolate(interpolate((loga,), x_H, Gridded(Linear())), Flat())
    T_m_i = extrapolate(interpolate((loga,), T_m, Gridded(Linear())), Line())
    x_e_i = extrapolate(interpolate((loga,), x_e, Gridded(Linear())), Flat())
    return x_He_i, x_H_i, T_m_i, x_e_i

end #func


function recombination_functions(c::Cosmology, z::Vector{T}, x_He::Vector{T}, x_H::Vector{T}, T_m::Vector{T}, x_e::Vector{T}; autodiff::Bool=false) where {T<:AbstractFloat}
    # https://www.uio.no/studier/emner/matnat/astro/AST5220/v20/forelesningsvideoer/project_mk2.pdf
    x_He_i, x_H_i, T_m_i, x_e_i = interpolate_recfast(z, x_He, x_H, T_m, x_e)
    loga = log.(1 ./ (1 .+ z))
    min_a = minimum(exp.(loga))
    if autodiff
        T_m′ = loga -> Zygote.gradient(x -> log(T_m_i(x)), loga)[1]
    else
        T′ = [(log(T_m[i + 1]) - log(T_m[i])) / (loga[i+1] - loga[i]) for i in 1:(length(loga)-1)]
        push!(T′, -2)
        loga_T = [0.5 * (loga[i+1] + loga[i]) for i in 1:(length(loga)-1)]
        push!(loga_T, 0.)
        T_m′ = extrapolate(interpolate((loga_T,), T′, Gridded(Linear())), Flat())
    end #if
    # x_e is taken to be x_H + x_He
    τ_dot_norm = σT * (c.ρ_crit * msun / mpc^3)  / (mp * 1e6 * evc²)
    τ_dot(c, a) = -τ_dot_norm * mpc / a^2 * x_e_i(log(a)) * c.Ω_b₀ * c.h * (1.0 - c.Yp) # [h Mpc^-1] (see Dodelson eq. 3.44, p73) + 1-Y factor from COSMICS) 

    Tm(c, a) = a < min_a ? c.T_cmb / a : T_m_i(log(a))
    μ(c, a) = 1.0 / (1.0 + x_e_i(log(a))) * (1.0 - c.Yp) + 0.25 * (1.0 + 2.0 * x_He_i(log(a))) * c.Yp
    cₛ(c, a) = sqrt(kb * Tm(c,a) / μ(c,a) / (mp * 1e6) * (1.0 - T_m′(log(a)) / 3.0))
    
    return cₛ, τ_dot
end #func