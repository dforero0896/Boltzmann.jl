using Boltzmann
using Test

@testset "Boltzmann.jl" begin


    @show c = Boltzmann.Cosmology(h = 0.7, ω_c = 0.25 * (0.7)^2, ω_b = 0.05 * (0.7)^2, Neff = 0, T_cmb = 2.7255)
    @show Boltzmann.H(c, 1)
    @show Boltzmann.comoving_distance(c, 1.)
    @show Boltzmann.ρx_over_ωx
    @show y, loga0 = Boltzmann.boltzmann_ics(c, 1e-2, 3)
    dy = similar(y)
    fill!(dy, 0.)
    Boltzmann.boltzmann_derivatives!(dy, y, (1e-2, c, 3), loga0)
    
    Boltzmann.solve_boltzmann(c, 1e-2, 20)



end
