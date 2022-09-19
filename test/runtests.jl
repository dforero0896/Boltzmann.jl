using Boltzmann
using Test

@testset "Boltzmann.jl" begin


    @show c = Boltzmann.Cosmology(h = 0.7, ω_c = 0.25 * (0.7)^2, ω_b = 0.05 * (0.7)^2, Neff = 0, T_cmb = 2.7255)

    sol = Boltzmann.solve_boltzmann(c, 1e-2, 3)
    @show sol.t
    @show [s[1] for s in sol.u]
    

    



end
