using Boltzmann
using Test

@testset "Boltzmann.jl" begin


    @show c = Boltzmann.Cosmology(h = 0.7, ω_c = (0.3 - 0.06) * (0.7)^2, ω_b = 0.06 * (0.7)^2, Neff = 3.044, T_cmb = 2.725)
    a = 0.5
    @test Boltzmann.comoving_distance(c, 1 / a - 1) ≈ 3303.5262800468654 rtol=1e-3
    @test Boltzmann.η(c, 1 / a - 1) ≈ 7413.48949805 rtol=1e-1
    @test Boltzmann.H(c, 1 / a - 1) ≈ 123.27314813645114 rtol=1e-4
    k = 0.2
    ℓ_max = 3
    
    #sol = Boltzmann.solve_boltzmann(c, 0.2, 3, 10 .^range(-6, stop=0, length=500))
    #@show sol.t[1:20:end]
    #@show exp.(sol.t[1:20:end])
    #@show [s[1] for s in sol.u][1:20:end]
    #@show [s[2] for s in sol.u][1:20:end]
    #k, δ = Boltzmann.solve_boltzmann(c, 10 .^(-4:0.1:0), 3)
    #@show δ
    
    

    



end
