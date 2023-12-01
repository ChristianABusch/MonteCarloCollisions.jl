using MonteCarloCollision
using Interpolations

function electron_collisions()

    # electron
    m₁  = 9.109e-31
    v⃗₁  = [0.0, 3e6, 0.0] 

    # background gas atom
    m₂  = 1.66e-27
    v⃗₂  = fill(100, 3)
    n   = 1e24   # gas density
    dt  = 1e-12  # timestep

    # cross sections and energy threshhholds
    E_exc = 10
    E_iz  = 12
    δ     = 0.1
    σ_ela = linear_interpolation([0.0, 1000], [1e-19, 1e-19])
    σ_exc = linear_interpolation([0.0, E_iz-δ, E_iz, 1000], [0.0, 0.0, 1e-19, 1e-19])
    σ_iz  = linear_interpolation([0.0, E_exc-δ, E_exc, 1000], [0.0, 0.0, 1e-19, 1e-19])
    σ_tot(E) = σ_ela(E) + σ_exc(E) + σ_iz(E)

    # check for collisions and execute them if one occurs    
    g⃗ = relative_velocity(v⃗₁, v⃗₂)
    w⃗ = center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)
    μ = reduced_mass(m₁, m₂)
    g = absolute_value(g⃗)
    E = to_eV(energy_in_collision(μ, g))

    P_tot = collision_probability(σ_tot(E), n, g, E, dt) # total probability of a collision

    if rand() < P_tot
        # select which process occurs
        r = rand()

        # segment the interval [0,1] according to the probability of a collision occuring
        P_ell = σ_ela(E)/σ_tot(E)
        P_exc = (σ_exc(E) + σ_ela(E))/σ_tot(E)
        P_iz  = (σ_exc(E) + σ_ela(E) + σ_iz(E))/σ_tot(E)

        if r < P_ell
            v⃗₁′ = isotropic_elastic_collision(g⃗, w⃗, m₁, m₂)

            E  = round(E, digits = 2)
            E′ = round(to_eV(kinetic_energy(m₁, absolute_value(v⃗₁′))), digits = 2)
            println("elastic collision occured. Initial energy : $(E), final energy: $(E′)")
        elseif r < P_exc
            v⃗₁′ = excitation_collision(g⃗, w⃗, m₁, m₂, E_exc)

            E  = round(E, digits = 2)
            E′ = round(to_eV(kinetic_energy(m₁, absolute_value(v⃗₁′))), digits = 2)
            println("excitation collision occured. Initial energy : $(E), final energy: $(E′), ΔE = $(E-E′)")
        elseif r < P_iz
            v⃗₁′, v⃗_new = ionization_collision(g⃗, w⃗, m₁, m₂, E_iz, Ē[:Ar])

            E  = round(E, digits = 2)
            E′ = round(to_eV(kinetic_energy(m₁, absolute_value(v⃗₁′))), digits = 2)
            println("ionization collision occured. Initial energy : $(E), final energy: $(E′), ΔE = $(E-E′)")
        end
    end
end