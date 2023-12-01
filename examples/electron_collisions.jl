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
    g = absolute_value(g⃗)
    μ = reduced_mass(m₁, m₂)   
    E = to_eV(kinetic_energy(μ, g⃗))

    P_tot = collision_probability(σ_tot(E), n, g, dt) # total probability of a collision

    if rand() < P_tot
        # select which process occurs
        i_process = select_collision([σ_ela(E), σ_exc(E), σ_iz(E)])

        if i_process == 1
            v⃗₁′ = isotropic_elastic_collision(v⃗₁, v⃗₂, m₁, m₂)
        elseif i_process == 2
            v⃗₁′ = excitation_collision(v⃗₁, v⃗₂, m₁, m₂, E_exc)
        elseif i_process == 3
            v⃗₁′, v⃗_new = ionization_collision(v⃗₁, v⃗₂, m₁, m₂, E_iz, Ē[:Ar])
        end
    end
end