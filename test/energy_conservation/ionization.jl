using MonteCarloCollisions

function energy_conservation_test_ionization()
    # electron
    m₁ = 9.109e-31
    v⃗₁ = [2e6, 2e6, 2e6] 
    E = kinetic_energy(m₁, v⃗₁)

    # background gas atom
    m₂ = 1.66e-27
    v⃗₂ = fill(100, 3)

    E_iz  = 12

    # ionization
    v⃗₁′, v⃗_new = ionization_collision(v⃗₁, v⃗₂, m₁, m₂, E_iz, Ē[:Ar])
    E₁′   = kinetic_energy(m₁, v⃗₁′)
    E_new = kinetic_energy(m₁, v⃗_new)
    ΔE_iz = E - (E₁′ + E_new + to_J(E_iz))

    return abs(ΔE_iz) /E <= 0.01
end 

energy_conservation_test_ionization()