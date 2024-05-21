using MonteCarloCollisions

function energy_conservation_test_excitation()
    # electron
    m₁ = 9.109e-31
    v⃗₁ = [2e6, 2e6, 2e6] 
    E = kinetic_energy(m₁, v⃗₁)

    # background gas atom
    m₂ = 1.66e-27
    v⃗₂ = fill(100, 3)

    E_exc = 10

    # excitation
    v⃗₁′ = excitation_collision(v⃗₁, v⃗₂, m₁, m₂, E_exc)
    E₁′ = kinetic_energy(m₁, v⃗₁′)
    ΔE_exc = E - (E₁′ + to_J(E_exc))

    return abs(ΔE_exc) /E <= 0.01
end 

energy_conservation_test_excitation()