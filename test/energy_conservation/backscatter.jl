using MonteCarloCollisions

function energy_conservation_test_backscatter()
    # electron
    m₁ = 9.109e-31
    v⃗₁ = [2e6, 2e6, 2e6] 
    E = kinetic_energy(m₁, v⃗₁)

    # background gas atom
    m₂ = 1.66e-27
    v⃗₂ = fill(100, 3)

    # backscatter elastic scattering
    v⃗₁′ = backscatter_elastic_collision(v⃗₁, v⃗₂, m₁, m₂)
    E₁′ = kinetic_energy(m₁, v⃗₁′)
    ΔE_el_back = E - E₁′

    return abs(ΔE_el_back) /E <= 0.01
end 

energy_conservation_test_backscatter()