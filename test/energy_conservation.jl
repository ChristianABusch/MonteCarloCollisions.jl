using MonteCarloCollisions

function energy_conservation_test()
    # electron
    m₁ = 9.109e-31
    v⃗₁ = [0.0, 3e6, 0.0] 
    E = kinetic_energy(m₁, v⃗₁)

    # background gas atom
    m₂ = 1.66e-27
    v⃗₂ = fill(100, 3)

    E_exc = 10
    E_iz  = 12

    # isotropic elastic scattering
    v⃗₁′ = isotropic_elastic_collision(v⃗₁, v⃗₂, m₁, m₂)
    E₁′ = kinetic_energy(m₁, v⃗₁′)
    ΔE_el_iso = E - E₁′

    # backscatter elastic scattering
    v⃗₁′ = backscatter_elastic_collision(v⃗₁, v⃗₂, m₁, m₂)
    E₁′ = kinetic_energy(m₁, v⃗₁′)
    ΔE_el_back = E - E₁′

    # excitation
    v⃗₁′ = excitation_collision(v⃗₁, v⃗₂, m₁, m₂, E_exc)
    E₁′ = kinetic_energy(m₁, v⃗₁′)
    ΔE_exc = E - (E₁′ + to_J(E_exc))

    # ionization
    v⃗₁′, v⃗_new = ionization_collision(v⃗₁, v⃗₂, m₁, m₂, E_iz, Ē[:Ar])
    E₁′   = kinetic_energy(m₁, v⃗₁′)
    E_new = kinetic_energy(m₁, v⃗_new)
    ΔE_iz = E - (E₁′ + E_new + to_J(E_iz))

    relative_energy_defects = abs.([ΔE_el_iso, ΔE_el_back, ΔE_exc, ΔE_iz]) ./E
    return all(relative_energy_defects .< 0.01)
end

energy_conservation_test()