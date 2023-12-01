using MonteCarloCollision

function energy_conservation_test()
    # electron
    m₁ = 9.109e-31
    v⃗₁ = [0.0, 3e6, 0.0] 

    # background gas atom
    m₂ = 1.66e-27
    v⃗₂ = fill(100, 3)

    E_exc = 10
    E_iz  = 12

    # calculate some relevant quantities 
    g⃗ = relative_velocity(v⃗₁, v⃗₂)
    w⃗ = center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)
    μ = reduced_mass(m₁, m₂)
    g = absolute_value(g⃗)
    E = to_eV(energy_in_collision(μ, g))

    # isotropic elastic scattering
    v⃗₁′ = isotropic_elastic_collision(g⃗, w⃗, m₁, m₂)
    E₁′ = to_eV(kinetic_energy(m₁, absolute_value(v⃗₁′)))
    ΔE_el_iso = E - E₁′

    # backscatter elastic scattering
    v⃗₁′ = backscatter_elastic_collision(g⃗, w⃗, m₁, m₂)
    E₁′ = to_eV(kinetic_energy(m₁, absolute_value(v⃗₁′)))
    ΔE_el_back = E - E₁′

    # excitation
    v⃗₁′ = excitation_collision(g⃗, w⃗, m₁, m₂, E_exc)
    E₁′ = to_eV(kinetic_energy(m₁, absolute_value(v⃗₁′)))
    ΔE_exc = E - (E₁′ + E_exc)

    # ionization
    v⃗₁′, v⃗_new = ionization_collision(g⃗, w⃗, m₁, m₂, E_iz, Ē[:Ar])
    E₁′   = to_eV(kinetic_energy(m₁, absolute_value(v⃗₁′)))
    E_new = to_eV(kinetic_energy(m₁, absolute_value(v⃗_new)))
    ΔE_iz = E - (E₁′ + E_new + E_iz)
    
    relative_energy_defects = abs.([ΔE_el_iso, ΔE_el_back, ΔE_exc, ΔE_iz]) ./E
    return all(relative_energy_defects .< 0.01)
end

energy_conservation_test()