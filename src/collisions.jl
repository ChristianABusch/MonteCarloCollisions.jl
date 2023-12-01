function Euler_angles(g⃗)
    gx, gy, gz = g⃗
    θ = atan(√(gy^2 + gz^2), gx)
    ϕ = atan(gz, gy)
    return θ, ϕ 
end

function post_collision_rotations(g ,θ, ϕ, χ, η)
    # θ, ϕ: Euler angles
    # χ: scattering angle
    # π: azimuthal angle

    st = sin(θ)
    ct = cos(θ)
    sp = sin(ϕ)
    cp = cos(ϕ)
    sc = sin(χ)
    cc = cos(χ)
    se = sin(η)
    ce = cos(η)

    # compute new relative velocity:
    gx = g * (ct * cc - st * sc * ce)
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se)
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se)
    g⃗′ = [gx, gy, gz]
    return g⃗′
end

function isotropic_elastic_collision(v⃗₁, v⃗₂, m₁, m₂)
    # m₁: mass of the scattered particle
    # m₂: mass of the target particle
    # g⃗: relative velocity vector
    # w⃗: center of mass velocity vector

    # relative and center of mass velocity
    g⃗ = relative_velocity(v⃗₁, v⃗₂)
    w⃗ = center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)

    # Euler angles
    θ, ϕ = Euler_angles(g⃗)

    # scattering angles
    χ = acos(1.0 - 2.0 * rand()) # isotropic scattering
    η = 2π * rand()              # azimuthal angle

    # compute new relative velocity:
    g  = absolute_value(g⃗)
    g⃗′ = post_collision_rotations(g ,θ, ϕ, χ, η)

    return w⃗ .+ m₂/(m₁+m₂) .* g⃗′
end

function backscatter_elastic_collision(v⃗₁, v⃗₂, m₁, m₂)
    # m₁: mass of the scattered particle
    # m₂: mass of the target particle
    # g⃗: relative velocity vector
    # w⃗: center of mass velocity vector

    # relative and center of mass velocity
    g⃗ = relative_velocity(v⃗₁, v⃗₂)
    w⃗ = center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)

    # Euler angles
    θ, ϕ = Euler_angles(g⃗)

    # scattering angles
    χ = π           # backscattering
    η = 2π * rand() # azimuthal angle

    # compute new relative velocity:
    g  = absolute_value(g⃗)
    g⃗′ = post_collision_rotations(g ,θ, ϕ, χ, η)

    return w⃗ .+ m₂/(m₁+m₂) .* g⃗′
end

function excitation_collision(v⃗₁, v⃗₂, m₁, m₂, E_exc)
    # m₁: electron mass
    # m₂: mass of the target particle
    # v⃗₁: velcoiity of electron
    # v⃗₂: velocity of collision partner
    # E_exc: excitation energy in eV

    # relative and center of mass velocity and reduced mass
    g⃗ = relative_velocity(v⃗₁, v⃗₂)
    w⃗ = center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)
    g = absolute_value(g⃗)
    μ = reduced_mass(m₁, m₂)

    # Euler angles
    θ, ϕ = Euler_angles(g⃗)

    # Calculate new energy
    E = kinetic_energy(μ, g)     # electron energy
    E = E - to_J(E_exc)          # subtract energy loss for excitation
    
    # scattering angles
    g = √(2.0 * E / μ)           # relative velocity after energy loss
    χ = acos(1.0 - 2.0 * rand()) # isotropic scattering
    η = 2π * rand()              # azimuthal angle

    # compute new velocities
    g⃗′ = post_collision_rotations(g ,θ, ϕ, χ, η)
    return w⃗ .+ m₂/(m₁+m₂) .* g⃗′
end

function ionization_collision(v⃗₁, v⃗₂, m₁, m₂, E_iz, Ē)
    # m₁: electron mass
    # m₂: mass of the target particle
    # v⃗₁: velcoiity of electron
    # v⃗₂: velocity of collision partner
    # E_iz: ionization energy in eV
    # Ē: electron spectrum shape parameter for ionzation (J. Chem. Phys. 55, 4100–4106 (1971))

    # relative and center of mass velocity and reduced mass
    g⃗ = relative_velocity(v⃗₁, v⃗₂)
    w⃗ = center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)
    g = absolute_value(g⃗)
    μ = reduced_mass(m₁, m₂)

    # Euler angles
    θ, ϕ = Euler_angles(g⃗)

    # Subtract ionization energy 
    E = kinetic_energy(μ, g)      # electron energy
    E = E - to_J(E_iz)            # subtract energy loss for ionization
    
    # split remaining energy up between the original electron and the newly freed one
    E_new  = to_J(10.0 * tan(rand() * atan(to_eV(E)/(2*Ē)))) # energy of new electron (J. Chem. Phys. 55, 4100–4106 (1971))
    E_orig = E - E_new                                       # energy of old electron
    
    # relative velocities of the two electrons
    g     = √(2.0 * E_orig / m₁)  # relative velocity of incoming (original) electron
    g_new = √(2.0 * E_new  / m₁)  # relative velocity of emitted (new) electron
    
    # scattering angles
    χ     = acos(√(E_orig / E))   # scattering angle for incoming electron
    χ_new = acos(√(E_new  / E))   # scattering angle for emitted electron
    η     = 2π * rand()           # azimuthal angle for incoming electron
    η_new = η + π                 # azimuthal angle for emitted electron

    # compute new velocity of original electron
    g⃗′ = post_collision_rotations(g ,θ, ϕ, χ, η)
    v⃗′ = w⃗ .+ m₂/(m₁+m₂) .* g⃗′

    # compute new velocity of new electron
    g⃗′_new = post_collision_rotations(g_new ,θ, ϕ, χ_new, η_new)
    v⃗′_new = w⃗ .+ m₂/(m₁+m₂) .* g⃗′_new
    return v⃗′, v⃗′_new
end