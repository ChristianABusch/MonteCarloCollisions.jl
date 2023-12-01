function collision_probability(σ, n, v, dt)
    # v : relative velocity in the center of mass system
    # E : E = 1/2 m v^2, where v is the relative velocity
    # n : is the density of the collision partner
    # dt: is the timestep taken

    ν = n * σ * v 
    p = 1 - exp(-ν * dt)
    return p
end

function reduced_mass(m₁, m₂)
    μ = m₁ * m₂ / (m₁ + m₂)
    return μ
end

function relative_velocity(v⃗₁, v⃗₂)
    g⃗ = v⃗₁ .- v⃗₂
    return g⃗
end

function center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)
    w⃗ = @. (m₁ * v⃗₁ + m₂ * v⃗₂) / (m₁ + m₂)
    return w⃗
end

function kinetic_energy(m, v)
    # m: mass
    # v: velocity
    return 1/2 * m * v^2
end

function kinetic_energy(m, v⃗::Vector{Float64})
    # m: mass
    # v: velocity
    return 1/2 * m * sum(v⃗.^2)
end

function absolute_value(v⃗)
    return √(sum(v⃗ .^2))
end

function relative_kinetic_energy(v⃗₁, m₁, v⃗₂, m₂)
    μ  = reduced_mass(m₁, m₂)
    g⃗  = relative_velocity(v⃗₁, v⃗₂)
    g  = absolute_value(g⃗′)
    E  = energy_in_collision(μ, g)
    return E
end

function partition_one_by_cross_sections(σ::Vector{Float64})
    # Partition the interval [0,1] into block of size Pₖ where Pₖ=σₖ/σₜₒₜ
    σₜₒₜ = sum(σ)
    partition = [sum(σ[1:i])/σₜₒₜ for i in 1:length(σ)]
    return partition
end

function select_collision(σ::Vector{Float64})
    # Partition the interval [0,1] into block of size Pₖ where Pₖ=σₖ/σₜₒₜ
    partition = partition_one_by_cross_sections(σ)
    # draw random number and select the collision process of the intervall it lands in
    r = rand()
    process_index = findfirst(Pₖ -> r < Pₖ, partition)
    return process_index
end

function to_eV(E)
    return E/1.602e-19
end

function to_J(E)
    return E*1.602e-19
end