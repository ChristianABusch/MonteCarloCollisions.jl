using Interpolations

function collision_probability(σ, n, v, E, dt)
    # v : relative velocity in the center of mass system
    # E : E = 1/2 m v^2, where v is the relative velocity
    # n : is the density of the collision partner
    # dt: is the timestep taken

    ν = n * σ(E) * v 
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
    w⃗ = (m₁ * v⃗₁ + m₂ * v⃗₂) / (m₁ + m₂)
    return w⃗
end

function energy_in_collision(μ, g)
    E = 0.5 * μ * g^2
end

function absolute_value(v⃗)
    return √(sum(v⃗ .^2))
end

function Euler_angles(g⃗)
    gx, gy, gz = g⃗
    θ = atan(√(gy^2 + gz^2), gx)
    ϕ = atan(gz, gy)
    return θ, ϕ 
end

function elastic_collision(g⃗, w⃗, m₁, m₂)
    θ, ϕ = Euler_angles(g⃗)

    χ = acos(1.0 - 2.0 * rand()) # isotropic scattering
    η = 2π * rand()              # azimuthal angle

    # rotations of vector
    st = sin(θ)
    ct = cos(θ)
    sp = sin(ϕ)
    cp = cos(ϕ)
    sc = sin(χ)
    cc = cos(χ)
    se = sin(η)
    ce = cos(η)

    # compute new relative velocity:
    g  = absolute_value(g⃗)
    gx = g * (ct * cc - st * sc * ce)
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se)
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se)
    g⃗′ = [gx, gy, gz]

    return w⃗ .+ m₂/(m₁+m₂) .* g⃗′
end

function particle_collision(σ, n, dt, v⃗₁, m₁, v⃗₂, m₂)   
    
    g⃗ = relative_velocity(v⃗₁, v⃗₂)
    w⃗ = center_of_mass_velocity(v⃗₁, m₁, v⃗₂, m₂)
    μ = reduced_mass(m₁, m₂)
    g = absolute_value(g⃗)
    E = energy_in_collision(μ, g) / 1.602e-19
    p = collision_probability(σ, n, g, E, dt)

    if rand() < p
        elastic_collision(g⃗, w⃗, m₁, m₂)
    else
        return v⃗₁
    end
end



v⃗₁  = [1.7e6, 0.0, 0.0]
v⃗₂  = fill(100, 3)
m₁  = 9.109e-31
m₂  = 1.66e-27
n   = 1e24
dt  = 1e-12


cs_test = linear_interpolation([0.0, 1000.0], [1e-20, 1e-20])

@time particle_collision(cs_test, n, dt, v⃗₁, m₁, v⃗₂, m₂)
@time particle_collision(cs_test, n, dt, v⃗₁, m₁, v⃗₂, m₂)


N = 1000000
x = zeros(N)
@time for i in 2:N
    global v⃗₁
    x[i] = x[i-1] + v⃗₁[1] * dt
    v⃗₁ = particle_collision(cs_test, n, dt, v⃗₁, m₁, v⃗₂, m₂)
end

lines(x)
