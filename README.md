# MonteCarloCollision
## Overview
This package implements Monte-Carlo collisions for 
- isotropic elastic scattering
- backward scattering (charge exchange)
- excitation of an atom/molecule
- ionization

Explanations for the calculations can be found in "Zoltán Donkó et al 2021 Plasma Sources Sci. Technol. 30 095017".

# API
Let $\vec{v}_1$ and $m_1$ be the velocity and mass of the scattering particle (electron here) and  $\vec{v}_2$ and $m_2$ the velocity and mass of the target particle (atom here). $\vec{g}=\vec{v}_1 - \vec{v}_2$ is the relative velocity and $\vec{w}$ the center of mass velocity. 
The velocities of the electron after te different types of collisions can then be calculated as:
```julia
v⃗₁′ = isotropic_elastic_collision(g⃗, w⃗, m₁, m₂)
v⃗₁′ = backscatter_elastic_collision(g⃗, w⃗, m₁, m₂)
v⃗₁′ = excitation_collision(g⃗, w⃗, m₁, m₂, E_exc)
v⃗₁′, v⃗_new_electron = ionization_collision(g⃗, w⃗, m₁, m₂, E_iz, Ē[:Ar])
```

The package also includes some helpful functions when working with collisions. The collision probability for a collision of a particle with velocity $v$ with a gas with density $n$ and cross section $\sigma(v)$ over a time-interval $dt$ is
```julia
P = collision_probability(σ, n, v, dt)
```
When it was determined that a collision occured by evaluating the total cross section, a collision process can be determined by partitioning the interval $[0,1]$ into blocks of size $P_k$, where $P_k=\sigma_k / \sigma_{tot}$. Then a random number is selected and the process into whose block the random number lands is selected.
```julia
σ = [σ₁, σ₂, σ₃]
i_process = select_collision(σ)
```

