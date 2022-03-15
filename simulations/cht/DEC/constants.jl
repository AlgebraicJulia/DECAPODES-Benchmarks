cₚ = 1004.703 # Specific Heat at constant pressure
kₜ = 0.0246295028571 #Thermal conductivity
k_cyl = kₜ * 4
kᵦ = 1.38064852e-23 # Boltzmann constant (m² kg/(s² K))

density = 0.000210322
R = kᵦ * 6.0221409e23 # kg⋅m²/(s²*K*mol)
mol_mass = 28.96 # g/mol
μ = 28.96 / 6.0221409e23 # mean molecular mass (g)
R₀ = R / (mol_mass / 1000)

# Heat diffusion constant in fluid
k₁ = kₜ / (density * cₚ)

# Heat diffusion constant in cylinder
k₂ = k_cyl / (density * cₚ)

kᵨ = 1e-3

ν = 1.716e-5 # 0.0005081150545826

kᵥ = ν# / density

γ = 1 + R₀ / (cₚ - R₀) # Mayer's formula + definition of adiabatic constant
println("γ: $γ")
# T = (γ - 1) * (mol_mass / 1000) / R₀ * e 
#e2t = (μ/1000) / kᵦ
#e2t = (γ - 1) * mol_mass / R₀
e2t = 1/(cₚ * density)
@show e2t
e2p = e2t * R₀
#e2p = (γ - 1)
t2e = 1/e2t
