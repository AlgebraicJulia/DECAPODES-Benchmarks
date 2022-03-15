##################
# Define Physics #
##################

@present Flow2DQuantities <: ExtendedFreeExtCalc2D begin
  X::Space

  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  kᵨ::Hom(Form1(X), Form1(X))   # REMOVE: Diffusion of density
  R₀::Hom(Form0(X), Form0(X))    # Ideal gas constant (usually scalar multiplication)
  kᵥ::Hom(Form1(X), Form1(X))    # viscosity (usually scalar multiplication)
  e2t::Hom(Form0(X), Form0(X))    # scale between internal energy and temperature (ideal gas)
  t2e::Hom(Form0(X), Form0(X))    # inverse of e2t (ideal gas)
  e2p::Hom(Form0(X), Form0(X))    # # scale between internal energy × density and pressure (ideal gas)
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
  L₁′::Hom(Form1(X)⊗Form1(X), Form1(X))
  ∂₀::Hom(Form0(X), Form0(X)) # all_wall boundary condition
  ∂₀₊::Hom(Form0(X), Form0(X)) # top/bottom wall boundary condition
  ∂₀ₛ::Hom(Form0(X), Form0(X)) # Sticking boundary condition
  ∂₁ₛ::Hom(Form1(X), Form1(X)) # Sticking boundary condition
  ∂₁ₗ₊::Hom(Form1(X), Form1(X)) # Sticking boundary condition
  ∂₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  ∂ₜₐ::Hom(Form0(X), Form0(X))  # Temperature advection boundary condition
  ∂ₜ::Hom(Form0(X), Form0(X))  # internal energy boundary condition
  ∂ₚ::Hom(Form0(X), Form0(X))  # Pressure boundary condition
  ∂ₚₐ::Hom(Form0(X), Form0(X))  # Pressure advection boundary condition
  ∂ᵥ::Hom(Form1(X), Form1(X))  # Velocity boundary condition
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # In/Out edge flow boundary condition
  neg₁::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  dneg₀::Hom(DualForm0(X), DualForm0(X)) # In/Out edge flow boundary condition
  neg₀::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  mask₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  mask₀ₗ::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  mask₀ₛ::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  half::Hom(Form0(X), Form0(X)) # half
  third::Hom(Form1(X), Form1(X))    # multiplication by 1/3
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  mult₀::Hom(Form0(X)⊗Form0(X), Form0(X)) # elementwise multiplication
  div₀::Hom(Form0(X)⊗Form0(X), Form0(X)) # elementwise division
  div₁::Hom(Form1(X)⊗Form1(X), Form1(X)) # elementwise division
  avg₀₁::Hom(Form0(X), Form1(X))
  sub₀::Hom(Form0(X)⊗Form0(X), Form0(X))
end

Diffusion = @decapode Flow2DQuantities begin
  (T, Ṫ)::Form0{X}
  ϕ::DualForm1{X}

  # Fick's first law
  ϕ ==  ⋆₁{X}(k(d₀{X}(T)))
  # Diffusion equation
  Ṫ == ⋆₀⁻¹{X}(dual_d₁{X}(ϕ))
end

Advection = @decapode Flow2DQuantities begin
  (T, Ṫ)::Form0{X}
  V::Form1{X}
  Ṫ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(T))))
end

Superposition = @decapode Flow2DQuantities begin
  (Ṫ₁, Ṫ₂, Ṫ, T)::Form0{X}
  Ṫ == Ṫ₁ + Ṫ₂
  ∂ₜ{Form0{X}}(T) == Ṫ
end

compose_diff_adv = @relation (T, Ṫ, V) begin
  diffusion(T, Ṫ₁)
  advection(T, Ṫ₂, V)
  superposition(Ṫ₁, Ṫ₂, Ṫ, T)
end
#=
DiffusionAdvection = oapply(compose_diff_adv,
                  [OpenDiagram(Diffusion, [:T, :Ṫ]),
                   OpenDiagram(Advection, [:T, :Ṫ, :V]),
                   OpenDiagram(Superposition, [:Ṫ₁, :Ṫ₂, :Ṫ, :T])])


dwd_res = diag2dwd(diff_adv.functor, in_vars = [:V, :T])
draw_equations(diff_adv.functor)
draw_dwd(dwd_res)
=#
NavierStokes = @decapode Flow2DQuantities begin
  (V, V̇, G)::Form1{X}
  (T, ρ, ṗ, p)::Form0{X}
  V̇ == neg₁(L₁′(V, V)) + 
        div₁(kᵥ(Δ₁{X}(V) + third(d₀{X}(δ₁{X}(V)))), avg₀₁(ρ)) +
        d₀{X}(half(i₁′(V, V))) +
        neg₁(div₁(d₀{X}(p),avg₀₁(ρ))) +
        G
  ∂ₜ{Form1{X}}(V) == V̇
  ṗ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(p))))# + ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(kᵨ(d₀{X}(ρ)))))
  ∂ₜ{Form0{X}}(p) == ṗ
end

Energy = @decapode Flow2DQuantities begin
  (V)::Form1{X}
  (ρ, p, T, Ṫ, Ṫₐ, Ṫ₁, bc₀)::Form0{X}

  ρ == div₀(p, R₀(T))
  Ṫₐ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(T))))
  ∂ₜₐ(Ṫₐ) == bc₀
  Ṫ == Ṫₐ + Ṫ₁
  ∂ₜ{Form0{X}}(T) == Ṫ
end

BoundaryConditions = @decapode Flow2DQuantities begin 
  (V, V̇, bc₁)::Form1{X}
  (Ṫ, ṗ, bc₀)::Form0{X} 
  # no-slip edges
  ∂ᵥ(V̇) == bc₁
  # No change on left/right boundaries
  ∂ₜ(Ṫ) == bc₀
  ∂ₚ(ṗ) == bc₀
end

dwd_res = diag2dwd(NavierStokes)
draw_equations(NavierStokes)
draw_dwd(dwd_res)

compose_heat_xfer = @relation (V, ρ) begin
  flow(V, V̇, T, ρ, ṗ, p)
  energy(Ṫ, V, ρ, p, T, Ṫ₁)
  diffusion(T, Ṫ₁)
  bcs(Ṫ, ṗ, V, V̇)
end

HeatXfer = oapply(compose_heat_xfer,
                  [OpenDiagram(NavierStokes, [:V, :V̇, :T, :ρ, :ṗ, :p]),
                   OpenDiagram(Energy, [:Ṫ, :V, :ρ, :p, :T, :Ṫ₁]),
                   OpenDiagram(Diffusion, [:T, :Ṫ]),
                   OpenDiagram(BoundaryConditions, [:Ṫ, :ṗ, :V, :V̇])])

draw_equations(HeatXfer.functor)

new_dwd = diag2dwd(HeatXfer.functor, in_vars=[:T, :p, :V, :G])
draw_dwd(new_dwd)
