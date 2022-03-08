using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using AlgebraicPetri

using Catlab
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Programs


using Decapodes.Simulations
using Decapodes.Examples
using Decapodes.Diagrams
using Decapodes.Schedules
using Decapodes.OpenDiagrams
using Decapodes.PetriNets
include("deca_schema.jl")

draw_uwd(d) = to_graphviz(d, box_labels=:name, junction_labels=:variable,
            graph_attrs=Dict(:start => "2"))

draw_equations(eq) = begin
  to_graphviz(eq,
  node_labels=true, prog="neato",
  node_attrs=Dict(:shape=>"oval"),
  graph_attrs=Dict(:nodesep=>"4.0"))
end

@present Flow2DQuantities <: ExtendedFreeExtCalc2D begin
  X::Space

  k::Hom(DualForm1(X), DualForm1(X))    # diffusivity (usually scalar multiplication)
  k⁻¹::Hom(DualForm1(X), DualForm1(X))  # inverse of diffusivity
  neg₀::Hom(Form0(X), Form0(X))
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
end

#####################
# Example DECAPODES #
#####################


Diffusion = @decapode Flow2DQuantities begin
  T::Form0{X}
  ϕ::DualForm1{X}

  # Fick's first law
  ϕ ==  k(⋆₁{X}(d₀{X}(T)))
  # Diffusion equation
  ∂ₜ{Form0{X}}(T) == ⋆₀⁻¹{X}(dual_d₁{X}(ϕ))
end

open("diffusion.svg", "w") do io
  run_graphviz(io, draw_equations(Diffusion), format="svg")
end

UncompilableDiffusion = @decapode Flow2DQuantities begin
  T::Form0{X}
  ϕ::DualForm1{X}

  # Fick's first law
  k⁻¹(ϕ) ==  ⋆₁{X}(d₀{X}(T))
  # Diffusion equation
  ∂ₜ{Form0{X}}(T) == ⋆₀⁻¹{X}(dual_d₁{X}(ϕ))
end

open("uncomp_diffusion.svg", "w") do io
  run_graphviz(io, draw_equations(UncompilableDiffusion), format="svg")
end

##############################
# Physics Theory Composition #
##############################

Diffusion = @decapode Flow2DQuantities begin
  (T, Ṫ)::Form0{X}
  ϕ::DualForm1{X}

  # Fick's first law
  ϕ ==  k(⋆₁{X}(d₀{X}(T)))
  Ṫ == ⋆₀⁻¹{X}(dual_d₁{X}(ϕ))
end

Advection = @decapode Flow2DQuantities begin
  (T, Ṫ)::Form0{X}
  V::Form1{X}
  Ṫ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(T))))
end

Dynamics = @decapode Flow2DQuantities begin
  (X, Ẋ)::Form0{X}
  ∂ₜ{Form0{X}}(X) == Ẋ
end

dynam_physics = @relation (X, Ẋ) begin
  physics(X, Ẋ)
  dynamics(X, Ẋ)
end

DiffusionDynam = oapply(dynam_physics,
                  [OpenDiagram(Diffusion, [:T, :Ṫ]),
                   OpenDiagram(Dynamics, [:X, :Ẋ])])

AdvectionDynam = oapply(dynam_physics,
                   [OpenDiagram(Advection, [:T, :Ṫ]),
                    OpenDiagram(Dynamics, [:X, :Ẋ])])

draw_equations(Diffusion)

mkpath("simple_composition")

open("simple_composition/diffusion.svg", "w") do io
  run_graphviz(io, draw_equations(Diffusion), format="svg")
end

open("simple_composition/advection.svg", "w") do io
  run_graphviz(io, draw_equations(Advection), format="svg")
end

open("simple_composition/dynamics.svg", "w") do io
  run_graphviz(io, draw_equations(Dynamics), format="svg")
end

open("simple_composition/comp_pattern.svg", "w") do io
  run_graphviz(io, draw_uwd(dynam_physics), format="svg")
end

open("simple_composition/diff_dynam.svg", "w") do io
  run_graphviz(io, draw_equations(DiffusionDynam), format="svg")
end

open("simple_composition/adv_dynam.svg", "w") do io
  run_graphviz(io, draw_equations(AdvectionDynam), format="svg")
end

mkpath("multi_composition")

Superposition = @decapode Flow2DQuantities begin
  (Ẋ, Ẋ₁, Ẋ₂)::Form0{X}
  Ẋ == Ẋ₁ + Ẋ₂
end

two_phys = @relation (X, Ẋ₁, Ẋ₂) begin
  physics1(X, Ẋ₁)
  physics2(X, Ẋ₂)
  superposition(Ẋ, Ẋ₁, Ẋ₂)
end

AdvectionDiffusion = oapply(two_phys,
                   [OpenDiagram(Advection, [:T, :Ṫ]),
                    OpenDiagram(Diffusion, [:T, :Ṫ]),
                    OpenDiagram(Superposition, [:Ẋ, :Ẋ₁, :Ẋ₂])])

AdvectionDiffusionDynam = oapply(dynam_physics,
  [AdvectionDiffusion,
  OpenDiagram(Dynamics, [:X, :Ẋ])])

open("multi_composition/adv_diff_pattern.svg", "w") do io
  run_graphviz(io, draw_uwd(two_phys), format="svg")
end

open("multi_composition/advection_diffusion.svg", "w") do io
  run_graphviz(io, draw_equations(AdvectionDiffusion), format="svg")
end

open("multi_composition/superposition.svg", "w") do io
  run_graphviz(io, draw_equations(Superposition), format="svg")
end

open("multi_composition/comp_pattern.svg", "w") do io
  run_graphviz(io, draw_uwd(dynam_physics), format="svg")
end

open("multi_composition/diff_adv_dynam.svg", "w") do io
  run_graphviz(io, draw_equations(AdvectionDiffusionDynam), format="svg")
end
#=
enz_deg = LabelledReactionNet{Float64, Float64}([:AB=>0.0, :E=>0.0, :A=>0.0, :B=>0.0], (:deg=>0.05)=>((:AB,:E)=>(:A, :B, :E)))
Graph(enz_deg)

expand_pres!(Flow2DQuantities, enz_deg)
EnzymesDyn = pn2dec(Flow2DQuantities, enz_deg)
draw_equations(EnzymesDyn)

enz_phys = @relation (A, B, E, AB) begin
  physicsA(A, Ȧₚ)
  physicsB(B, Ḃₚ)
  physicsAB(AB, ∂ₜABₚ)
  enzyme(A, B, E, AB, Ȧₑ, Ḃₑ, ∂ₜABₑ)
  superposition(Ȧ, Ȧₑ, Ȧₚ)
  superposition(Ḃ, Ḃₑ, Ḃₚ)
  superposition(∂ₜAB, ∂ₜABₑ, ∂ₜABₚ)
end

draw_uwd(enz_phys)
AdvectionDiffusion = oapply(two_phys,
                   [AdvectionDiffusion,
                    AdvectionDiffusion,
                    AdvectionDiffusion,
                    OpenDiagram(EnzymesDyn, [:A, :B, :E, :AB, :Ḃ, :Ḃ, :∂ₜAB]),
                    OpenDiagram(Superposition, [:Ẋ, :Ẋ₁, :Ẋ₂]),
                    OpenDiagram(Superposition, [:Ẋ, :Ẋ₁, :Ẋ₂]),
                    OpenDiagram(Superposition, [:Ẋ, :Ẋ₁, :Ẋ₂])])

mkpath("chem_composition")=#

