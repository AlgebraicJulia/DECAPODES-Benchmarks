using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus

using Catlab
using Catlab.Graphics
using Catlab.Graphics.Graphviz


using Decapodes.Simulations
using Decapodes.Examples
using Decapodes.Diagrams
using Decapodes.Schedules
using Decapodes.OpenDiagrams


draw_equations(eq) = begin
  to_graphviz(eq,
  node_labels=true, prog="neato",
  node_attrs=Dict(:shape=>"oval"),
  graph_attrs=Dict(:nodesep=>"4.0"))
end

@present Flow2DQuantities(FreeExtCalc2D) begin
  X::Space

  k::Hom(DualForm1(X), DualForm1(X))    # diffusivity (usually scalar multiplication)
  k⁻¹::Hom(DualForm1(X), DualForm1(X))  # inverse of diffusivity
end

Diffusion = @decapode Flow2DQuantities begin
  (T, Ṫ)::Form0{X}
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
  (T, Ṫ)::Form0{X}
  ϕ::DualForm1{X}

  # Fick's first law
  k⁻¹(ϕ) ==  ⋆₁{X}(d₀{X}(T))
  # Diffusion equation
  ∂ₜ{Form0{X}}(T) == ⋆₀⁻¹{X}(dual_d₁{X}(ϕ))
end

open("uncomp_diffusion.svg", "w") do io
  run_graphviz(io, draw_equations(UncompilableDiffusion), format="svg")
end