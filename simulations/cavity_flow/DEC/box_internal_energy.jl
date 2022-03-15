# Our developed libraries
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus

using Catlab
using Catlab.Present
using Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.Programs
using LinearAlgebra

# Julia community libraries
using MeshIO
using CairoMakie
using DifferentialEquations
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using JSON
global_logger(TerminalLogger())

using Decapodes.Simulations
using Decapodes.Examples
using Decapodes.Diagrams
using Decapodes.Schedules
using Decapodes.Debug
using Decapodes.OpenDiagrams

using CombinatorialSpaces: volume
using SparseArrays

######################
# Functions for Curl #
######################
include("./special_ops.jl")
include("./deca_schema.jl")

draw_dwd(dwd) = begin
  new_res = deepcopy(dwd)
  new_res.diagram[:outer_out_port_type] .= nothing
  new_res.diagram[:outer_in_port_type] .= nothing
  new_res.diagram[:out_port_type] .= nothing
  new_res.diagram[:in_port_type] .= nothing
  to_graphviz(new_res)
end
draw_equations(eq) = begin
  to_graphviz(eq,
  node_labels=true, prog="neato",
  node_attrs=Dict(:shape=>"oval"),
  graph_attrs=Dict(:nodesep=>"4.0"))
end

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
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # In/Out edge flow boundary condition
  neg₁::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  dneg₀::Hom(DualForm0(X), DualForm0(X)) # In/Out edge flow boundary condition
  neg₀::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  mask₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  mask₀ₗ::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  mask₀ₛ::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  half::Hom(Form0(X), Form0(X)) # half
  R₀⁻¹::Hom(Form0(X), Form0(X)) # division by R₀
  third::Hom(Form1(X), Form1(X))    # multiplication by 1/3
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  mult₀::Hom(Form0(X)⊗Form0(X), Form0(X)) # elementwise multiplication
  div₁::Hom(Form1(X)⊗Form1(X), Form1(X)) # elementwise division
  div₀::Hom(Form0(X)⊗Form0(X), Form0(X)) # elementwise division
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
  (V, V̇, G, τ)::Form1{X}
  (T, ρ, ρ̇, p, δV)::Form0{X}
  δV == δ₁{X}(V)
  τ == div₁(kᵥ(Δ₁{X}(V) + third(d₀{X}(δV))), avg₀₁(ρ))
  V̇ == neg₁(L₁′(V, V)) + 
        d₀{X}(half(i₁′(V, V))) +
        τ +
        neg₁(div₁(d₀{X}(p),avg₀₁(ρ))) +
        G
  ∂ₜ{Form1{X}}(V) == V̇
  ρ̇ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(ρ))))# + ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(kᵨ(d₀{X}(ρ)))))
  ∂ₜ{Form0{X}}(ρ) == ρ̇
end

Energy = @decapode Flow2DQuantities begin
  (V, τ)::Form1{X}
  (e, ėₐ, ėₜ, ė, ρ, p, δV, T, Ṫ)::Form0{X}

  T == R₀⁻¹(div₀(p, ρ))
  ėₜ == Ṫ + i₁′(τ, V)

  p == e2p(mult₀((e + neg₀(half(i₁′(V, V)))) , ρ))

  ėₐ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(mult₀(e,ρ) + p))))
  ė == div₀(ėₜ + ėₐ, ρ)
  ∂ₜ{Form0{X}}(e) == ė
end

BoundaryConditions = @decapode Flow2DQuantities begin 
  (V, V̇, bc₁)::Form1{X}
  (ė, bc₀)::Form0{X} 
  # no-slip edges
  ∂₁ₛ(V) == bc₁
  ∂₁ₛ(V̇) == bc₁
  # No change on left/right boundaries
  ∂₀ₛ(ė) == bc₀
end

dwd_res = diag2dwd(NavierStokes)
draw_equations(NavierStokes)
draw_dwd(dwd_res)

compose_heat_xfer = @relation (V, ρ) begin
  flow(V, V̇, δV, T, ρ, p, τ)
  energy(ė, V, ρ, p, T, Ṫ, δV, τ)
  diffusion(T, Ṫ)
  bcs(ė, V, V̇)
end

HeatXfer = oapply(compose_heat_xfer,
                  [OpenDiagram(NavierStokes, [:V, :V̇, :δV, :T, :ρ, :p, :τ]),
                   OpenDiagram(Energy, [:ė, :V, :ρ, :p, :T, :Ṫ, :δV, :τ]),
                   OpenDiagram(Diffusion, [:T, :Ṫ]),
                   OpenDiagram(BoundaryConditions, [:ė, :V, :V̇])])

draw_equations(HeatXfer.functor)

new_dwd = diag2dwd(HeatXfer.functor, in_vars=[:e, :ρ, :V, :G])
draw_dwd(new_dwd)

                   

################################
# Complex Operator Definitions #
################################
@present ExtendedOperators <: ExtendedFreeExtCalc2D begin
  X::Space
  neg::Hom(DualForm1(X), DualForm1(X)) # negative
  neg₁::Hom(Form1(X), Form1(X)) # negates 1-forms
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # negates 1-forms
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
  i₀::Hom(Form1(X)⊗DualForm2(X), DualForm1(X))
  i₁::Hom(Form1(X)⊗DualForm1(X), DualForm0(X))
    
  L₁′::Hom(Form1(X)⊗Form1(X), Form1(X))
  i₀′::Hom(Form1(X)⊗Form2(X), Form1(X))
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  ∧₁₁′::Hom(Form1(X)⊗DualForm1(X), DualForm2(X))
  ∧₁₀′::Hom(Form1(X)⊗DualForm0(X), Form1(X))
end

rules = gen_dec_rules()

# Lie derivative between two primal 1-forms
Lie1Imp′ = @decapode ExtendedOperators begin
  (F1, F1′, Ḟ1)::Form1{X}
  Ḟ1 == i₀′(F1, d₁{X}(F1′)) + d₀{X}(i₁′(F1, F1′))
end
lie1_imp′ = diag2dwd(Lie1Imp′, in_vars = [:F1, :F1′], out_vars = [:Ḟ1])
rules[:L₁′] = lie1_imp′

# Internal product between a primal 1-form and a primal 2-form
I0Imp′ = @decapode ExtendedOperators begin
  (F1, Ḟ1)::Form1{X}
  F2::Form2{X}
  Ḟ1 == neg₁(∧₁₀′(F1, ⋆₂{X}(F2)))
end
i0_imp′ = diag2dwd(I0Imp′, in_vars = [:F1, :F2], out_vars = [:Ḟ1])
rules[:i₀′] = i0_imp′

# Internal product between two primal 1-forms
I1Imp′ = @decapode ExtendedOperators begin
  (F1, F1′)::Form1{X}
  F0::Form0{X}
  F0 == ⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
end
i1_imp′ = diag2dwd(I1Imp′, in_vars = [:F1, :F1′], out_vars = [:F0])
rules[:i₁′] = i1_imp′

##########################
# Simulation Computation #
##########################

exp_dwd = Examples.expand_dwd(new_dwd, rules)
draw_dwd(exp_dwd)

s1 = EmbeddedDeltaSet2D("../cavity_unstructured_64.stl")
s = EmbeddedDeltaSet2D{Bool, Point{3, Float64}}()
copy_parts!(s, s1)

mod_p = findall(p -> (0.71 <= p[1] <= 0.72) && (0.99 <= p[2] <= 0.998), s[:point])
s[mod_p[1], :point] = s[mod_p[1], :point] .- (0.0,0.001,0.0)


sd = dual(s);
if ⋆(Val{1}, sd)[1,1] < 0.0
  orient_component!(s, 1, false)
  sd = dual(s); 
end

cₚ = 1004.703 # Specific Heat at constant pressure
kₜ = 10 * 0.0246295028571 #Thermal conductivity
kᵦ = 1.38064852e-23 # Boltzmann constant (m² kg/(s² K))

density = 0.00597782417156
R = kᵦ * 6.0221409e23 # kg⋅m²/(s²*K*mol)
mol_mass = 28.96 # g/mol
μ = 28.96 / 6.0221409e23 # mean molecular mass (g)
R₀ = R / (mol_mass / 1000)

# Heat diffusion constant in fluid
k₁ = kₜ
k_col = fill(k₁, ne(s))
k = diagm(k_col)

kᵨ = 1e-3

ν = 1.716e-5 # 0.0005081150545826

kᵥ = ν# / density

γ = 1 + R₀ / (cₚ - R₀) # Mayer's formula + definition of adiabatic constant
println("γ: $γ")
# T = (γ - 1) * (mol_mass / 1000) / R₀ * e 
#e2t = (μ/1000) / kᵦ
#e2t = (γ - 1) * mol_mass / R₀
e2t = (γ - 1) / R₀  #at 0 velocity
@show e2t
e2p = (γ - 1)
#e2p = (γ - 1)
t2e = 1/e2t

############################
# Define BCs and Operators #
############################

funcs = sym2func(sd)
∂₀ = Examples.boundary_inds(Val{0}, s)
∂₁ = Examples.boundary_inds(Val{1}, s)

∂ₗ₀ = ∂₀[findall(p-> p[1] <= 1e-4, s[∂₀, :point])]
∂ᵣ₀ = ∂₀[findall(p-> p[1] >= 1.0 - 1e-4, s[∂₀, :point])]
∂ₜ₀ = ∂₀[findall(p-> p[2] >= 1.0 - 1e-4, s[∂₀, :point])]
∂ᵦ₀ = ∂₀[findall(p-> p[2] <= 1e-4, s[∂₀, :point])]
∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)

∂ₗ₁ = Examples.bound_edges(s, ∂ₗ₀)
∂ᵣ₁ = Examples.bound_edges(s, ∂ᵣ₀)
∂ₑ₁ = Examples.bound_edges(s, ∂ₑ₀)

# Just induced edges
∂₁₊ = Examples.bound_edges(s, ∂₀)

# Full adjacent boundary
#∂₁₊ = Examples.adj_edges(s, ∂₀)
# Just edges and corners
#=
∂_corners = Examples.adj_edges(s,
                vcat(intersect(∂ₗ₀, ∂ₜ₀), intersect(∂ₗ₀, ∂ᵦ₀),
                     intersect(∂ᵣ₀, ∂ₜ₀), intersect(∂ᵣ₀, ∂ᵦ₀)))
@show ∂_corners
∂₁₊ = collect(unique(vcat(∂ₑ₁, ∂_corners)))
=#
# Full triangle on boundary

#∂₁₊ = Examples.adj_edges(s, ∂₀)
#∂_points = unique(vcat(s[∂₁₊, :∂v0], s[∂₁₊, :∂v1]))
#∂₁₊ = Examples.bound_edges(s, ∂_points)


∂ₗ₁₊ = Examples.adj_edges(s, ∂ₗ₀)
∂ᵣ₁₊ = Examples.adj_edges(s, ∂ᵣ₀)
∂ₑ₁₊ = Examples.adj_edges(s, ∂ₑ₀)


∂ₗ₀₊ = unique(vcat(s[∂ₗ₁₊, :∂v1], s[∂ₗ₁₊, :∂v0]))
∂ᵣ₀₊ = unique(vcat(s[∂ᵣ₁₊, :∂v1], s[∂ᵣ₁₊, :∂v0]))
∂ₑ₀₊ = unique(vcat(s[∂ₑ₁₊, :∂v1], s[∂ₑ₁₊, :∂v0]))

funcs[:mcopy] = Dict(:operator => I(ne(sd)), :type => MatrixFunc())
funcs[:k] = Dict(:operator => k, :type => MatrixFunc())
funcs[:e2p] = Dict(:operator => e2p * I(nv(sd)), :type => MatrixFunc())
funcs[:e2t] = Dict(:operator => e2t * I(nv(sd)), :type => MatrixFunc())
funcs[:t2e] = Dict(:operator => (1/e2t) * I(nv(sd)), :type => MatrixFunc())
funcs[:half] = Dict(:operator => 0.5 * I(nv(sd)), :type => MatrixFunc())
funcs[:R₀⁻¹] = Dict(:operator => 1/R₀ * I(nv(sd)), :type => MatrixFunc())
funcs[:third] = Dict(:operator => (1/3) * I(ne(sd)), :type => MatrixFunc())
funcs[:R₀] = Dict(:operator => R₀ * I(nv(sd)), :type => MatrixFunc())
funcs[:kᵥ] = Dict(:operator => kᵥ * I(ne(sd)), :type => MatrixFunc())
funcs[:kᵨ] = Dict(:operator => kᵨ * I(ne(sd)), :type => MatrixFunc())
funcs[:dneg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:mult₀] = Dict(:operator => (x′,x,y) -> (x′ .= x .* y), :type => InPlaceFunc())
funcs[:plus] = Dict(:operator => (x′,x,y) -> (x′ .= x .+ y), :type => InPlaceFunc())
funcs[:sum₀] = Dict(:operator => (x′,x,y) -> (x′ .= x .+ y), :type => InPlaceFunc())
funcs[:sum₁] = Dict(:operator => (x′,x,y) -> (x′ .= x .+ y), :type => InPlaceFunc())
funcs[:sum₂] = Dict(:operator => (x′,x,y) -> (x′ .= x .+ y), :type => InPlaceFunc())
funcs[:sum₀̃] = Dict(:operator => (x′,x,y) -> (x′ .= x .+ y), :type => InPlaceFunc())
funcs[:sum₁̃] = Dict(:operator => (x′,x,y) -> (x′ .= x .+ y), :type => InPlaceFunc())
funcs[:sum₂̃] = Dict(:operator => (x′,x,y) -> (x′ .= x .+ y), :type => InPlaceFunc())
funcs[:sub₀] = Dict(:operator => (x′,x,y) -> (x′ .= x .- y), :type => InPlaceFunc())
funcs[:div₁] = Dict(:operator => (x′,x,y) -> (x′ .= x ./ y), :type => InPlaceFunc())
funcs[:div₀] = Dict(:operator => (x′,x,y) -> (x′ .= x ./ y), :type => InPlaceFunc())
funcs[:avg₀₁] = Dict(:operator => avg_mat(Val{(0,1)}, s), :type => MatrixFunc())
funcs[:neg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₀] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())
funcs[:avg] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())

funcs[:∂ₗ₀₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())

const ramp_up = 0.1 # seconds
# Final temperatures should be
temperature = 0.0 #t2e * 172.89

# Pre-allocate boundary condition values
const ∂ₗ₀_z = zeros(length(∂ₗ₀))
const ∂ᵣ₀_z = zeros(length(∂ᵣ₀))
const ∂ₗ₀_t = fill(temperature / ramp_up, length(∂ₗ₀))
const ∂ᵣ₀_t = fill(-(temperature / ramp_up), length(∂ᵣ₀))
ramp_up_temp(x′, x, t) = begin
  x′ .= x
  if t < ramp_up
    setindex!(x′, ∂ₗ₀_t, ∂ₗ₀)
    setindex!(x′, ∂ᵣ₀_t, ∂ᵣ₀)
  else
    setindex!(x′, ∂ₗ₀_z, ∂ₗ₀)
    setindex!(x′, ∂ᵣ₀_z, ∂ᵣ₀)
  end
end

const ∂₁₊_z = zeros(length(∂₁₊))
const ∂ₗ₁₊_z = zeros(length(∂ₗ₁₊))
const ∂ᵣ₁₊_z = zeros(length(∂ᵣ₁₊))
const ∂ₑ₀₊_z = zeros(length(∂ₑ₀₊))
funcs[:∂₀ₛ] = Dict(:operator => ramp_up_temp, :type => TDInPlaceFunc())
funcs[:∂₁ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; setindex!(x′, ∂₁₊_z ,  ∂₁₊)), :type => InPlaceFunc())
funcs[:∂₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; setindex!(x′, ∂ₗ₁₊_z, ∂ₗ₁₊); setindex!(x′, ∂ᵣ₁₊_z, ∂ᵣ₁₊)), :type => InPlaceFunc())
funcs[:mask₀ₗ] = Dict(:operator => (x′,x) -> (x′ .= x; setindex!(x′, ∂ₑ₀₊_z, ∂ₑ₀₊)), :type => InPlaceFunc())

wedge_cache = init_wedge_ops(sd)
v2comp = comp_support(sd);
cache_mat = Dict(:t2c => tri2comp(s, v2comp), :e2c => edge2comp(s, v2comp), :cross => changes(sd, v2comp),
                 :α_cache => zeros(ntriangles(sd)*3), :β_cache => zeros(ntriangles(sd)*3))

funcs[:∧₁₀′] = Dict(:operator => (x′, α, β) -> (cp_2_1!(x′, β, α, cache_mat)), :type => InPlaceFunc())
funcs[:∧₁₁′] = Dict(:operator => (x′, α, β) -> (pd_wedge!(x′, Val{(1,1)},s, α, β; wedge_cache...)), :type => InPlaceFunc())

const hodge_LU = lu(⋆(Val{1}, sd))
funcs[:⋆₁⁻¹] = Dict(:operator => (x′,x) -> (x′ .= -1 * (hodge_LU \ x)), :type => InPlaceFunc())

t2e = 1/e2t
c_objs = fill(t2e * 288.15, nv(s))
#c_objs[∂ₗ₀] .= t2e * 461.04
#c_objs[∂ᵣ₀] .= t2e * 115.26
velocity(p) = [0.0, 0.0, 0.0]
gravity(p) = [0.0,-9.81,0.0]
v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;
g = ♭(sd, DualVectorField(gravity.(sd[triangle_center(sd),:dual_point]))).data;
ρ = [density for p in s[:point]]

  
new_dwd = deepcopy(exp_dwd)
to_graphviz(new_dwd, orientation=LeftToRight)
Examples.zip_dwd!(new_dwd)
to_graphviz(new_dwd, orientation=LeftToRight)
Examples.contract_matrices!(new_dwd, funcs)
to_graphviz(new_dwd, orientation=LeftToRight)

func, code = gen_sim(new_dwd, funcs, sd; autodiff=false, params=[:G]);
dt = 0.01
tseg = 1.0
for i in 1:10
  prev_loc = "res_$(i-1)"
  loc = "res_$i"
  mkpath(loc)
  u0 = vcat(c_objs, ρ, v)
  if i > 1
    pre_cond = open("$(prev_loc)/sim_res.json", "r") do f
      JSON.parse(f)
    end
    u0 .= pre_cond["$tseg"]
  end
  GC.gc()
  
  ##############
  # Simulation #
  ##############
  prob = ODEProblem(func, u0, (0.0, tseg))
  sol1 = solve(prob, Tsit5(), progress=true, progress_steps=1, dtmax=1.0e-5, saveat=dt, p=g)
  
  #################
  # Write Results #
  #################
  t_range = 1:nv(s)
  ρ_range = (1:nv(s)) .+ nv(s)
  v_range = (1:ne(s)) .+ (2 * nv(s))
   
  using JSON
  res = Dict{Float64, Vector}()
  for i in 0.0:dt:sol1.t[end]
      res[i] = sol1(i)
  end
  open("$(loc)/sim_res.json", "w") do f
      JSON.print(f, res)
  end
  
  using WriteVTK
  for i in 1:length(dt:dt:sol1.t[end])
    tris = triangle_vertices(s)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [tris[1][t],tris[2][t],tris[3][t] ]) for t in 1:ntriangles(s)];
    x =  [p[1] for p in s[:point]]
    y =  [p[2] for p in s[:point]]
    z =  [p[3] for p in s[:point]]
    vtkfile = vtk_grid("$(loc)/su2_coupled_v1_$(lpad(i, 4, "0"))", x, y, z, cells)
    vel_form = ♯(sd, sol1(dt*i)[v_range], CombinatorialSpaces.DiscreteExteriorCalculus.PPSharp())
  
    vtkfile["temperature", VTKPointData()] = e2t * sol1(dt*i)[t_range]
    vtkfile["density", VTKPointData()] = sol1(dt*i)[ρ_range]
    vtkfile["vel", VTKPointData()] = vel_form
    vtkfile["v_mag", VTKPointData()] = sqrt.(abs.(1e-7 .+ inv_hodge_star(Val{0}, sd)*pd_wedge(Val{(1,1)}, sd, sol1(dt*i)[v_range], ⋆(Val{1}, sd) * sol1(dt*i)[v_range])))
    vtk_save(vtkfile)
  end
  
  fig = Figure()
  ax, ob = mesh(fig[1,1], s, color = e2t * sol1[end][t_range], colormap=:seismic)#, colorrange=(115 + 100,462 - 100))
  
  xlims!(ax, (0,1))
  ylims!(ax, (0,1))
  ax.aspect = AxisAspect(1)
  Colorbar(fig[1,2], ob)
  fig
  save("$(loc)/heat_prof.png", fig)
  
  times = range(0.01, sol1.t[end]/2, length=150)
  colors = [e2t * sol1(t)[t_range] for t in times]
  fig = Figure()
  axis, ob = mesh(fig[1,1], s, color=colors[1], colormap=:seismic)
  axis.aspect = AxisAspect(1.0)
  Colorbar(fig[1,2], ob)
  framerate = 30
  fig
  
  record(fig, "$(loc)/conc_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
  
  
  times = range(0.01, sol1.t[end], length=150)
  colors = [sol1(t)[ρ_range] for t in times]
  
  fig = Figure()
  axis, ob = mesh(fig[1,1], s, color=colors[1],
                                     colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
  
  axis.aspect = AxisAspect(1.0)
  Colorbar(fig[1,2], ob)
  framerate = 30
  fig
  
  record(fig, "$(loc)/density_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
  
  times = range(0, sol1.t[end], length=150)
  colors = [sol1(t)[t_range] .* sol1(t)[ρ_range] for t in times]
  
  figure, axis, ob = mesh(s, color=colors[1],
                                     colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
  
  axis.aspect = AxisAspect(1.0)
  framerate = 30
  
  record(figure, "$(loc)/press_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
  
  magnitudes(t) = sqrt.(abs.(1e-4 .+ inv_hodge_star(Val{0}, sd)*pd_wedge(Val{(1,1)}, sd, sol1(t)[v_range], ⋆(Val{1}, sd) * sol1(t)[v_range])))
  times = range(0.01, sol1.t[end], length=150)
  colors = [magnitudes(t) for t in times]
  fig = Figure()
  axis, ob = mesh(fig[1,1], s, color=colors[3], colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
  
  axis.aspect = AxisAspect(1.0)
  Colorbar(fig[1,2], ob)
  framerate = 30
  fig
  record(fig, "$(loc)/vel_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
end
