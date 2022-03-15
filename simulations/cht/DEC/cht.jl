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
global_logger(TerminalLogger())

using Revise
using Decapodes.Simulations
using Decapodes.Examples
using Decapodes.Diagrams
using Decapodes.Schedules
using Decapodes.Debug
using Decapodes.OpenDiagrams

using CombinatorialSpaces: volume
using SparseArrays


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


######################
# Functions for Curl #
######################
include("special_ops.jl")
include("deca_schema.jl")
include("physics.jl")
include("complex_operators.jl")
include("constants.jl")

###############
# Import Mesh #
###############

s1 = EmbeddedDeltaSet2D("../su2_mesh_square_small_51.stl")
s = EmbeddedDeltaSet2D{Bool, Point{3, Float64}}()
copy_parts!(s, s1)
sd = dual(s);
if ⋆(Val{1}, sd)[1,1] < 0.0
  orient_component!(s, 1, false)
  sd = dual(s); 
end

include("boundaries.jl")
include("functions.jl")

c_objs = fill(288.15, nv(s))
c_objs[∂ₒ₀] .= 350.0
#c_objs[∂ₗ₀] .= t2e * 461.04
#c_objs[∂ᵣ₀] .= t2e * 115.26
velocity(p) = [3.402, 0.0, 0.0]
gravity(p) = [0.0,0.0,0.0]
v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;
g = ♭(sd, DualVectorField(gravity.(sd[triangle_center(sd),:dual_point]))).data;
p = [density for p in s[:point]] * (288.15 * R₀)



@show ⋆(Val{1},sd)[1,1]

exp_dwd = Examples.expand_dwd(new_dwd, rules)
new_dwd = deepcopy(exp_dwd)
Examples.zip_dwd!(new_dwd)
Examples.contract_matrices!(new_dwd, funcs)

func, _ = gen_sim(new_dwd, funcs, sd; autodiff=false, params=[:G]);

tstep = 2.0

using JSON
for i in 1:10
  GC.gc()
  prev_loc = "res_$(i-1)"
  loc = "res_$i"
  mkpath(loc)
  u0 = vcat(c_objs, p, v)
  if i == 1
    pre_cond = open("../cht_sims/final_cond_v2.json", "r") do f
        JSON.parse(f)
    end
    u0[(1:ne(s)) .+ 2*nv(s)] .= pre_cond[(nv(s)+1):(nv(s)+ne(s))]
  else
    pre_cond = open("$prev_loc/sim_res.json", "r") do f
        JSON.parse(f)
    end
    u0 .= pre_cond["$tstep"]
    u0[∂ₑ₀₊] .= c_objs[∂ₑ₀₊]
    u0[∂ₑ₀₊ .+ nv(s)] .= p[∂ₑ₀₊]
    u0[∂ₑ₁₊ .+ (2*nv(s))] .= v[∂ₑ₁₊]
  end
  ##############
  # Simulation #
  ##############
  dt = 0.01
  prob = ODEProblem(func, u0, (0.0, tstep))
  sol1 = solve(prob, Tsit5(), progress=true, progress_steps=1, dtmax=#=4e-5=#8e-4, saveat=dt, p=g)
  
  #=@show extrema(sol1.t[2:(end-1)] .- sol1.t[1:(end-2)])
  fig, ax, ob = plot(sol1.t[2:(end-1)] .- sol1.t[1:(end-2)])
  save("$(loc)/timesteps.png", fig)=#

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
  open("$loc/sim_res.json", "w") do f
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
    vtkfile["pressure", VTKPointData()] = sol1(dt*i)[ρ_range]
    vtkfile["vel", VTKPointData()] = vel_form
    vtkfile["v_mag", VTKPointData()] = sqrt.(abs.(1e-7 .+ inv_hodge_star(Val{0}, sd)*pd_wedge(Val{(1,1)}, sd, sol1(dt*i)[v_range], ⋆(Val{1}, sd) * sol1(dt*i)[v_range])))
    vtk_save(vtkfile)
  end
  
  fig = Figure()
  ax, ob = mesh(fig[1,1], s, color = e2t * sol1[end][t_range], colormap=:seismic)#, colorrange=(115 + 100,462 - 100))
  xl, xr = (0.0, 11.0)
  yl, yr = (-1, 1)
  asp = 11 / 2.0
  xlims!(ax, (xl,xr))
  ylims!(ax, (yl,yr))
  ax.aspect = AxisAspect(asp)
  Colorbar(fig[1,2], ob)
  fig
  save("$(loc)/heat_prof.png", fig)
  
  times = range(0.01, sol1.t[end]/2, length=150)
  colors = [sol1(t)[t_range] for t in times]
  fig = Figure()
  axis, ob = mesh(fig[1,1], s, color=colors[1], colormap=:seismic)
  axis.aspect = AxisAspect(asp)
  Colorbar(fig[1,2], ob)
  framerate = 30
  fig
  
  xlims!(axis, (xl, xr))
  ylims!(axis, (yl, yr))
  record(fig, "$(loc)/conc_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
  
  
  times = range(0.01, sol1.t[end], length=150)
  colors = [sol1(t)[ρ_range] for t in times]
  
  fig = Figure()
  axis, ob = mesh(fig[1,1], s, color=colors[1],
                                     colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
  
  xlims!(axis, (xl, xr))
  ylims!(axis, (yl, yr))
  axis.aspect = AxisAspect(asp)
  Colorbar(fig[1,2], ob)
  framerate = 30
  fig
  
  record(fig, "$(loc)/press_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
  
  times = range(0, sol1.t[end], length=150)
  colors = [ sol1(t)[ρ_range] ./ (R₀ * sol1(t)[t_range]) for t in times]
  
  figure, axis, ob = mesh(s, color=colors[1],
                                     colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
  
  xlims!(axis, (xl, xr))
  ylims!(axis, (yl, yr))
  axis.aspect = AxisAspect(asp)
  framerate = 30
  
  record(figure, "$(loc)/density_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
  
  magnitudes(t) = sqrt.(abs.(1e-4 .+ inv_hodge_star(Val{0}, sd)*pd_wedge(Val{(1,1)}, sd, sol1(t)[v_range], ⋆(Val{1}, sd) * sol1(t)[v_range])))
  times = range(0.01, sol1.t[end], length=150)
  colors = [magnitudes(t) for t in times]
  fig = Figure()
  axis, ob = mesh(fig[1,1], s, color=colors[3], colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
  
  xlims!(axis, (xl, xr))
  ylims!(axis, (yl, yr))
  axis.aspect = AxisAspect(asp)
  Colorbar(fig[1,2], ob)
  framerate = 30
  fig
  record(fig, "$(loc)/vel_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
  end
end
