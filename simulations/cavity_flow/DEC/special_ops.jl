function edge_to_support(s)
  vals = Dict{Tuple{Int64, Int64}, Float64}()
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(s)
      de = elementary_duals(Val{1},s, e)
      ev = volume(Val{1}, s, e)
      dv = sum([dual_volume(Val{1}, s, d) for d in de])
      for d in de
          dt = incident(s, d, :D_∂e0)
          append!(I, dt)
          append!(J, fill(e, length(dt)))
          append!(V, fill(1/(dv*ev), length(dt)))
      end
  end
  sparse(I,J,V)
end

function tri_to_support(s)
  vals = Dict{Tuple{Int64, Int64}, Float64}()
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for t in 1:ntriangles(s)
      dt = vcat(incident(s, incident(s, triangle_center(s, t), :D_∂v0), :D_∂e1)...)
      tv = volume(Val{2}, s, t)
      append!(I, dt)
      append!(J, fill(t, length(dt)))
      append!(V, fill(1/tv#= * sign(Val{2}, s, t)=#, length(dt)))
  end
  sparse(I,J,V)
end

function support_to_tri(s)
  vals = Dict{Tuple{Int64, Int64}, Float64}()
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for t in 1:nv(s)
      dt = elementary_duals(Val{0},s, t)
      for d in dt
          push!(I, t)
          push!(J, d)
          push!(V, 1)
      end
  end
  sparse(I,J,V)
end

diag_vols(s) = spdiagm([dual_volume(Val{2}, s, dt) for dt in 1:nparts(s, :DualTri)])

wedge_mat(::Type{Val{(1,1)}}, s) = 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))
wedge_mat(::Type{Val{(2,0)}}, s) = support_to_tri(s)*diag_vols(s)*tri_to_support(s)

wedge_edge(::Type{Val{(1,1)}}, s) = dual_boundary(Val{2}, s) * 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))
wedge_edge(::Type{Val{(2,0)}}, s) = dual_boundary(Val{2}, s) * support_to_tri(s)*diag_vols(s)*tri_to_support(s)

function pd_wedge!(x, ::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat(Val{(1,1)}, s)), caches=[zeros(nv(s)), zeros(ne(s))], kw...)
  broadcast!(*, caches[2], α, β)
  mul!(x, wedge_t[(1,1)], caches[2])
end

function pd_wedge(::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat(Val{(1,1)}, s)), caches=[zeros(nv(s)), zeros(ne(s))], kw...)
  broadcast!(*, caches[2], α, β)
  wedge_t[(1,1)] * caches[2]
end

function pd_wedge(::Type{Val{(2,0)}}, s, α, β; wedge_t = Dict((2,0)=>wedge_mat(Val{(2,0)},s)), kw...)
  wedge_t[(2,0)] * (α .* β)
end

function pd_wedge(::Type{Val{(0,2)}}, s, α, β; wedge_t = nothing, kw...)
  α .* β
end

function init_wedge_ops(s)
  (d_mat=Dict(:d₁=>d(Val{1}, s), :dual_d₀=>dual_derivative(Val{0}, s)),
   wedge_e=Dict((1,1)=>wedge_edge(Val{(1,1)},s), (2,0)=>wedge_edge(Val{(2,0)},s)),
   wedge_t=Dict((1,1)=>wedge_mat(Val{(1,1)}, s), (2,0)=>wedge_mat(Val{(2,0)},s)),
   caches=[zeros(nv(s)), zeros(ne(s)), zeros(ntriangles(s))])
end

vect(s, e) = (s[s[e,:∂v1], :point] - s[s[e,:∂v0], :point]) * sign(1, s, e)
vect(s, e::AbstractVector) = [vect(s, el) for el in e]
t_vects(s,t) = vect(s, triangle_edges(s,t)) .* ([1,-1,1] * sign(2,s,t))
function comp_support(sd)
  vects = []
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(vects, (t, inds[i]))
      end
  end
  v2comp = Dict{Tuple{Int64, Int64}, Int64}()
  for (i, v) in enumerate(vects)
      v2comp[v] = i
  end
  v2comp
end
function changes(s, v2comp)
  orient_vals = [1,-1,1]
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for t in 1:ntriangles(s)
      inds = triangle_edges(s, t)
      e_vects = t_vects(s,t)
      vals = zeros(1:3)
      for i in 1:3
          ns = [(i+1)%3 + 1, i%3+1]
          ort = e_vects[i] × (e_vects[i] × e_vects[ns[1]])
          n_ort = ort / norm(ort)
          append!(J, v2comp[(t,inds[i])])
          append!(I, inds[ns[1]])
          append!(V, dot(n_ort, e_vects[ns[1]]) * orient_vals[ns[1]] * sign(1, s, ns[1])* orient_vals[i]* sign(2,s,t) / 3.0)
          append!(J, v2comp[(t,inds[i])])
          append!(I, inds[ns[2]])
          append!(V, dot(n_ort, e_vects[ns[2]]) * orient_vals[ns[2]] * sign(1, s, ns[2])* orient_vals[i]* sign(2,s,t) / 3.0)
      end
  end
  sparse(I,J,V, ne(s), ntriangles(s)*3)
end
function edge2comp(s, v2comp)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(I, v2comp[(t,inds[i])])
          push!(J, inds[i])
          push!(V, 1 / volume(Val{1}, s, inds[i]))
      end
  end
  sparse(I,J,V)
end
function tri2comp(s, v2comp)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(I, v2comp[(t,inds[i])])
          push!(J, t)
          push!(V, 1)
      end
  end
  sparse(I,J,V)
end

function cp_2_1!(x, α, β, matrices)
  mul!(matrices[:α_cache], matrices[:t2c], α)
  mul!(matrices[:β_cache], matrices[:e2c], β)
  broadcast!(*, matrices[:β_cache], matrices[:α_cache], matrices[:β_cache])
  mul!(x, matrices[:cross], matrices[:β_cache])
  #x .= matrices[:cross] * ((matrices[:t2c]*α).*(matrices[:e2c]*β))
end


function avg_mat(::Type{Val{(0,1)}},s)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(s)
      append!(J, [s[e,:∂v0],s[e,:∂v1]])
      append!(I, [e,e])
      append!(V, [0.5, 0.5])
  end
  sparse(I,J,V)
end
