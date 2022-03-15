############################
# Define BCs and Operators #
############################

# Get boundaries of cylinders

locs = [(10.5, 0.0), (5.5, 0.0), (0.5, 0.0)]
cyl = vcat(map(locs) do loc
    findall(p -> ( ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
end...)
fuzzy_bound = unique(vcat(incident(s, cyl, :∂v1)..., incident(s, cyl, :∂v0)...))
cyl_edge = filter(e -> (s[e, :∂v1] ∈ cyl)&&(s[e, :∂v0] ∈ cyl), fuzzy_bound)

cyl_bound = vcat(map(locs) do loc
    findall(p -> (0.5^2 - 1e-4 <= ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
end...)
fuzzy_boundary = unique(vcat(incident(s, cyl_bound, :∂v1)..., incident(s, cyl_bound, :∂v0)...))
cyl_bound_edge = filter(e -> (s[e, :∂v1] ∈ cyl_bound)&&(s[e, :∂v0] ∈ cyl_bound), fuzzy_boundary)
cyl_inner = filter(p -> !(p ∈ cyl_bound), cyl)
slip_edge = filter!(p -> !(p ∈ cyl_bound_edge), cyl_edge)


k_col = fill(k₁, ne(s))
k_col[cyl_edge] .= k₂
k = diagm(k_col)

# Get other boundaries
funcs = sym2func(sd)
∂₀ = Examples.boundary_inds(Val{0}, s)
∂₁ = Examples.boundary_inds(Val{1}, s)

∂ₒ₀ = ∂₀[findall(p-> -9 <= p[1] <= 20 && -9 <= p[2] <= 9, s[∂₀, :point])]

lx = -10.0
rx = 21.0
ty = 15.5
by = -15.5

∂ₗ₀ = ∂₀[findall(p-> p[1] <= lx + 1e-4, s[∂₀, :point])]
∂ᵣ₀ = ∂₀[findall(p-> p[1] >= rx - 1e-4, s[∂₀, :point])]
∂ₜ₀ = ∂₀[findall(p-> p[2] >= ty - 1e-4, s[∂₀, :point])]
∂ᵦ₀ = ∂₀[findall(p-> p[2] <= by + 1e-4, s[∂₀, :point])]
∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)

∂ₗ₁ = Examples.bound_edges(s, ∂ₗ₀)
∂ᵣ₁ = Examples.bound_edges(s, ∂ᵣ₀)
∂ₑ₁ = Examples.bound_edges(s, ∂ₑ₀)

# Just induced edges
#∂₁₊ = Examples.bound_edges(s, ∂₀)

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

∂₁₊ = Examples.adj_edges(s, ∂₀)
#∂_points = unique(vcat(s[∂₁₊, :∂v0], s[∂₁₊, :∂v1]))
#∂₁₊ = Examples.bound_edges(s, ∂_points)


∂ₗ₁₊ = Examples.adj_edges(s, ∂ₗ₀)
∂ᵣ₁₊ = Examples.adj_edges(s, ∂ᵣ₀)
∂ₑ₁₊ = Examples.adj_edges(s, ∂ₑ₀)
∂_points = unique(vcat(s[∂ₑ₁₊, :∂v0], s[∂ₑ₁₊, :∂v1]))
∂ₑ₁₊ = Examples.bound_edges(s, ∂_points)


∂ₗ₀₊ = unique(vcat(s[∂ₗ₁₊, :∂v1], s[∂ₗ₁₊, :∂v0]))
∂ᵣ₀₊ = unique(vcat(s[∂ᵣ₁₊, :∂v1], s[∂ᵣ₁₊, :∂v0]))
∂ₑ₀₊ = unique(vcat(s[∂ₑ₁₊, :∂v1], s[∂ₑ₁₊, :∂v0]))


