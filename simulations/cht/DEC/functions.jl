k_col = fill(k₁, ne(s))
k_col[cyl] .= k₂
k = diagm(k_col)

funcs[:mcopy] = Dict(:operator => I(ne(sd)), :type => MatrixFunc())
funcs[:k] = Dict(:operator => k, :type => MatrixFunc())
funcs[:e2p] = Dict(:operator => e2p * I(nv(sd)), :type => MatrixFunc())
funcs[:e2t] = Dict(:operator => e2t * I(nv(sd)), :type => MatrixFunc())
funcs[:t2e] = Dict(:operator => (1/e2t) * I(nv(sd)), :type => MatrixFunc())
funcs[:half] = Dict(:operator => 0.5 * I(nv(sd)), :type => MatrixFunc())
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
funcs[:div₀] = Dict(:operator => (x′,x,y) -> (x′ .= x ./ y), :type => InPlaceFunc())
funcs[:div₁] = Dict(:operator => (x′,x,y) -> (x′ .= x ./ y), :type => InPlaceFunc())
funcs[:avg₀₁] = Dict(:operator => avg_mat(Val{(0,1)}, s), :type => MatrixFunc())
funcs[:neg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₀] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())
funcs[:avg] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())

funcs[:∂ₗ₀₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())
funcs[:∂ₜₐ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[cyl_inner] .= 0), :type => InPlaceFunc())
funcs[:∂ₜ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀₊] .= 0; x′[∂ₒ₀] .= 0), :type => InPlaceFunc())
funcs[:∂ₚ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀₊] .= 0; x′[cyl_inner] .= 0), :type => InPlaceFunc())
funcs[:∂ᵥ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[cyl_edge] .= 0; x′[∂ₑ₁₊] .= 0), :type => InPlaceFunc())

const ramp_up = 0.1 # seconds
# Final temperatures should be
temperature = t2e * 172.89

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


