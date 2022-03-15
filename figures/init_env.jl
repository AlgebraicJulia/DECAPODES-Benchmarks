using Pkg
Pkg.activate(".")
Pkg.add(PackageSpec(url="https://github.com/AlgebraicJulia/Decapods.jl.git"))
Pkg.instantiate()
