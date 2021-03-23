### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f6ad04d4-8b2d-11eb-377a-8d687e3d41a0
begin
	# instantiate project
	using Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	
	# load packages
	using GeoStats
	using PlutoUI
	using Plots
	gr(format=:png)
end;

# ╔═╡ 5da15a84-8b31-11eb-0fbf-f104e95a7818
md"""
# Automatic variography
"""

# ╔═╡ 25a7181e-8b6b-11eb-1f42-11cff1be2b52
md"""
Author: [Júlio Hoffimann](https://juliohm.github.io) @ julio.hoffimann@impa.br

[VI Workshop of Mathematical Solutions](http://www.cemeai.icmc.usp.br/6WSMPI)

This notebook provides a baseline solution for the Deere challenge.
"""

# ╔═╡ cc3fce5a-8b6b-11eb-3ef7-c1d30acb1e5a
md"""
## Project setup
"""

# ╔═╡ 372046e2-8b31-11eb-20b3-c10d04e905be
md"""
## Load and plot data
"""

# ╔═╡ 5083df82-8b6a-11eb-1bce-fbcf900ddd9d
md"""
We load the data from a CSV file and specify which columns represent geospatial coordinates. We also specify the types of some of the columns so that they are interpreted correctly as continuous variables:
"""

# ╔═╡ 1e0bc97a-8b2e-11eb-1e7d-99639a945454
geotable = readgeotable("data/data.csv", coordnames = (:POINT_X, :POINT_Y),
	                    types = Dict(:Sand => Float64, :Silt => Float64,
		                             :Clay => Float64, :Ca => Float64));

# ╔═╡ 011d10ee-8b6a-11eb-17a3-dbaad6028cb7
md"""
After loading the data, we make sure that we don't have repeated coordinates:
"""

# ╔═╡ 9ef6a654-8b65-11eb-160d-f5bc9cc14bf2
𝒮 = uniquecoords(geotable)

# ╔═╡ 393dc7a2-8b6a-11eb-2f8e-8d082ccc722c
md"""
We can plot the data directly given that it is georeferenced already:
"""

# ╔═╡ 49c1439c-8b2e-11eb-0bbf-2fd6ad8b6fbf
plot(𝒮, (:Sand,:Silt,:Clay,:Ca), ms = 4, size = (700,500))

# ╔═╡ 7c469242-8b31-11eb-28ee-0b722bd46c25
md"""
## Fit variogram model
"""

# ╔═╡ 49354c34-8b6a-11eb-188f-eb898e2eaa31
md"""
In order to estimate an empirical variogram, we need to specify the maximum lag. We will use half of the minimum side of the bounding box as a rule of thumb:
"""

# ╔═╡ c31433ec-8b62-11eb-10c9-c740d7bfcd8e
ℬ = boundingbox(𝒮)

# ╔═╡ 9eecf284-8b2f-11eb-16cd-1109694718e9
maxlag = minimum(sides(ℬ)) / 2

# ╔═╡ 7c138c10-8b6a-11eb-20c4-7f704301e18a
md"""
Here is how we can fit variogram models automatically:
"""

# ╔═╡ 3735889a-8b2f-11eb-367f-cfd40581d114
begin
	Zs = [:Sand,:Silt,:Clay,:Ca]
	
	γs = [EmpiricalVariogram(𝒮, Z, maxlag = maxlag) for Z in Zs]
	
	γt = [fit(Variogram, γ) for γ in γs]
end

# ╔═╡ 73c32e60-8b30-11eb-016a-4d18dd46a7b5
begin
	plots = []
	for (i, Z) in enumerate(Zs)
		p = plot(γs[i], legend=:bottomright, label="empirical", size=(700,500))
		plot!(γt[i], 0, maxlag, label="$Z model")
		push!(plots, p)
	end
	plot(plots...)
end

# ╔═╡ 8b014104-8b6a-11eb-2d35-1fe531c5d208
md"""
Let's save these models in a dictionary mapping the variable name to the model:
"""

# ╔═╡ 035274aa-8b63-11eb-157e-e5923d85dc3c
model = Dict(Zs .=> γt)

# ╔═╡ e8741980-8b68-11eb-0246-118378c3ed1e
md"""
## Interpolation problem

Select a variable of interest in the drop down tp define an interpolation problem:
"""

# ╔═╡ 0bd3799e-8b68-11eb-02a7-03869dfeac09
md"""
Z: $(@bind Z Select(string.(Zs)))
"""

# ╔═╡ cbae4984-8b62-11eb-070a-2fb52c773230
𝒟 = CartesianGrid(minimum(ℬ), maximum(ℬ), dims=(100,100))

# ╔═╡ 847dbfe4-8b63-11eb-3fc2-59921190f2cd
problem = EstimationProblem(𝒮, 𝒟, Symbol(Z))

# ╔═╡ 2a60a3ee-8b65-11eb-2214-039e7762db1d
plot(problem)

# ╔═╡ a3c095e6-8b6a-11eb-00e5-5d7ac2f4c13f
md"""
Now we can solve the problem with different solvers as shown below.
"""

# ╔═╡ 8c1f2554-8b62-11eb-373f-194e459fcff3
md"""
### Kriging
"""

# ╔═╡ d971483a-8b67-11eb-2d28-2da42fa4818a
begin
	krig = Kriging(Symbol(Z) => (variogram=model[Symbol(Z)],))
	
	krigsol = solve(problem, krig)
	
	plot(krigsol, size = (900,300))
end

# ╔═╡ 7fa83a52-8b69-11eb-067b-d768376028fb
md"""
### IDW
"""

# ╔═╡ 5aef4462-8b69-11eb-2505-7752c5c4e19c
begin
	idw = IDW()
	
	idwsol = solve(problem, idw)
	
	plot(idwsol, size = (900,300))
end

# ╔═╡ Cell order:
# ╟─5da15a84-8b31-11eb-0fbf-f104e95a7818
# ╟─25a7181e-8b6b-11eb-1f42-11cff1be2b52
# ╟─cc3fce5a-8b6b-11eb-3ef7-c1d30acb1e5a
# ╠═f6ad04d4-8b2d-11eb-377a-8d687e3d41a0
# ╟─372046e2-8b31-11eb-20b3-c10d04e905be
# ╟─5083df82-8b6a-11eb-1bce-fbcf900ddd9d
# ╠═1e0bc97a-8b2e-11eb-1e7d-99639a945454
# ╟─011d10ee-8b6a-11eb-17a3-dbaad6028cb7
# ╠═9ef6a654-8b65-11eb-160d-f5bc9cc14bf2
# ╟─393dc7a2-8b6a-11eb-2f8e-8d082ccc722c
# ╠═49c1439c-8b2e-11eb-0bbf-2fd6ad8b6fbf
# ╟─7c469242-8b31-11eb-28ee-0b722bd46c25
# ╟─49354c34-8b6a-11eb-188f-eb898e2eaa31
# ╠═c31433ec-8b62-11eb-10c9-c740d7bfcd8e
# ╠═9eecf284-8b2f-11eb-16cd-1109694718e9
# ╟─7c138c10-8b6a-11eb-20c4-7f704301e18a
# ╠═3735889a-8b2f-11eb-367f-cfd40581d114
# ╟─73c32e60-8b30-11eb-016a-4d18dd46a7b5
# ╟─8b014104-8b6a-11eb-2d35-1fe531c5d208
# ╠═035274aa-8b63-11eb-157e-e5923d85dc3c
# ╟─e8741980-8b68-11eb-0246-118378c3ed1e
# ╟─0bd3799e-8b68-11eb-02a7-03869dfeac09
# ╠═cbae4984-8b62-11eb-070a-2fb52c773230
# ╠═847dbfe4-8b63-11eb-3fc2-59921190f2cd
# ╠═2a60a3ee-8b65-11eb-2214-039e7762db1d
# ╟─a3c095e6-8b6a-11eb-00e5-5d7ac2f4c13f
# ╟─8c1f2554-8b62-11eb-373f-194e459fcff3
# ╠═d971483a-8b67-11eb-2d28-2da42fa4818a
# ╟─7fa83a52-8b69-11eb-067b-d768376028fb
# ╠═5aef4462-8b69-11eb-2505-7752c5c4e19c
