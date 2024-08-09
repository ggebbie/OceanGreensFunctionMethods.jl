### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 8f520c8b-19d7-48a8-be9f-3f167f07d188
import Pkg

# ╔═╡ c536e9f3-0457-499e-958c-384d6e388ef9
Pkg.activate(".")

# ╔═╡ de3c6443-5ca1-4e97-82c8-5c4c9f204480
using Revise

# ╔═╡ 69147ae0-1c89-48a6-831b-ff325a984817
using OceanGreensFunctionMethods

# ╔═╡ 97df5706-1829-419a-b96b-5dbb4d704434
using DimensionalData

# ╔═╡ 157462cf-b6f9-4de3-80f6-a3f846b5ea1a
using DimensionalData: @dim

# ╔═╡ 0a9a45e2-a561-4a21-afb9-b96ec884de4a
using Unitful

# ╔═╡ 5d30be92-5266-4700-ba7b-ac88a7f066e3
# define our own unit: the sverdrup
module UnitfulOcean; using Unitful; 
@unit sverdrup "Sv" Sverdrup (10^6)u"m^3/s" false;
end

# ╔═╡ 10b07d8a-aee4-4b64-b9eb-f22f408877ba
md"""
# Pedagogical Box Model 

Julia package to complement "A Review of Green's Function Methods in Ocean Circulation Models," by Haine et al. This package goes toward one of the stated goals of the manuscript, namely to make Green's Function methods accessible for learning purposes. Here, we also aim to make a Julia package that is useful and computationally efficient for research purposes."""


# ╔═╡ bac2fd5b-70b5-4afb-8580-a13552f9c482
html"""<style>
main {
    max-width: 96%;
    margin-left: 1%;
    margin-right: 2% !important;
}
"""

# ╔═╡ 17a2dd87-04b7-45a2-975e-f0e670a07eac


# ╔═╡ 27b7af71-e396-45b3-8723-8b2fc804a77f
md"""## Activate a reproducible project environment

Uses Julia's built-in package manager `Pkg.jl` """

# ╔═╡ 07f01269-cfd8-4d3d-8d85-0b1132ff2005
md""" ## Load some helpful packages """


# ╔═╡ 6d11d809-9902-4a6a-b85e-18aed70e352f
md""" ## Define and label the boxes """

# ╔═╡ 55a42a84-587b-41ec-8b18-96f83245ee7d
@dim Meridional "meridional location"; @dim Vertical "vertical location"

# ╔═╡ 543f8a23-2e43-426e-87f0-e7750bcadd2b
# labels for the three latitudes
meridional_locs = ["1 High latitudes", "2 Mid-latitudes", "3 Low latitudes"]

# ╔═╡ 8fb8a936-9f06-4944-8c18-02eaa32f2dd0
# labels for the three depths
vertical_locs = ["1 Thermocline", "2 Deep", "3 Abyssal"]

# ╔═╡ 931c8f1a-d97f-4543-8e7b-a24304651c0b
# 3 x 3 box model
Ny = length(meridional_locs); Nz = length(vertical_locs)

# ╔═╡ d936f769-aa79-41a6-ba56-4b99eb8738bc
# define "dimensions" to be used with `DimensionalData.jl` 
# permits numerical quantities to be bundled with their meta-data
model_dims = (Meridional(meridional_locs),Vertical(vertical_locs))

# ╔═╡ fdd4e823-bdb5-4f02-8de3-aade687a94c6
md""" ## Point of emphasis: physical units

We aren't just solving mathematical constructs, but we are interpreting physical scenarios. Use `Unitful.jl` to put units on all quantities. """

# ╔═╡ ae1f5365-78c5-4ae5-8aaf-c0818fa8c474
# define units using a convenient shorthand
const kg = u"kg" # kilograms

# ╔═╡ f1fff88c-0357-4f56-bea5-75b9a63807c0
const m = u"m" # meters

# ╔═╡ 5e51cd6d-11db-4d1d-9a22-38d5a7efeea1
const yr = u"yr" # years

# ╔═╡ 09f35864-f687-4513-8c3d-3d14961c27bc
const km = u"km" # kilometers

# ╔═╡ 8134c137-fcd4-47b9-8ffd-8d524c004ced
const Tg = u"Tg" # teragrams

# ╔═╡ 2be17af4-9dc5-4179-ac30-0f8ca6da64e9
const s = u"s" # seconds

# ╔═╡ ec6602e6-deeb-4358-98b5-6bdf69bafd35
# display units in a nicer fashion
ENV["UNITFUL_FANCY_EXPONENTS"] = true

# ╔═╡ 00f69a96-d8d7-4e7c-905c-461f8132c565
# some boilerplate to register the new unit
Unitful.register(UnitfulOcean)

# ╔═╡ 8306d2c4-8d50-4309-add1-6d1eef56cd4a
# set the units of a quantity, then use a pipe to convert units
Vol0 = 1e16m^3 |> km^3 # uniform value of volume for all boxes

# ╔═╡ 1e91d3e1-c26f-4118-9630-d654d352da76
Vol = DimArray(fill(Vol0, Ny, Nz), model_dims)


# ╔═╡ 51d0e115-5859-4eab-8a91-b8193afd52b5
Vol' # take transpose or complex conjugate transpose to view more intuitively

# ╔═╡ e17220c4-d45d-4d36-a05f-c245393b05ef


# ╔═╡ 75f344b4-e273-4369-89dc-5ebfdb675d21


# ╔═╡ 1a5b40f1-74de-4693-b838-a3e4cb0a3f77


# ╔═╡ c685a1d0-00fc-4e4b-bbe7-72174f53b919


# ╔═╡ b47e501c-86a0-41b5-99ba-06aef5522384


# ╔═╡ ba360ff5-8ad9-44ba-8af3-c1e8d25e9040


# ╔═╡ cbd95539-4ac8-49ce-b2c7-89549dfc911d


# ╔═╡ c2c7c17e-6dee-4423-87ee-4db359a5ac5a


# ╔═╡ 8472db81-951e-4e81-9949-4183fa91c338


# ╔═╡ 33bc9c53-ab2b-4406-81d0-d316ac1862c0


# ╔═╡ Cell order:
# ╟─10b07d8a-aee4-4b64-b9eb-f22f408877ba
# ╠═bac2fd5b-70b5-4afb-8580-a13552f9c482
# ╠═17a2dd87-04b7-45a2-975e-f0e670a07eac
# ╟─27b7af71-e396-45b3-8723-8b2fc804a77f
# ╠═8f520c8b-19d7-48a8-be9f-3f167f07d188
# ╠═c536e9f3-0457-499e-958c-384d6e388ef9
# ╠═07f01269-cfd8-4d3d-8d85-0b1132ff2005
# ╠═de3c6443-5ca1-4e97-82c8-5c4c9f204480
# ╠═69147ae0-1c89-48a6-831b-ff325a984817
# ╠═97df5706-1829-419a-b96b-5dbb4d704434
# ╠═157462cf-b6f9-4de3-80f6-a3f846b5ea1a
# ╠═0a9a45e2-a561-4a21-afb9-b96ec884de4a
# ╟─6d11d809-9902-4a6a-b85e-18aed70e352f
# ╠═55a42a84-587b-41ec-8b18-96f83245ee7d
# ╠═543f8a23-2e43-426e-87f0-e7750bcadd2b
# ╠═8fb8a936-9f06-4944-8c18-02eaa32f2dd0
# ╠═931c8f1a-d97f-4543-8e7b-a24304651c0b
# ╠═d936f769-aa79-41a6-ba56-4b99eb8738bc
# ╟─fdd4e823-bdb5-4f02-8de3-aade687a94c6
# ╠═ae1f5365-78c5-4ae5-8aaf-c0818fa8c474
# ╠═f1fff88c-0357-4f56-bea5-75b9a63807c0
# ╠═5e51cd6d-11db-4d1d-9a22-38d5a7efeea1
# ╠═09f35864-f687-4513-8c3d-3d14961c27bc
# ╠═8134c137-fcd4-47b9-8ffd-8d524c004ced
# ╠═2be17af4-9dc5-4179-ac30-0f8ca6da64e9
# ╠═ec6602e6-deeb-4358-98b5-6bdf69bafd35
# ╠═5d30be92-5266-4700-ba7b-ac88a7f066e3
# ╠═00f69a96-d8d7-4e7c-905c-461f8132c565
# ╠═8306d2c4-8d50-4309-add1-6d1eef56cd4a
# ╠═1e91d3e1-c26f-4118-9630-d654d352da76
# ╠═51d0e115-5859-4eab-8a91-b8193afd52b5
# ╠═e17220c4-d45d-4d36-a05f-c245393b05ef
# ╠═75f344b4-e273-4369-89dc-5ebfdb675d21
# ╠═1a5b40f1-74de-4693-b838-a3e4cb0a3f77
# ╠═c685a1d0-00fc-4e4b-bbe7-72174f53b919
# ╠═b47e501c-86a0-41b5-99ba-06aef5522384
# ╠═ba360ff5-8ad9-44ba-8af3-c1e8d25e9040
# ╠═cbd95539-4ac8-49ce-b2c7-89549dfc911d
# ╠═c2c7c17e-6dee-4423-87ee-4db359a5ac5a
# ╠═8472db81-951e-4e81-9949-4183fa91c338
# ╠═33bc9c53-ab2b-4406-81d0-d316ac1862c0
