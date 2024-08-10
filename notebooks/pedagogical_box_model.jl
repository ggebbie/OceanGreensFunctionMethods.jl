### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8f520c8b-19d7-48a8-be9f-3f167f07d188
import Pkg

# ╔═╡ c536e9f3-0457-499e-958c-384d6e388ef9
Pkg.activate(".")

# ╔═╡ de3c6443-5ca1-4e97-82c8-5c4c9f204480
using Revise

# ╔═╡ 69147ae0-1c89-48a6-831b-ff325a984817
using OceanGreensFunctionMethods

# ╔═╡ b85c6513-5a1a-4fdf-b4be-2efc3c1db830
using PlutoUI

# ╔═╡ bbc2b198-ca1f-461d-a72c-e37695de357c
using LinearAlgebra

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

# ╔═╡ 852f36b7-170b-4d20-bd31-af6ae5c716a5
Sv = u"sverdrup" # a convenient shortcut

# ╔═╡ 8306d2c4-8d50-4309-add1-6d1eef56cd4a
# set the units of a quantity, then use a pipe to convert units
Vol0 = 1e16m^3 |> km^3 # uniform value of volume for all boxes

# ╔═╡ 1e91d3e1-c26f-4118-9630-d654d352da76
# If your screen is big enough, you should see a labeled, 3 x 3 table of volume values
Vol = DimArray(fill(Vol0, Ny, Nz), model_dims)


# ╔═╡ 51d0e115-5859-4eab-8a91-b8193afd52b5
Vol' # take transpose or complex conjugate transpose to view more intuitively

# ╔═╡ b9f2165e-2d18-4179-a69f-ab0fc6ceb8b6
md""" abyssal overturning rate $(@bind Ψ_abyssal Slider((2:40)Sv,show_value = true, default = 20Sv)) """

# ╔═╡ 2d21fdae-8d7d-4ef5-a447-9a2f37e695a4
md""" intermediate overturning rate $(@bind Ψ_intermediate Slider((2:40)Sv,show_value = true, default = 10Sv)) """

# ╔═╡ 246b677a-24bc-41da-96d5-1bff248657b8
md""" vertical diffusion (exchange flux) $(@bind Fv_exchange Slider((1:30)Sv,show_value = true, default = 5Sv)) """

# ╔═╡ cd6dc878-4442-430b-a263-3651719f2f11
# abyssal volume flux
Fv_abyssal = abyssal_overturning(Ψ_abyssal, model_dims) # volume fluxes

# ╔═╡ ff24a30f-56ee-4095-8bb3-0c7e4a72fe87
# volume flux in intermediate overturning
Fv_intermediate = intermediate_overturning(Ψ_intermediate, model_dims) # volume fluxes

# ╔═╡ 4bc2c6b0-ad10-4035-be82-d02060d1b3d7
# vertical diffusive-like flux
Fv_diffusion = vertical_diffusion(Fv_exchange, model_dims) # volume fluxes

# ╔═╡ e17220c4-d45d-4d36-a05f-c245393b05ef
# combine all of the volume fluxes
Fv = Fv_abyssal + Fv_intermediate + Fv_diffusion

# ╔═╡ 75f344b4-e273-4369-89dc-5ebfdb675d21
# do the volume fluxes conserve mass?
deldotFm = mass_convergence(Fv)

# ╔═╡ 095bd0d6-0249-4e26-b91a-ff488980e119
# interrogate the mass convergence by box number
deldotFm[2,2]

# ╔═╡ fa9f454b-b7ca-4e37-a67c-a28ff91a5e11
# or use the box names using the `At` notation
deldotFm[At("2 Mid-latitudes"),At("2 Deep")]

# ╔═╡ 17bf7f50-78a7-4c7c-bc2f-4d9086dd2181
# it's ok if you don't remember the order of the dimensions
# here, get the three boxes in the Mid-latitudes
deldotFm[Meridional=At("2 Mid-latitudes")]

# ╔═╡ 08464a04-7652-4f24-978d-cd329e7fe0a7
# Given a tracer distribution C and volume fluxes Fv, find the tracer fluxes
C = rand(model_dims) # first generate a uniform U(0,1) tracer distribution

# ╔═╡ 5dfddc9c-6313-4679-a994-15a771ee4a90
# ╠═╡ disabled = true
#=╠═╡
# then solve for tracer fluxes
J = advective_diffusive_flux(C, Fv)
  ╠═╡ =#

# ╔═╡ 2da2069f-4035-448f-9833-6fcd1156b2e0
J = advective_diffusive_flux(C, Fv)

# ╔═╡ e325d781-ae5c-4f64-a608-170b4df77882
# fluxes are stored in a structure that is organized by directionality
# again, use transpose to see the screen output in a reasonable order
J.poleward' 

# ╔═╡ 1e92642c-396f-4353-aa5c-8849cf26af1d
# greatest poleward fluxes in Thermocline
J.poleward[Meridional=At("2 Mid-latitudes")] # hit rightward arrow above to see flux values

# ╔═╡ 86076566-a96b-4faf-bdef-93b95733dcff
deldotJ = convergence(J) # tracer flux convergence

# ╔═╡ c9abc24c-d3f2-4d64-8dfc-b0fdf42d1502
md"""## Boundary conditions """

# ╔═╡ 6ed0c4d7-8fbb-439e-85fd-382e6b6e030f
# boundary exchange: define the locations affected by boundary fluxes
meridional_boundary = ["1 High latitudes", "2 Mid-latitudes"]; vertical_boundary = ["1 Thermocline"]

# ╔═╡ c5381568-9e1a-4763-8da3-86638b468ec4
# define dimensions that label the numerical arrays
boundary_dims = (Meridional(meridional_boundary), Vertical(vertical_boundary))

# ╔═╡ 378e4e6c-d399-458d-85a9-23c8ceda2b43
# prescribe boundary volume fluxes
Fb = DimArray(hcat([10Sv, 5Sv]), boundary_dims) # boundary flux

# ╔═╡ d344750e-e335-4e3c-baaa-a2937c2497df
# example: B_D (Dirichlet boundary conditions) to 1
f = ones(boundary_dims) # boundary tracer values

# ╔═╡ ba4789f3-7576-423b-9e94-abf4c3259eb4
C0 = zeros(model_dims) # zero interior tracer, will help identify boundary source in example

# ╔═╡ f51bfbed-0c3b-415c-9aef-80574b905b17
# boundary flux (already posed as a convergence or net effect in each box) 
Jb = boundary_flux( f, C0, Fb)

# ╔═╡ f44447e3-5e8e-4fbc-b2cf-83176fb93c9f
md""" ## Construct transport matrix """

# ╔═╡ 6897f4af-ca8e-43a7-b741-5f2dd48c97cb
# example: find the tracer tendency for a given box-model state
dCdt = tracer_tendency(C, f, Fv, Fb, Vol)

# ╔═╡ 9656a0f3-59ff-4bf2-85aa-33a5b29fd7d9
# find A matrix.
# If f = 0, q = 0, then dC/dt  = Ac
A =  linear_probe(tracer_tendency, C, f, Fv, Fb, Vol)

# ╔═╡ 15e6cead-7de1-4cdd-ae84-f7537e789900
# A is stored with box labels for both the rows and columns
# to view matrix in usual mathematical form
Matrix(A)

# ╔═╡ 6a25a144-0ccc-4604-85a7-b724eaa4cfed
# select of column of A corresponding to a tracer location
A[Vertical=At("3 Abyssal"),Meridional=At("2 Mid-latitudes")] # still displayed with info about spatial-locations

# ╔═╡ 382db56a-d39b-4835-bf13-6dd0088b0b39
# select an entry of A, caution: first index=column, second index=row
A[5][5]

# ╔═╡ 1312b135-a267-4736-8e56-ff44bc7be59b
# or get the same information using labels, but it gets long
A[Vertical=At("2 Deep"),Meridional=At("2 Mid-latitudes")][Vertical=At("2 Deep"),Meridional=At("2 Mid-latitudes")] # still displayed with info about spatial-locations

# ╔═╡ 59b47e9c-784a-4ed5-aeb6-79b5c756fff6
md""" ## Construct boundary matrix """

# ╔═╡ e74b7cee-01d6-4984-be3d-16a71b350c99
# probe for B (boundary matrix)
B =  linear_probe(tracer_tendency, f, C, Fb, Vol)

# ╔═╡ 769d63db-36d1-437a-a3da-0bc9f6e14b69
# view in the usual mathematical way where order information is lost
Matrix(B)

# ╔═╡ 34e0f62a-9e14-4b9d-bad3-e6b23eb86c59
md""" ## Eigenstructure """

# ╔═╡ 3a777fc4-4770-4fb3-8074-2f66881a78ee
# Find eigenvalues of A. 
# destructuring via iteration
μ, V = eigen(A)  # type \mu + TAB

# ╔═╡ e63dfd51-6d85-47ad-9e07-d5164506ea91
# for stability, all eigenvalues must be non-positive

# ╔═╡ 2376d66b-8a50-46cf-88cb-1721c986c248


# ╔═╡ 5e05d040-475f-4f8b-b094-f52d2fd1511b


# ╔═╡ 5ed4eadc-396a-4e56-96f4-b9bbe605796a


# ╔═╡ d9f77a6e-dada-476c-9e7a-25676c34518a
diag(μ)

# ╔═╡ 1e3f4bd2-94cf-43a1-af98-11373a4d8561
# maximum timescale is related to the smallest negative eigenvalue
Tmax = maximum_timescale(μ)

# ╔═╡ 3628ccd7-38d8-45bc-a0b6-4d74c1cb7bd9
# water-mass fractions
a = watermass_fraction(μ, V, B)

# ╔═╡ 2739e74f-6c5d-4d2c-ac29-9b929c2fadeb
# see the water-mass fraction related to the first boundary of interest

# ╔═╡ 2175673e-5232-4804-84cb-0d5b11f31413
first(a)'

# ╔═╡ 01484ca5-ed33-4b94-b188-780e9e3ef8c7
# water-mass fraction from second source
last(a)'

# ╔═╡ c33d09fb-fbf8-43c9-8d4b-345d90e7b40f
Matrix(a) # all water-mass information concatenated

# ╔═╡ bcac98bd-a250-4497-9757-7d357128ba4a
# ╠═╡ disabled = true
#=╠═╡
        deldotJ = convergence(J)

  ╠═╡ =#

# ╔═╡ d0649e4c-6b78-43b5-aec6-3be449d03c09


# ╔═╡ d0c8e5d6-733f-44c3-a2f4-a776f3137222


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
# ╠═b85c6513-5a1a-4fdf-b4be-2efc3c1db830
# ╠═bbc2b198-ca1f-461d-a72c-e37695de357c
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
# ╠═852f36b7-170b-4d20-bd31-af6ae5c716a5
# ╠═8306d2c4-8d50-4309-add1-6d1eef56cd4a
# ╠═1e91d3e1-c26f-4118-9630-d654d352da76
# ╠═51d0e115-5859-4eab-8a91-b8193afd52b5
# ╟─b9f2165e-2d18-4179-a69f-ab0fc6ceb8b6
# ╟─2d21fdae-8d7d-4ef5-a447-9a2f37e695a4
# ╠═246b677a-24bc-41da-96d5-1bff248657b8
# ╠═cd6dc878-4442-430b-a263-3651719f2f11
# ╠═ff24a30f-56ee-4095-8bb3-0c7e4a72fe87
# ╠═4bc2c6b0-ad10-4035-be82-d02060d1b3d7
# ╠═e17220c4-d45d-4d36-a05f-c245393b05ef
# ╠═75f344b4-e273-4369-89dc-5ebfdb675d21
# ╠═095bd0d6-0249-4e26-b91a-ff488980e119
# ╠═fa9f454b-b7ca-4e37-a67c-a28ff91a5e11
# ╠═17bf7f50-78a7-4c7c-bc2f-4d9086dd2181
# ╠═08464a04-7652-4f24-978d-cd329e7fe0a7
# ╠═5dfddc9c-6313-4679-a994-15a771ee4a90
# ╠═e325d781-ae5c-4f64-a608-170b4df77882
# ╠═1e92642c-396f-4353-aa5c-8849cf26af1d
# ╠═2da2069f-4035-448f-9833-6fcd1156b2e0
# ╠═86076566-a96b-4faf-bdef-93b95733dcff
# ╟─c9abc24c-d3f2-4d64-8dfc-b0fdf42d1502
# ╠═6ed0c4d7-8fbb-439e-85fd-382e6b6e030f
# ╠═c5381568-9e1a-4763-8da3-86638b468ec4
# ╠═378e4e6c-d399-458d-85a9-23c8ceda2b43
# ╠═d344750e-e335-4e3c-baaa-a2937c2497df
# ╠═ba4789f3-7576-423b-9e94-abf4c3259eb4
# ╠═f51bfbed-0c3b-415c-9aef-80574b905b17
# ╟─f44447e3-5e8e-4fbc-b2cf-83176fb93c9f
# ╠═6897f4af-ca8e-43a7-b741-5f2dd48c97cb
# ╠═9656a0f3-59ff-4bf2-85aa-33a5b29fd7d9
# ╠═15e6cead-7de1-4cdd-ae84-f7537e789900
# ╠═6a25a144-0ccc-4604-85a7-b724eaa4cfed
# ╠═382db56a-d39b-4835-bf13-6dd0088b0b39
# ╠═1312b135-a267-4736-8e56-ff44bc7be59b
# ╟─59b47e9c-784a-4ed5-aeb6-79b5c756fff6
# ╠═e74b7cee-01d6-4984-be3d-16a71b350c99
# ╠═769d63db-36d1-437a-a3da-0bc9f6e14b69
# ╟─34e0f62a-9e14-4b9d-bad3-e6b23eb86c59
# ╠═3a777fc4-4770-4fb3-8074-2f66881a78ee
# ╠═e63dfd51-6d85-47ad-9e07-d5164506ea91
# ╠═2376d66b-8a50-46cf-88cb-1721c986c248
# ╠═5e05d040-475f-4f8b-b094-f52d2fd1511b
# ╠═5ed4eadc-396a-4e56-96f4-b9bbe605796a
# ╠═d9f77a6e-dada-476c-9e7a-25676c34518a
# ╠═1e3f4bd2-94cf-43a1-af98-11373a4d8561
# ╠═3628ccd7-38d8-45bc-a0b6-4d74c1cb7bd9
# ╠═2739e74f-6c5d-4d2c-ac29-9b929c2fadeb
# ╠═2175673e-5232-4804-84cb-0d5b11f31413
# ╠═01484ca5-ed33-4b94-b188-780e9e3ef8c7
# ╠═c33d09fb-fbf8-43c9-8d4b-345d90e7b40f
# ╠═bcac98bd-a250-4497-9757-7d357128ba4a
# ╠═d0649e4c-6b78-43b5-aec6-3be449d03c09
# ╠═d0c8e5d6-733f-44c3-a2f4-a776f3137222
