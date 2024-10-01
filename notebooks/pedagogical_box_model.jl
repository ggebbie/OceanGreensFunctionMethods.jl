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

# ╔═╡ 2fe46717-3f77-4afa-9e74-1ddb594e40ea
using Plots

# ╔═╡ cc363185-cdc4-47be-a926-5178e1535f0d
using Distributions

# ╔═╡ 5d30be92-5266-4700-ba7b-ac88a7f066e3
# define our own unit: the sverdrup
module UnitfulOcean; using Unitful; 
@unit sverdrup "Sv" Sverdrup (10^6)u"m^3/s" false;
end

# ╔═╡ 10b07d8a-aee4-4b64-b9eb-f22f408877ba
md"""
# Pedagogical Box Model 

Julia package to complement "A Review of Green's Function Methods in Ocean Circulation Models," by Haine et al. This package goes toward one of the stated goals of the manuscript, namely to make Green's Function methods accessible for learning purposes. Here, we also aim to make a Julia package that is useful and computationally efficient for research purposes."""


# ╔═╡ 27b7af71-e396-45b3-8723-8b2fc804a77f
md"""## Activate a reproducible project environment

Uses Julia's built-in package manager `Pkg.jl`  """

# ╔═╡ 28e6a6c1-4bdf-49aa-afdd-b27f1b88661b
Pkg.instantiate()

# ╔═╡ 07f01269-cfd8-4d3d-8d85-0b1132ff2005
md""" ## Load some helpful packages """


# ╔═╡ 01f84c5f-8881-401f-a0a8-8ae69385f9fe
#using MultipliableDimArrays

# ╔═╡ 39045ccd-fd9a-4d87-a2d9-79171a3366dc
plotly()

# ╔═╡ abe2697f-3bcd-49ae-bbcb-dd0a04c3f147
md""" ## Suggested modifiable circulation inputs"""

# ╔═╡ 6d11d809-9902-4a6a-b85e-18aed70e352f
md""" ## Define and label the boxes 

Using the outstanding [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl)"""

# ╔═╡ ccc6b783-6cca-4d03-8fca-f3c312316c34
# define "dimensions" to be used with `DimensionalData.jl` 
# permits numerical quantities to be bundled with their meta-data
model_dims = model_dimensions()

# ╔═╡ 82a174f2-1849-4d67-85dd-944c6e445d53
Nb = prod(size(model_dims)) # number of boxes

# ╔═╡ f3a7040d-2177-4693-a423-48e833718b43
Ny, Nz = size(model_dims) # size in each dimension

# ╔═╡ fdd4e823-bdb5-4f02-8de3-aade687a94c6
md""" ## Embed physical units with numerical values

We aren't just solving mathematical constructs, but we are interpreting physical scenarios. Use [Unitful.jl]("https://github.com/PainterQubits/Unitful.jl") to put units on all quantities. """

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

# ╔═╡ b9f2165e-2d18-4179-a69f-ab0fc6ceb8b6
md""" abyssal overturning rate $(@bind Ψ_abyssal Slider((2:40)Sv,show_value = true, default = 20Sv)) """

# ╔═╡ 2d21fdae-8d7d-4ef5-a447-9a2f37e695a4
md""" intermediate overturning rate $(@bind Ψ_intermediate Slider((2:40)Sv,show_value = true, default = 10Sv)) """

# ╔═╡ 246b677a-24bc-41da-96d5-1bff248657b8
md""" vertical diffusion (exchange flux) $(@bind Fv_exchange Slider((1:30)Sv,show_value = true, default = 5Sv)) """

# ╔═╡ c2a29255-e95a-4dd6-b97c-03a09337136e
md""" high latitude boundary exchange $(@bind Fb_high Slider((1:40)Sv,show_value = true, default = 20Sv)) """

# ╔═╡ b96b1c34-ef03-4874-a3f3-d5ade9a62c70
md""" mid-latitude boundary exchange $(@bind Fb_mid Slider((1:40)Sv,show_value = true, default = 10Sv)) """

# ╔═╡ 8306d2c4-8d50-4309-add1-6d1eef56cd4a
# set the units of a quantity, then use a pipe to convert units
#Vol0 = 1e16m^3 |> km^3 # uniform value of volume for all boxes
Vol0 = 300.0Sv*yr |> km^3 # use MATLAB value not manuscript value (5% difference)

# ╔═╡ c9c96f53-3fab-4591-91cd-911ba4c26329
Vol = fill(Vol0, model_dims)

# ╔═╡ 51d0e115-5859-4eab-8a91-b8193afd52b5
Vol' # take transpose or complex conjugate transpose to view more intuitively

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
deldotFm[At("Mid-latitudes"),At("Deep")]

# ╔═╡ 17bf7f50-78a7-4c7c-bc2f-4d9086dd2181
# it's ok if you don't remember the order of the dimensions
# here, get the three boxes in the Mid-latitudes
deldotFm[Meridional=At("Mid-latitudes")]

# ╔═╡ 08464a04-7652-4f24-978d-cd329e7fe0a7
# Given a tracer distribution C and volume fluxes Fv, find the tracer fluxes
# As an example, consider a randomly generated tracer field
C = rand(model_dims) # first generate a uniform U(0,1) tracer distribution

# ╔═╡ 100e928b-679b-48a3-b817-ac0dae73476b
 # extract tracer value by using geographic indices
 C[2,2]

# ╔═╡ 18f8bdce-9ce9-4e27-bf8e-53be37dc3fd0
# or extracer tracer using dimensional labels
C[Meridional=At("Mid-latitudes")]

# ╔═╡ df1cc59e-9e5f-48ea-b82f-65ab89b3e80a
Plots.heatmap(transpose(C),xflip=true)

# ╔═╡ 5dfddc9c-6313-4679-a994-15a771ee4a90
# then solve for tracer fluxes
J = advective_diffusive_flux(C, Fv)

# ╔═╡ e325d781-ae5c-4f64-a608-170b4df77882
# fluxes are stored in a structure that is organized by directionality
# again, use transpose to see the screen output in a reasonable order
J.poleward' 

# ╔═╡ 1e92642c-396f-4353-aa5c-8849cf26af1d
# greatest poleward fluxes in Thermocline
J.poleward[Meridional=At("Mid-latitudes")] # hit rightward arrow above to see flux values

# ╔═╡ 86076566-a96b-4faf-bdef-93b95733dcff
deldotJ = convergence(J) # tracer flux convergence

# ╔═╡ c9abc24c-d3f2-4d64-8dfc-b0fdf42d1502
md"""## Boundary conditions """

# ╔═╡ 20d23e54-8eec-4c6e-894a-d7d90d82ce54
# boundary exchange: define the locations affected by boundary fluxes
boundary_dims = boundary_dimensions()

# ╔═╡ 378e4e6c-d399-458d-85a9-23c8ceda2b43
# prescribe boundary volume fluxes
Fb = DimArray(hcat([Fb_high, Fb_mid]), boundary_dims) # boundary flux

# ╔═╡ d344750e-e335-4e3c-baaa-a2937c2497df
# example: set B_D (Dirichlet boundary conditions) to 1
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
# Instead, to view matrix in usual mathematical form, use `Matrix`
Matrix(A)

# ╔═╡ 6a25a144-0ccc-4604-85a7-b724eaa4cfed
# select of column of A corresponding to a tracer location
A[Vertical=At("Abyssal"),Meridional=At("Mid-latitudes")] # still displayed with info about spatial-locations even if out of order

# ╔═╡ 382db56a-d39b-4835-bf13-6dd0088b0b39
# select an entry of A, caution: first index=column, second index=row
A[5][5]

# ╔═╡ 1312b135-a267-4736-8e56-ff44bc7be59b
# or get the same information about one element using labels, but it gets long
A[Vertical=At("Deep"),Meridional=At("Mid-latitudes")][Vertical=At("Deep"),Meridional=At("Mid-latitudes")] 

# ╔═╡ 59b47e9c-784a-4ed5-aeb6-79b5c756fff6
md""" ## Construct boundary matrix """

# ╔═╡ e74b7cee-01d6-4984-be3d-16a71b350c99
# probe for B (boundary matrix)
B =  linear_probe(tracer_tendency, f, C, Fb, Vol)

# ╔═╡ 769d63db-36d1-437a-a3da-0bc9f6e14b69
# view in the usual mathematical way where order information is obscured
Matrix(B)

# ╔═╡ 34e0f62a-9e14-4b9d-bad3-e6b23eb86c59
md""" ## Eigenstructure """

# ╔═╡ 3a777fc4-4770-4fb3-8074-2f66881a78ee
# Find eigenvalues of A. 
# destructuring via iteration
μ, V = eigen(A)  # type \mu + TAB

# ╔═╡ e63dfd51-6d85-47ad-9e07-d5164506ea91
# for stability, all eigenvalues must be non-positive

# ╔═╡ d9f77a6e-dada-476c-9e7a-25676c34518a
diag(μ)

# ╔═╡ c191889e-b3eb-4839-b494-8fad1f0ed9ce
# real part of all eigenvalues is negative
plot(real.(diag(μ)),xlabel="eigenvalue i",ylabel="μᵢ",legend=false)

# ╔═╡ 1e3f4bd2-94cf-43a1-af98-11373a4d8561
# maximum timescale is related to the smallest negative eigenvalue
Tmax = maximum_timescale(μ)

# ╔═╡ 6eac27ef-647c-4884-aaf3-69f6705da3a8
md"""## Tracer histories """

# ╔═╡ e9653ffb-4da5-4b45-948d-3656dbb6df66
projectdir()

# ╔═╡ a45c8594-9fc7-46c2-833d-c44ece6648e5
BD = read_transient_tracer_histories() # Dirichlet boundary conditions

# ╔═╡ 96ab52bc-1863-4056-ba32-894747f6ae0e
BD_iodine129 = read_iodine129_history() # iodine-129

# ╔═╡ 6f979bb9-733d-4981-9a53-d75162cbd372
md""" Choose tracers """

# ╔═╡ f53b4b2f-cda2-45a2-96f8-2dd348bc3c1f
md""" $(@bind use_CFC11 CheckBox(default=true)) CFC-11 $(@bind use_CFC12 CheckBox(default=true)) CFC-12  $(@bind use_SF6 CheckBox(default=true)) SF₆ $(@bind use_argon39 CheckBox(default=true)) ³⁹Ar $(@bind use_iodine129 CheckBox(default=true)) ¹²⁹I"""

# ╔═╡ cf38b164-4414-4344-824e-68a09cc38f6b
md""" Source history """

# ╔═╡ ea2e8fe1-94cc-4cb6-bb88-05c708eac5a3
md""" $(@bind mbound Select(meridional_names()[1:2])) $(@bind vbound Select([vertical_names()[1]]))""" 


# ╔═╡ e34ae847-d82e-49f4-aa22-6753596c4ea0
begin
	tlims = (1930yr,2015yr)
	tvals = tlims[1]:1yr:tlims[2]
	source_plot = plot(xlims=(1930yr,2015yr),
		yscale=:log10,
		ylims = (1e-1,1e3),
		legend = :outerbottomright,
		titlelabel="")	

	if use_CFC11 
		plot!(tvals,
			[tracer_source_history(t,:CFC11NH, 0.75, BD)[At(mbound),At(vbound)] for t in tvals],
			label="CFC-11")
	end
	
	if use_CFC12 
		plot!(tvals,
			[tracer_source_history(t,:CFC12NH, 0.75, BD)[At(mbound),At(vbound)] for t in tvals],
			label="CFC-12")
	end
	
	if use_SF6
		plot!(tvals,
			[tracer_source_history(t,:SF6NH, 0.75, BD)[At(mbound),At(vbound)] for t in tvals],
			label="SF₆")
	end
	
	if use_argon39
		plot!(tvals,
			[tracer_source_history(t,:argon39, 1)[At(mbound),At(vbound)] for t in tvals],
			label="³⁹Ar")
	end

	if use_iodine129
		plot!(tvals,
			[tracer_source_history(t,:iodine129, 0.25, BD_iodine129)[At(mbound),At(vbound)] for t in tvals],
			label="¹²⁹I")
	end
	
	title!("")
	source_plot
end

# ╔═╡ 897deef3-d754-4ca4-8c6f-00b67313a5a0
md""" Interior evolution """

# ╔═╡ 1ecd9ce2-cea7-417e-b965-24784cd0f563
md""" $(@bind mbox1 Select(meridional_names())) $(@bind vbox1 Select(vertical_names())) """

# ╔═╡ ec1439f9-7f02-439a-ac70-d67869cdae35
begin 
	

	tracer_plot = plot(xlims=(1930yr,2015yr),
		yscale=:log10,
		ylims = (1e-1,1e3),
		legend = :outerbottomright,
		title = mbox1*", "*vbox1,
		titlefontsize=6)

	if use_CFC11 
		tlist = (1900.25:0.25:2015.0)yr
		ct = tracer_timeseries(:CFC11NH, A, B, tlist, mbox1, vbox1, BD = BD)
		plot!(tlist, ct, label="CFC-11")
	end
	
	if use_CFC12 
		tlist = (1900.25:0.25:2015.0)yr
		ct = tracer_timeseries(:CFC12NH, A, B, tlist, mbox1, vbox1, BD = BD)
		plot!(tlist, ct, label="CFC-12")
	end
	
	if use_SF6 
		tlist = (1900.25:0.25:2015.0)yr
		ct = tracer_timeseries(:SF6NH, A, B, tlist, mbox1, vbox1, BD = BD)
		plot!(tlist, ct, label="SF₆")
	end
	
	if use_argon39 
		tlist = (1750.0:0.25:2015.0)yr # start early enough to get equilibrium
		ct = tracer_timeseries(:argon39, A, B, tlist, mbox1, vbox1, halflife = 269yr)
		plot!(tlist, ct, label="³⁹Ar")
	end

	if use_iodine129 
		tlist = (1957.5:0.25:2015.0)yr
		ct = tracer_timeseries(:iodine129, A, B, tlist, mbox1, vbox1, BD = BD_iodine129, halflife = 1.57e7yr)
		plot!(tlist, ct, label="¹²⁹I")
	end
	
	tracer_plot
end


# ╔═╡ 11eb59cf-de62-4fb4-9963-defe594e6b92
md""" ## Transport matrix diagnostics """

# ╔═╡ 3628ccd7-38d8-45bc-a0b6-4d74c1cb7bd9
# water-mass fractions
a = watermass_fraction(μ, V, B)

# ╔═╡ 2175673e-5232-4804-84cb-0d5b11f31413
# see the water-mass fraction related to the first boundary of interest
first(a)'

# ╔═╡ 01484ca5-ed33-4b94-b188-780e9e3ef8c7
# water-mass fraction from second source
last(a)'

# ╔═╡ c33d09fb-fbf8-43c9-8d4b-345d90e7b40f
Matrix(a) # all water-mass information concatenated

# ╔═╡ e5841ad8-dfb9-47d5-bcb0-0f7448f43645
a

# ╔═╡ 0071aa97-27c3-469f-b1bb-e07337489f0e
begin
	msource1 = "High latitudes"
	vsource1 = "Thermocline"
	Plots.heatmap(transpose(a[At(msource1),At(vsource1)]),
		title="Water mass fraction: "*msource1*" "*vsource1,
		titlefontsize=6,
		xflip=true,
		color=:heat,
		clims=(0.25,0.75))
end

# ╔═╡ 0b804941-fed3-4830-980b-8d383d473858
 a[At(msource1),At(vsource1)]

# ╔═╡ 5786b2d4-d049-4119-8e1c-5ecf8e8c683e
begin
	msource2 = "Mid-latitudes"
	vsource2 = "Thermocline"
	Plots.heatmap(transpose(a[At(msource2),At(vsource2)]),
		title="Water mass fraction: "*msource2*" "*vsource2,
		titlefontsize=6,
		xflip=true,
		color=:heat,
		clims=(0.25,0.75))
end

# ╔═╡ cf5bb364-5336-4dd1-8bb6-6e3f944673bf
begin
	Γ = mean_age(μ, V, B)
	
	Plots.heatmap(transpose(Γ),
		title="Mean Age ["*string(unit(first(Γ)))*"]",
		titlefontsize=6,
		xflip=true,
		color=:heat,
		clims=(0yr,200yr))
end

# ╔═╡ 4021feb1-36ac-42f6-a5f6-391c0f064dc7
# very similar values - matches with MATLAB results
Δ = ttd_width(μ, V, B)

# ╔═╡ 93c9614e-70a1-49ef-933b-b86fec342597
md"""### Green's functions """

# ╔═╡ 96240170-eacb-4d5a-9316-eb6615a78f0a
md""" Select interior box for diagnostics """

# ╔═╡ 7a71a95a-8523-4cb8-9f69-00bf374acf67
md""" $(@bind mbox Select(meridional_names())) $(@bind vbox Select(vertical_names())) """

# ╔═╡ cd492316-d6b2-4645-80ba-c5817ec5877c
Δτ = 0.25yr # time resolution

# ╔═╡ 4c258084-da30-4393-b844-c379c9e79efd
τ = 0yr:Δτ:2000yr # list of time lags

# ╔═╡ 589ab455-2e9c-47d6-abd7-f89f367a5ed5
G(t) = greens_function(t,A) # a closure that captures A

# ╔═╡ c122abb6-185c-4894-a2c4-8ab6224e83d2
G′(t) = boundary_propagator(t,A,B,alg=:forward) # type G + \prime + TAB

# ╔═╡ 595fba3f-65ec-461f-a257-92456d4f94a0
# global (or total) TTD
𝒢(t) = global_ttd(t,A,B) # type \scr + G + TAB

# ╔═╡ 00902450-ceb7-4c33-be7e-906502990813
# a list comprehension
ttd1 = [G′(τ[i])[Meridional=At("High latitudes"),Vertical=At("Thermocline")][Meridional=At(mbox),Vertical=At(vbox)] for i in eachindex(τ)]

# ╔═╡ c2a38bc2-ef10-4fe5-8642-857f9acdadd7
# could be written as a for loop instead
ttd2 = [G′(τ[i])[Meridional=At("Mid-latitudes"),Vertical=At("Thermocline")][Meridional=At(mbox),Vertical=At(vbox)] for i in eachindex(τ)] 

# ╔═╡ 8d69c375-6a0c-400e-af85-3013a364fa1d
ttd_global = [𝒢(τ[i])[Meridional=At(mbox),Vertical=At(vbox)] for i in eachindex(τ)] 

# ╔═╡ 09a85965-d1dc-47a3-9eba-dd1dc129db36
Γ_ = Γ[Meridional=At(mbox),Vertical=At(vbox)] 

# ╔═╡ 19ef1da1-9b1a-4300-83aa-bb503027122b
Δ_ = Δ[Meridional=At(mbox),Vertical=At(vbox)]

# ╔═╡ a183e31d-8bab-46e0-a6b1-0a181c5f0f69
a1 = a[Meridional=At("High latitudes"),Vertical=At("Thermocline")][Meridional=At(mbox),Vertical=At(vbox)]

# ╔═╡ 9537166f-054f-441e-a001-3ba59a4b59e0
a2 = a[Meridional=At("Mid-latitudes"),Vertical=At("Thermocline")][Meridional=At(mbox),Vertical=At(vbox)]

# ╔═╡ fd907198-8e2e-4296-b640-c0aebbd0a796
G_inversegaussian = TracerInverseGaussian(Γ_, Δ_)

# ╔═╡ 1bb59934-17be-40d3-b227-b73bb1b9c4df
ttd_inversegaussian = pdf.(G_inversegaussian,τ)

# ╔═╡ 1f163cd9-5db4-4b1b-a8df-a967156e6b59
# label names for two boundary boxes
boxname = ["("*meridional_names()[i][1:6]*","*vertical_names()[1][1:6]*")" for i in 1:2]

# ╔═╡ e8fabe44-3a7d-47fc-84af-02baebf5f45a
begin 
	# to do: put plotting into functions
	p = plot(τ,
		normalized_exponential_decay.(τ,Tmax),
		linestyle = :dash,
		yscale = :log10,
		ylabel = "G′(τ)",
		xlabel = "τ",
		label = "Tmax",
		legend = :topright,
		titlefontsize = 6,
		title = mbox*", "*vbox,
		xlims = (0yr,400yr),
		ylims = (1e-4/yr,1e-1/yr))
	
	plot!([Γ_,Γ_],
		[1e-4,1e-2]/yr,
		label="Γ")	
	
	plot!([Γ_ + Δ_/2,
		Γ_ - Δ_/2],
		[1e-4,1e-4]/yr,
		width=4,
		color=:grey,
		label="Δ")
	
	plot!(τ,ttd1,label="TTD "*boxname[1],width=4*a1)
	plot!(τ,ttd2,label="TTD "*boxname[2],width=4*a2)
	plot!(τ,ttd_global,label="Total TTD",width=4*(a1+a2),color=:black)
	plot!(τ,ttd_inversegaussian,label="Fitted inverse Gaussian")
end

# ╔═╡ 7c725552-883e-4fb3-b22e-292518913dfd
md""" ## Adjoint Green's functions """

# ╔═╡ 4bd0734f-d3f9-49e5-a7cb-ef719acb23f4
md""" $(@bind mbox_adj Select(meridional_names())) $(@bind vbox_adj Select(vertical_names())) """

# ╔═╡ ab31341c-ff59-41bc-8a7f-752931bb8e9d
# † is invalid in Julia as an identifier 
G′dagger(t) = boundary_propagator(t,A,B,alg=:adjoint) # type G + \prime + TAB

# ╔═╡ 1df15962-dd41-4f07-82c8-37d2d60511fb
ttd1_adj = [G′dagger(τ[i])[At(mbox_adj),At(vbox_adj)][At("High latitudes"),At("Thermocline")] for i in eachindex(τ)]

# ╔═╡ 48449ccf-df3f-4b71-a160-53d39baa9a90
ttd2_adj = [G′dagger(τ[i])[At(mbox_adj),At(vbox_adj)][At("Mid-latitudes"),At("Thermocline")] for i in eachindex(τ)]

# ╔═╡ d7822a1b-a780-4df8-9f71-cc3de36ed4c0
a1_adj = watermass_fraction(μ, V, B, alg=:adjoint)[At(mbox_adj),At(vbox_adj)][At("High latitudes"),At("Thermocline")]

# ╔═╡ c6fd37b2-39cc-4c34-b021-6483285e6d56
a2_adj = watermass_fraction(μ, V, B, alg=:adjoint)[At(mbox_adj),At(vbox_adj)][At("Mid-latitudes"),At("Thermocline")]

# ╔═╡ b3522980-6beb-4e05-901d-0859c7a8cb58
# global adjoint TTD
𝒢dagger(t) = global_ttd(t,A,B,alg=:adjoint)

# ╔═╡ b719ab41-4226-40c7-9682-5385d076dc7a
𝒢dagger(1yr)

# ╔═╡ 257c6649-d003-42bc-9e17-0c33b7cd304c
ttd_global_adjoint = [𝒢dagger(τ[i])[At(mbox_adj),At(vbox_adj)] for i in eachindex(τ)] 

# ╔═╡ cf82fade-07ac-4aa9-bd06-7a10820a724f
#Γ_adjoint = adjoint_mean_age(A,B)[At(mbox_adj),At(vbox_adj)]
Γ_adjoint = mean_age(μ, V, B, alg=:adjoint)[At(mbox_adj),At(vbox_adj)]

# ╔═╡ c6460013-d800-4280-97db-50c5aa84e709
Δ_adjoint = ttd_width(μ, V, B,alg=:adjoint)[At(mbox_adj),At(vbox_adj)]

# ╔═╡ 4e0ce7d3-a1fd-4995-83d0-bdc74bc5e339
G_inversegaussian_adjoint = TracerInverseGaussian(Γ_adjoint, Δ_adjoint)

# ╔═╡ f861d37b-427b-4c12-b0ff-c55be4d82523
ttd_inversegaussian_adjoint = pdf.(G_inversegaussian_adjoint,τ)

# ╔═╡ c7a4d285-25e3-42eb-8e5b-7967aad1a366
begin 
	# to do: put plotting into functions
	p_adj = plot(τ,
		normalized_exponential_decay.(τ,Tmax),
		linestyle = :dash,
		yscale = :log10,
		ylabel = "G′†",
		xlabel = "τ",
		label = "Tmax",
		legend = :right,
		titlefontsize = 6,
		title = mbox_adj*", "*vbox_adj,
		xlims = (0yr,400yr),
		ylims = (1e-4/yr,1e-1/yr))
	
	plot!([Γ_adjoint,Γ_adjoint],
		[1e-4,1e-2]/yr,
		label="Γ")	
	
	plot!([Γ_adjoint + Δ_adjoint/2,
		Γ_adjoint - Δ_adjoint/2],
		[1e-4,1e-4]/yr,
		width=4,
		color=:grey,
		label="Δ")
	
	plot!(τ,ttd1_adj,label="TTD "*boxname[1],width=4*a1_adj)
	plot!(τ,ttd2_adj,label="TTD "*boxname[2],width=4*a2_adj)
	plot!(τ,ttd_global_adjoint,label="Total TTD",width=4*(a1_adj+a2_adj),color=:black)
	plot!(τ,ttd_inversegaussian_adjoint,label="Fitted inverse Gaussian")
end

# ╔═╡ 13d659ac-d820-404e-bdcb-c66b05381309
md""" ## Residence time distributions """

# ╔═╡ 42ca866d-9c14-4761-9d0f-131870e25d9e
md""" Destination box $(@bind mbox_destination Select(meridional_names()[1:2])) $(@bind vbox_destination Select([vertical_names()[1]]))"""

# ╔═╡ 0a62e096-f375-4053-bc88-7ef89ce1173a
RTD(t) = residence_time(t,A,B)

# ╔═╡ f6f550a5-d04d-4d2a-89e7-484734370416
rtd1 = [RTD(τ[i])[At("High latitudes"),At("Thermocline")][At(mbox_destination),At(vbox_destination)] for i in eachindex(τ)]

# ╔═╡ 6d1b4753-aabb-4274-a5fb-de26270c4378
rtd2 = [RTD(τ[i])[At("Mid-latitudes"),At("Thermocline")][At(mbox_destination),At(vbox_destination)] for i in eachindex(τ)]

# ╔═╡ 29c38299-d49f-422c-9065-2faa9d2db491
a_residence = watermass_fraction(μ, V, B, alg=:residence)

# ╔═╡ 31f0f55f-5e2c-42da-9598-9b0bc1ce262f
a_residence1 = a_residence[At("High latitudes"),At("Thermocline")][At(mbox_destination),At(vbox_destination)]

# ╔═╡ 9f82d814-f4fb-4184-8177-13f8fe7eceef
a_residence2 = a_residence[At("Mid-latitudes"),At("Thermocline")][At(mbox_destination),At(vbox_destination)]

# ╔═╡ 7e8cd9ef-60cf-4909-93ef-643be08e4bc2
Γ_residence = mean_age(μ, V, B, alg=:residence)

# ╔═╡ 989f966a-2078-408d-b3ca-5f4fa332f8b6
Δ_residence = ttd_width(μ, V, B, alg=:residence)

# ╔═╡ 025e7a9d-d587-44d6-ba0c-1343ad18121a
begin 
	p_source = plot(τ,
		normalized_exponential_decay.(τ,Tmax),
		linestyle = :dash,
		yscale = :log10,
		ylabel = "R(τ)",
		xlabel = "τ",
		label = "Tmax",
		legend = :topright,
		titlefontsize = 6,
		title = "Destination: "*mbox_destination*", "*" Thermocline",
		xlims = (0yr,400yr),
		ylims = (1e-4/yr,1e-1/yr)) 

	plot!([Γ_residence,Γ_residence],
		[1e-4,1e-2]/yr,
		label="Γ")	
	
	plot!([Γ_residence + Δ_residence/2,
		Γ_residence - Δ_residence/2],
		[1e-4,1e-4]/yr,
		width=4,
		color=:grey,
		label="Δ")

	plot!(τ,rtd1,label="R: src "*boxname[1],
		legendfontsize=6,
		width=8*a_residence1)
	plot!(τ,rtd2,label="R: src "*boxname[2],
		legendfontsize=6,
		width=8*a_residence2)
end

# ╔═╡ 58701b47-1669-484c-ab88-904f31fedb97
sum(Matrix(a_residence)[:]) # a test that all mass is taken into account

# ╔═╡ 55325cc7-0b71-44c9-ac41-837f6451647b
md""" ## Path density """

# ╔═╡ 7550d6cf-9bfb-470c-b74f-880a298ab19e
md""" $(@bind mbox_pdensity Select(meridional_names())) $(@bind vbox_pdensity Select(vertical_names())) """

# ╔═╡ 175e67c2-c458-4bde-87a1-1fdf49e923c5
path_density_11 = [path_density(μ, V, B, τ[i], mbox_pdensity, vbox_pdensity)[At("High latitudes"),At("Thermocline")][At("High latitudes"),At("Thermocline")] for i in eachindex(τ)]

# ╔═╡ 00ded914-2b67-42e8-a913-d936f4442da8
path_density_12 = [path_density(μ, V, B, τ[i], mbox_pdensity, vbox_pdensity)[At("High latitudes"),At("Thermocline")][At("Mid-latitudes"),At("Thermocline")] for i in eachindex(τ)]

# ╔═╡ 696608f7-49c1-45e9-8522-ab222b0da243
path_density_21 = [path_density(μ, V, B, τ[i], mbox_pdensity, vbox_pdensity)[At("Mid-latitudes"),At("Thermocline")][At("High latitudes"),At("Thermocline")] for i in eachindex(τ)]

# ╔═╡ 7465a888-bdf4-4a9f-a111-d53fb23a63bc
path_density_22 = [path_density(μ, V, B, τ[i], mbox_pdensity, vbox_pdensity)[At("Mid-latitudes"),At("Thermocline")][At("Mid-latitudes"),At("Thermocline")] for i in eachindex(τ)]

# ╔═╡ d4504bca-c858-4434-ba95-df5bdae8d61c
begin 
	plot_pdensity = plot(τ,
		normalized_exponential_decay.(τ,Tmax),
		linestyle = :dash,
		yscale = :log10,
		ylabel = "η(τ)",
		xlabel = "τ",
		label = "Tmax",
		legend = :right,
		titlefontsize = 6,
		title = mbox_pdensity*", "*vbox_pdensity,
		xlims = (0yr,400yr),
		ylims = (1e-4/yr,1e-1/yr)) 

	plot!(τ,path_density_11,label=boxname[1]*" to "*boxname[1],width=2)
	plot!(τ,path_density_12,label=boxname[1]*" to "*boxname[2],width=2)
	plot!(τ,path_density_21,label=boxname[2]*" to "*boxname[1],width=2,linestyle= :dash)
	plot!(τ,path_density_22,label=boxname[2]*" to "*boxname[2],width=2, linestyle= :dash)
end

# ╔═╡ Cell order:
# ╟─10b07d8a-aee4-4b64-b9eb-f22f408877ba
# ╠═27b7af71-e396-45b3-8723-8b2fc804a77f
# ╠═8f520c8b-19d7-48a8-be9f-3f167f07d188
# ╠═c536e9f3-0457-499e-958c-384d6e388ef9
# ╠═28e6a6c1-4bdf-49aa-afdd-b27f1b88661b
# ╟─07f01269-cfd8-4d3d-8d85-0b1132ff2005
# ╠═de3c6443-5ca1-4e97-82c8-5c4c9f204480
# ╠═69147ae0-1c89-48a6-831b-ff325a984817
# ╠═b85c6513-5a1a-4fdf-b4be-2efc3c1db830
# ╠═bbc2b198-ca1f-461d-a72c-e37695de357c
# ╠═97df5706-1829-419a-b96b-5dbb4d704434
# ╠═157462cf-b6f9-4de3-80f6-a3f846b5ea1a
# ╠═0a9a45e2-a561-4a21-afb9-b96ec884de4a
# ╠═2fe46717-3f77-4afa-9e74-1ddb594e40ea
# ╠═cc363185-cdc4-47be-a926-5178e1535f0d
# ╠═01f84c5f-8881-401f-a0a8-8ae69385f9fe
# ╠═39045ccd-fd9a-4d87-a2d9-79171a3366dc
# ╟─abe2697f-3bcd-49ae-bbcb-dd0a04c3f147
# ╟─b9f2165e-2d18-4179-a69f-ab0fc6ceb8b6
# ╟─2d21fdae-8d7d-4ef5-a447-9a2f37e695a4
# ╟─246b677a-24bc-41da-96d5-1bff248657b8
# ╟─c2a29255-e95a-4dd6-b97c-03a09337136e
# ╟─b96b1c34-ef03-4874-a3f3-d5ade9a62c70
# ╟─6d11d809-9902-4a6a-b85e-18aed70e352f
# ╠═ccc6b783-6cca-4d03-8fca-f3c312316c34
# ╠═82a174f2-1849-4d67-85dd-944c6e445d53
# ╠═f3a7040d-2177-4693-a423-48e833718b43
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
# ╠═c9c96f53-3fab-4591-91cd-911ba4c26329
# ╠═51d0e115-5859-4eab-8a91-b8193afd52b5
# ╠═cd6dc878-4442-430b-a263-3651719f2f11
# ╠═ff24a30f-56ee-4095-8bb3-0c7e4a72fe87
# ╠═4bc2c6b0-ad10-4035-be82-d02060d1b3d7
# ╠═e17220c4-d45d-4d36-a05f-c245393b05ef
# ╠═75f344b4-e273-4369-89dc-5ebfdb675d21
# ╠═095bd0d6-0249-4e26-b91a-ff488980e119
# ╠═fa9f454b-b7ca-4e37-a67c-a28ff91a5e11
# ╠═17bf7f50-78a7-4c7c-bc2f-4d9086dd2181
# ╠═08464a04-7652-4f24-978d-cd329e7fe0a7
# ╠═100e928b-679b-48a3-b817-ac0dae73476b
# ╠═18f8bdce-9ce9-4e27-bf8e-53be37dc3fd0
# ╠═df1cc59e-9e5f-48ea-b82f-65ab89b3e80a
# ╠═5dfddc9c-6313-4679-a994-15a771ee4a90
# ╠═e325d781-ae5c-4f64-a608-170b4df77882
# ╠═1e92642c-396f-4353-aa5c-8849cf26af1d
# ╠═86076566-a96b-4faf-bdef-93b95733dcff
# ╟─c9abc24c-d3f2-4d64-8dfc-b0fdf42d1502
# ╠═20d23e54-8eec-4c6e-894a-d7d90d82ce54
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
# ╠═d9f77a6e-dada-476c-9e7a-25676c34518a
# ╠═c191889e-b3eb-4839-b494-8fad1f0ed9ce
# ╠═1e3f4bd2-94cf-43a1-af98-11373a4d8561
# ╟─6eac27ef-647c-4884-aaf3-69f6705da3a8
# ╠═e9653ffb-4da5-4b45-948d-3656dbb6df66
# ╠═a45c8594-9fc7-46c2-833d-c44ece6648e5
# ╠═96ab52bc-1863-4056-ba32-894747f6ae0e
# ╟─6f979bb9-733d-4981-9a53-d75162cbd372
# ╟─f53b4b2f-cda2-45a2-96f8-2dd348bc3c1f
# ╟─cf38b164-4414-4344-824e-68a09cc38f6b
# ╟─ea2e8fe1-94cc-4cb6-bb88-05c708eac5a3
# ╠═e34ae847-d82e-49f4-aa22-6753596c4ea0
# ╟─897deef3-d754-4ca4-8c6f-00b67313a5a0
# ╟─1ecd9ce2-cea7-417e-b965-24784cd0f563
# ╠═ec1439f9-7f02-439a-ac70-d67869cdae35
# ╟─11eb59cf-de62-4fb4-9963-defe594e6b92
# ╠═3628ccd7-38d8-45bc-a0b6-4d74c1cb7bd9
# ╠═2175673e-5232-4804-84cb-0d5b11f31413
# ╠═01484ca5-ed33-4b94-b188-780e9e3ef8c7
# ╠═c33d09fb-fbf8-43c9-8d4b-345d90e7b40f
# ╠═e5841ad8-dfb9-47d5-bcb0-0f7448f43645
# ╠═0071aa97-27c3-469f-b1bb-e07337489f0e
# ╠═0b804941-fed3-4830-980b-8d383d473858
# ╠═5786b2d4-d049-4119-8e1c-5ecf8e8c683e
# ╠═cf5bb364-5336-4dd1-8bb6-6e3f944673bf
# ╟─4021feb1-36ac-42f6-a5f6-391c0f064dc7
# ╟─93c9614e-70a1-49ef-933b-b86fec342597
# ╟─96240170-eacb-4d5a-9316-eb6615a78f0a
# ╟─7a71a95a-8523-4cb8-9f69-00bf374acf67
# ╠═e8fabe44-3a7d-47fc-84af-02baebf5f45a
# ╠═cd492316-d6b2-4645-80ba-c5817ec5877c
# ╠═4c258084-da30-4393-b844-c379c9e79efd
# ╠═589ab455-2e9c-47d6-abd7-f89f367a5ed5
# ╠═c122abb6-185c-4894-a2c4-8ab6224e83d2
# ╠═595fba3f-65ec-461f-a257-92456d4f94a0
# ╠═00902450-ceb7-4c33-be7e-906502990813
# ╠═c2a38bc2-ef10-4fe5-8642-857f9acdadd7
# ╠═8d69c375-6a0c-400e-af85-3013a364fa1d
# ╠═09a85965-d1dc-47a3-9eba-dd1dc129db36
# ╠═19ef1da1-9b1a-4300-83aa-bb503027122b
# ╠═a183e31d-8bab-46e0-a6b1-0a181c5f0f69
# ╠═9537166f-054f-441e-a001-3ba59a4b59e0
# ╠═fd907198-8e2e-4296-b640-c0aebbd0a796
# ╠═1bb59934-17be-40d3-b227-b73bb1b9c4df
# ╠═1f163cd9-5db4-4b1b-a8df-a967156e6b59
# ╟─7c725552-883e-4fb3-b22e-292518913dfd
# ╟─4bd0734f-d3f9-49e5-a7cb-ef719acb23f4
# ╠═c7a4d285-25e3-42eb-8e5b-7967aad1a366
# ╠═ab31341c-ff59-41bc-8a7f-752931bb8e9d
# ╠═1df15962-dd41-4f07-82c8-37d2d60511fb
# ╠═48449ccf-df3f-4b71-a160-53d39baa9a90
# ╠═d7822a1b-a780-4df8-9f71-cc3de36ed4c0
# ╠═c6fd37b2-39cc-4c34-b021-6483285e6d56
# ╠═b3522980-6beb-4e05-901d-0859c7a8cb58
# ╠═b719ab41-4226-40c7-9682-5385d076dc7a
# ╠═257c6649-d003-42bc-9e17-0c33b7cd304c
# ╠═cf82fade-07ac-4aa9-bd06-7a10820a724f
# ╠═c6460013-d800-4280-97db-50c5aa84e709
# ╠═4e0ce7d3-a1fd-4995-83d0-bdc74bc5e339
# ╠═f861d37b-427b-4c12-b0ff-c55be4d82523
# ╟─13d659ac-d820-404e-bdcb-c66b05381309
# ╟─42ca866d-9c14-4761-9d0f-131870e25d9e
# ╠═025e7a9d-d587-44d6-ba0c-1343ad18121a
# ╠═0a62e096-f375-4053-bc88-7ef89ce1173a
# ╠═f6f550a5-d04d-4d2a-89e7-484734370416
# ╠═6d1b4753-aabb-4274-a5fb-de26270c4378
# ╠═29c38299-d49f-422c-9065-2faa9d2db491
# ╠═31f0f55f-5e2c-42da-9598-9b0bc1ce262f
# ╠═9f82d814-f4fb-4184-8177-13f8fe7eceef
# ╠═7e8cd9ef-60cf-4909-93ef-643be08e4bc2
# ╠═989f966a-2078-408d-b3ca-5f4fa332f8b6
# ╠═58701b47-1669-484c-ab88-904f31fedb97
# ╟─55325cc7-0b71-44c9-ac41-837f6451647b
# ╟─7550d6cf-9bfb-470c-b74f-880a298ab19e
# ╠═d4504bca-c858-4434-ba95-df5bdae8d61c
# ╠═175e67c2-c458-4bde-87a1-1fdf49e923c5
# ╠═00ded914-2b67-42e8-a913-d936f4442da8
# ╠═696608f7-49c1-45e9-8522-ab222b0da243
# ╠═7465a888-bdf4-4a9f-a111-d53fb23a63bc
