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

Uses Julia's built-in package manager `Pkg.jl` """

# ╔═╡ 28e6a6c1-4bdf-49aa-afdd-b27f1b88661b
Pkg.instantiate()

# ╔═╡ 07f01269-cfd8-4d3d-8d85-0b1132ff2005
md""" ## Load some helpful packages """


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
Vol0 = 1e16m^3 |> km^3 # uniform value of volume for all boxes

# ╔═╡ 1e91d3e1-c26f-4118-9630-d654d352da76
# If your screen is big enough, you should see a labeled, 3 x 3 table of volume values
Vol = DimArray(fill(Vol0, Ny, Nz), model_dims)

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
deldotFm[At("2 Mid-latitudes"),At("2 Deep")]

# ╔═╡ 17bf7f50-78a7-4c7c-bc2f-4d9086dd2181
# it's ok if you don't remember the order of the dimensions
# here, get the three boxes in the Mid-latitudes
deldotFm[Meridional=At("2 Mid-latitudes")]

# ╔═╡ 08464a04-7652-4f24-978d-cd329e7fe0a7
# Given a tracer distribution C and volume fluxes Fv, find the tracer fluxes
# As an example, consider a randomly generated tracer field
C = rand(model_dims) # first generate a uniform U(0,1) tracer distribution

# ╔═╡ 100e928b-679b-48a3-b817-ac0dae73476b
 # extract tracer value by using geographic indices
 C[2,2]

# ╔═╡ 18f8bdce-9ce9-4e27-bf8e-53be37dc3fd0
# or extracer tracer using dimensional labels
C[Meridional=At("2 Mid-latitudes")]

# ╔═╡ df1cc59e-9e5f-48ea-b82f-65ab89b3e80a
Plots.heatmap(transpose(C),yflip=true)

# ╔═╡ 5dfddc9c-6313-4679-a994-15a771ee4a90
# then solve for tracer fluxes
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

# ╔═╡ 20d23e54-8eec-4c6e-894a-d7d90d82ce54
# boundary exchange: define the locations affected by boundary fluxes
boundary_dims = boundary_dimensions()

# ╔═╡ 378e4e6c-d399-458d-85a9-23c8ceda2b43
# prescribe boundary volume fluxes
Fb = DimArray(hcat([Fb_high, Fb_mid]), boundary_dims) # boundary flux

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
# Instead, to view matrix in usual mathematical form, use `Matrix`
Matrix(A)

# ╔═╡ 6a25a144-0ccc-4604-85a7-b724eaa4cfed
# select of column of A corresponding to a tracer location
A[Vertical=At("3 Abyssal"),Meridional=At("2 Mid-latitudes")] # still displayed with info about spatial-locations

# ╔═╡ 382db56a-d39b-4835-bf13-6dd0088b0b39
# select an entry of A, caution: first index=column, second index=row
A[5][5]

# ╔═╡ 1312b135-a267-4736-8e56-ff44bc7be59b
# or get the same information using labels, but it gets long
A[Vertical=At("2 Deep"),Meridional=At("2 Mid-latitudes")][Vertical=At("2 Deep"),Meridional=At("2 Mid-latitudes")] 

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

# ╔═╡ a45c8594-9fc7-46c2-833d-c44ece6648e5
BD = read_tracer_histories() # Dirichlet boundary conditions

# ╔═╡ 6f979bb9-733d-4981-9a53-d75162cbd372
md""" Choose tracers """

# ╔═╡ f53b4b2f-cda2-45a2-96f8-2dd348bc3c1f
md""" $(@bind use_CFC11 CheckBox(default=true)) CFC-11"""

# ╔═╡ ef360b00-8f37-45b4-9d95-4992922ee03d
md""" $(@bind use_CFC12 CheckBox(default=true)) CFC-12"""

# ╔═╡ fc1d1240-e96a-4d5e-a5f4-773d84074e22
md""" $(@bind use_SF6 CheckBox(default=true)) SF₆"""

# ╔═╡ e34ae847-d82e-49f4-aa22-6753596c4ea0
begin
	plot()
	use_CFC11 && plot!(BD[Tracer=At(:CFC11NH)])
	use_CFC12 && plot!(BD[Tracer=At(:CFC12NH)])
	use_SF6 && plot!(BD[Tracer=At(:SF6NH)])
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

# ╔═╡ cf5bb364-5336-4dd1-8bb6-6e3f944673bf
Γ = mean_age(μ, V, B)

# ╔═╡ 4021feb1-36ac-42f6-a5f6-391c0f064dc7
# very similar values; is this correct?
Δ = ttd_width(μ, V, B)

# ╔═╡ 93c9614e-70a1-49ef-933b-b86fec342597
md"""### Green's functions """

# ╔═╡ cd492316-d6b2-4645-80ba-c5817ec5877c
Δτ = 0.25yr # time resolution

# ╔═╡ 4c258084-da30-4393-b844-c379c9e79efd
τ = 0yr:Δτ:2000yr # list of time lags

# ╔═╡ 589ab455-2e9c-47d6-abd7-f89f367a5ed5
G(t) = greens_function(t,A) # a closure that captures A

# ╔═╡ c122abb6-185c-4894-a2c4-8ab6224e83d2
G′(t) = forward_boundary_propagator(t,A,B) # type G + \prime + TAB

# ╔═╡ 595fba3f-65ec-461f-a257-92456d4f94a0
# global (or total) TTD
𝒢(t) = global_ttd(t,A,B) # type \scr + G + TAB

# ╔═╡ 96240170-eacb-4d5a-9316-eb6615a78f0a
md"""## Select interior box for diagnostics """

# ╔═╡ 07eccdb4-894d-4cd5-a639-0c01a70a84ec
@bind mbox Select(meridional_locs)

# ╔═╡ 6e1fe604-4c47-4967-bcc0-fa80fbe5bfa5
@bind vbox Select(vertical_locs)

# ╔═╡ 00902450-ceb7-4c33-be7e-906502990813
# a list comprehension
ttd1 = [G′(τ[i])[Meridional=At("1 High latitudes"),Vertical=At("1 Thermocline")][Meridional=At(mbox),Vertical=At(vbox)] for i in eachindex(τ)]

# ╔═╡ c2a38bc2-ef10-4fe5-8642-857f9acdadd7
# could be written as a for loop instead
ttd2 = [G′(τ[i])[Meridional=At("2 Mid-latitudes"),Vertical=At("1 Thermocline")][Meridional=At(mbox),Vertical=At(vbox)] for i in eachindex(τ)] 

# ╔═╡ 8d69c375-6a0c-400e-af85-3013a364fa1d
ttd_global = [𝒢(τ[i])[Meridional=At(mbox),Vertical=At(vbox)] for i in eachindex(τ)] 

# ╔═╡ 09a85965-d1dc-47a3-9eba-dd1dc129db36
Γ_ = Γ[Meridional=At(mbox),Vertical=At(vbox)] 

# ╔═╡ 19ef1da1-9b1a-4300-83aa-bb503027122b
Δ_ = Δ[Meridional=At(mbox),Vertical=At(vbox)]

# ╔═╡ fd907198-8e2e-4296-b640-c0aebbd0a796
G_inversegaussian = TracerInverseGaussian(Γ_, Δ_)

# ╔═╡ 1bb59934-17be-40d3-b227-b73bb1b9c4df
ttd_inversegaussian = pdf.(G_inversegaussian,τ)


# ╔═╡ a183e31d-8bab-46e0-a6b1-0a181c5f0f69
a1 = a[Meridional=At("1 High latitudes"),Vertical=At("1 Thermocline")][Meridional=At(mbox),Vertical=At(vbox)]

# ╔═╡ 9537166f-054f-441e-a001-3ba59a4b59e0
a2 = a[Meridional=At("2 Mid-latitudes"),Vertical=At("1 Thermocline")][Meridional=At(mbox),Vertical=At(vbox)]

# ╔═╡ e8fabe44-3a7d-47fc-84af-02baebf5f45a
begin 

	#boxloc = (Meridional=At(meridional_box),Vertical=At(vertical_box))
	# to do: put plotting into functions
	p = plot(τ,
		normalized_exponential_decay.(τ,Tmax),
		linestyle = :dash,
		yscale = :log10,
		ylabel = "Density",
		xlabel = "τ",
		label = "Tmax",
		legend = :topright,
		titlefontsize = 8,
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
	
	plot!(τ,ttd1,label="TTD 1",width=4*a1)
	plot!(τ,ttd2,label="TTD 2",width=4*a2)
	plot!(τ,ttd_global,label="Total TTD",width=4*a2,color=:black)
	plot!(τ,ttd_inversegaussian,label="Fitted inverse Gaussian")
end

# ╔═╡ Cell order:
# ╟─10b07d8a-aee4-4b64-b9eb-f22f408877ba
# ╟─27b7af71-e396-45b3-8723-8b2fc804a77f
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
# ╠═1e91d3e1-c26f-4118-9630-d654d352da76
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
# ╠═a45c8594-9fc7-46c2-833d-c44ece6648e5
# ╟─6f979bb9-733d-4981-9a53-d75162cbd372
# ╟─f53b4b2f-cda2-45a2-96f8-2dd348bc3c1f
# ╟─ef360b00-8f37-45b4-9d95-4992922ee03d
# ╟─fc1d1240-e96a-4d5e-a5f4-773d84074e22
# ╠═e34ae847-d82e-49f4-aa22-6753596c4ea0
# ╟─11eb59cf-de62-4fb4-9963-defe594e6b92
# ╠═3628ccd7-38d8-45bc-a0b6-4d74c1cb7bd9
# ╠═2175673e-5232-4804-84cb-0d5b11f31413
# ╠═01484ca5-ed33-4b94-b188-780e9e3ef8c7
# ╠═c33d09fb-fbf8-43c9-8d4b-345d90e7b40f
# ╠═cf5bb364-5336-4dd1-8bb6-6e3f944673bf
# ╠═4021feb1-36ac-42f6-a5f6-391c0f064dc7
# ╟─93c9614e-70a1-49ef-933b-b86fec342597
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
# ╟─96240170-eacb-4d5a-9316-eb6615a78f0a
# ╟─07eccdb4-894d-4cd5-a639-0c01a70a84ec
# ╟─6e1fe604-4c47-4967-bcc0-fa80fbe5bfa5
# ╟─e8fabe44-3a7d-47fc-84af-02baebf5f45a
