var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = OceanGreensFunctionMethods","category":"page"},{"location":"#OceanGreensFunctionMethods","page":"Home","title":"OceanGreensFunctionMethods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for OceanGreensFunctionMethods.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [OceanGreensFunctionMethods]","category":"page"},{"location":"#OceanGreensFunctionMethods.TracerInverseGaussian","page":"Home","title":"OceanGreensFunctionMethods.TracerInverseGaussian","text":"TracerInverseGaussian(Γ,Δ)\n\nusing LinearAlgebra: NumberArray\n\nThe tracer inverse Gaussian distribution with mean Γ and width Δ has probability density function\n\nG(𝐱 tau) = sqrtfracGamma^3 4 pi Delta^2 tau^3  exp left( - fracGamma (tau - Gamma)^24 Delta ^2 tauright) \n\nTracerInverseGaussian()              # Tracer Inverse Gaussian distribution with unit mean and unit width, i.e. TracerInverseGaussian(1, 1)\nTracerInverseGaussian(Γ, Δ)          # Tracer Inverse Gaussian distribution with mean Γ and width Δ\n\nparams(d)           # Get the parameters, i.e. (Γ, Δ)\nmean(d)             # Get the mean parameter, i.e. Γ\nshape(d)            # Get the shape parameter, i.e. Δ\n\nExternal links\n\nCompare to Inverse Gaussian distribution on Wikipedia\n\n\n\n\n\n","category":"type"},{"location":"#OceanGreensFunctionMethods.abyssal_overturning-Tuple{Any, Any}","page":"Home","title":"OceanGreensFunctionMethods.abyssal_overturning","text":"abyssal_overturning(Ψ,model_dims)\n\nSet volume flux, Ψ, in an abyssal overturning loop that satisfies the conservation of volume. Return a structure of Fluxes.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.advective_diffusive_flux-Tuple{DimensionalData.DimArray, DimensionalData.DimArray}","page":"Home","title":"OceanGreensFunctionMethods.advective_diffusive_flux","text":"advective_diffusive_flux(C, Fv; ρ)\n\nAdvective-diffusive flux of tracer C given volume fluxes Fv and optional density ρ.\n\nArguments\n\nC::DimArray: tracer distribution\nFv::DimArray: volume fluxes\nρ::Number=1035kg/m^3: uniform density\n\nReturns\n\nFc::DimArray: tracer flux\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.advective_diffusive_flux-Tuple{DimensionalData.DimArray, Fluxes}","page":"Home","title":"OceanGreensFunctionMethods.advective_diffusive_flux","text":"advective_diffusive_flux(C, Fv; ρ)\n\nAdvective-diffusive flux of tracer C given volume fluxes Fv and optional density ρ.\n\nArguments\n\nC::DimArray: tracer distribution\nFv::Fluxes: volume fluxes\nρ::Number=1035kg/m^3: uniform density\n\nReturns\n\nFc::Fluxes: tracer flux\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.boundary_dimensions-Tuple{}","page":"Home","title":"OceanGreensFunctionMethods.boundary_dimensions","text":"boundary_dimensions()\n\nDefine labels for the boundary's physical dimensions, as well as labels for the box names, consistently with the model dimensions. Use the format of DimensionalData.jl. Permits numerical quantities to be bundled with their meta-data. Dimensions are Unordered to avoid issues related to the alphabetical order.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.boundary_flux-Tuple{DimensionalData.DimArray, DimensionalData.DimArray, DimensionalData.DimArray}","page":"Home","title":"OceanGreensFunctionMethods.boundary_flux","text":"boundary_flux(f::DimArray, C::DimArray, Fb::DimArray)\n\nConvergence or net effect of boundary fluxes.\n\nArguments\n\nf::DimArray: Dirichlet boundary condition\nC::DimArray: tracer distribution \nFb::DimArray: boundary exchange volume flux\n\nReturns\n\nJb::Fluxes: boundary tracer flux\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.convergence-Tuple{Fluxes}","page":"Home","title":"OceanGreensFunctionMethods.convergence","text":"convergence(J)\n\nConvergence of fluxes J of type Fluxes.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.evolve_concentration-NTuple{5, Any}","page":"Home","title":"OceanGreensFunctionMethods.evolve_concentration","text":"evolve_concentration(C₀, A, B, tlist, source_history; halflife = nothing)\n\nIntegrate forcing vector over time to compute the concentration history. Find propagator by analytical expression using eigen-methods.\n\nArguments\n\nC₀: initial tracer concentration\nA: tracer transport information used in matrix calculations\nB: boundary condition information used in matrix calculations\ntlist: list of times to save tracer concentration\nsource_history::Function: returns Dirichlet boundary condition at a given time\nhalflife=nothing: radioactive half life (optional)\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.forcing_integrand-NTuple{6, Any}","page":"Home","title":"OceanGreensFunctionMethods.forcing_integrand","text":"forcing_integrand(t, tf, μ, V, B, source_history)\n\nIntegrand for boundary condition term in equation 10 (Haine et al., 2024).\n\nArguments\n\nt: time\ntf: final time \nμ: eigenvalue diagonal matrix\nV: eigenvector matrix\nB: boundary condition matrix\nsource_history::Function: returns Dirichlet boundary condition at a given time\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.greens_function-Union{Tuple{DM}, Tuple{Q}, Tuple{Any, DimensionalData.DimMatrix{DM, D, R, A} where {D<:Tuple, R<:Tuple, A<:AbstractMatrix{DM}}}} where {Q<:Unitful.Quantity, DM<:(DimensionalData.DimMatrix{Q, D, R, A} where {D<:Tuple, R<:Tuple, A<:AbstractMatrix{Q}})}","page":"Home","title":"OceanGreensFunctionMethods.greens_function","text":"greens_function(τ,A)\n\nGreen's function for a box model (for steady transport given by the matrix 𝐀 for response at time t to a source at time t′ where τ = t - t′): the matrix exponential function of the elapsed time between the source time and field time.\n\n𝐆(t) = eᴬᵗ\n\nwhere 𝐆(t) is a  N × N matrix with the spatial locations of field points (boxes) down its N rows and source points (boxes) along its N columns. Thus, the element 𝐆{i,j}(τ) quantifies transfer from a source at time t′ in box j to receiver at time t in box i.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.integrate_forcing-NTuple{6, Any}","page":"Home","title":"OceanGreensFunctionMethods.integrate_forcing","text":"integrate_forcing(t0, tf, μ, V, B, source_history)\n\nIntegrate boundary condition term in equation 10 (Haine et al., 2024).\n\nArguments\n\nt0: initial time\ntf: final time \nμ: eigenvalue diagonal matrix\nV: eigenvector matrix\nB: boundary condition matrix\nsource_history::Function: returns Dirichlet boundary condition at a given time\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.intermediate_overturning-Tuple{Any, Any}","page":"Home","title":"OceanGreensFunctionMethods.intermediate_overturning","text":"intermediate_overturning(Ψ,model_dims)\n\nSet the volume flux, Ψ, in an intermediate overturning loop that satisfies the conservation of volume. Return a structure of Fluxes.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.linear_probe-Union{Tuple{T}, Tuple{N}, Tuple{Function, DimensionalData.DimArray{T, N, D, R, A} where {D<:Tuple, R<:Tuple, A<:AbstractArray{T, N}}, Vararg{Any}}} where {N, T<:Number}","page":"Home","title":"OceanGreensFunctionMethods.linear_probe","text":"linear_probe(funk, x, args...)\n\nProbe a function to determine its linear response in matrix form. Assumes units are needed and available. A simpler function to handle cases without units would be nice.\n\nArguments\n\nfunk: function to be probed\nx: input (independent) variable\nhalflife::Number: radioactive half life\nargs: the arguments that follow x in funk\n\nReturns\n\nA::DimArray{DimArray}: labeled transport information used in matrix operations \n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.location_iodine129_history-Tuple{}","page":"Home","title":"OceanGreensFunctionMethods.location_iodine129_history","text":"location_transient_tracer_histories()\n\nURL of iodine-129 source history file.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.location_transient_tracer_histories-Tuple{}","page":"Home","title":"OceanGreensFunctionMethods.location_transient_tracer_histories","text":"location_transient_tracer_histories()\n\nURL of tracer source history file.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.mass-Tuple{Any}","page":"Home","title":"OceanGreensFunctionMethods.mass","text":"mass(V; ρ)\n\nSeawater mass derived from the volume V and an optional input of density ρ.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.mass_convergence-Tuple{Any}","page":"Home","title":"OceanGreensFunctionMethods.mass_convergence","text":"mass_convergence(Fv)\n\nConvergence of volume derived from a field of volume fluxes Fv, translated into a mass flux convergence with the assumption of uniform density.  \n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.model_dimensions-Tuple{}","page":"Home","title":"OceanGreensFunctionMethods.model_dimensions","text":"model_dimensions()\n\nDefine labels for the model's physical dimensions, as well as labels for the box names. Use the format of DimensionalData.jl. Permits numerical quantities to be bundled with their meta-data. Dimensions are Unordered to avoid issues related to the alphabetical order.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.radioactive_decay-Tuple{DimensionalData.DimArray, Number}","page":"Home","title":"OceanGreensFunctionMethods.radioactive_decay","text":"radioactive_decay(C, halflife)\n\nRadioactive decay rate of tracer C with half life of halflife.\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.read_iodine129_history-Tuple{}","page":"Home","title":"OceanGreensFunctionMethods.read_iodine129_history","text":"read_iodine129_history()\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.read_transient_tracer_histories-Tuple{}","page":"Home","title":"OceanGreensFunctionMethods.read_transient_tracer_histories","text":"read_transient_tracer_histories()\n\nRead transient tracer source histories and save as a DimArray. \n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.steady_tracer_timeseries-NTuple{7, Any}","page":"Home","title":"OceanGreensFunctionMethods.steady_tracer_timeseries","text":"steady_tracer_timeseries(tracername, A, B, halflife, tlist, mbox1, vbox1)\n\nSimulate non-transient tracers and return tracer timeseries from one box.\n\nArguments\n\ntracername: name of tracer\nA: tracer transport matrix\nB: boundary condition matrix\nhalflife: radioactive half life\nBD: Dirichlet boundary condition\ntlist: list of times to save tracer concentration\nmbox: name of meridional box of interest\nvbox: name of vertical box of interest\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.timestep_initial_condition-NTuple{5, Any}","page":"Home","title":"OceanGreensFunctionMethods.timestep_initial_condition","text":"timestepinitialcondition(C, μ, V, ti, tf)\n\nArguments\n\nC::DimArray: tracer distribution at ti\nμ: eigenvalue diagonal matrix\nV: eigenvector matrix\nti: initial time\ntf: final time\n\nReturns\n\nCf::DimArray: tracer distribution at tf\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.tracer_point_source_history-Tuple{Any, Any}","page":"Home","title":"OceanGreensFunctionMethods.tracer_point_source_history","text":"tracer_point_source_history(tracername, BD)\n\nReturn a function that yields transient tracer source history (such as CFCs) at given time.\n\nArguments\n\ntracername: name of tracer in source history file\nBD::DimArray: Dirichlet boundary condition compendium for many tracers\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.tracer_source_history","page":"Home","title":"OceanGreensFunctionMethods.tracer_source_history","text":"tracer_source_history(t, tracername, box2_box1_ratio, BD = nothing)\n\nReturn source history values for all boundary points.\n\nArguments\n\nt: time\ntracername: name of tracer in source history file\nbox2_box1_ratio: ratio of boundary condition value in Mid-latitudes to High Latitudes\nBD::DimArray=nothing: Dirichlet boundary condition compendium (optional)\n\n\n\n\n\n","category":"function"},{"location":"#OceanGreensFunctionMethods.tracer_tendency-Union{Tuple{N}, Tuple{DimensionalData.DimArray{var\"#s11\", N, D, R, A} where {var\"#s11\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s11\", N}}, Number}} where N","page":"Home","title":"OceanGreensFunctionMethods.tracer_tendency","text":"tracer_tendency(C)\n\nTracer tendency ∂C/∂t for the radioactive decay of a tracer C with half life halflife, for use with finding the radioactive contribution to a tracer transport matrix.\n\nArguments\n\nC::DimArray: tracer distribution\nhalflife::Number: radioactive half life\n\nReturns\n\ndCdt::DimArray: tracer tendency\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.tracer_tendency-Union{Tuple{N}, Tuple{DimensionalData.DimArray{var\"#s7\", N, D, R, A} where {var\"#s7\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s7\", N}}, DimensionalData.DimArray{var\"#s6\", N, D, R, A} where {var\"#s6\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s6\", N}}, Fluxes{<:Number, N}, DimensionalData.DimArray{var\"#s3\", N, D, R, A} where {var\"#s3\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s3\", N}}, DimensionalData.DimArray{var\"#s2\", N, D, R, A} where {var\"#s2\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s2\", N}}}} where N","page":"Home","title":"OceanGreensFunctionMethods.tracer_tendency","text":"tracer_tendency(C, f, Fv, Fb, V)\n\nTracer tendency ∂C/∂t for a tracer C, especially useful for finding a tracer transport matrix. \n\nArguments\n\nC::DimArray: tracer distribution\nf::DimArray: Dirichlet boundary condition\nFv::Fluxes: volume fluxes\nFb::Fluxes: volume fluxes\nV::DimArray: box volume\n\nReturns\n\ndCdt::DimArray: tracer tendency\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.tracer_tendency-Union{Tuple{N}, Tuple{DimensionalData.DimArray{var\"#s8\", N, D, R, A} where {var\"#s8\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s8\", N}}, DimensionalData.DimArray{var\"#s7\", N, D, R, A} where {var\"#s7\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s7\", N}}, DimensionalData.DimArray{var\"#s6\", N, D, R, A} where {var\"#s6\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s6\", N}}, DimensionalData.DimArray{var\"#s5\", N, D, R, A} where {var\"#s5\"<:Number, D<:Tuple, R<:Tuple, A<:AbstractArray{var\"#s5\", N}}}} where N","page":"Home","title":"OceanGreensFunctionMethods.tracer_tendency","text":"tracer_tendency(f, C, Fv, Fb, V)\n\nTracer tendency ∂C/∂t for a boundary flux f, for use with finding B boundary matrix.\n\nArguments\n\nf::DimArray: Dirichlet boundary condition\nC::DimArray: tracer distribution\nFv::Fluxes: volume fluxes\nFb::Fluxes: volume fluxes\nV::DimArray: box volume\n\nReturns\n\ndCdt::DimArray: tracer tendency\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.tracer_timeseries-NTuple{6, Any}","page":"Home","title":"OceanGreensFunctionMethods.tracer_timeseries","text":"tracer_timeseries(tracername, A, B, tlist, mbox1, vbox1; BD=nothing, halflife=nothing)\n\nSimulate tracers and return tracer timeseries from one box.\n\nArguments\n\ntracername: name of tracer\nA: tracer transport matrix\nB: boundary condition matrix\ntlist: list of times to save tracer concentration\nmbox: name of meridional box of interest\nvbox: name of vertical box of interest\nBD=nothing: Dirichlet boundary condition\nhalflife=nothing: radioactive half life\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.transient_tracer_timeseries-NTuple{7, Any}","page":"Home","title":"OceanGreensFunctionMethods.transient_tracer_timeseries","text":"transient_tracer_timeseries(tracername, A, B, BD, tlist, mbox1, vbox1; halflife = nothing)\n\nSimulate transient tracers and return tracer timeseries from one box.\n\nArguments\n\ntracername: name of tracer\nA: tracer transport matrix\nB: boundary condition matrix\nBD: Dirichlet boundary condition\ntlist: list of times to save tracer concentration\nmbox: name of meridional box of interest\nvbox: name of vertical box of interest\nhalflife=nothing: radioactive half life\n\n\n\n\n\n","category":"method"},{"location":"#OceanGreensFunctionMethods.vertical_diffusion-Tuple{Any, Any}","page":"Home","title":"OceanGreensFunctionMethods.vertical_diffusion","text":"vertical_diffusion(Fv_exchange,model_dims)\n\nSet vertical diffusive-like exchange flux Fv_exchange. Return a structure of Fluxes.\n\n\n\n\n\n","category":"method"}]
}
