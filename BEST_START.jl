#= BEST method - Beerkan Estimation of Soil Transfer params####
Jesus Fernandez and Joseph =#

#= =============== INITIALIZATION ============== =#
cd("C:\\JOE\\Main\\MODELS\\SOIL\\BEST")
DIR_Working = pwd() # Saving the current working directory
println( "Working path : $DIR_Working")
push!(LOAD_PATH, DIR_Working) # All the subroutines are in the current working directory
 
 # Do not change order
include("Path.jl")
include("Package.jl")
include("Cst.jl")
include("Param.jl")
include("Tools.jl")
include("Option.jl")
include("Wrc.jl")
include("Kunsat.jl")
include("Diffusivity.jl")
include("Stat.jl")
include("HydroRelationship.jl")
include("Best.jl")
include("ThetaTime.jl")
include("QuasiExact.jl")
include("ConvertHydroModel.jl")
include("Reading.jl")
include("Sorptivity.jl")
include("KunsatModel.jl")
include("TestBest.jl")
include("Err.jl")
include("Plots.jl")
# include("Plots.jl")

using PGFPlots
using BlackBoxOptim, DataFrames, CSV, SpecialFunctions, Optim, Suppressor, plots
export reading, option, cst, wrc, best, sorptivity, diffusivity, kunsat, param, array, quasiExact, PACKAGES, testbest, stats, relationship

function START_BEST()
	@suppress_err begin
	 println("START running:      $(path.FileName) ...\n")
	 #Import packages if not yet imported
	 if option.DownloadPackage == true
		 PACKAGES(Option_PackageUpdate = false)
	 end

	 # ====================BEFORE THE LOOPS==================================== #
	 
	 # Preparing the plots
	 	Plot_CharacUnsat = GroupPlot(2, 8, groupStyle = "horizontal sep = 2cm, vertical sep = 1.5cm")
		Path_CharacUnsat = string(DIR_Working, "//OUTPUT//Figure//CharacUnsat//CharacUnsat_", path.FileName, ".svg") #~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println(Path_CharacUnsat)
		if isfile(Path_CharacUnsat) # If the file is there than delete it before saving figure
			rm(Path_CharacUnsat)
		end

		Plot_InfHydro = GroupPlot(2, 8, groupStyle =  "horizontal sep = 2cm, vertical sep = 1.5cm")
		Path_InfHydro = string(DIR_Working, "//OUTPUT//Figure//InfiltrationHydro//InfHydro_", path.FileName, ".svg") #~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if isfile(Path_InfHydro) # If the file is there than delete it before saving figure
			rm(Path_InfHydro)
		end

		Plot_Infiltration = GroupPlot(2, 8, groupStyle =  "horizontal sep = 2cm, vertical sep = 1.5cm")
		Path_Infiltration = string(DIR_Working, "//OUTPUT//Figure//Infiltration//Inf_", path.FileName, ".svg") #~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println(Path_Infiltration)
		if isfile(Path_Infiltration) # If the file is there than delete it before saving figure
			rm(Path_Infiltration)
		end

		#= =============== printing header of table =============== =# 
		# We will remove table at the very beginning since we are appending
		Path_Table =  string(DIR_Working, path.Table, path.FileName, ".csv")
		if isfile(Path_Table) # If the file is there than delete it before saving figure
			rm(Path_Table)
		end
		Data = [  ]
		CSV.write(Path_Table, DataFrame(Data); header=true, colnames=param.Header, append=false)	 
	
	
		#= =============== reading infiltration =============== =#
	   println("START to read infiltration data...")
		  Path_Inf = string(DIR_Working, path.Read_Inf)
		
		  Time, Inf_3D_Obs, Time_Ni, Sample_N = reading.INFILTRATION(Path_Inf)
		
		  Id = linspace(1,Sample_N, Sample_N) # This is just an array 1 ... Sample_N
	  println("END of reading infiltration data \n")


	 #= =============== reading hydraulic =============== =#
	 println("START to read hydraulic params...")
		 Path_Hydraulic = string(DIR_Working, path.Read_Hydraulic)

		 Se_Ini, θs, θr, σ, Hkg, Ks, Hvg, N, Km, RingRadius, SampleTrue, Name = reading.HYDRAULIC(Path_Hydraulic)
		 
		 # We measure effective porosity
		 θs = θs - θr 
		 θr = θr * 0.
	 println("END to read hydraulic params \n")


	 Flag_Pσ_1_2 = false
	 if Flag_Pσ_1_2
	 	Pσ_1, Pσ_2 = relationship.OPTIMIZATION_σ_2_Hkg(SampleTrue[1:Sample_N], Hkg[1:Sample_N], σ[1:Sample_N])
	 end
	 
	 #= =============== DECLARING ARRAYS ============== =#
	Time_Ni_Max = 2000

	DθDh_Accur_Kg = Array{Float64}(Sample_N)
	DθDh_Accur_Vg = Array{Float64}(Sample_N)
	Err_Ks_BestG_Kg = Array{Float64}(Sample_N)
	Err_Ks_BestG_Vg = Array{Float64}(Sample_N)
	Err_Ks_BestGi_Kg = Array{Float64}(Sample_N)
	Err_Ks_BestGi_Vg = Array{Float64}(Sample_N)
	Err_Ks_Qe_Kg = Array{Float64}(Sample_N)
	Err_Ks_Qe_Vg = Array{Float64}(Sample_N)
	Err_Qe = Array{Float64}(Sample_N)
	Err_Sorpt_BestG_Kg = Array{Float64}(Sample_N)
	Err_Sorpt_BestG_Vg = Array{Float64}(Sample_N)
	Err_Sorpt_BestGi_Kg = Array{Float64}(Sample_N)
	Err_Sorpt_BestGi_Vg = Array{Float64}(Sample_N)
	Err_Sorpt_Qe_Kg = Array{Float64}(Sample_N)
	Err_Sorpt_Qe_Vg = Array{Float64}(Sample_N)
	Hkg_BestG = Array{Float64}(Sample_N)
	Hkg_BestGi = Array{Float64}(Sample_N)
	Hkg_BestImp = Array{Float64}(Sample_N)
	Hkg_Qe = Array{Float64}(Sample_N)
	Hkg_Qe_Imp = Array{Float64}(Sample_N)
	Hkg_Qe_Select = Array{Float64}(Sample_N)
	Hvg_BestG= Array{Float64}(Sample_N)
	Hvg_BestGi= Array{Float64}(Sample_N)
	Hvg_BestImp= Array{Float64}(Sample_N)
	Hvg_Inf = Array{Float64}(Sample_N, Time_Ni_Max)
	Hvg_Qe = Array{Float64}(Sample_N)
	Hvg_Qe_Imp = Array{Float64}(Sample_N)
	Hvg_Qe_Select = Array{Float64}(Sample_N)
	Inf_BestG_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_BestG_Vg=Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_BestGi_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_BestGi_Vg=Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_BestImp_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_BestImp_Vg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Hydro_BestG_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Hydro_BestG_Vg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Hydro_BestGi_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Hydro_BestGi_Vg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Hydro_Qe_Kg= Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Hydro_Qe_Vg = Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Qe_Imp_Kg =  Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Qe_Imp_Vg =  Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Qe_Kg =  Array{Float64}(Sample_N, Time_Ni_Max)
	Inf_Qe_Vg =  Array{Float64}(Sample_N, Time_Ni_Max)
	iT_TransStead_BestG_Kg= Array{Int}(Sample_N)
	iT_TransStead_BestGi_Kg= Array{Int}(Sample_N)
	iT_TransStead_Hydro_Kg= Array{Int}(Sample_N)
	iT_TransStead_Sim = Array{Int}(Sample_N)
	iT_TransStead_Small= Array{Int}(Sample_N)
	iT_TransStead_Vg = Array{Int}(Sample_N)
	Kr_θini_Best_Kg = Array{Float64}(Sample_N)
	Kr_θini_BestG_Kg = Array{Float64}(Sample_N)
	Kr_θini_BestG_Vg= Array{Float64}(Sample_N)
	Kr_θini_BestGi_Kg = Array{Float64}(Sample_N)
	Kr_θini_BestGi_Vg= Array{Float64}(Sample_N)
	Kr_θini_Inf_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	Kr_θini_Inf_Vg = Array{Float64}(Sample_N, Time_Ni_Max)
	Kr_θini_Kg = Array{Float64}(Sample_N)
	Kr_θini_Qe_Imp_Kg= Array{Float64}(Sample_N)
	Kr_θini_Qe_Imp_Vg = Array{Float64}(Sample_N)
	Kr_θini_Qe_Kg= Array{Float64}(Sample_N)
	Kr_θini_Qe_Vg = Array{Float64}(Sample_N)
	Kr_θini_Vg = Array{Float64}(Sample_N)
	Ks_BestG_Kg = Array{Float64}(Sample_N)
	Ks_BestG_Vg= Array{Float64}(Sample_N)
	Ks_BestGi_Kg = Array{Float64}(Sample_N)
	Ks_BestGi_Vg= Array{Float64}(Sample_N)
	Ks_BestImp_Kg = Array{Float64}(Sample_N)
	Ks_BestImp_Vg = Array{Float64}(Sample_N)
	Ks_Qe_Imp_Kg = Array{Float64}(Sample_N)
	Ks_Qe_Imp_Vg = Array{Float64}(Sample_N)
	Ks_Qe_Kg = Array{Float64}(Sample_N)
	Ks_Qe_Vg = Array{Float64}(Sample_N)
	N_Best = Array{Float64}(Sample_N)
	N_BestG= Array{Float64}(Sample_N)
	N_BestGi= Array{Float64}(Sample_N)
	N_BestImp= Array{Float64}(Sample_N) 
	N_Inf = Array{Float64}(Sample_N, Time_Ni_Max)
	N_Qe = Array{Float64}(Sample_N)
	N_Qe_Imp = Array{Float64}(Sample_N)
	N_Qe_Select = Array{Float64}(Sample_N)
	NSE_BestG_Kg = Array{Float64}(Sample_N)
	NSE_BestG_Vg = Array{Float64}(Sample_N)
	NSE_BestGi_Kg = Array{Float64}(Sample_N)
	NSE_BestGi_Vg = Array{Float64}(Sample_N)
	NSE_Hydro_Inf_BestG_Kg = Array{Float64}(Sample_N)
	NSE_Hydro_Inf_BestG_Vg = Array{Float64}(Sample_N)
	NSE_Hydro_Inf_BestGi_Kg = Array{Float64}(Sample_N)
	NSE_Hydro_Inf_BestGi_Vg = Array{Float64}(Sample_N)
	NSE_Hydro_Inf_Qe_Kg = Array{Float64}(Sample_N)
	NSE_Hydro_Inf_Qe_Vg = Array{Float64}(Sample_N)
	NSE_Inf_BestG_Kg = Array{Float64}(Sample_N)
	NSE_Inf_BestG_Vg = Array{Float64}(Sample_N)
	NSE_Inf_BestGi_Kg = Array{Float64}(Sample_N)
	NSE_Inf_BestGi_Vg = Array{Float64}(Sample_N)
	NSE_Inf_Qe_Kg = Array{Float64}(Sample_N)
	NSE_Inf_Qe_Vg = Array{Float64}(Sample_N)
	NSE_Qe_Kg = Array{Float64}(Sample_N)
	NSE_Qe_Vg = Array{Float64}(Sample_N)
	Sorpt_BestG_Kg = Array{Float64}(Sample_N)
	Sorpt_BestG_Vg= Array{Float64}(Sample_N)
	Sorpt_BestGi_Kg = Array{Float64}(Sample_N)
	Sorpt_BestGi_Vg= Array{Float64}(Sample_N)
	Sorpt_BestImp_Kg = Array{Float64}(Sample_N)
	Sorpt_BestImp_Vg = Array{Float64}(Sample_N)
	Sorpt_Hydro_Kg = Array{Float64}(Sample_N)
	Sorpt_Hydro_Vg = Array{Float64}(Sample_N)
	Sorpt_Inf_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	Sorpt_Inf_Vg = Array{Float64}(Sample_N, Time_Ni_Max)
	Sorpt_Qe_Imp_Kg= Array{Float64}(Sample_N)
	Sorpt_Qe_Imp_Vg = Array{Float64}(Sample_N)
	Sorpt_Qe_Kg= Array{Float64}(Sample_N)
	Sorpt_Qe_Vg = Array{Float64}(Sample_N)
	Time_TransStead_Best  = Array{Float64}(Sample_N)
	Time_TransStead_Best_Vg = Array{Float64}(Sample_N) 
	Time_TransStead_BestG_Kg = Array{Float64}(Sample_N)
	Time_TransStead_BestGi_Kg = Array{Float64}(Sample_N)
	Time_TransStead_BestGi_Vg = Array{Float64}(Sample_N)
	Time_TransStead_BestG_Vg= Array{Float64}(Sample_N)
	Time_TransStead_Hydro_Kg = Array{Float64}(Sample_N)
	Time_TransStead_Hydro_Vg = Array{Float64}(Sample_N)
	Time_TransStead_Sim = Array{Float64}(Sample_N)	
	Time_TransStead_Small= Array{Float64}(Sample_N)
	θ_Inf_Kg = Array{Float64}(Sample_N, Time_Ni_Max)
	θ_Inf_Vg = Array{Float64}(Sample_N, Time_Ni_Max)
	θ_Ini = Array{Float64}(Sample_N)
	σ_Best = Array{Float64}(Sample_N)
	σ_BestG = Array{Float64}(Sample_N)
	σ_BestGi = Array{Float64}(Sample_N)
	σ_BestImp = Array{Float64}(Sample_N)  
	σ_Inf = Array{Float64}(Sample_N, Time_Ni_Max)
	σ_Qe = Array{Float64}(Sample_N)
	σ_Qe_Imp = Array{Float64}(Sample_N)
	σ_Qe_Select = Array{Float64}(Sample_N)

	for iS in param.i_Sample_Start:min(param.i_Sample_End,Sample_N)
	if SampleTrue[iS] == 1
		println("iS: $iS -------------------------------------\n\n")

		iT_N_Max = Time_Ni[iS] # Maximum End of iT for each iS

		# Deriving Time_TransStead
		iT_TransStead_Sim[iS], Time_TransStead_Sim[iS] = stats.INFILTRATION_2_iSTEADYTRANSIT(param.ΔInf_SteadyTransit/1., Time[iS, 1:iT_N_Max], Inf_3D_Obs[iS,1:iT_N_Max])

		iT_TransStead_Small[iS], Time_TransStead_Small[iS] = stats.INFILTRATION_2_iSTEADYTRANSIT(param.ΔInf_SteadyTransit/10., Time[iS, 1:iT_N_Max], Inf_3D_Obs[iS,1:iT_N_Max])

		# The time of simulation is too long so we reduduced to iT_TransStead_Obs x Param
		iT_N = Int(min(array.SEARCH_INDEX(Time[iS, 1:iT_N_Max], Time_TransStead_Sim[iS]*param.TransStead_Multiply), iT_N_Max))
		println("iT_N $iT_N")


	#= =============== GETTING H_Ini, θ_Ini, Kr =============== =# 
		θ_Ini[iS] = wrc.se.Se_2_θ(Se_Ini[iS], θs[iS], θr[iS])

		Hini = wrc.kg.Se_2_H(Se_Ini[iS], Hkg[iS], σ[iS])

		Kr_θini_Vg[iS] = kunsat.vg.KUNSAT(Se_Ini[iS], N[iS], 1., Km[iS])

		Kr_θini_Kg[iS] = kunsat.kg.KUNSAT(Se_Ini[iS], θs[iS], θr[iS], σ[iS], 1., θs[iS], σ[iS])
			

	 #= =============== DERIVING EXPECTED ERROR OF SIMULATIONS =============== =#	
		Err_Qe[iS] = err.kg.ERR_QE(θs[iS], θr[iS], Hkg[iS], σ[iS], Ks[iS])
		println("Err_Qe $(Err_Qe[iS]) \n")

		if option.DataAvailable.HydraulicParam
			println("MEASURED HYDRAULIC.....")


	#= =============== DERIVE SORPTIVITY =============== =#	
		Flag_PlotSorptivity = true

		if Flag_PlotSorptivity
			plots.kg.DIFFUSIVITY_PLOT( iS, θ_Ini[iS], θs[iS], θr[iS], Hkg[iS], σ[iS], Ks[iS])
		end


	 #= =============== DERIVE FROM HYDRAULIC PARAMETERS =============== =#
			Time_TransStead_Hydro_Kg[iS] = best.kg.STEADY_TRANSIT_BESTG(θ_Ini[iS], θs[iS], θr[iS], Hkg[iS], σ[iS], Ks[iS], RingRadius[iS])
			println("Time_TransStead_Hydro_Kg =  ", Time_TransStead_Hydro_Kg[iS])

			Sorpt_Hydro_Kg[iS] = sorptivity.kg.SORPTIVITY(θ_Ini[iS], θs[iS], θs[iS], θr[iS], Hkg[iS], σ[iS], Ks[iS])

			for iT in 1: iT_N  # Looping for every infiltration point
				Inf_Hydro_BestG_Kg[iS,iT] = best.BESTG(Time[iS, iT], θ_Ini[iS], θs[iS], θr[iS], Ks[iS], RingRadius[iS], Sorpt_Hydro_Kg[iS], Kr_θini_Kg[iS]; Flag_Best="Best_G")

				Inf_Hydro_BestGi_Kg[iS,iT] = best.BESTG(Time[iS, iT], θ_Ini[iS], θs[iS], θr[iS], Ks[iS], RingRadius[iS], Sorpt_Hydro_Kg[iS], Kr_θini_Kg[iS]; Flag_Best="Best_GI")
			end

			# Statistics on infiltration HYDRO
			NSE_Hydro_Inf_BestG_Kg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Hydro_BestG_Kg[iS, 2:iT_N])) # Remark iT=1=>Inf=0

			NSE_Hydro_Inf_BestGi_Kg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Hydro_BestGi_Kg[iS,2:iT_N])) # Remark iT=1=>Inf=0

			# println("Kosugi Hydro Qe")
			Inf_Hydro_Qe_Kg[iS, 1:iT_N] = quasiExact.HYDRO_2_INFILTRATION3D(Time[iS,1:iT_N], Sorpt_Hydro_Kg[iS], Ks[iS], Ks[iS]*Kr_θini_Kg[iS], θs[iS], θ_Ini[iS], iT_N, RingRadius[iS], iT_TransStead_Sim[iS])

			# Statistics on infiltration HYDRO
			NSE_Hydro_Inf_Qe_Kg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Hydro_Qe_Kg[iS, 2:iT_N]))

			println( iS, " ,Sorpt_Hydro_Kg=", Sorpt_Hydro_Kg[iS]," ,σ=", σ[iS], " ,Hkg= ", Hkg[iS], " ,Ks=", 360.*Ks[iS], " ,NSE_Hydro_Inf_Best_Kg=", ",",NSE_Hydro_Inf_BestG_Kg[iS], " ,NSE_Hydro_Inf_Qe_Kg=", ",",NSE_Hydro_Inf_Qe_Kg[iS], "\n")

			
			
			#-----------------------------------------------------------------------------------------------
			# iT_TransStead = Int(min(array.SEARCH_INDEX(Time[iS, 1:iT_N_Max], Time_TransStead_Hydro_Kg[iS] , iT_N_Max)))
			# σ_Mod = thetaTime.kg.θTIME_σ(Time[iS, 1:iT_N], Inf_Hydro_Qe_Kg[iS, 1:iT_N], RingRadius[iS],  Sorpt_Hydro_Kg[iS], θ_Ini[iS], θs[iS], θr[iS], Ks[iS], 15,  Time[iS, 15])

			if option.HydroModel.HydraulicParam == "vangenuchten"
				# println("Vangenuchten Hydro BestSorpt")
				Sorpt_Hydro_Vg[iS] = sorptivity.vg.SORPTIVITY(θ_Ini[iS], θs[iS], θs[iS], θr[iS], Hvg[iS], N[iS], Ks[iS], Km[iS])

				Time_TransStead_Hydro_Vg[iS] = best.vg.STEADY_TRANSIT_BESTG(θ_Ini[iS], θs[iS], θr[iS], Hvg[iS], N[iS], Ks[iS], Km[iS], RingRadius[iS])

				for iT in 1: iT_N
					Inf_Hydro_BestG_Vg[iS, iT] = best.BESTG(Time[iS, iT], θ_Ini[iS], θs[iS], θr[iS], Ks[iS], RingRadius[iS], Sorpt_Hydro_Vg[iS], Kr_θini_Vg[iS]; Flag_Best="Best_G")

					Inf_Hydro_BestGi_Vg[iS, iT] = best.BESTG(Time[iS, iT], θ_Ini[iS], θs[iS], θr[iS], Ks[iS], RingRadius[iS], Sorpt_Hydro_Vg[iS], Kr_θini_Vg[iS]; Flag_Best="Best_Gi")
				end

				# Statistics on infiltration HYDRO
				NSE_Hydro_Inf_BestG_Vg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Hydro_BestG_Vg[iS,2:iT_N]))

				NSE_Hydro_Inf_BestGi_Vg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Hydro_BestGi_Vg[iS,2:iT_N]))

				# println("vanGenuchten Hydro Qe")
				Inf_Hydro_Qe_Vg[iS, 1:iT_N] = quasiExact.HYDRO_2_INFILTRATION3D(Time[iS, 1:iT_N], Sorpt_Hydro_Vg[iS], Ks[iS], Ks[iS]*Kr_θini_Vg[iS], θs[iS], θ_Ini[iS], iT_N, RingRadius[iS], iT_TransStead_Sim[iS])

				NSE_Hydro_Inf_Qe_Vg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Hydro_Qe_Vg[iS,2:iT_N]))

				println( iS, " ,Sorpt_Hydro_Vg=", Sorpt_Hydro_Vg[iS], "  N=", N[iS]," ,Hvg= ", Hvg[iS], " ,Ks=", 360.*Ks[iS], " ,NSE_Hydro_Inf_Best_Vg=", ",", NSE_Hydro_Inf_BestG_Vg[iS]," ,NSE_Hydro_Qe_Vg=", NSE_Hydro_Inf_Qe_Vg[iS], "\n")
			end
		end


	 #= =============== BEST SORPTIVITY =============== =#
			println("BEST_Gi_Kg.......................................................... \n")

			Sorpt_BestGi_Kg[iS], Hkg_BestGi[iS], σ_BestGi[iS], Ks_BestGi_Kg[iS], Kr_θini_BestGi_Kg[iS], Inf_BestGi_Kg[iS,1:iT_N] = best.kg.BESTG_INVERSE(Time[iS,1:iT_N], iT_N, Inf_3D_Obs[iS,1:iT_N], θ_Ini[iS], Se_Ini[iS], θs[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS],Option_Opt_σ=true, Flag_Best="Best_Gi")

			σ_BestGi[iS], Hkg_BestGi[iS] = best.kg.BESTG_INVERSE_σmodel(Sorpt_BestGi_Kg[iS], Ks_BestGi_Kg[iS], Time[iS,1:iT_N], iT_N, Inf_3D_Obs[iS,1:iT_N], θ_Ini[iS], Se_Ini[iS], θs[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS]; Flag_Best="Best_Gi")

			# σ_BestGi[iS], Hkg_BestGi[iS] = relationship.SORPTIVITY_2_σ_Hkg(Sorpt_BestGi_Kg[iS], θ_Ini[iS], θs[iS], θr[iS], Ks_BestGi_Kg[iS])

			Time_TransStead_BestGi_Kg[iS] = best.kg.STEADY_TRANSIT_BESTG(θ_Ini[iS], θs[iS], θr[iS], Hkg_BestGi[iS], σ_BestGi[iS], Ks_BestGi_Kg[iS], RingRadius[iS])

			iT_TransStead_BestGi_Kg[iS] = array.SEARCH_INDEX(Time[iS,1:iT_N], Time_TransStead_BestGi_Kg[iS])

			θ_Time_BestGi_Kg = Array{Float64}(iT_N)
			θ_Time_BestGi_Kg[1:iT_N], ~ = thetaTime.kg.θTIME(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius[iS], Sorpt_BestGi_Kg[iS], θ_Ini[iS], θs[iS], θr[iS],  σ_BestGi[iS], Ks_BestGi_Kg[iS], iT_TransStead_BestGi_Kg[iS], iT_N)

			NSE_Inf_BestGi_Kg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS, 2:iT_N]), log10.(Inf_BestGi_Kg[iS, 2:iT_N]))

			Err_Sorpt_BestGi_Kg[iS] = stats.RELATIVEerr(Sorpt_Hydro_Kg[iS], Sorpt_BestGi_Kg[iS])

			Err_Ks_BestGi_Kg[iS] = stats.RELATIVEerr((Ks[iS]), (Ks_BestGi_Kg[iS]))

			println( iS, " ,Sorpt_BestGi_Kg=", Sorpt_BestGi_Kg[iS], " ,σ_BestGi=", σ_BestGi[iS], " ,Hkg_BestG_Kg= ", Hkg_BestGi[iS], " ,Ks_BestGi_Kg=", 360.*Ks_BestGi_Kg[iS], " ,NSE_Inf_BestGi_Kg=", NSE_Inf_BestGi_Kg[iS], "\n ")

			println("BEST_G_Kg..........................................................\n")
			Sorpt_BestG_Kg[iS], Hkg_BestG[iS], σ_BestG[iS], Ks_BestG_Kg[iS], Kr_θini_BestG_Kg[iS], Inf_BestG_Kg[iS,1:iT_N] = best.kg.BESTG_INVERSE(Time[iS,1:iT_N], iT_N, Inf_3D_Obs[iS,1:iT_N], θ_Ini[iS], Se_Ini[iS], θs[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS],Option_Opt_σ=true, Flag_Best="Best_G")

			σ_BestG[iS], Hkg_BestG[iS] = best.kg.BESTG_INVERSE_σmodel(Sorpt_BestG_Kg[iS], Ks_BestG_Kg[iS], Time[iS,1:iT_N], iT_N, Inf_3D_Obs[iS,1:iT_N], θ_Ini[iS], Se_Ini[iS], θs[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS]; Flag_Best="Best_G")

			Time_TransStead_BestG_Kg[iS] = best.kg.STEADY_TRANSIT_BESTG(θ_Ini[iS], θs[iS], θr[iS], Hkg_BestG[iS], σ_BestG[iS], Ks_BestG_Kg[iS], RingRadius[iS])

			iT_TransStead_BestG_Kg[iS] = array.SEARCH_INDEX(Time[iS,1:iT_N], Time_TransStead_BestG_Kg[iS])

			θ_Time_BestG_Kg = Array{Float64}(iT_N)
			θ_Time_BestG_Kg[1:iT_N], ~ = thetaTime.kg.θTIME(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius[iS], Sorpt_BestG_Kg[iS], θ_Ini[iS], θs[iS], θr[iS],  σ_BestG[iS], Ks_BestG_Kg[iS], iT_TransStead_BestG_Kg[iS], iT_N)

			NSE_Inf_BestG_Kg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS, 2:iT_N]), log10.(Inf_BestG_Kg[iS, 2:iT_N]))

			Err_Sorpt_BestG_Kg[iS] = stats.RELATIVEerr(Sorpt_Hydro_Kg[iS], Sorpt_BestG_Kg[iS])

			Err_Ks_BestG_Kg[iS] = stats.RELATIVEerr((Ks[iS]), (Ks_BestG_Kg[iS]))

			

			println("BEST_Gi_Vg..........................................................\n")

			Sorpt_BestGi_Vg[iS], Hvg_BestGi[iS], N_BestGi[iS], Ks_BestGi_Vg[iS], Kr_θini_BestGi_Vg[iS], Inf_BestGi_Vg[iS, 1:iT_N] = best.vg.BESTG_INVERSE(Time[iS,1:iT_N], iT_N, Inf_3D_Obs[iS, 1:iT_N], θ_Ini[iS], Se_Ini[iS], θs[iS], Km[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS]; Option_Opt_N=true, Flag_Best="Best_Gi")

			Time_TransStead_BestGi_Vg[iS] = best.vg.STEADY_TRANSIT_BESTG(θ_Ini[iS], θs[iS], θr[iS], Hvg_BestGi[iS], N_BestGi[iS], Ks_BestGi_Vg[iS], Km[iS], RingRadius[iS])

			NSE_Inf_BestGi_Vg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_BestGi_Vg[iS,2:iT_N]))

			Err_Sorpt_BestGi_Vg[iS] = stats.RELATIVEerr(Sorpt_Hydro_Vg[iS], Sorpt_BestGi_Vg[iS])

			Err_Ks_BestGi_Vg[iS] = stats.RELATIVEerr((Ks[iS]), (Ks_BestGi_Vg[iS]))

			println( iS, " ,Sorpt_BestGi_Vg=", Sorpt_BestGi_Vg[iS], " ,N_Best=",N_BestGi[iS], " ,Hvg_Best= ", Hvg_BestGi[iS], " ,Ks_Best_Vg=", 360.*Ks_BestGi_Vg[iS], " ,NSE_Inf_Best_Vg=",  NSE_Inf_BestGi_Vg[iS], "\n ")


			println("BEST_G_Vg..........................................................\n")
			Sorpt_BestG_Vg[iS], Hvg_BestG[iS], N_BestG[iS], Ks_BestG_Vg[iS], Kr_θini_BestG_Vg[iS], Inf_BestG_Vg[iS, 1:iT_N] = best.vg.BESTG_INVERSE(Time[iS,1:iT_N], iT_N, Inf_3D_Obs[iS, 1:iT_N], θ_Ini[iS], Se_Ini[iS], θs[iS], Km[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS]; Option_Opt_N=true, Flag_Best="Best_G")

			NSE_Inf_BestG_Vg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_BestG_Vg[iS,2:iT_N]) )

			Err_Sorpt_BestG_Vg[iS] = stats.RELATIVEerr(Sorpt_Hydro_Vg[iS], Sorpt_BestG_Vg[iS])

			Err_Ks_BestG_Vg[iS] = stats.RELATIVEerr((Ks[iS]), (Ks_BestG_Vg[iS]))


	 #= =============== QUASI EXACT SOLUTION =============== =#
		if option.output.QuasiExact == true
		println("QUASI_EXACT KG..........................................................\n")
		
		Ks_Qe_Kg[iS], Kr_θini_Qe_Kg[iS], Sorpt_Qe_Kg[iS], σ_Qe[iS], Hkg_Qe[iS] = quasiExact.kg.INFILTRATION3D_2_HYDRO(Time[iS,1:iT_N], Inf_3D_Obs[iS,1:iT_N], iT_N, θs[iS], θ_Ini[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS], true)	
		
		σ_Qe[iS], Hkg_Qe[iS] = quasiExact.kg.INFILTRATION3D_2_HYDRO_σMOD(Time[iS,1:iT_N], Inf_3D_Obs[iS,1:iT_N], iT_N, θs[iS], θ_Ini[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS], Ks_Qe_Kg[iS], Sorpt_Qe_Kg[iS])

			# σ_Qe[iS], Hkg_Qe[iS] = quasiExact.kg.INFILTRATION3D_2_HYDRO_σmod(Sorpt_Qe_Kg[iS], Ks_Qe_Kg[iS], Time[iS,1:iT_N], Inf_3D_Obs[iS,1:iT_N], iT_N, θs[iS], θ_Ini[iS], RingRadius[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS],)

			Inf_Qe_Kg[iS,1:iT_N] = quasiExact.HYDRO_2_INFILTRATION3D(Time[iS,1:iT_N], Sorpt_Qe_Kg[iS], Ks_Qe_Kg[iS], Ks_Qe_Kg[iS]*Kr_θini_Qe_Kg[iS], θs[iS], θ_Ini[iS], iT_N, RingRadius[iS], iT_TransStead_Sim[iS])

			NSE_Inf_Qe_Kg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Qe_Kg[iS,2:iT_N]))

			Err_Sorpt_Qe_Kg[iS] = stats.RELATIVEerr(Sorpt_Hydro_Kg[iS], Sorpt_Qe_Kg[iS])

			Err_Ks_Qe_Kg[iS] = stats.RELATIVEerr((Ks[iS]), (Ks_Qe_Kg[iS]))

			println("Sorpt_Qe_Kg=", Sorpt_Qe_Kg[iS]," ,σ_Qe=", σ_Qe[iS], " ,Hkg_Qe=", Hkg_Qe[iS], " ,Ks=",360.*Ks_Qe_Kg[iS]," ,NSE_Inf_Qe_Kg=", NSE_Inf_Qe_Kg[iS], "\n")

			θ_Time_Qe_Kg = Array{Float64}(iT_N)
			θ_Time_Qe_Kg[1:iT_N], ~ = thetaTime.kg.θTIME(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius[iS], Sorpt_Qe_Kg[iS], θ_Ini[iS], θs[iS], θr[iS],  σ_Qe[iS], Ks_Qe_Kg[iS], iT_TransStead_Small[iS], iT_N)


			# VAN GENUCHTEN MODEL ===
			Option_Vg = true
			if Option_Vg

			# Compute Ks ====
			println("QUASI_EXACT VG..........................................................\n")
			Ks_Qe_Vg[iS], Kr_θini_Qe_Vg[iS], Sorpt_Qe_Vg[iS], N_Qe[iS], Hvg_Qe[iS] =quasiExact.vg.INFILTRATION3D_2_HYDRO(Time[iS, 1:iT_N], Inf_3D_Obs[iS,1:iT_N], iT_N, θs[iS], θ_Ini[iS], RingRadius[iS], Km[iS], iT_TransStead_Sim[iS], Time_TransStead_Sim[iS], true)

			Inf_Qe_Vg[iS, 1:iT_N] = quasiExact.HYDRO_2_INFILTRATION3D(Time[iS,1:iT_N], Sorpt_Qe_Vg[iS], Ks_Qe_Vg[iS], Ks_Qe_Vg[iS]*Kr_θini_Qe_Vg[iS], θs[iS], θ_Ini[iS], iT_N, RingRadius[iS], iT_TransStead_Sim[iS])

			NSE_Inf_Qe_Vg[iS] = stats.NASH_SUTCLIFFE(log10.(Inf_3D_Obs[iS,2:iT_N]), log10.(Inf_Qe_Vg[iS,2:iT_N]))

			Err_Sorpt_Qe_Vg[iS] = stats.RELATIVEerr(Sorpt_Hydro_Vg[iS], Sorpt_Qe_Vg[iS])

			Err_Ks_Qe_Vg[iS] = stats.RELATIVEerr((Ks[iS]), (Ks_Qe_Vg[iS]))

			println("Sorpt_Qe_Vg=", Sorpt_Qe_Vg[iS]," ,N_Qe=", N_Qe[iS], " ,Hvg_Qe=", Hvg_Qe[iS], " ,Ks=",360.*Ks_Qe_Vg[iS], "\n")
		end
	end


	#= =============== TESTING =============== =#
		 # Testing derivatives
		 if option.output.TestingDθDh == true
			  println("START to test DθDh... ")
					testbest.DθDH_ACCURACY(iS, θs[iS], θr[iS], Hkg[iS], σ[iS], Hvg[iS], N[iS], Km[iS])
			  println("END to test DθDh \n")
		 end

		 # Testing finding Hvg and Hm from soprtivity
		 if option.output.TestingSorpt == true
			  println("START to test sorptivity... ")
				Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorpt_Hydro_Kg[iS], θ_Ini[iS], θs[iS], θr[iS],  σ[iS], Ks[iS])
				Hkg_2 = Hkg[iS]
				println("$iS Hkg : $Hkg_2 Hkg_Sorpt : $Hkg_Sorpt")
			println("\n")

			if option.HydroModel.HydraulicParam == "vangenuchten"
				Hvg_Sorpt= sorptivity.vg.SORPTIVITY_2_Hvg(Sorpt_Hydro_Vg[iS],  θ_Ini[iS], θs[iS], θr[iS], N[iS], Ks[iS], Km[iS])
				Hvg_2 = Hvg[iS]
				println("$iS Hvg_Sorpt : $Hvg_Sorpt $Hvg_2 ")
				println("END to test sorptivity \n")
			end
		end



	#= =============== PREPERARING PLOTTING & TABLES =============== =#
	if option.Plot.WRC_Kunsat
		iMax = 5000

		Diffusivity_Kg= Array{Float64}(iMax)
		DθDh_Kg=Array{Float64}(iMax)
		DθDr_Vg=Array{Float64}(iMax)
		H_Best_Kg = Array{Float64}(iMax)
		H_Best_Vg =  Array{Float64}(iMax)
		H_Best_Vg = Array{Float64}(iMax)
		# H_BestImp_Kg =  Array{Float64}(iMax)
		# H_BestImp_Vg =  Array{Float64}(iMax)
		H_BestG_Kg = Array{Float64}(iMax)
		H_BestG_Vg = Array{Float64}(iMax)
		H_BestGi_Kg = Array{Float64}(iMax)
		H_BestGi_Vg = Array{Float64}(iMax)
		H_Kg = Array{Float64}(iMax)
		# H_Qe_Imp_Kg = Array{Float64}(iMax)
		# H_Qe_Imp_Vg = Array{Float64}(iMax)
		H_Qe_Kg = Array{Float64}(iMax)
		H_Qe_Select_Kg = Array{Float64}(iMax)
		H_Qe_Select_Vg = Array{Float64}(iMax)
		H_Qe_Vg = Array{Float64}(iMax)
		H_Vg = Array{Float64}(iMax)
		Kunsat_Best_Kg=Array{Float64}(iMax)
		Kunsat_Best_Vg = Array{Float64}(iMax)
		Kunsat_Best_Vg = Array{Float64}(iMax)
		# Kunsat_BestImp_Kg=Array{Float64}(iMax)
		# Kunsat_BestImp_Vg=Array{Float64}(iMax)
		Kunsat_BestG_Kg=Array{Float64}(iMax)
		Kunsat_BestG_Vg = Array{Float64}(iMax)
		Kunsat_BestGi_Kg=Array{Float64}(iMax)
		Kunsat_BestGi_Vg = Array{Float64}(iMax)
		Kunsat_Kg=Array{Float64}(iMax)
		# Kunsat_Qe_Imp_Kg = Array{Float64}(iMax)
		# Kunsat_Qe_Imp_Vg = Array{Float64}(iMax)
		Kunsat_Qe_Kg = Array{Float64}(iMax)
		Kunsat_Qe_Select_Kg = Array{Float64}(iMax)
		Kunsat_Qe_Select_Vg = Array{Float64}(iMax)
		Kunsat_Qe_Vg = Array{Float64}(iMax)
		Kunsat_Vg=Array{Float64}(iMax)

		# Plotting characteristic and unsaturated curves

		Se = linspace(0.001,1., iMax) # Range of Se
		for i in 1:iMax
			θ = wrc.se.Se_2_θ(Se[i], θs[iS], θr[iS])

			# VG H
			H_BestG_Vg[i] = 10.*wrc.vg.Se_2_H(Se[i], Hvg_BestG[iS], N_BestG[iS], Km[iS]) #cm
			H_BestGi_Vg[i] = 10.*wrc.vg.Se_2_H(Se[i], Hvg_BestGi[iS], N_BestGi[iS], Km[iS]) #cm
			H_Qe_Vg[i] = 10.*wrc.vg.Se_2_H(Se[i], Hvg_Qe[iS], N_Qe[iS], Km[iS])
			H_Vg[i] = 10.*wrc.vg.Se_2_H(Se[i], Hvg[iS], N[iS], Km[iS])

			#KG H
			H_BestG_Kg[i] = 10.* wrc.kg.Se_2_H(Se[i], Hkg_BestG[iS], σ_BestG[iS]) #cm	
			H_BestGi_Kg[i] = 10.* wrc.kg.Se_2_H(Se[i], Hkg_BestGi[iS], σ_BestGi[iS]) #cm
			H_Kg[i] = 10.* wrc.kg.Se_2_H(Se[i], Hkg[iS], σ[iS]) #cm
			H_Qe_Kg[i] =  10.* wrc.kg.Se_2_H(Se[i], Hkg_Qe[iS], σ_Qe[iS]) #cm

			# VG Kunsat
			Kunsat_BestG_Vg[i] = 360.*kunsat.vg.KUNSAT(Se[i], N_BestG[iS], Ks_BestG_Vg[iS], Km[iS])
			Kunsat_BestGi_Vg[i] = 360.*kunsat.vg.KUNSAT(Se[i], N_BestGi[iS], Ks_BestGi_Vg[iS], Km[iS])
			Kunsat_Qe_Vg[i] = 360.*kunsat.vg.KUNSAT(Se[i], N_Qe[iS], Ks_Qe_Vg[iS], Km[iS])
			Kunsat_Vg[i] = 360.*kunsat.vg.KUNSAT(Se[i], N[iS], Ks[iS], Km[iS])
			
			#KG Kunsat
			Kunsat_BestG_Kg[i] = 360.*kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ_BestG[iS], Ks_BestG_Kg[iS], θs[iS], σ_BestG[iS])
			Kunsat_BestGi_Kg[i] = 360.*kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ_BestGi[iS], Ks_BestGi_Kg[iS], θs[iS], σ_BestGi[iS])
			Kunsat_Kg[i] = 360.*kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ[iS], Ks[iS], θs[iS], σ[iS]) #cm/h
			Kunsat_Qe_Kg[i] = 360.*kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ_Qe[iS], Ks_Qe_Kg[iS], θs[iS], σ_Qe[iS])
		end

		# Statistics of the hydraulic parameters
		NSE_BestG_Kg[iS] = (stats.NASH_SUTCLIFFE(log10.(H_Kg[1:iMax]), log10.(H_BestG_Kg[1:iMax])) + stats.NASH_SUTCLIFFE(log10.(Kunsat_Kg[1:iMax]), log10.(Kunsat_BestG_Kg[1:iMax]))) / 2.
		NSE_BestG_Vg[iS] =  (stats.NASH_SUTCLIFFE(log10.(H_Vg[1:iMax]), log10.(H_BestG_Vg[1:iMax])) + stats.NASH_SUTCLIFFE(log10.(Kunsat_Vg[1:iMax]), log10.(Kunsat_BestG_Vg[1:iMax]))) / 2.
		NSE_BestGi_Kg[iS] = (stats.NASH_SUTCLIFFE(log10.(H_Kg[1:iMax]), log10.(H_BestGi_Kg[1:iMax])) + stats.NASH_SUTCLIFFE(log10.(Kunsat_Kg[1:iMax]), log10.(Kunsat_BestGi_Kg[1:iMax]))) / 2.
		NSE_BestGi_Vg[iS] =  (stats.NASH_SUTCLIFFE(log10.(H_Vg[1:iMax]), log10.(H_BestGi_Vg[1:iMax])) + stats.NASH_SUTCLIFFE(log10.(Kunsat_Vg[1:iMax]), log10.(Kunsat_BestGi_Vg[1:iMax]))) / 2.
		NSE_Qe_Kg[iS] = (stats.NASH_SUTCLIFFE(log10.(H_Kg[1:iMax]), log10.(H_Qe_Kg[1:iMax])) + stats.NASH_SUTCLIFFE(log10.(Kunsat_Kg[1:iMax]), log10.(Kunsat_Qe_Kg[1:iMax]))) / 2.
		NSE_Qe_Vg[iS] = (stats.NASH_SUTCLIFFE(log10.(H_Vg[1:iMax]), log10.(H_Qe_Vg[1:iMax])) + stats.NASH_SUTCLIFFE(log10.(Kunsat_Vg[1:iMax]), log10.(Kunsat_Qe_Vg[1:iMax]))) / 2.

		println("NSE_BestG_Kg = ", NSE_BestG_Kg[iS],", ", "NSE_Qe_Kg= ",NSE_Qe_Kg[iS])


		#------------------WRITING TABLES----------------------------=====================================#
		println("START to write output... ")
		Data = [Id[iS], Se_Ini[iS], θs[iS], θr[iS], Ks[iS], N[iS], Hvg[iS], Km[iS], σ[iS], Hkg[iS], Err_Ks_BestG_Kg[iS], Err_Ks_BestG_Vg[iS], Err_Ks_BestGi_Kg[iS], Err_Ks_BestGi_Vg[iS], Err_Ks_Qe_Kg[iS], Err_Ks_Qe_Vg[iS], Err_Qe[iS], Err_Sorpt_BestG_Kg[iS], Err_Sorpt_BestG_Vg[iS], Err_Sorpt_BestGi_Kg[iS], Err_Sorpt_BestGi_Vg[iS], Err_Sorpt_Qe_Kg[iS], Err_Sorpt_Qe_Vg[iS], Hkg_Qe[iS], Hvg_BestG[iS], Hvg_BestGi[iS], Hvg_Qe[iS], Hkg_BestG[iS], Kr_θini_Vg[iS], Ks_BestG_Kg[iS], Ks_BestG_Vg[iS], Ks_BestGi_Kg[iS], Ks_BestGi_Vg[iS], Ks_Qe_Kg[iS], Ks_Qe_Vg[iS], N_BestG[iS], N_BestGi[iS], N_Qe[iS], NSE_Hydro_Inf_BestG_Kg[iS], NSE_Hydro_Inf_BestG_Vg[iS], NSE_Hydro_Inf_BestG_Vg[iS], NSE_Hydro_Inf_BestGi_Kg[iS], NSE_Hydro_Inf_BestGi_Vg[iS], NSE_Hydro_Inf_Qe_Kg[iS], NSE_Hydro_Inf_Qe_Vg[iS], NSE_Inf_BestG_Vg[iS],	NSE_Inf_BestGi_Kg[iS],	NSE_Inf_BestGi_Kg[iS],	NSE_Inf_BestGi_Vg[iS],	NSE_Inf_Qe_Kg[iS], NSE_Inf_Qe_Vg[iS], NSE_BestGi_Kg[iS], NSE_BestGi_Vg[iS],NSE_BestG_Kg[iS], NSE_BestG_Vg[iS], NSE_Qe_Kg[iS],	NSE_Qe_Vg[iS],	Sorpt_BestG_Kg[iS],	Sorpt_BestGi_Kg[iS],	Sorpt_BestGi_Vg[iS],	Sorpt_Hydro_Kg[iS], Sorpt_Hydro_Vg[iS],	Sorpt_Qe_Kg[iS],	Time_TransStead_BestG_Vg[iS],	Time_TransStead_Hydro_Kg[iS], Time_TransStead_Hydro_Vg[iS], Hkg_BestGi[iS], NSE_Hydro_Inf_BestGi_Vg[iS], NSE_Inf_BestG_Kg[iS], Sorpt_BestG_Kg[iS], σ_BestG[iS], σ_BestGi[iS], σ_Qe[iS]]

		# We are appending the table
		CSV.write(Path_Table, DataFrame(Data); append=true)
		println("END of writing output \n ")

			plots.kg.UNSAT_CHARAC_PLOT(Name[:])

			# Boundaries for plotting
			Se_Min = 0
			Se_Max = 1.1
			H_Min = 0.1
			H_Max = 10^8 # Equivalent to 2000kpa, we are working in cm
			Kmin = 10^-5.
			Kmax = 360. * Ks[iS]
			
			Title = string("\$\\sigma =", string(round.(σ[iS],2)),", " ,Name[iS], " \$")
			println(Title)

			# Plotting the model characteristic curves and hydraulic conductivity ===================================================================================================
			if iS <= 10
				push!(Plot_CharacUnsat, Axis([
					Plots.Linear(Se, H_Kg, mark="none", style="red, very thick"),
					Plots.Linear(Se, H_BestGi_Kg, mark="none", style="loosely dashed, blue, very thick"),
					# Plots.Linear(Se, H_BestGi_Kg, mark="none", style="densely dashed, teal, very thick"),
					Plots.Linear(Se, H_Qe_Kg, mark="none", style="densely dashdotdotted, orange, very thick"),
				], ylabel=L"$h \ [cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=H_Min, ymax=H_Max, ymode="log", style="width=10.4cm, height=6.4cm"))

				push!(Plot_CharacUnsat, Axis([
					Plots.Linear(Se, Kunsat_Kg, mark="none", style="red, very thick"),
					# Plots.Linear(Se, Kunsat_Qe_Vg, mark="none",  style="densely dashed, teal, very thick"),
					Plots.Linear(Se, Kunsat_BestGi_Kg, mark="none", style="densely dashed, blue, very thick"),
					Plots.Linear(Se, Kunsat_Qe_Kg, mark="none", style="densely dashdotdotted, orange, very thick"),
					Plots.Node(eval(Title),0.05, Kmax - Kmax / 10., style="right"),
				] , ylabel=L"$K(S_e) \ [cm \ h^{-1}]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin, ymax=Kmax, style="width=10.4cm, height=6.4cm"))
			else
				push!(Plot_CharacUnsat, Axis([
					Plots.Linear(Se, H_Kg, mark="none", style="red, very thick"),
					Plots.Linear(Se, H_BestGi_Kg, mark="none", style="loosely dashed, blue, very thick"),
					# Plots.Linear(Se, H_BestGi_Kg, mark="none", style="densely dashed, teal, very thick"),
					Plots.Linear(Se, H_Qe_Kg, mark="none", style="densely dashdotdotted, orange, very thick"),
				], xlabel=L"$S_e \ [-]$", ylabel=L"$h \ [cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=H_Min, ymax=H_Max, ymode="log", style="width=10.4cm, height=6.4cm", legendStyle = "{at={(-0.2,-0.4)}, anchor=south east, legend columns=-1}"))

				push!(Plot_CharacUnsat, Axis([
					Plots.Linear(Se, Kunsat_Kg, mark="none", style="red, very thick", legendentry=L"$_{HYDRUS}$"),
					Plots.Linear(Se, Kunsat_BestGi_Kg, mark="none", style="densely dashed, blue, very thick", legendentry=L"$_{BEST_{SA\_KG}}$"),
					Plots.Linear(Se, Kunsat_Qe_Kg, mark="none",  style="densely dashdotdotted, orange, very thick", legendentry=L"$_{BEST_{QEI\_KG}}$"),
					# Plots.Linear(Se, Kunsat_Qe_Vg, mark="none", style="densely dashed, teal, very thick", legendentry=L"$_{BEST_{QEI\_VG}}$"),
					Plots.Node(eval(Title),0.05, Kmax - Kmax / 10., style="right"),
				],xlabel=L"$S_e \ [-]$", ylabel=L"$K(S_e) \ [cm \ h^{-1}]$",
				style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin, ymax=Kmax, style="width=10.4cm, height=6.4cm", legendStyle = "{at={(-0.2,-0.4)}, anchor=south west, legend columns=-1}"))
		
				save(Path_CharacUnsat, Plot_CharacUnsat)
			end
			save(Path_CharacUnsat, Plot_CharacUnsat)
			
			# # van Genuchten
			# push!(Plot1, Axis([
			# 	# Plots.Linear(Se, H_Best_Vg, mark="none", style="cyan"),
			# 	Plots.Linear(Se, H_Vg, mark="none", style="red, very thick"),
			# 	Plots.Linear(Se, H_BestG_Vg, mark="none", style="dashed, blue, very thick"),
			# 	Plots.Linear(Se, H_BestGi_Vg, mark="none", style="dashed, orange, very thick"),
			# 	Plots.Linear(Se, H_Qe_Vg, mark="none", style="dashdotdotted, brown, very thick"),

			# ], title="VANGENUCHTEN H(Se)", xlabel=L"$Se[-]$", ylabel=L"$H[cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymode="log")
			# )

			# push!(Plot1, Axis([
			# 	Plots.Linear(Se, Kunsat_Vg, mark="none", style="red, very thick", legendentry=L"$Obs$"),
			# 	Plots.Linear(Se, Kunsat_BestG_Vg, mark="none", style="dashed, blue, very thick", legendentry=L"$BestG$"),
			# 	Plots.Linear(Se, Kunsat_BestGi_Vg, mark="none", style="dashed, orange, very thick", legendentry=L"$BestGi$"),
			# 	Plots.Linear(Se, Kunsat_Qe_Vg, mark="none", style="dashdotdotted, brown, very thick", legendentry=L"$QE$"),
			# ], title="VANGENUCHTEN K(Se)", xlabel=L"$Se[-]$", ylabel=L"$K(Se)[cm/h]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin, legendStyle = "{at={(-0.3,-0.4)},anchor=south west, legend columns=-1}"))
		

			Flag_Plot_InfHydro = true #--------------------------------------------------------------------------------------------------------------
			if Flag_Plot_InfHydro
				
	
			# Inf curves
			Time_2 = Time[iS, 1:iT_N] / 60. # To minutes
			Inf_2 = Inf_3D_Obs[iS, 1: iT_N] / 10. #To cm 
			Inf_Hydro_BestG_Kg2 = Inf_Hydro_BestG_Kg[iS, 1:iT_N]/10.
			Inf_Hydro_BestGi_Kg2 = Inf_Hydro_BestGi_Kg[iS, 1:iT_N]/10.
			Inf_Hydro_Qe_Kg2 = Inf_Hydro_Qe_Kg[iS, 1:iT_N]/10.
			Inf_Hydro_BestG_Vg2 = Inf_Hydro_BestG_Vg[iS, 1:iT_N]/10.
			Inf_Hydro_BestGi_Vg2 = Inf_Hydro_BestGi_Vg[iS, 1:iT_N]/10.
			Inf_Hydro_Qe_Vg2 = Inf_Hydro_Qe_Vg[iS, 1:iT_N]/10

			iT_TransStead_Hydro_Kg = array.SEARCH_INDEX(Time[iS, 1:iT_N], Time_TransStead_Hydro_Kg[iS])
			TransStead_X_Best_Kg = Time[iS, iT_TransStead_Hydro_Kg] / 60. # Points of where TransStead occures
			TransStead_Y_Best_Kg = Inf_Hydro_BestG_Kg[iS, iT_TransStead_Hydro_Kg ]/ 10.

			iT_TransStead_Hydro_Vg = array.SEARCH_INDEX(Time[iS, 1:iT_N], Time_TransStead_Hydro_Vg[iS])
			TransStead_X_Best_Vg = Time[iS, iT_TransStead_Hydro_Vg] / 60. # Points of where TransStead occures
			TransStead_Y_Best_Vg = Inf_Hydro_BestG_Vg[iS, iT_TransStead_Hydro_Vg]/ 10.

			Inf_Max = max(maximum(Inf_2[:]), maximum(Inf_Hydro_BestGi_Kg2[:]), maximum(Inf_Hydro_Qe_Vg2[:]), maximum(Inf_Hydro_Qe_Kg2[:]))

			if iS <= 10
				push!(Plot_InfHydro, Axis([
					Plots.Linear(Time_2, Inf_2, mark="none", style="red, thick"),
					# Plots.Linear(Time_2, Inf_Hydro_BestG_Kg2, mark="none", style="loosely dashed, blue, very thick"),
					Plots.Linear(Time_2, Inf_Hydro_BestGi_Kg2, mark="none", style="loosely dashed, blue, very thick"),
					Plots.Linear(Time_2, Inf_Hydro_Qe_Kg2, mark="none", style="densely dashdotdotted, orange, very thick"),
# 					Plots.Scatter(TransStead_X_Best_Kg, TransStead_Y_Best_Kg, style="magenta, very thick", onlyMarks=true, mark="o", markSize=5),
				], ylabel=L"$Infiltration \ [cm]$", style="width=10.4cm, height=6.4cm"))

				push!(Plot_InfHydro, Axis([
					Plots.Linear(Time_2, Inf_2, mark="none", style="red, thick"),
					# Plots.Linear(Time_2, Inf_Hydro_BestG_Vg2, mark="none", style="loosely dashed, blue, very thick"),
					# Plots.Linear(Time_2, Inf_Hydro_BestGi_Vg2, mark="none", style="densely dashed, teal, very thick"),
					Plots.Linear(Time_2, Inf_Hydro_Qe_Kg2, mark="none", style="densely dashdotdotted, orange, very thick"),
					Plots.Linear(Time_2, Inf_Hydro_Qe_Vg2, mark="none", style="densely dashed, teal, very thick"),
					Plots.Node(eval(Title),maximum(Time_2)/100., Inf_Max, style="right"),
# 					Plots.Scatter(TransStead_X_Best_Vg, TransStead_Y_Best_Vg, style="magenta, very thick", onlyMarks=true, mark="o", markSize=5),
				], ylabel=L"$Infiltration \ [cm]$", style="smooth", style="width=10.4cm,
				height=6.4cm"))
				
			else
				push!(Plot_InfHydro, Axis([
					Plots.Linear(Time_2, Inf_2, mark="none", style="red, thick", legendentry=L"$_{HYDRUS}$"),
					# Plots.Linear(Time_2, Inf_Hydro_BestG_Kg2, mark="none", style="loosely dashed, blue, very thick"),
					Plots.Linear(Time_2, Inf_Hydro_BestGi_Kg2, mark="none", style="loosely dashed, blue, very thick",  legendentry=L"$_{BEST_{SA\_KG}}$"),
					Plots.Linear(Time_2, Inf_Hydro_Qe_Kg2, mark="none", style="densely dashdotdotted, orange, very thick",  legendentry=L"$_{BEST_{QEI\_KG}}$"),
# 					Plots.Scatter(TransStead_X_Best_Kg, TransStead_Y_Best_Kg, style="magenta, very thick", onlyMarks=true, mark="o", markSize=5),
				], xlabel=L"$Time \ [min]$", ylabel=L"$Infiltration \ [cm]$", style="smooth", style="width=10.4cm, height=6.4cm", legendStyle = "{at={(-0.5,-0.4)}, anchor=south west, legend columns=-1}",))

				push!(Plot_InfHydro, Axis([
					Plots.Linear(Time_2, Inf_2, mark="none", style="red, thick", legendentry=L"$_{HYDRUS}$"),
					Plots.Linear(Time_2, Inf_Hydro_Qe_Kg2, mark="none",  style="densely dashdotdotted, orange, very thick", legendentry=L"$_{BEST_{QEI\_KG}}$"),
					# Plots.Linear(Time_2, Inf_Hydro_BestGi_Vg2, mark="none", style="densely dashed, teal, very thick", legendentry=L"$_{BestGi}$"),
					Plots.Linear(Time_2, Inf_Hydro_Qe_Vg2, mark="none", style="densely dashed, teal, very thick", legendentry=L"$_{BEST_{QEI\_VG}}$"),
					Plots.Node(eval(Title),maximum(Time_2)/100., Inf_Max, style="right"),
# 					Plots.Scatter(TransStead_X_Best_Vg, TransStead_Y_Best_Vg, style="magenta, very thick", onlyMarks=true, mark="o", markSize=5, legendentry=L"$_{Tmod}$"),
				],  xlabel=L"$Time \ [min]$", ylabel=L"$Infiltration \ [cm]$", style="smooth", style="width=10.4cm, height=6.4cm", legendStyle = "{at={(0.2,-0.4)}, anchor=south east, legend columns=-1}",))
				
				save(Path_InfHydro, Plot_InfHydro)
			end
			save(Path_InfHydro, Plot_InfHydro)
			end

			# PLOTING INFILTRATION #--------------------------------------------------------------------------------------------------------------
			# iT_TransStead = iT_TransStead_Sim[iS]
			# TimeTransit = Time[iS, iT_TransStead_Sim[iS]]
				
			# Inf_Qe_Vg2 = Inf_Qe_Vg[iS, 1:iT_N]/10.
			# iT_TransStead_Sim[iS], Time_TransStead_Sim[iS]
			# Time_TransStead_BestG_Kg Time_TransStead_BestG_Vg

			Inf_BestGi_Kg2 = Inf_BestGi_Kg[iS, 1:iT_N] / 10.
			Inf_BestG_Kg2 = Inf_BestG_Kg[iS, 1:iT_N] / 10.
			Inf_Qe_Kg2 = Inf_Qe_Kg[iS, 1:iT_N] / 10.
			
			TransStead_X_Mod_Kg = Time[iS, iT_TransStead_Sim[iS]] / 60. # Points of where TransStead occures
			TransStead_Y_Mod_Kg = Inf_3D_Obs[iS, iT_TransStead_Sim[iS]] / 10.

			# iT_TransStead_BestG_Kg = array.SEARCH_INDEX(Time[iS, 1:iT_N], Time_TransStead_BestG_Kg[iS])
			TransStead_X_Best_Kg = Time[iS, iT_TransStead_BestG_Kg[iS]] / 60. # Points of where TransStead occures
			TransStead_Y_Best_Kg = Inf_3D_Obs[iS, iT_TransStead_BestG_Kg[iS]] / 10.

			TransStead_X_Mod_Vg = Time[iS, iT_TransStead_Sim[iS]] / 60. # Points of where TransStead occures
			TransStead_Y_Mod_Vg = Inf_3D_Obs[iS, iT_TransStead_Sim[iS]] / 10.

			iT_TransStead_BestG_Vg = array.SEARCH_INDEX(Time[iS, 1:iT_N], Time_TransStead_BestG_Vg[iS])
			TransStead_X_Best_Vg = Time[iS, iT_TransStead_BestG_Vg] / 60. # Points of where TransStead occures
			TransStead_Y_Best_Vg = Inf_3D_Obs[iS, iT_TransStead_BestG_Vg] / 10.

			Inf_Max = max(maximum(Inf_2[:]), maximum(Inf_BestG_Kg2[:]), maximum( Inf_Qe_Kg2[:]))
			
			if iS <= 10
				push!(Plot_Infiltration, Axis([
				Plots.Linear(Time_2, Inf_2, mark="none", style="red, thick"),
				Plots.Linear(Time_2, Inf_BestG_Kg2, mark="none", style="loosely dashed, blue, very thick"),
				# Plots.Linear(Time_2, Inf_BestGi_Kg2, mark="none", style="densely dashed, teal, very thick"),
				Plots.Linear(Time_2, Inf_Qe_Kg2, mark="none", style="densely dashdotdotted, orange, very thick"),
				Plots.Scatter(TransStead_X_Mod_Kg, TransStead_Y_Mod_Kg, style="cyan, very thick",  mark="square", onlyMarks=true, markSize=4),
				Plots.Scatter(TransStead_X_Best_Kg, TransStead_Y_Best_Kg, style=" brown, very thick",  mark="o", onlyMarks=true, markSize=4),
				Plots.Node(eval(Title),maximum(Time_2)/100., Inf_Max, style="right"),
			], ylabel=L"$Infiltration \ [cm]$", style="smooth, width=10.4cm, height=6.4cm"))
			else
				push!(Plot_Infiltration, Axis([
				Plots.Linear(Time_2, Inf_2, mark="none", style="red, thick", legendentry=L"$_{HYDRUS}$"),
				Plots.Linear(Time_2, Inf_BestG_Kg2, mark="none", style="loosely dashed, blue, very thick", legendentry=L"$_{BEST_{SA}}$"),
				# Plots.Linear(Time_2, Inf_BestGi_Kg2, mark="none", style="densely dashed, teal, very thick", legendentry=L"$_{BestGi}$"),
				Plots.Linear(Time_2, Inf_Qe_Kg2, mark="none", style="densely dashdotdotted, orange, very thick", legendentry=L"$_{BEST_{QEI}}$"),
				Plots.Scatter(TransStead_X_Mod_Kg, TransStead_Y_Mod_Kg, style="cyan, very thick",
				mark="square", onlyMarks=true, markSize=4, legendentry=L"$_{T\_trans\_steady\_QEI}$"),
				Plots.Scatter(TransStead_X_Best_Kg, TransStead_Y_Best_Kg, style=" brown, very thick",
				  mark="o", onlyMarks=true, markSize=4, legendentry=L"$_{T\_trans\_steady\_SA}$"),
				Plots.Node(eval(Title),maximum(Time_2)/100.,Inf_Max, style="right"),
			],  xlabel=L"$Time \ [min]$", ylabel=L"$Infiltration \ [cm]$", style="smooth, width=10
			.4cm, height=6.4cm", legendStyle = "{at={(-0.2,-0.4)}, anchor=south west, legend
			columns=-1}",))

			save(Path_Infiltration, Plot_Infiltration)
			end
			save(Path_Infiltration, Plot_Infiltration)
			
			# push!(Plot2b, Axis([
			# 	Plots.Linear(Time_2, Inf_2, mark="none", style="red, thick", legendentry=L"$Obs$"),
			# 	Plots.Linear(Time_2, Inf_BestG_Vg2, mark="none", style="dashed, blue, very thick", legendentry=L"$GBEST$"),
			# 	Plots.Linear(Time_2, Inf_Qe_Vg2 , mark="none", style="dashdotdotted, brown, very thick", legendentry=L"$Qe$"),
			# 	Plots.Scatter(TransStead_X_Mod_Vg, TransStead_Y_Mod_Vg, style="magenta, very thick", onlyMarks=true, markSize=2, legendentry=L"$Tmod$"),
			# 	Plots.Scatter(TransStead_X_Best_Vg, TransStead_Y_Best_Vg, style="brown, very thick", onlyMarks=true, markSize=3,  legendentry=L"$Tbest$"),
			# ],  xlabel=L"$Time \ [min]$", ylabel=L"$Infiltration \ [cm]$", style="smooth, width=10.4cm, height=6.4cm", legendStyle = "{at={(-0.4,-0.4)}, anchor=south west, legend columns=-1}",))


			# PLOTING SOIL MOISTURE-------------------------------------------------------------------------------------------------------------
			Path_Sm = string(DIR_Working, "//OUTPUT//Figure//SoilMoisture//SoilMoisture_", path.FileName, "_", string(iS), ".svg")
			println(Path_Sm)
			if isfile(Path_Sm) # If the file is there than delete it before saving figure
				rm(Path_Sm)
			end
			Time_3 = Time[iS, 1:iT_N] / 60.
			θ_Time_BestG_Kg_2 =  θ_Time_BestG_Kg[1:iT_N]
			θ_Time_BestGi_Kg_2 =  θ_Time_BestGi_Kg[1:iT_N]
			θ_Time_Qe_Kg_2 = θ_Time_Qe_Kg[1:iT_N]

		
			Plot_Sm = GroupPlot(2, 8, groupStyle = "horizontal sep = 2cm, vertical sep = 2cm")
			push!(Plot_Sm,Axis([
				Plots.Linear(Time_2 , θ_Time_BestG_Kg_2, mark="none", style="loosely dashed, blue, very thick", legendentry=L"$_{BEST_{SA}}$")
				# Plots.Linear(Time_2 , θ_Time_BestGi_Kg_2, mark="none", style="densely dashed, teal, very thick", legendentry=L"$_{BestGi}$")
				Plots.Linear(Time_2, θ_Time_Qe_Kg_2, mark="none", style="densely dashdotdotted, orange, very thick", legendentry=L"$_{BEST_{QEI}}$")
				# Plots.Linear(Time_3, θ_Time_Vg_2, mark="none", style="cyan", legendentry=L"$\theta_{vg}$")
			],  xlabel=L"$Time \ [min]$", ylabel=L"$\theta(t)$", style="smooth, width=10.4cm,
			height=6.4cm", legendStyle = "{at={(-0.2,-0.4)}, anchor=south west, legend columns=-1}",))
			
			save(Path_Sm, Plot_Sm) 		
			
		end	
		println("END plotting \n")	
	end

end

end # Core true 

		#End of plotting
		plots.kg.SORPTIVITY_PLOT()
		plots.kg.KSAT_PLOT()
end

START_BEST()
#     tic()
#     t1 = round(toq(), 3)
#     println("Time: serial = $(t1)s"
