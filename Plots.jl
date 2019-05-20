cd("C:\\JOE\\Main\\MODELS\\SOIL\\BEST")
DIR_Working = pwd() # Saving the current working directory
println( "Working path : $DIR_Working")
push!(LOAD_PATH, DIR_Working) # All the subroutines are in the current working directory
 
 # Do not change order
# include("Path.jl")
# include("Cst.jl")
# include("Param.jl")
# include("Tools.jl")
# include("Option.jl")
# include("Wrc.jl")
# include("Kunsat.jl")
# include("Diffusivity.jl")
# include("Stat.jl")
# include("HydroRelationship.jl")
# include("Best.jl")
# include("ThetaTime.jl")
# include("QuasiExact.jl")
# include("ConvertHydroModel.jl")
# include("Reading.jl")
# include("Sorptivity.jl")
# include("KunsatModel.jl")
# include("TestBest.jl")
# include("Err.jl")


module plots
	 module kg
			using wrc, kunsat, diffusivity, PGFPlots, path, LaTeXStrings, reading, param, array, relationship
			export DIFFUSIVITY_PLOT, KSAT_PLOT, SORPTIVITY_PLOT, UNSAT_CHARAC_PLOT, σMODEL_PLOT

			function DIFFUSIVITY_PLOT(iS, θ_Ini, θs, θr, Hkg, σ, Ks)           
				 # Axis(Plots.Linear([1,2], [1,2]), xlabel=L"x_\text{awesome}", ylabel="\$\\text{mathOperation}(x)\$")			
				 N = 1000
				 Se2 = linspace(0.0001,0.999999999999, N)
				 Diff =  Array{Float64}(N)
				 for i in 1:N
						Diff[i] = diffusivity.kg.DIFFUSIVITY(wrc.se.Se_2_θ(Se2[i], θs, θ_Ini) , θs, θr, Hkg, σ, Ks) * 3600 /  100.
				 end
				 
				 DIR_Working = pwd() 
				 Path = string(DIR_Working, "//OUTPUT//Figure//Diffusivity//Diff_", path.FileName, "_", string(iS), ".svg")
				 println(Path)
				 if isfile(Path) # If the file is there than delete it before saving figure
						rm(Path)
				 end

				 Legend = string("\$\\sigma =", string( round(σ,2))," \$")

				 pushPGFPlotsPreamble("\\usepackage{amsmath}")
				 Plot_Sorpt = Axis([	
						Plots.Linear(Se2, Diff, mark="none", style="blue, very thick"),
						Plots.Node(eval(Legend),0., maximum(Diff)/5., style="right"),
				 ],
				 style="width=16.18cm, height=10cm", title =" ", xlabel=L"$S_e[-]$", ylabel=L"$Diffusivity \ [cm^2 \ h^{-1}]$", style="smooth", xmin=0., xmax=1.05, ymode="log", legendPos="north west")
				 save(Path, Plot_Sorpt)
			end



			function KSAT_PLOT()
				 DIR_Working = pwd() 

				 # path.FileName
				 Path_Table =  string(DIR_Working, path.Table, path.FileName, ".csv")
				 Data = Array(readdlm(Path_Table, ',', header=true, skipstart=0)[1]) # skip first line with header\
				 
				 # Converting from mm/s -> cm/h
				 Ks_Obs = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Ks")] * 360.
				 Ks_Qe = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Ks_Qe_Kg")] * 360.
				 Ks_BestG = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Ks_BestG_Kg")]  * 360.
				 Ks_BestGi = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Ks_BestGi_Kg")]  * 360.

				 Ks_Min = min(minimum(Ks_Obs), minimum(Ks_Qe), minimum(Ks_BestG), minimum(Ks_BestGi))
				 Ks_Min -= Ks_Min / 10.

				 Ks_Max = max(maximum(Ks_Obs), maximum(Ks_Qe), maximum(Ks_BestG), maximum(Ks_BestGi))
				 Ks_Max += Ks_Max / 10

				 Path = string(DIR_Working, "//OUTPUT//Figure//Ksat//Ksat_", path.FileName,  ".svg")
				 println(Path)
				 if isfile(Path) # If the file is there than delete it before saving figure
						rm(Path)
				 end

				 pushPGFPlotsPreamble("\\usepackage{amsmath}")
				 Plot_Ks = Axis([	
						Plots.Scatter(Ks_Obs, Ks_Qe, style="orange, very thick", onlyMarks=true, mark="o", markSize=4, legendentry=L"$_{BEST_{QEI}}$"),
						Plots.Scatter(Ks_Obs, Ks_BestG, style="blue, very thick", onlyMarks=true, mark="square", markSize=4, legendentry=L"$_{BEST_{SA}}$"),
# 						Plots.Scatter(Ks_Obs, Ks_BestGi, style="teal, very thick", onlyMarks=true, mark="diamond", markSize=4, legendentry=L"$_{BestGi}$"),
						Plots.Linear(x-> x, (Ks_Min,Ks_Max), xbins=50, style="dashdotdotted")
				 ],
				 style="width=8cm, height=8cm", title =" ", xlabel=L"$K_s \ Obs \ [cm \ h^{-1}]$", ylabel=L"$K_s \ Sim \ [cm \ h^{-1}]$", xmin=Ks_Min, xmax=Ks_Max, ymin=Ks_Min, ymax=Ks_Max, xmode="log", ymode="log", legendStyle = "{at={(0.14,-0.26)}, anchor=south west, legend columns=-1}")
				 save(Path, Plot_Ks)
			end



			function SORPTIVITY_PLOT()
				 DIR_Working = pwd() 

				 # path.FileName
				 Path_Table =  string(DIR_Working, path.Table, path.FileName, ".csv")
				 Data = Array(readdlm(Path_Table, ',', header=true, skipstart=0)[1]) # skip first line with header\
				 S_Obs = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Sorpt_Hydro_Kg")] * 3600^0.5 / 10.
				 S_Qe = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Sorpt_Qe_Kg")] * 3600^0.5 / 10.
				 S_BestG = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Sorpt_BestG_Kg")] * 3600^0.5 / 10.
				 S_BestGi = Data[:,array.SEARCH_INDEX_STRING(param.Header, "Sorpt_BestGi_Kg")] * 3600^0.5 / 10.

				 S_Min = min(minimum(S_Obs), minimum(S_Qe), minimum(S_BestG), minimum(S_BestGi)) 
				 S_Max = max(maximum(S_Obs), maximum(S_Qe), maximum(S_BestG), maximum(S_BestGi)) + 1.

				 Path = string(DIR_Working, "//OUTPUT//Figure//Sorptivity//Sorpty_", path.FileName,  ".svg")
				 println(Path)
				 if isfile(Path) # If the file is there than delete it before saving figure
						rm(Path)
				 end

				 Plot_S = Axis([	
						Plots.Scatter(S_Obs, S_Qe, style="orange, very thick", onlyMarks=true, mark="o", markSize=4, legendentry=L"$_{BEST_{QEI}}$"),
						Plots.Scatter(S_Obs, S_BestG, style="blue, very thick", onlyMarks=true, mark="square", markSize=4, legendentry=L"$_{BEST_{SA}}$"),
# 						Plots.Scatter(S_Obs, S_BestGi, style="teal, very thick", onlyMarks=true, mark="diamond", markSize=4, legendentry=L"$_{BestGi}$"),
						Plots.Linear(x-> x, (0., S_Max), xbins=50, style="dashdotdotted")
				 ],
				 style="width=8cm, height=8cm", title =" ", xlabel=L"$S \ Obs \ [cm \ h^{-0.5}]$", ylabel=L"$S \ Sim \ [cm \ h^{-0.5}]$", xmin=0., xmax=S_Max, ymin=0., ymax=S_Max, legendStyle = "{at={(0.14,-0.26)}, anchor=south west, legend columns=-1}")
				 save(Path, Plot_S)
			end



			function σMODEL_PLOT()
				DIR_Working = pwd() 

				Path = string(DIR_Working, "//OUTPUT//Figure//SigmaRelationship//SigmaRelationship.svg") 

				Path_Hydraulic = string(DIR_Working, path.Read_Hydraulic)
				Se_Ini, θs, θr, σ, Hkg, Ks, Hvg, N, Km, RingRadius, SampleTrue = reading.HYDRAULIC(Path_Hydraulic)

				σ_Obs = [0.98, 1.60,1.65,1.79, 2.00, 2.74]


				N_Data = length(σ_Obs)
				σ_Mod =  Array{Float64}(N_Data)
				# σ_Obs =  Array{Float64}(N_Data)

				
				iCount = 1
				for i in 1:N_Data
					# if SampleTrue[i] == 1
						σ_Mod[iCount] = relationship.Hkg_2_σ(Hkg[iCount])
						# σ_Obs[iCount] = σ[iCount]
						iCount += 1
					# end
				end

				iCount = iCount - 1

				σ_Max = max(maximum(σ_Mod[1:iCount]), maximum(σ_Obs[1:iCount])) + .2

				Plot_σ = Axis([	
					Plots.Scatter(σ_Obs[1:iCount], σ_Mod[1:iCount], style="blue, very thick", onlyMarks=true, mark="o", markSize=4),
					Plots.Linear(x-> x, (0., σ_Max), xbins=50, style="dashdotdotted")
			 ],
			 style="width=8cm, height=8cm", title =" ", xlabel=L"$\sigma \ Obs \ [-]$", ylabel=L"$\sigma \ Sim \ [-]$", xmin=0.65, xmax=σ_Max, ymin=0.65, ymax=σ_Max)
			 save(Path, Plot_σ)
		
			end



			function UNSAT_CHARAC_PLOT(Name)
				 DIR_Working = pwd() 

				 Path_Hydraulic = string(DIR_Working, path.Read_Hydraulic)
				 Se_Ini, θs, θr, σ, Hkg, Ks, Hvg, N, Km, RingRadius, SampleTrue = reading.HYDRAULIC(Path_Hydraulic)

				 Ratio_Width_Height = 0.8
 
				 N_Plots = 1000
				 H_Vg = Array{Float64}(N_Plots)
				 H_Kg = Array{Float64}(N_Plots)
				 Kunsat_Kg = Array{Float64}(N_Plots)
				 Kunsat_Vg = Array{Float64}(N_Plots)

				 DIR_Working = pwd()
				 Path = string(DIR_Working, "//OUTPUT//Figure//KosVang//KosVang.svg") #~~~~~~~~~~~~~~~~~~~~~~~~~~~
				 if isfile(Path) # If the file is there than delete it before saving figure
						rm(Path)
				 end

				 println(Path)
				 Plot1b = GroupPlot(2, 8, groupStyle = "horizontal sep = 2cm, vertical sep = 1.5cm")
 
				 for iS in 1:11
						if SampleTrue[iS] == 1
						
						 # Plotting characteristic and unsaturated curves
							 Se = linspace(0.0001, 1., N_Plots) # Range of Se
							 for i in 1:1000
									θ = wrc.se.Se_2_θ(Se[i], θs[iS], θr[iS])
							 
									H_Kg[i] = 10.0 * wrc.kg.Se_2_H(Se[i], Hkg[iS], σ[iS])
									H_Vg[i] = 10.0 *wrc.vg.Se_2_H(Se[i], Hvg[iS], N[iS], Km[iS])
									
									Kunsat_Kg[i] = 360. * kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ[iS], Ks[iS], θs[iS], σ[iS])
									Kunsat_Vg[i] = 360. * kunsat.vg.KUNSAT(Se[i], N[iS], Ks[iS], Km[iS])
							 end

							 # Boundaries for plotting
							 Se_Min = 0.0
							 Se_Max = 1.01
							 H_Min = 0.1
							 H_Max = 10^8 # Equivalent to 2000kpa, we are working in cm
							 Kmin = 360. * 10. ^ -5.
							 Kmax =  360. * Ks[iS] 

							 Title = string("\$\\sigma =", string(round.(σ[iS],2)), ", " ,Name[iS], " \$")

							 if iS <= 10
									push!(Plot1b, Axis([
										 Plots.Linear(Se, H_Kg, mark="none", style="red, very thick"),
										 Plots.Linear(Se, H_Vg, mark="none", style="dashed, blue, very thick"),
									], ylabel=L"$h \ [cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=H_Min, ymax=H_Max, ymode="log", style="width=10.4cm, height=6.4cm"))  
					
									push!(Plot1b, Axis([
										 Plots.Linear(Se, Kunsat_Kg, mark="none", style="red, very thick"),
										 Plots.Linear(Se, Kunsat_Vg, mark="none", style="dashed, blue, very thick"),
										 Plots.Node(eval(Title),0.1,  Kmax - Kmax / 10., style="right"),
									], ylabel=L"$K(S_e) \ [cm \ h^{-1}]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin, ymax=Kmax, style="width=10.4cm, height=6.4cm"))
							 else
									push!(Plot1b, Axis([
										 Plots.Linear(Se, H_Kg, mark="none", style="red, very thick"),
										 Plots.Linear(Se, H_Vg, mark="none", style="dashed, blue, very thick"),
									], xlabel=L"$S_e \ [-]$", ylabel=L"$h \ [cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=H_Min, ymax=H_Max, ymode="log",  style="width=10.4cm, height=6.4cm"))  
		
									push!(Plot1b, Axis([
										 Plots.Linear(Se, Kunsat_Kg, mark="none", style="red, very thick", legendentry=L"$KG$"),
										 Plots.Linear(Se, Kunsat_Vg, mark="none", style="dashed, blue, very thick", legendentry=L"$VG$"),
										 Plots.Node(eval(Title),0.1, Kmax - Kmax / 10., style="right"),
									], xlabel=L"$S_e \ [-]$", ylabel=L"$K(S_e) \ [cm \ h^{-1}]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin, ymax=Kmax, legendStyle = "{at={(-0.3,-0.4)}, anchor=south west, legend columns=-1}",  style="width=10.4cm, height=6.4cm"))
							 end
						end
				 end
				 save(Path, Plot1b)
		 end

	 end #KG
end # MODULE

plots.kg.σMODEL_PLOT()

