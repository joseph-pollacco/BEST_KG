cd("C:\\JOE\\Main\\MODELS\\SOIL\\BEST")
DIR_Working = pwd() # Saving the current working directory
println( "Working path : $DIR_Working")
push!(LOAD_PATH, DIR_Working) # All the subroutines are in the current working directory

# include("Path.jl")
# include("Package.jl")
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
# include("TestBest.jl")
# include("Err.jl")

module kunsatmodel
# For scientific purpose
    using PGFPlots
    using SpecialFunctions, QuadGK, cst, kunsat, wrc, convertHydroModel, best, Sorptivity, path, reading, stats
    export θh_2_KunsatModel, RUN_KUNSAT_MODEL, RUN_CONVERT_KG_2_VG_MODEL, RUN_CONVERT_VG_2_KG_MODEL
    Ks_Correction= 15.14328956

    function  READ_HYDRO()
        DIR_Working = pwd()
        Path_Kg = string(DIR_Working, path.Read_Hydraulic)
        Data = Array(readdlm(Path_Kg, ',', header=true, skipstart=0)[1])
        ID = Data[:,1]
        Se =  Data[:,2]
        θs =  Data[:,3]
        θr = Data[:,4]
        σ = Data[:,5]
        Hkg = Data[:,6]

        σ_Mac = σ
        Hm_Mac = Hkg
        θs_Mac = θs

        N_Data = length(ID)

        return N_Data, θs, θr, σ, Hkg, σ_Mac, Hm_Mac, θs_Mac, Se
    end


    function θh_2_KunsatModel(Se_Ini, Se_Final, θs, θs_Mac, θr, Hm_Mac, σ_Mac, Hkg, σ)
        # UNIMODAL top soils
		T1 = 5.859
        T2 = 0.967
        T3 = 0.530
        
        #UNIMODAL bottom soils		
        # T1 = 6.484
		# T2 = 0.854
		# T3 = 0.316
		T1_Mac = T1
		T2_Mac = T2
		T3_Mac = T3

    	T1_Mac_Transf = 10. ^ -T1_Mac
    	T2_Mac_Transf = 2. * (1. - T2_Mac)
    	T3_Mac_Transf = 1. / (1. - T3_Mac)

    	T1_Transf =  10.^ -T1
    	T2_Transf = (2. * (1. - T2))
    	T3_Transf = (1. / (1. - T3))

        # Model Kunsat_Mac
        function KUNSAT_MAC(Se, θs, θs_Mac, Hm_Mac, σ_Mac, T1_Mac_Transf, T2_Mac_Transf, T3_Mac_Transf)
            Kunsat_Mac =   ((cst.Y / Hm_Mac) / (exp(erfcinv(2. * Se) * σ_Mac * sqrt(2.)))) ^ T2_Mac_Transf
            return Kunsat_Mac
        end

        # Model Kunsat_Mat
        function KUNSAT_MAT(Se, θs_Mac, θr, Hkg, σ, T1_Transf, T2_Transf, T3_Transf)
        	Kunsat_Mat =  ((cst.Y / Hkg) / (exp(erfcinv(2. * Se) * σ * sqrt(2.)))) ^ T2_Transf
        end

        Kunsat_Mac =  cst.Kconst * T1_Mac_Transf * (max(θs - θs_Mac, 0.) ^ T3_Mac_Transf) * quadgk( Se->  KUNSAT_MAC(Se, θs, θs_Mac, Hm_Mac, σ_Mac, T1_Mac_Transf, T2_Mac_Transf, T3_Mac_Transf), Se_Ini , Se_Final, abstol=10^-8.)[1]

        Kunsat_Mat =  cst.Kconst * T1_Transf * ((θs_Mac - θr) ^ T3_Transf) * quadgk( Se-> KUNSAT_MAT(Se, θs_Mac, θr, Hkg, σ, T1_Transf, T2_Transf, T3_Transf), Se_Ini, Se_Final, abstol=10^-8.)[1]

    	Kunsat = Kunsat_Mac + Kunsat_Mat

    	return Kunsat
    end


    function RUN_KUNSAT_MODEL()
        N_Data, θs, θr, σ, Hkg, σ_Mac, Hm_Mac, θs_Mac, Se = READ_HYDRO()

        H_Vg = Array{Float64}(10000)
        H_Kg = Array{Float64}(10000)
        Kunsat_Kg = Array{Float64}(10000)
        Kunsat_Vg = Array{Float64}(10000)
        Ks = Array{Float64}(N_Data)
        Hvg = Array{Float64}(N_Data)
        N = Array{Float64}(N_Data)
        Km = Array{Float64}(N_Data)

        for iS in 1:45
            Km[iS] = 1
            Ks[iS] =   Ks_Correction * θh_2_KunsatModel(0., 1., θs[ iS], θs_Mac[iS], θr[ iS], Hm_Mac[ iS], σ_Mac[ iS], Hkg[ iS], σ[ iS])
            Hvg[iS], N[iS], NSE, NSE_Hse, NSE_Kunsat =convertHydroModel.KOSUGI_2_VANGENUCHTEN(Hkg[iS], σ[iS])

            Kr_θini = kunsat.kg.KUNSAT(Se[iS], θs[iS], θr[iS], σ[iS], 1., θs[iS], σ[iS])
            B = best.B(Kr_θini)
            θ_Ini = wrc.se.Se_2_θ(Se[iS],θs[iS], θr[iS])
            Sorptivity_Kg = Sorptivity.kg.Sorptivity(θ_Ini[iS], θs[iS], θs[iS], θ_Ini, θs[iS], θr[iS], Hkg[iS], σ[iS], Ks[iS])
            Sorptivity_Vg = Sorptivity.vg.Sorptivity(θ_Ini[iS], θs[iS], θs[iS], θr[iS], Hvg[iS], N[iS], Ks[iS], Km[iS])
            TimeSteadyTrans = best.TIME_STEADYTRANSIT(Sorptivity_Kg, B, Ks[iS])

            println("Hvg= ,",Hvg[iS], " ,N, ", N[iS],  " ,Ks, ", Ks[iS], " ,Tmax, ", TimeSteadyTrans, " ,NSE, ", NSE," ,NSE_Hse, ", NSE_Hse," ,NSE_Kunsat, ", NSE_Kunsat, " ,S_Kg, ", Sorptivity_Kg ,  " ,S_Vg, ", Sorptivity_Vg)

        Plotting = true
        if Plotting
            # Plotting characteristic and unsaturated curves
                Se = linspace(0.001, 1., 10000) # Range of Se
                for i in 1:10000
                    θ = wrc.se.Se_2_θ(Se[i], θs[iS], θr[iS])
                
                    H_Kg[i] = 10.* wrc.kg.Se_2_H(Se[i], Hkg[iS], σ[iS])
                    H_Vg[i] = 10.*wrc.vg.Se_2_H(Se[i], Hvg[iS], N[iS], Km[iS])
                    
                    Kunsat_Kg[i] = 360.*kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ[iS], Ks[iS], θs[iS], σ[iS])
                    Kunsat_Vg[i] = 360.*kunsat.vg.KUNSAT(Se[i], N[iS], Ks[iS], Km[iS])
                end

                # Boundaries for plotting
                Se_Min = 0
                Se_Max = 1.01
                H_Min = 0.1
                H_Max = 10^8 # Equivalent to 2000kpa, we are working in cm
                Kmin = 0.

                DIR_Working = pwd()
                Path = string(DIR_Working, "//OUTPUT//Figure//KosVang//KosVang_", string(iS), ".pdf") #~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # println(Path)
                if isfile(Path) # If the file is there than delete it before saving figure
                    rm(Path)
                end

                pushPGFPlotsOptions("scale=1.5")
                Plot1b = GroupPlot(2, 2, groupStyle = "horizontal sep = 1.8cm, vertical sep = 2.cm")
                # Kosugi

                Title = string("Se(H) Sigma = ", σ[iS])
                push!(Plot1b, Axis([
                    Plots.Linear(Se, H_Kg, mark="none", style="red, very thick"),
                    Plots.Linear(Se, H_Vg, mark="none", style="dashed, blue, very thick"),
                ], title= Title, xlabel=L"$Se[-]$", ylabel=L"$H[cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=H_Min, ymax=H_Max, ymode="log")
                )  

                push!(Plot1b, Axis([
                    Plots.Linear(Se, Kunsat_Kg, mark="none", style="red, very thick", legendentry=L"$KG$"),
                    Plots.Linear(Se, Kunsat_Vg, mark="none", style="dashed, blue, very thick", legendentry=L"$VG$"),
                ], title="K(Se)", xlabel=L"$Se[-]$", ylabel=L"$K(Se)[cm/h]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin,  ymode="log", legendStyle = "{at={(-0.3,-0.4)}, anchor=south west, legend columns=-1}")
                )
                save(Path, Plot1b)
            end
        end  # Plotting true
    end



    function RUN_CONVERT_KG_2_VG_MODEL()
        println("START to read hydraulic params...")
        DIR_Working = pwd()
        Path_Hydraulic = string(DIR_Working, path.Read_Hydraulic)
        Se_Ini, θs, θr, σ, Hkg, Ks, ~, ~, Km, RingRadius = reading.HYDRAULIC(Path_Hydraulic)
        println("END to read hydraulic params \n")

        σ_Mac = σ
        Hm_Mac = Hkg
        θs_Mac = θs
        Se = Se_Ini
        
        N_Data = length(θs)

        H_Vg = Array{Float64}(10000)
        H_Kg = Array{Float64}(10000)
        Kunsat_Kg = Array{Float64}(10000)
        Kunsat_Vg = Array{Float64}(10000)

        Hvg = Array{Float64}(N_Data)
        N = Array{Float64}(N_Data)
        # Km = Array{Float64}(N_Data)

		  pushPGFPlotsOptions("scale=1.5")
        for iS in 1:45
       
            Hvg[iS], N[iS], NSE, NSE_Hse, NSE_Kunsat =convertHydroModel.KOSUGI_2_VANGENUCHTEN(Hkg[iS], σ[iS])

            Kr_θini = kunsat.kg.KUNSAT(Se[iS], θs[iS], θr[iS], σ[iS], 1., θs[iS], σ[iS])
            B = best.B(Kr_θini)
            θ_Ini = wrc.se.Se_2_θ(Se[iS],θs[iS], θr[iS])
            Sorptivity_Kg = Sorptivity.kg.Sorptivity(θ_Ini, θs[iS], θs[iS], θr[iS], Hkg[iS], σ[iS], Ks[iS])
            Sorptivity_Vg = sorptivity.vg.SORPTIVITY(θ_Ini, θs[iS], θs[iS], θr[iS], Hvg[iS], N[iS], Ks[iS], Km[iS])
            TimeSteadyTrans = best.TIME_STEADYTRANSIT(Sorptivity_Kg, B, Ks[iS])

            println(iS, " , ", "Hvg= ,",Hvg[iS], " ,N, ", N[iS],  " ,Ks, ", Ks[iS], " ,Tmax, ", TimeSteadyTrans, " ,NSE, ", NSE," ,NSE_Hse, ", NSE_Hse," ,NSE_Kunsat, ", NSE_Kunsat, " ,S_Kg, ", Sorptivity_Kg ,  " ,S_Vg, ", Sorptivity_Vg)

        Plotting = true
        if Plotting
            # Plotting characteristic and unsaturated curves
                Se = linspace(0.001, 1., 10000) # Range of Se
                for i in 1:10000
                    θ = wrc.se.Se_2_θ(Se[i], θs[iS], θr[iS])
                
                    H_Kg[i] = 10.* wrc.kg.Se_2_H(Se[i], Hkg[iS], σ[iS])
                    H_Vg[i] = 10.*wrc.vg.Se_2_H(Se[i], Hvg[iS], N[iS], Km[iS])
                    
                    Kunsat_Kg[i] = 360.*kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ[iS], Ks[iS], θs[iS], σ[iS])
                    Kunsat_Vg[i] = 360.*kunsat.vg.KUNSAT(Se[i], N[iS], Ks[iS], Km[iS])
                end

                # Boundaries for plotting
                Se_Min = 0
                Se_Max = 1.01
                H_Min = 0.1
                H_Max = 10^8 # Equivalent to 2000kpa, we are working in cm
                Kmin = 0.

                DIR_Working = pwd()
                Path = string(DIR_Working, "//OUTPUT//Figure//KosVang//KosVang_", string(iS), ".svg") #~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # println(Path)
                if isfile(Path) # If the file is there than delete it before saving figure
                    rm(Path)
                end

                Plot1b = GroupPlot(2, 2, groupStyle = "horizontal sep = 1.8cm, vertical sep = 2.cm")
                # Kosugi

                Title = string("Se(H) Sigma = ", σ[iS])
                push!(Plot1b, Axis([
                    Plots.Linear(Se, H_Kg, mark="none", style="red, very thick"),
                    Plots.Linear(Se, H_Vg, mark="none", style="dashed, blue, very thick"),
                ], title= Title, xlabel=L"$Se[-]$", ylabel=L"$H[cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=H_Min, ymax=H_Max, ymode="log")
                )  

                push!(Plot1b, Axis([
                    Plots.Linear(Se, Kunsat_Kg, mark="none", style="red, very thick", legendentry=L"$KG$"),
                    Plots.Linear(Se, Kunsat_Vg, mark="none", style="dashed, blue, very thick", legendentry=L"$VG$"),
                ], title="K(Se)", xlabel=L"$Se[-]$", ylabel=L"$K(Se)[cm/h]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin, legendStyle = "{at={(-0.3,-0.4)}, anchor=south west, legend columns=-1}")
                )
                save(Path, Plot1b)
            end
        end  # Plotting true
    end # RUN_CONVERT_KG_2_VG_MODEL



    function RUN_CONVERT_VG_2_KG_MODEL()
        println("START to read hydraulic params...")
        DIR_Working = pwd()
        Path_Hydraulic = string(DIR_Working, path.Read_Hydraulic)
        Se_Ini, θs, θr, ~, ~, Ks, Hvg, N, Km, RingRadius = reading.HYDRAULIC(Path_Hydraulic)
        println("END to read hydraulic params \n")
        
        N_Data = length(θs)
        H_Vg = Array{Float64}(1000)
        H_Kg = Array{Float64}(1000)
        Kunsat_Kg = Array{Float64}(1000)
        Kunsat_Vg = Array{Float64}(1000)
        Hkg = Array{Float64}(N_Data)
        σ = Array{Float64}(N_Data)
       
        for iS in 1:N_Data
       
            Hkg[iS], σ[iS], NSE, NSE_Hse, NSE_Kunsat = convertHydroModel.VANGENUCHTEN_2_KOSUGI(Hvg[iS], N[iS], Km[iS])

            println(iS, " , ", "Hkg= ,",Hkg[iS], " ,σ, ", σ[iS],  " ,Ks, ", Ks[iS], " ,NSE, ", NSE," ,NSE_Hse, ", NSE_Hse," ,NSE_Kunsat, ", NSE_Kunsat)

        Plotting = true
        if Plotting
            # Plotting characteristic and unsaturated curves
                Se = linspace(0.0001, 1., 1000) # Range of Se
                for i in 1:1000
                    θ = wrc.se.Se_2_θ(Se[i], θs[iS], θr[iS])
                
                    H_Kg[i] = 10.* wrc.kg.Se_2_H(Se[i], Hkg[iS], σ[iS])
                    H_Vg[i] = 10.*wrc.vg.Se_2_H(Se[i], Hvg[iS], N[iS], Km[iS])
                    
                    Kunsat_Kg[i] = 360.*kunsat.kg.KUNSAT(Se[i], θs[iS], θr[iS], σ[iS], Ks[iS], θs[iS], σ[iS])
                    Kunsat_Vg[i] = 360.*kunsat.vg.KUNSAT(Se[i], N[iS], Ks[iS], Km[iS])
                end

                # Boundaries for plotting
                Se_Min = 0
                Se_Max = 1.01
                H_Min = 0.1
                H_Max = 10^8 # Equivalent to 2000kpa, we are working in cm
					 Kmin = 360. * 10. ^ -5.
					 Kmax =  360. * Ks[iS] 

                DIR_Working = pwd()
                Path = string(DIR_Working, "//OUTPUT//Figure//KosVang//KosVang_", string(iS), ".svg") #~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # println(Path)
                if isfile(Path) # If the file is there than delete it before saving figure
                    rm(Path)
                end

                Plot1b = GroupPlot(2, 2, groupStyle = "horizontal sep = 1.8cm, vertical sep = 2.cm")
                # Kosugi

                Title = string("Se(H) Sigma = ", σ[iS])
                push!(Plot1b, Axis([
                    Plots.Linear(Se, H_Kg, mark="none", style="red, very thick"),
                    Plots.Linear(Se, H_Vg, mark="none", style="dashed, blue, very thick"),
                ], title= Title, xlabel=L"$S_e[-]$", ylabel=L"$H[cm]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=H_Min, ymax=H_Max, ymode="log")
                )  

                push!(Plot1b, Axis([
                    Plots.Linear(Se, Kunsat_Kg, mark="none", style="red, very thick", legendentry=L"$KG$"),
                    Plots.Linear(Se, Kunsat_Vg, mark="none", style="dashed, blue, very thick", legendentry=L"$VG$"),
                ], title="K(Se)", xlabel=L"$S_e[-]$", ylabel=L"$K(S_e)[cm/h]$", style="smooth", xmin=Se_Min, xmax=Se_Max, ymin=Kmin, ymax=Kmax, ymode="log", legendStyle = "{at={(-0.3,-0.4)}, anchor=south west, legend columns=-1}")
                )
                save(Path, Plot1b)
            end
        end  # Plotting true
    end # RUN_CONVERT_VG_2_KG_MODEL
end

#  kunsatmodel.RUN_CONVERT_KG_2_VG_MODEL()
 # or
#  kunsatmodel.RUN_CONVERT_VG_2_KG_MODEL()
