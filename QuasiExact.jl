module quasiExact # quasi-exact objective function
    using cst, BlackBoxOptim, param, wrc, kunsat, Optim, sorptivity
    export INFILTRATION3D_2_INFILTRATION1D, OF_INF_2_HYDRO

    #= COMPUTE NORMALISED TIME: Tη =#
    function TIME_2_TIMEη(iT, Sorpt, ΔK)
        return Time_η = iT * 2. * (ΔK/ Sorpt)^2.
    end



    # TRANSFORMS INFILTRATION_3D TO INFILTRATION_1D
    function INFILTRATION3D_2_INFILTRATION1D(iT, Inf_3D_Obs, Sorpt, RingRadius, θs, θr)
		Δθ = θs - θr
        return Inf_1D = max(Inf_3D_Obs - (iT * cst.γ * Sorpt^2.) / (RingRadius * Δθ), cst.ϵ)
    end



    # TRANSFORMS NORMALIZED INFILTRATION TO INFILTRATION-3D
    function INFILTTRATIONη_2_INFILTRATION3D(iT, Time_η, Inf_η, Sorpt, ΔK,  K_θini, Δθ, RingRadius)
        # Function compute infiltration-1d from normalized infiltration
        function INFILTTRATIONη_2_INFILTRATION1D(iT, Inf_η, Sorpt, ΔK,  K_θini)
            return INFILTRATION_1D = K_θini*iT + Inf_η * (Sorpt^2.) / (2.*ΔK)
        end
        ΔI_η = Time_η * cst.γ
       Inf_3D_Obs = INFILTTRATIONη_2_INFILTRATION1D(iT, Inf_η, Sorpt, ΔK, K_θini) + ΔI_η * (Sorpt^4.) / (2. *RingRadius*Δθ*(ΔK^2.))
        return Inf_3D_Obs
    end



    # TRANSFORMS INFILTRATION-3D TO NORMALIZED INFILTRATION
    function INFILTRATION3D_2_INFILTTRATIONη(iT,Inf_3D_Obs, Sorpt, ΔK, K_θini, Δθ, RingRadius)
          Inf_η = max((2.*ΔK / Sorpt^2.) * (Inf_3D_Obs  - K_θini*iT -  cst.γ * TIME_2_TIMEη(iT, Sorpt, ΔK) * (Sorpt^4.) /  (RingRadius * Δθ * 2. * (ΔK^2.)) ), cst.ϵ)
         return Inf_η
    end



    # NORMALISED QUASI-EXACT 3D CUMULATIVE INFILTRATION EQUATION
    function OF_QUASIEXACTη(Time_η, Inf_η)
		Left_Term = Time_η
		Right_Term =  (1./(1.-cst.β)) * (Inf_η - log((exp(cst.β*Inf_η) + cst.β-1.) / cst.β))
		if  Right_Term  < 0.
            # OF = 100000.*(abs(exp(cst.β*Inf_η) + cst.β-1.))
            OF = 10000.* exp(Inf_η)
		else
            OF =  abs(Left_Term - Right_Term)
        end
		return OF
    end



    function OF_INF_2_HYDRO(iT_N, iT_TransStead, Time, Inf_3D_Obs, Sorpt, ΔK, K_θini, Δθ, RingRadius)
        Left_Term = zeros(iT_N)
        Right_Term = zeros(iT_N)

        Of_Penalty = 0.
        Of_Stead = 0.
        Of_Trans = 0.
        W = 0.2
        for iT in 2:iT_N
            Time_η = quasiExact.TIME_2_TIMEη(Time[iT], Sorpt, ΔK)
            Inf_η = quasiExact.INFILTRATION3D_2_INFILTTRATIONη(Time[iT], Inf_3D_Obs[iT], Sorpt, ΔK, K_θini, Δθ, RingRadius)

            Left_Term[iT] = Time_η
            Right_Term[iT] =  (1./(1.-cst.β)) * (Inf_η - log((exp(cst.β*Inf_η) + cst.β-1.) / cst.β))
            
            if Right_Term[iT] > 0.
                if iT <= iT_TransStead
                    Of_Trans = Of_Trans + ((Left_Term[iT]) - (Right_Term[iT]))^2.
                else
                    Of_Stead = Of_Stead + (log10(Left_Term[iT]) - log10(Right_Term[iT]))^2.
                end
            else
                Of_Penalty = Of_Penalty + 1000. * exp(Inf_η)
                Right_Term[iT] = 0.
            end
        end
        Of = W * Of_Trans / Float64(iT_TransStead-1)  + (1.-W) * Of_Stead / Float64(iT_N - iT_TransStead + 1) + Of_Penalty

        return Of
    end

  
    # COMPUTE INFILTRATION_3D FROM OPTIMIZED HYDRAULIC PARAMETERS
    # Solving quasiexact solution
    function HYDRO_2_INFILTRATION3D(Time, Sorpt, Ks, K_θini, θs, θ_Ini, iT_N, RingRadius, iT_TransStead)
        Inf_3D_Sim= Array{Float64}(iT_N) # preparing the matrix
        Inf_η= Array{Float64}(iT_N) # preparing the matrix
		Δθ = θs - θ_Ini
		ΔK = Ks - K_θini
       
        # At t=1
        Inf_3D_Sim[1] = 0.
        Inf_η[1] = 0.
        Inf_η_Min = 0.0001
        Inf_η_Max = 0.05 * Time[2] #Since Time[1] = 0

        for iT in 2:iT_N # Looping for every time step
            Time_η= TIME_2_TIMEη(Time[iT], Sorpt, ΔK)

            # We are solving for Inf_η
            Optimization =  Optim.optimize(Inf_η ->  OF_QUASIEXACTη(Time_η, Inf_η), Inf_η_Min, Inf_η_Max, GoldenSection())
            Inf_η[iT] = Optim.minimizer(Optimization)[1]

            # Deriving the new bounds such that infiltration increases with time & the slope decreases with time
            Inf_η_Min = Inf_η[iT] + 0.0001

            # Maximum infiltration rate for T+1: (Inf[T2] - Inf[T1]) / (T2 - T1) which is 1 seconds
            if iT <= iT_N - 1
                Inf_η_Max = Inf_η[iT] + (Time[iT+1]- Time[iT]) * (Inf_η[iT] - Inf_η[iT-1]) / (Time[iT] - Time[iT-1])
            else
                Inf_η_Max = Inf_η[iT] + (Inf_η[iT] - Inf_η[iT-1])
            end

            # Transforming INFILTRATION3D 2 INFILTTRATIONη
            Inf_3D_Sim[iT] =  INFILTTRATIONη_2_INFILTRATION3D(Time[iT], Time_η, Inf_η[iT], Sorpt, ΔK, K_θini, Δθ, RingRadius)
        end
        return Inf_3D_Sim
    end



    #= =============== KOSUGI =============== =#
    module kg
        using cst, BlackBoxOptim, param, wrc, kunsat, Optim, sorptivity, quasiExact, stats, thetaTime, array, relationship
        export INFILTRATION3D_2_HYDRO, INFILTRATION3D_2_HYDRO_SORPTIVITY, INFILTRATION3D_2_HYDRO_σMOD

        # OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
        function INFILTRATION3D_2_HYDRO_σMOD(Time, Inf_3D_Obs, iT_N, θs, θ_Ini, RingRadius, iT_TransStead, Time_TransStead, Ks, Sorpt)
            θr = 0.
            Δθ = θs - θ_Ini
            Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

            # Objective function which matches observed with simulated infiltration
            function OF_Fast_η(Time, Inf_3D_Obs, Δθ, iT_N, RingRadius, θr, θ_Ini, σ, iT_TransStead)

                Hkg_σ = relationship.σ_2_Hkg(σ)

                Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorpt, θ_Ini, θs, θr, σ, Ks)

                K_θini = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)

                ΔK = Ks - K_θini

                Of_Hkg = abs(log10(Hkg_Sorpt) - log10(Hkg_σ))

                Wof = quasiExact.OF_INF_2_HYDRO(iT_N, iT_TransStead, Time[1:iT_N], Inf_3D_Obs[1:iT_N], Sorpt, ΔK, K_θini, Δθ, RingRadius) + Of_Hkg / 10.

                return Wof
            end

            # OPTIMIZATION
        
                # Optimization = BlackBoxOptim.bboptimize(σ ->  OF_Fast_η(Time[1:iT_N], Inf_3D_Obs[1:iT_N ], Δθ, iT_N , RingRadius, θr, θ_Ini, σ[1], iT_TransStead) ; SearchRange =[ (param.σ_Min, param.σ_Max)], NumDimensions=1, TraceMode=:silent)

                Optimization =  Optim.optimize(σ ->  OF_Fast_η(Time[1:iT_N], Inf_3D_Obs[1:iT_N ], Δθ, iT_N , RingRadius, θr, θ_Ini, σ, iT_TransStead), param.σ_Min, param.σ_Max, GoldenSection() )

                # Values of the optimal hydraulic params
                # σ = BlackBoxOptim.best_candidate(Optimization)[1]
                σ = Optim.minimizer(Optimization)[1]
            
          

            Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorpt, θ_Ini, θs, θr, σ, Ks)
           
            return σ, Hkg_Sorpt
        end


            # OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
            function INFILTRATION3D_2_HYDRO(Time, Inf_3D_Obs, iT_N, θs, θ_Ini, RingRadius, iT_TransStead, Time_TransStead, Option_Opt_σ, σ=1.)
                θr = 0.
                Δθ = θs - θ_Ini
    
                # Objective function which matches observed with simulated infiltration
                function OF_Fast_η(Time, Inf_3D_Obs, Δθ, iT_N, RingRadius, θr, θ_Ini, Hkg, Ks, σ, iT_TransStead)
    
                    Sorpt = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg, σ, Ks)
    
                    Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)
    
                    K_θini = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)
    
                    ΔK = Ks - K_θini
    
                    Wof = quasiExact.OF_INF_2_HYDRO(iT_N, iT_TransStead, Time[1:iT_N], Inf_3D_Obs[1:iT_N], Sorpt, ΔK, K_θini, Δθ, RingRadius)
    
                    Of_Sorpt = 0.
                    # Of_Sorpt = abs(log10(Hkg) - log10(relationship.σ_2_Hkg(σ)))
    
                    return Wof =  Wof +  Of_Sorpt
                end
    
                # OPTIMIZATION
                if Option_Opt_σ # If Ks is not known
                    Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:iT_N], Inf_3D_Obs[1:iT_N ], Δθ, iT_N , RingRadius, θr, θ_Ini, 10.^Param[1], 10.^Param[2], Param[3], iT_TransStead) ; SearchRange =[ (log10(param.Hkg_Min), log10(param.Hkg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max), (param.σ_Min, param.σ_Max)], NumDimensions=3, TraceMode=:silent)
    
                    # Values of the optimal hydraulic params
                    Hkg = 10.^(BlackBoxOptim.best_candidate(Optimization)[1])
                    Ks = 10.^(BlackBoxOptim.best_candidate(Optimization)[2])
                    σ = BlackBoxOptim.best_candidate(Optimization)[3]
                
                else # If σ is known
                    Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:iT_N], Inf_3D_Obs[1:iT_N], Δθ, iT_N, RingRadius, θr, θ_Ini, 10.^Param[1], 10.^Param[2], σ, iT_TransStead) ; SearchRange = [ (log10(param.Hkg_Min), log10(param.Hkg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=2, TraceMode=:silent)
                    # Values of the optimal hydraulic params
                    Hkg = 10.^ (BlackBoxOptim.best_candidate(Optimization)[1])
                    Ks = 10.^(BlackBoxOptim.best_candidate(Optimization)[2])
                end
    
                Sorpt = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg, σ, Ks)
    
                # Converting θ_Ini to Kr_θini
                Se_Ini= wrc.se.θ_2_Se(θ_Ini, θs, θr)
                K_θini= kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)
                Kr_θini = K_θini / Ks
                
                return Ks, Kr_θini, Sorpt, σ, Hkg
            end



        # # OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
        # function INFILTRATION3D_2_HYDRO_σmod(Sorpt, Ks, Time, Inf_3D_Obs, iT_N, θs, θ_Ini, RingRadius, iT_TransStead, Time_TransStead)
        #     θr = 0.
        #     Δθ = θs - θ_Ini
 
        #     Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

        #     # Objective function which matches observed with simulated infiltration
        #     function OF_Fast_η_σmod(Sorpt, Time, Inf_3D_Obs, Δθ, iT_N, RingRadius, θr, θ_Ini, Ks, σ, iT_TransStead, Se_Ini)
        #         K_θini = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)

        #         ΔK = Ks - K_θini

        #         Wof = quasiExact.OF_INF_2_HYDRO(iT_N, iT_TransStead, Time[1:iT_N], Inf_3D_Obs[1:iT_N], Sorpt, ΔK, K_θini, Δθ, RingRadius)

        #         Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorpt, θ_Ini, θs, θr, σ, Ks)

		# 		Hkg_σ = relationship.σ_2_Hkg(σ)
				
		# 		Of_Hkg = abs(log10(Hkg_Sorpt) - log10(Hkg_σ))

        #         return Wof =  Wof +  Of_Hkg 
        #     end

        #     # OPTIMIZATION
        
        #         # Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η_σmod(Sorpt, Time[1:iT_N], Inf_3D_Obs[1:iT_N ], Δθ, iT_N , RingRadius, θr, θ_Ini, Ks, Param[1], iT_TransStead, Se_Ini) ; SearchRange =[ (param.σ_Min, param.σ_Max)], NumDimensions=1, TraceMode=:silent)
        #         # σ = BlackBoxOptim.best_candidate(Optimization)[1]

        #         Optimization =  Optim.optimize(σ -> OF_Fast_η_σmod(Sorpt, Time[1:iT_N], Inf_3D_Obs[1:iT_N ], Δθ, iT_N , RingRadius, θr, θ_Ini, Ks, σ, iT_TransStead, Se_Ini), param.σ_Min, param.σ_Max, GoldenSection() )
        
		# 	σ = Optim.minimizer(Optimization)[1]

        #         Hkg = sorptivity.kg.SORPTIVITY_2_Hkg(Sorpt, θ_Ini, θs, θr, σ, Ks)
        #     return σ, Hkg
        # end




        function INFILTRATION3D_2_HYDRO_SORPTIVITY(Time, Inf_3D_Obs, Se_Ini, θ_Ini, θs, θr, Ks, RingRadius, Time_TransStead, iT_TransStead, iT_N, Sorptivity_1)

            θr = 0.

            Se_Ini_2 = 0.3
            Se_Ini_3 = 0.6
            Se_Ini_4 = 0.6

            θ_Ini_2 = wrc.se.Se_2_θ(Se_Ini_2, θs, θr)
            θ_Ini_3 = wrc.se.Se_2_θ(Se_Ini_3, θs, θr)
            θ_Ini_4 = wrc.se.Se_2_θ(Se_Ini_4, θs, θr)
            
            function OF_TIME(Time, Inf_3D_Obs, RingRadius, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N, Sorptivity_1, Se_Ini_2, θ_Ini_2,  Se_Ini_3 ,θ_Ini_3)

                Of_θ, θ_Time = thetaTime.kg.θTIME(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius, Sorptivity_1, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N)
                Wof = Of_θ

                return Wof, θ_Time
            end #OF_BEST
            
            function OF_θ_Ini(Time, Inf_3D_Obs, RingRadius, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N, Sorptivity_1, Se_Ini_2, θ_Ini_2,  Se_Ini_3 ,θ_Ini_3)

					 Of_θ, θ_Time = thetaTime.kg.θTIME(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius, Sorptivity_1, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N)
					 

                # iθ_Ini_2 = array.SEARCH_INDEX(θ_Time[1:iT_N], θ_Ini_2) # getting the closest index
                # iθ_Ini_3 = array.SEARCH_INDEX(θ_Time[1:iT_N], θ_Ini_3) # getting the closest index
                iθ_Ini_4 = array.SEARCH_INDEX(θ_Time[1:iT_N], θ_Ini_4) # getting the closest index

                Hkg = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity_1, θ_Ini, θs, θr, σ, Ks)

                # K_θini_2 = kunsat.kg.KUNSAT(Se_Ini_2, θs, θr, σ, Ks, θs, σ)
                # K_θini_3 = kunsat.kg.KUNSAT(Se_Ini_3, θs, θr, σ, Ks, θs, σ)
                K_θini_4 = kunsat.kg.KUNSAT(Se_Ini_4, θs, θr, σ, Ks, θs, σ)
                
                # Sorptivity_2 = sorptivity.kg.SORPTIVITY(θ_Ini_2, θs, θs, θr, Hkg, σ, Ks)
                # Sorptivity_3 = sorptivity.kg.SORPTIVITY(θ_Ini_3, θs, θs, θr, Hkg, σ, Ks)
                Sorptivity_4 = sorptivity.kg.SORPTIVITY(θ_Ini_4, θs, θs, θr, Hkg, σ, Ks)
    
                # ΔK_2 = Ks - K_θini_2
                # Δθ_2 = θs - θ_Ini_2
                # Wof_2 = quasiExact.OF_INF_2_HYDRO(iT_N-iθ_Ini_2+1, iT_TransStead-iθ_Ini_2+1, Time[iθ_Ini_2:iT_N]-Time[iθ_Ini_2-1], Inf_3D_Obs[iθ_Ini_2:iT_N]-Inf_3D_Obs[iθ_Ini_2-1], Sorptivity_2, ΔK_2, K_θini_2, Δθ_2, RingRadius)

                # ΔK_3 = Ks - K_θini_3
                # Δθ_3 = θs - θ_Ini_3
                # Wof_3 = quasiExact.OF_INF_2_HYDRO(iT_N-iθ_Ini_3+1, iT_TransStead-iθ_Ini_3+1, Time[iθ_Ini_3:iT_N]-Time[iθ_Ini_3-1], Inf_3D_Obs[iθ_Ini_3:iT_N]-Inf_3D_Obs[iθ_Ini_3-1], Sorptivity_3, ΔK_3, K_θini_3, Δθ_3, RingRadius)
                
                ΔK_4 = Ks - K_θini_4
                Δθ_4 = θs - θ_Ini_4
                Wof_4 = quasiExact.OF_INF_2_HYDRO(iT_N-iθ_Ini_4+1, iT_TransStead-iθ_Ini_4+1, Time[iθ_Ini_4:iT_N]-Time[iθ_Ini_4-1], Inf_3D_Obs[iθ_Ini_4:iT_N]-Inf_3D_Obs[iθ_Ini_4-1], Sorptivity_4, ΔK_4, K_θini_4, Δθ_4, RingRadius)

					 Wof = Wof_4 / (iT_TransStead-iθ_Ini_4+1)
					#  println("$Of_θ  $Wof_4")
					 return Wof
            end
                         


            Optimization = BlackBoxOptim.bboptimize(σ  ->  OF_TIME(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius, θ_Ini, θs, θr, σ[1], Ks, iT_TransStead, iT_N, Sorptivity_1, Se_Ini_2, θ_Ini_2, Se_Ini_3, θ_Ini_3)[1]; SearchRange = [(param.σ_Min, param.σ_Max)], NumDimensions=1, TraceMode=:silent)
            
            σ = BlackBoxOptim.best_candidate(Optimization)[1]
        
            Hkg = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity_1, θ_Ini, θs, θr, σ, Ks)

            Wof, θ_Time = OF_TIME(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N, Sorptivity_1, Se_Ini_2, θ_Ini_2, Se_Ini_3, θ_Ini_3)

            println( "Se_Ini_2=  ", Se_Ini_2, "  ,σ=" , σ, "  ,Hkg=", Hkg,  "\n ")
            println("OF= $Wof")

            return σ, Hkg, Ks, θ_Time
        end #BESTG_INVERSE_SORPTIVITY

    end


    #= =============== VAN GENUCHTEN =============== =#
    module vg
    using cst, BlackBoxOptim, param, wrc, kunsat, Optim, sorptivity, quasiExact
        export INFILTRATION3D_2_HYDRO

        # OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
        function INFILTRATION3D_2_HYDRO(Time, Inf_3D_Obs, iT_N, θs, θ_Ini, RingRadius, Km, iT_TransStead, Time_TransStead, Option_Opt_N, N=1.)
            θr = 0.
            Δθ = θs - θ_Ini
            
            # Objective function which matches observed with simulated infiltration
            function OF_Fast_η(Time, Inf_3D_Obs, Δθ, iT_N, RingRadius, θr, θ_Ini, Hvg, Ks, N, Km)
            
                Sorpt = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs,  θr, Hvg, N, Ks, Km)

                Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)
            
                K_θini = kunsat.vg.KUNSAT(Se_Ini, N, Ks, Km)
            
                ΔK = Ks - K_θini

                OF_Cumul = quasiExact.OF_INF_2_HYDRO(iT_N, iT_TransStead, Time[1:iT_N], Inf_3D_Obs[1:iT_N], Sorpt, ΔK, K_θini, Δθ, RingRadius)
                return OF_Cumul
            end

            # OPTIMIZATION
            # N_Min and N_Max depends on Km
            if Km == 1
                N_Min = param.N_Km1_Min
                N_Max = param.N_Km1_Max
            elseif  Km == 2
                N_Min = param.N_Km2_Min
                N_Max = param.N_Km2_Max
            end

            if Option_Opt_N # If Ks is not known
                Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:iT_N], Inf_3D_Obs[1:iT_N], Δθ, iT_N, RingRadius,  θr, θ_Ini, 10.^Param[1], 10.^Param[2], Param[3], Km) ; SearchRange =[ (log10(param.Hvg_Min), log10(param.Hvg_Max)), (log10(param.Ks_Min), log10(param.Ks_Max)), (N_Min, N_Max)], NumDimensions=3, TraceMode=:silent)
                # Values of the optimal hydraulic params
                Hvg = 10.^(BlackBoxOptim.best_candidate(Optimization)[1])
                Ks = 10.^(BlackBoxOptim.best_candidate(Optimization)[2])
                N = BlackBoxOptim.best_candidate(Optimization)[3]
            else # if N is known
                Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:iT_N], Inf_3D_Obs[1:iT_N], Δθ, iT_N, RingRadius,  θr, θ_Ini, 10.^Param[1], 10.^Param[2], N, Km) ; SearchRange = [ (log10(param.Hvg_Min), log10(param.Hvg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=2, TraceMode=:silent)
                # Values of the optimal hydraulic params
                Hvg = 10.^ (BlackBoxOptim.best_candidate(Optimization)[1])
                Ks = 10.^(BlackBoxOptim.best_candidate(Optimization)[2])      
            end

            Sorpt = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg, N, Ks, Km)

            Se_Ini= wrc.se.θ_2_Se(θ_Ini, θs, θr)
 
            Kr_θini=  kunsat.vg.KUNSAT(Se_Ini, N, 1., Km)
            
            return Ks, Kr_θini, Sorpt, N, Hvg
        end
    end

end