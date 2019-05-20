#= Compute time series θ(Time),kr_Θ with the following steps
1) Compute Infiltration_1D from Inf_3D_Obs
2) Compute  θ_Time
3) Compute from the quasi-exact solution
4) Inverse σ from the \
(Time, T_TransSteady)=#

module thetaTime
	using best, quasiExact, wrc, kunsat, quasiExact, BlackBoxOptim, sorptivity, param
	export TIME_INFILTRATION_SHIFT
	
    
    # # Shifting the time and Infiltration if numerically chaning θ_Ini
    # function TIME_INFILTRATION_SHIFT(iT_Start, iT_N, Time, Inf_3D_Obs)
    #     # This is when we are starting at different θ_Ini where we are numerically shifting the starting point       
    #     iT_N_Shift = iT_N - iT_Start + 1

    #     Time_Shift = Array{Float64}( iT_N_Shift)
    #     Inf_3D_Obs_Shift = Array{Float64}( iT_N_Shift)

    #     if iT_Start >= 2

    #         iT_Shift = 1
    #         for iT in iT_Start:iT_N
    #             Time_Shift[iT_Shift] = Time[iT] - Time[iT_Start-1] # at time 1 infiltration=/= 0 this is why iT_Start-1
    #             Inf_3D_Obs_Shift[iT_Shift] =  Inf_3D_Obs[iT] - Inf_3D_Obs[iT_Start-1]
    #             iT_Shift +=  1
    #         end         
    #     else
    #         Time_Shift = Time[:]
    #         Inf_3D_Obs_Shift =  Inf_3D_Obs[:]    
    #     end

    #     return Inf_3D_Obs_Shift, Time_Shift, iT_N_Shift
    # end


    #= =============== KOSUGI =============== =#
    module kg
        using best, quasiExact, wrc, kunsat, quasiExact, BlackBoxOptim, sorptivity, param, thetaTime, Optim
        export HYDRAU_VARY_INFILTRATION,  SELECT_Hvg_σ, θTIME

                #= =============== θTIME =============== =#
        # We are only interested in the transit state as in the steady steady state we assum that θ= θs and k = Ks
        function θTIME(Time, Inf_3D_Obs, RingRadius, Sorptivity, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N)
            # Dedicating memory
            Infiltration_1D = Array{Float64}(iT_N)
            ΔInfiltration_1D = Array{Float64}(iT_N)
            ΔTime = Array{Float64}(iT_N)
     
            #= Transforming Infiltration-3D to Infiltration-1D =#
            Infiltration_1D[1] = 0.
            ΔInfiltration_1D[1] = 0.
            ΔTime[1] = 0.
            for iT in 2:iT_N
				 Infiltration_1D[iT] = quasiExact.INFILTRATION3D_2_INFILTRATION1D(Time[iT], Inf_3D_Obs[iT], Sorptivity, RingRadius, θs, θr)
                
                # Compute the infiltration rate for every time step =#
                ΔInfiltration_1D[iT] = Infiltration_1D[iT] - Infiltration_1D[iT-1]

                ΔTime[iT] = Time[iT] - Time[iT-1]
            end

				Power = 1.

            # Compute θ_TIME
            function θ_TIME_OF(iT_TransStead, ΔInfiltration_1D, ΔTime, θ_Ini, θs, θr, σ, Ks, ZinfMax, iT_N, Power)
					θ_Time = Array{Float64}(iT_N)
					Se = Array{Float64}(iT_N)

					θ_Time[1]  = θ_Ini

					Se[1] = wrc.se.θ_2_Se(θ_Ini, θs, θr ) 

					Of_Penalty = 0.
					Of = 0.
                for iT in 2:iT_N
                    ΔBalance = ΔInfiltration_1D[iT] - ΔTime[iT] * kunsat.kg.KUNSAT(Se[iT-1], θs, θr, σ, Ks, θs, σ)
                    
						  θ_Time[iT]  = max(θ_Time[iT-1] + ΔBalance /  ZinfMax, 0.)

						  θ_Time[iT]  =  max(θ_Time[iT] ,  θ_Time[iT-1])

                    Se[iT] = wrc.se.θ_2_Se(min(θ_Time[iT], θs), θs, θr )
                    
                    if θ_Time[iT-1] > θ_Time[iT]
                        Of_Penalty += θ_Time[iT-1] - θ_Time[iT]
                    end 
						  
                    if  iT >= iT_TransStead
                        Of += abs(θ_Time[iT] - θs)
                    end
					end
					 
					Of_θs = abs(θ_Time[iT_TransStead] - θs)

					Wof = 10.0 *((10.0 * Of  + Of_Penalty) / (iT_N - iT_TransStead + 1))^(1./Power)
     
                return Of_θs, Wof, θ_Time
            end            
            
            Optimization =  Optim.optimize(ZinfMax ->  θ_TIME_OF(iT_TransStead, ΔInfiltration_1D[1:iT_N], ΔTime[1:iT_N], θ_Ini, θs, θr, σ, Ks, ZinfMax, iT_N, Power)[1], param.ZinfMax_Min, param.ZinfMax_Max, Optim.GoldenSection())
            
            ZinfMax = Optim.minimizer(Optimization)[1]
           
            Of_θs, Wof, θ_Time = θ_TIME_OF(iT_TransStead, ΔInfiltration_1D[1:iT_N], ΔTime[1:iT_N], θ_Ini, θs, θr, σ, Ks, ZinfMax, iT_N, Power)
           
            # Just in case
            # θ_Time = min.(θ_Time[1:iT_N] *  θs / maximum(θ_Time[1:iT_N]), θs)

            return θ_Time, Wof       
        end


     
    #= =============== HYDRAU_VARY_INFILTRATION =============== =# 
    # function HYDRAU_VARY_INFILTRATION(θ_Time, Ks, Time, Inf_3D_Obs, RingRadius, Sorptivity, θs, θr, iT_TransStead, iT_N)
    #     # Dedicating memory
    #     Kr_Inf = Array{Float64}(iT_TransStead)
    #     Kr_Inf[1:iT_TransStead] = 1. # When the soil is saturated it is equal to 1
    #     Sorptivity_Inf = Array{Float64}(iT_TransStead)
    #     σ_Inf =  Array{Float64}(iT_TransStead)
    #     Hkg_Inf =  Array{Float64}(iT_TransStead)
    
    #     θ_Time_Previous = θ_Time[1]
    #     # = Compute Kr at every θ_Time(iT), but to speed up we compute for every  ΔSe_Inf_Max =#
    #     for i_Inf in 1:iT_TransStead-1
    #         # To save time just computing for equaly intervalled Δθ_Time_Max
    #         if ((θ_Time[i_Inf] - θ_Time_Previous) / θs >= param.ΔSe_Inf_Max) && (param.Se_Inf_Max >= θ_Time[i_Inf] / θs) || ( i_Inf == 1)

    #             Infiltration_3D_Shift, Time_Shift, iT_N_Shift = thetaTime.TIME_INFILTRATION_SHIFT(i_Inf, iT_N, Time[1:iT_N], Inf_3D_Obs[1:iT_N])
            
    #             ~, Kr_Inf[i_Inf], Sorptivity_Inf[i_Inf], σ_Inf[i_Inf], Hkg_Inf[i_Inf] = quasiExact.kg.INFILTRATION3D_2_HYDRO(Time_Shift[1:iT_N_Shift], Infiltration_3D_Shift[1:iT_N_Shift], iT_N_Shift, θs, θ_Time[i_Inf], RingRadius, false, Ks)
            
    #             println( i_Inf," , ", "Se_Inf ,", θ_Time[i_Inf] / θs, ",  Kr_Inf ,", Kr_Inf[i_Inf], " , Sorptivity ,", Sorptivity_Inf[i_Inf], ",σ ," ,σ_Inf[i_Inf], " , Hkg , ", Hkg_Inf[i_Inf], "\n" )

    #             θ_Time_Previous = θ_Time[i_Inf] # Making a copy of θ_Time 
    #         else
    #             # Just copy the values of the previous timestep since there are too many time steps
    #             Kr_Inf[i_Inf] =  Kr_Inf[i_Inf-1]
    #             Sorptivity_Inf[i_Inf] =  Sorptivity_Inf[i_Inf-1]
    #             σ_Inf[i_Inf] =  σ_Inf[i_Inf-1]
    #             Hkg_Inf[i_Inf] = Hkg_Inf[i_Inf-1]
    #         end 
    #     end
    #     println("End HYDRAU_VARY_INFILTRATION /n")
    #     #  Sorptivity_Inf[iT_TransStead-1: iT_N] = Sorptivity_Inf[iT_TransStead-1]
    #     return Kr_Inf, Sorptivity_Inf, σ_Inf, Hkg_Inf
    # end

    #      #= =============== SELECTING BEST HYDRAULIC parameters =============== =#  
    # #= There are different options of Hvg and σ and we are looking for the best hydraulic param set which would give the best infiltration for any initilial soil moisture =#
    #     function SELECT_Hvg_σ(Time, Inf_3D_Obs, iT_N, iT_TransStead, θs, θr, σ_Inf, Hkg_Inf, θ_Time, Kr_Inf, Ks, Sorptivity, RingRadius)

    #         println("START of selecting the best hydraulic parameter sets")

    #         Of_Cum = Array{Float64}(iT_N)
    #         θ_Time_Previous = θ_Time[1]  
    #         for i_Inf in 1:iT_TransStead-1
    #             # To speed up
    #             if ((θ_Time[i_Inf] - θ_Time_Previous) / θs >= param.ΔSe_Inf_Max) && (param.Se_Inf_Max >= θ_Time[i_Inf] / θs) || ( i_Inf == 1)

    #                 # Infiltration_3D_Shift, Time_Shift, iT_N_Shift = thetaTime.TIME_INFILTRATION_SHIFT(i_Inf, iT_N, Time[1:iT_N], Inf_3D_Obs[1:iT_N])

    #                 Sorptivity = sorptivity.kg.SORPTIVITY(θ_Time[1], θs, θ_Time[1], θs, θr, Hkg_Inf[i_Inf], σ_Inf[i_Inf], Ks)

    #                 Infiltration_Sim = quasiExact.HYDRO_2_INFILTRATION3D(Time[1: iT_N], Sorptivity[1], Ks, Ks*Kr_Inf[1], θs, θ_Time[1],  iT_N, RingRadius)

    #                 Of_Cum[i_Inf] = 0. 
    #                 for iT in 1:iT_N
    #                     Of = abs(Infiltration_Sim[iT] - Inf_3D_Obs[iT])
    #                     Of_Cum[i_Inf] = Of_Cum[i_Inf] + Of
    #                     # println(Infiltration_Sim[iT], "," ,Inf_3D_Obs[iT])
    #                 end
    #                 θ_Time_Previous = θ_Time[i_Inf] # Making a copy of θ_Time 
    #                 println( i_Inf," , ", "Se_Inf: ",  θ_Time[i_Inf]/ θs,  " Of_Cum: ", Of_Cum[i_Inf], "\n")
    #             else
    #                 # Just provide a high value of OF for those that we skeeped
    #                 Of_Cum[i_Inf] = 1000000000000000000.
    #             end        
    #         end
    #         # The position i of the hydraulic param sets which gives the best predictions for all initial θ
    #         Of_Cum_Min = minimum(Of_Cum[1:iT_TransStead-1])
    #         iOF_Min = searchsortedfirst( Of_Cum[1:iT_TransStead-1], Of_Cum_Min)
    #         println("iOF_Min, $iOF_Min")

    #         Hkg_Qe = Hkg_Inf[iOF_Min]
    #         σ_Qe = σ_Inf[iOF_Min]

    #         println("σ_Best $σ_Qe Hkg_Best $Hkg_Qe")

    #         # Derive the optimal infiltration
    #         Sorptivity_Qe_Kg= sorptivity.kg.SORPTIVITY(θ_Time[1], θs, 0, Hkg_Qe, σ_Qe, Ks)
    #         Se_Ini = wrc.se.θ_2_Se(θs, θ_Time[1], 0)
	# 		K_θini = kunsat.kg.KUNSAT(Se_Ini, θs , 0., σ_Qe, Ks, θs, σ_Qe ) 
	# 		Infiltration_Qe_Kg = quasiExact.HYDRO_2_INFILTRATION3D(Time[1:iT_N], Sorptivity_Qe_Kg, Ks, K_θini, θs, θ_Time[1], iT_N, RingRadius)

    #         println("END of selecting the best hydraulic parameter sets")
        
    #         return Infiltration_Qe_Kg, σ_Qe, Hkg_Qe  
    #     end 

    end 
    

     #= =============== VAN GENUCHTEN =============== =#
    module vg
        using best, quasiExact, wrc, kunsat, quasiExact, BlackBoxOptim, sorptivity, param, thetaTime
        export HYDRAU_VARY_INFILTRATION, SELECT_Hkg_N, θTIME, θTIME_N


        function θTIME_N(Time, Inf_3D_Obs, RingRadius, Sorptivity, θ_Ini, θs, θr, Ks, Km, iT_TransStead)

            # Compute θ_INF
            function θ_INF(iT_TransStead, ΔInfiltration_1D, ΔTime, θ_Ini, θs, θr, N, Ks, Km, ZinfMax )
                θ_Time = Array{Float64}(iT_TransStead)
                Se = Array{Float64}(iT_TransStead)

                Of_Penalty = 0.
                θ_Time[1] = θ_Ini
                for iT in 2:iT_TransStead
                    Se[iT-1] = wrc.se.θ_2_Se(θ_Time[iT-1], θs, θr )
                    ΔBalance = ΔInfiltration_1D[iT] - ΔTime[iT] * kunsat.vg.KUNSAT(Se[iT-1], N, Ks, Km)
                    if ΔBalance < 0.
                        Of_Penalty = -ΔBalance + Of_Penalty
                        ΔBalance = 0.
                    end
                    θ_Time[iT]  = θ_Time[iT-1] +  ΔBalance /  ZinfMax
                end
                return θ_Time, Of_Penalty
            end
        
            # Searching for θ_Time = 0
            function OF_θInf(iT_TransStead, ΔInfiltration_1D, ΔTime, θ_Ini, θs, θr, N, Ks, Km, ZinfMax)
                θ_Time, Of_Penalty = θ_INF(iT_TransStead, ΔInfiltration_1D, ΔTime, θ_Ini, θs, θr, N, Ks, Km, ZinfMax)
                OF = 1000*abs(θ_Time[iT_TransStead] - θs) + 1000*Of_Penalty + 10.0 * (1. / N)
                return OF
            end
            
            # Dedicating memory
            Infiltration_1D = Array{Float64}(iT_TransStead)
            ΔInfiltration_1D = Array{Float64}(iT_TransStead)
            ΔTime = Array{Float64}(iT_TransStead)
            
            #= STEP A: Transforming Infiltration-3D to Infiltration-1D =#
            ΔInfiltration_1D[1] = Infiltration_1D[1] - 0.
            for iT in 1:iT_TransStead
                Infiltration_1D[iT] = quasiExact.INFILTRATION3D_2_INFILTRATION1D(Time[iT], Inf_3D_Obs[iT], Sorptivity, RingRadius, θs, θr)
                #=  STEP B: Compute the infiltration rate for every time step =#
                if iT>=2
                    ΔInfiltration_1D[iT] = Infiltration_1D[iT] - Infiltration_1D[iT-1]
                    ΔTime[iT] = Time[iT] - Time[iT-1]
                end
            end

            # N_Min and N_Max depends on Km
			if Km == 1
				N_Min = param.N_Km1_Min
				N_Max = param.N_Km1_Max
			elseif  Km == 2
				N_Min = param.N_Km2_Min
				N_Max = param.N_Km2_Max
			end

            Optimization = bboptimize( Param -> OF_θInf(iT_TransStead, ΔInfiltration_1D[1:iT_TransStead], ΔTime[1:iT_TransStead], θ_Ini, θs, θr, Param[1], Ks, Km, Param[2]) ; SearchRange = [(N_Min, N_Max), (0.0001, 100000.)], NumDimensions = 2, TraceMode = :silent)
            N = best_candidate(Optimization)[1]
            ZinfMax = best_candidate(Optimization)[2]

            # println( "N_Opt= $N ZinfMax= $ZinfMax \n")

            θ_Time, Of_Penalty  = θ_INF(iT_TransStead, ΔInfiltration_1D[1:iT_TransStead], ΔTime[1:iT_TransStead], θ_Ini, θs, θr, N, Ks, Km, ZinfMax )

            return N         
        end

		# We are only interested in the transit state as in the steady steady state we assum that θ= θs and k = Ks
		function θTIME(Time, Inf_3D_Obs, RingRadius, Sorptivity, θ_Ini, θs, θr, n, Ks, iT_TransStead, iT_N)
			# Dedicating memory
			Infiltration_1D = Array{Float64}(iT_N)
			ΔInfiltration_1D = Array{Float64}(iT_N)
			ΔTime = Array{Float64}(iT_N)

			#= Transforming Infiltration-3D to Infiltration-1D =#
			Infiltration_1D[1] = 0.
			ΔInfiltration_1D[1] = 0.
			ΔTime[1] = 0.
			for iT in 2:iT_N
			Infiltration_1D[iT] = quasiExact.INFILTRATION3D_2_INFILTRATION1D(Time[iT], Inf_3D_Obs[iT], Sorptivity, RingRadius, θs, θr)
				
				#=  STEP B: Compute the infiltration rate for every time step =#
				ΔInfiltration_1D[iT] = Infiltration_1D[iT] - Infiltration_1D[iT-1]

				ΔTime[iT] = Time[iT] - Time[iT-1]
			end

			if iT_TransStead < iT_N - 4
				Power = 1.
			else
				Power = 4.
			end

			# Compute θ_TIME
			function θ_TIME_OF(iT_TransStead, ΔInfiltration_1D, ΔTime, θ_Ini, θs, θr, n, Ks, ZinfMax, iT_N, Power)
				θ_Time = Array{Float64}(iT_N)
				Se = Array{Float64}(iT_N)

				θ_Time[1]  = θ_Ini

				Se[1] = wrc.se.θ_2_Se(θ_Ini, θs, θr ) 

				Of_Penalty = 0.
				Of = 0.
				for iT in 2:iT_N
					ΔBalance = ΔInfiltration_1D[iT] - ΔTime[iT] * kunsat.kg.KUNSAT(Se[iT-1], θs, θr, n, Ks, θs, n)
					
					θ_Time[iT]  = θ_Time[iT-1] + ΔBalance /  ZinfMax

					Se[iT] = wrc.se.θ_2_Se( max(min(θ_Time[iT], θs), θ_Time[iT-1]), θs, θr )
					
					if θ_Time[iT-1] > θ_Time[iT]
							Of_Penalty += θ_Time[iT-1] - θ_Time[iT]
					end 
					
				if  iT >= iT_TransStead
					Of += abs(θ_Time[iT] - θs)^1.
				end
				end

			Of = 10.0 *((Of  + 10.0 *Of_Penalty) / (iT_N - iT_TransStead + 1))^(1./ 1.)

				return Of, θ_Time
			end            
			
			Optimization =  Optim.optimize(ZinfMax ->  θ_TIME_OF(iT_TransStead, ΔInfiltration_1D[1:iT_N], ΔTime[1:iT_N], θ_Ini, θs, θr, n, Ks, ZinfMax, iT_N, Power)[1], param.ZinfMax_Min, param.ZinfMax_Max, Optim.GoldenSection())
			
			ZinfMax = Optim.minimizer(Optimization)[1]
		
			Of, θ_Time = θ_TIME_OF(iT_TransStead, ΔInfiltration_1D[1:iT_N], ΔTime[1:iT_N], θ_Ini, θs, θr, n, Ks, ZinfMax, iT_N, Power)
		
			# Just in case
			# θ_Time = min.(θ_Time[1:iT_N] *  θs / maximum(θ_Time[1:iT_N]), θs)

			return Of, θ_Time       
		end
        
        
        # #= =============== HYDRAU_VARY_INFILTRATION =============== =#
        #  function HYDRAU_VARY_INFILTRATION(θ_Time, Ks, Time, Inf_3D_Obs, RingRadius, Sorptivity, θs, θr, iT_TransStead, Km)
        #     # Dedicating memory
        #     Kr_Inf = Array{Float64}(iT_TransStead)
        #     Kr_Inf[1:iT_TransStead] = 1. # When the soil is saturated it is equal to 1
        #     Sorptivity_Inf = Array{Float64}(iT_TransStead)
        #     N_Inf =  Array{Float64}(iT_TransStead)
        #     Hvg_Inf =  Array{Float64}(iT_TransStead)
        
        #     θ_Time_Previous = θ_Time[1]
        #     # = Compute Kr at every θ_Time(iT), but to speed up we compute for every  ΔSe_Inf_Max =#
        #     for i_Inf in 1:iT_TransStead-1
        #         # To save time just computing for equaly intervalled Δθ_Time_Max
        #         if ((θ_Time[i_Inf] - θ_Time_Previous) / θs >= param.ΔSe_Inf_Max) && (param.Se_Inf_Max >= θ_Time[i_Inf] / θs) || ( i_Inf == 1)

        #         Infiltration_3D_Shift, Time_Shift, iT_N_Shift = thetaTime.TIME_INFILTRATION_SHIFT(i_Inf, iT_TransStead, Time[1:iT_TransStead], Inf_3D_Obs[1:iT_TransStead])


        #             ~, Kr_Inf[i_Inf], Sorptivity_Inf[i_Inf], N_Inf[i_Inf], Hvg_Inf[i_Inf] = quasiExact.vg.INFILTRATION3D_2_HYDRO(Time_Shift[1:iT_N_Shift], Infiltration_3D_Shift[1:iT_N_Shift], iT_N_Shift, θs, θ_Time[i_Inf], RingRadius, Km, false, Ks)
                
        #             println( i_Inf," , ", "Se_Inf ,", θ_Time[i_Inf] / θs, ",  Kr_Inf ,", Kr_Inf[i_Inf], " , Sorptivity ,", Sorptivity_Inf[i_Inf], ",N ," ,N_Inf[i_Inf], " , Hvg , ", Hvg_Inf[i_Inf], "\n" )

        #             θ_Time_Previous = θ_Time[i_Inf] # Making a copy of θ_Time
        #         else
        #             # Just copy the values of the previous timestep since there are too many time steps
        #             Kr_Inf[i_Inf] =  Kr_Inf[i_Inf-1]
        #             Sorptivity_Inf[i_Inf] =  Sorptivity_Inf[i_Inf-1]
        #             N_Inf[i_Inf] =  N_Inf[i_Inf-1]
        #             Hvg_Inf[i_Inf] = Hvg_Inf[i_Inf-1]
        #         end 
        #     end
        #     println("HYDRAU_VARY_INFILTRATION \n")
        #     #  Sorptivity_Inf[iT_TransStead-1: iT_N] = Sorptivity_Inf[iT_TransStead-1]
        #     return Kr_Inf, Sorptivity_Inf, N_Inf, Hvg_Inf
        # end

        #     #= =============== SELECTING BEST HYDRAULIC parameters =============== =#  
        # #= There are different options of Hvg and σ and we are looking for the best hydraulic param set which would give the best infiltration for any initilial soil moisture =#
        # function SELECT_Hkg_N(Time, Inf_3D_Obs, iT_N, iT_TransStead, θs, θr, N_Inf, Hvg_Inf, Km, θ_Time, Kr_Inf, Ks, Sorptivity, RingRadius)

        #     println("START of selecting the best hydraulic parameter sets")

        #     Of_Cum = Array{Float64}(iT_N)
        #     θ_Time_Previous = θ_Time[1]  
        #     for i_Inf in 1:iT_TransStead-1
        #         # To speed up
        #         if ((θ_Time[i_Inf] - θ_Time_Previous) / θs >= param.ΔSe_Inf_Max) && (param.Se_Inf_Max >= θ_Time[i_Inf] / θs) || ( i_Inf == 1)

        #             Infiltration_3D_Shift, Time_Shift, iT_N_Shift = thetaTime.TIME_INFILTRATION_SHIFT(i_Inf, iT_N, Time[1:iT_N], Inf_3D_Obs[1:iT_N])

        #             Infiltration_Sim = quasiExact.HYDRO_2_INFILTRATION3D(Time[1:iT_N], Sorptivity[1], Ks, Ks*Kr_Inf[1], θs, θ_Time[1], iT_N, RingRadius)
                    
        #             Of_Cum[i_Inf] = 0. 
        #             for iT in 1:iT_N
        #                 Of = (Infiltration_Sim[iT] - Inf_3D_Obs[iT])^2.
        #                 Of_Cum[i_Inf] = Of_Cum[i_Inf] + Of
        #             end
        #             θ_Time_Previous = θ_Time[i_Inf] # Making a copy of θ_Time 
        #             println( i_Inf," , ", "Se_Inf: ",  θ_Time[i_Inf]/ θs,  " Of_Cum: ", Of_Cum[i_Inf])
        #         else
        #             # Just provide a high value of OF for those that we skeeped
        #             Of_Cum[i_Inf] = 1000000000000000000.
        #         end        
        #     end
        #     # The position i of the hydraulic param sets which gives the best predictions for all initial θ
        #     Of_Cum_Min = minimum(Of_Cum[1:iT_TransStead-1])
        #     iOF_Min = searchsortedfirst( Of_Cum[1:iT_TransStead-1], Of_Cum_Min)
        #     println("iOF_Min, $iOF_Min")

        #     Hvg_Qe = Hvg_Inf[iOF_Min]
        #     N_Qe = N_Inf[iOF_Min]

        #      # Derive the optimal infiltration
        #     Sorptivity = sorptivity.vg.SORPTIVITY(θ_Time[1], θs, θ_Time[1], θs, 0., Hvg_Qe, N_Qe, Ks, Km)
        #     Se_Ini= wrc.se.θ_2_Se(θ_Time[1], θs, 0.)
        #     K_θini=  kunsat.vg.KUNSAT(Se_Ini, N_Qe, Ks, Km)
        #     Kr_θini = K_θini / Ks

        #      Infiltration_Qe_Vg = quasiExact.HYDRO_2_INFILTRATION3D(Time[1:iT_N], Sorptivity, Ks, K_θini, θs, θ_Time[1], iT_End, RingRadius)
 
 
        #     println("END of selecting the best hydraulic parameter sets \n")
        
        #     return Infiltration_Qe_Vg, N_Qe, Hvg_Qe  
        # end 

    end


end # Module thetaTime``