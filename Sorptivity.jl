#= SORPTIVITY =#
module sorptivity

	#= =============== KOSUGI =============== =#
    module kg
        using wrc, kunsat, cst, diffusivity, QuadGK, SpecialFunctions, Optim, param, BlackBoxOptim, Roots,  SymEngine, Cuba, path
        export  SORPTIVITY, DθDh_CUMUL, SORPTIVITY_2_Hkg, SORPTIVITY_IniFinal, SORPTIVITY_TIME
     

        function SORPTIVITY(θ_Start, θ_End, θs, θr, Hkg, σ, Ks)
            Se_Max = 0.999999999999

            θ_Start = θ_Start + 10.^-15.
            Se_Start = wrc.se.θ_2_Se(θ_Start, θs, θr)

            θ_End = wrc.se.Se_2_θ(Se_Max, θs, θr)
            Se_End = wrc.se.θ_2_Se(θ_End, θs, θr)
         
            function SORPTIVITY_FUNC(Se, θ_End, θ_Start, θs, θr, Hkg, σ, Ks)
                Sorptivity_Func = (θ_End + wrc.se.Se_2_θ(Se, θs, θr) - 2.0 * θ_Start) * diffusivity.kg.DIFFUSIVITY(wrc.se.Se_2_θ(Se, θs, θr), θs, θr, Hkg, σ, Ks)
                return Sorptivity_Func
            end

            Sorptivity =  (θ_End - θ_Start) *QuadGK.quadgk(Se -> SORPTIVITY_FUNC(Se, θ_End, θ_Start, θs, θr, Hkg, σ, Ks), Se_Start, Se_End, abstol=10^-20.)[1] 

            return Sorptivity = (Sorptivity)^0.5
        end
        
        

        function SORPTIVITY2(θ_Start, θ_End, θ_Ini, θs, θr, Hkg, σ, Ks)
            Se_Max = 0.999999999999

            θ_Start = θ_Start + 10.^-15.

            θ_End = wrc.se.Se_2_θ(Se_Max, θs, θr)

            # θ_End = min(wrc.se.Se_2_θ(Se_Max, θs, θr), θ_End)

            function SEtran_2_θ(Se_Tran, θs, θr)
                Se = min(1. - (10.^-Se_Tran), Se_Max)
                θ = wrc.se.Se_2_θ(Se, θs, θr)
                return θ
            end

            function SEtran_2_Se(Se_Tran)
                Se = (1. - (10.^-Se_Tran))
                return Se
            end

            function SE_2_SEtran(Se)
                Se_Tran = -log10(1.0 -Se)
                return Se_Tran
            end

         
            function SORPTIVITY_FUNC(Se, θ_End, θ_Start, θs, θr, Hkg, σ, Ks)
                Sorptivity_Func1 = (θ_End + wrc.se.Se_2_θ(Se, θs, θr) - 2.0 * θ_Start) * diffusivity.kg.DIFFUSIVITY(wrc.se.Se_2_θ(Se, θs, θr), θs, θr, Hkg, σ, Ks)
                return Sorptivity_Func1
            end

            function SORPTIVITY_FUNC_2(Se_Tran, θ_End, θ_Start, θs, θr, Hkg, σ, Ks)
                Sorptivity_Func2 =  (Se_Tran * 10.^(-Se_Tran-1.))  * (θ_End + wrc.se.Se_2_θ(SEtran_2_Se(Se_Tran), θs, θr) - 2.0 * θ_Start) * diffusivity.kg.DIFFUSIVITY(wrc.se.Se_2_θ(SEtran_2_Se(Se_Tran), θs, θr), θs, θr, Hkg, σ, Ks)
                return Sorptivity_Func2
            end

            Sorptivity_1 =  (θs - θr) * QuadGK.quadgk(Se -> SORPTIVITY_FUNC(Se, θ_End, θ_Start, θs, θr, Hkg, σ, Ks), 0., 0.9, abstol=10^-15.)[1] 

            Sorptivity_2 =  (θs - θr) * QuadGK.quadgk(Se -> SORPTIVITY_FUNC_2(Se, θ_End, θ_Start, θs, θr, Hkg, σ, Ks), SE_2_SEtran(0.9), SE_2_SEtran(Se_Max), abstol=10^-15.)[1] 

            return Sorptivity = (Sorptivity_1 + Sorptivity_2)^0.5
        end
    


        function SORPTIVITY_CUBA(θ_Ini, θs, θr, Hkg, σ, Ks)
            θ_Start = θ_Ini + 10.^-15.

            Se_Max = wrc.kg.H_2_Se(param.H_Min_Kg, Hkg, σ)
            θ_End = min(wrc.se.Se_2_θ(Se_Max, θs, θr), θs-10.^-10.)
            Se_Max = min(Se_Max, Se_Max-10.^-10.)


            function SORPTIVITY_FUNC(Se, θ_Ini, θs, θr, Hkg, σ, Ks)
                Sorptivity_Func = (θs + wrc.se.Se_2_θ(Se, θ_End,  θ_Start) - 2.0 * θ_Ini) * diffusivity.kg.DIFFUSIVITY(wrc.se.Se_2_θ(Se, θ_End, θ_Start) , θs, θr, Hkg, σ, Ks) * ( θ_End -   θ_Start)
                return Sorptivity_Func
            end
            # try
            #vegas ; suave ;  cuhre
            # Integral = suave((Se, f) -> f[1] = SORPTIVITY_FUNC(Se[1], θ_Ini, θs, θr, Hkg, σ, Ks), θ_Start, θ_End)
            Integral = cuhre((Se, f) -> f[1] = SORPTIVITY_FUNC(Se[1], θ_Ini, θs, θr, Hkg, σ, Ks))
            Sorptivity = Integral[1][1]
            Sorptivity =  Sorptivity^0.5
           
            return Sorptivity
        end

        

        function SORPTIVITY_2_Hkg(Sorptivity, θ_Ini, θs, θr, σ, Ks)
        
            function SORPTIVITY_ROOT(Sorptivity, θ_Ini, θs, θr, Hkg, σ, Ks)

                 Sorptivity_Root = (log10(1.+Sorptivity) - log10(1.+ SORPTIVITY( θ_Ini, θs, θs, θr, Hkg, σ, Ks)))^2.
                 return  Sorptivity_Root
            end

            Optimization =  Optim.optimize(Hkg -> SORPTIVITY_ROOT(Sorptivity, θ_Ini, θs, θr, 10.^Hkg, σ, Ks), log10(param.Hkg_Min), log10(param.Hkg_Max), GoldenSection() )
            
            Hkg = 10.^Optim.minimizer(Optimization)[1]
            
            return Hkg
        end



        # CUMULATIVE DERIVATIVE Dθ/Dh to check the code
        function DθDh_CUMUL(θs, θr, Hkg, σ)
            CumulDerivative= quadgk( H -> wrc.kg.DθDh(H, θs, θr, Hkg, σ), 0., Inf, abstol=10^-8.)[1]
            return CumulDerivative
        end



        function SORPTIVITY_TIME(iT_TransStead, θ_Time, θ_Ini, θs, θr, Hkg, σ, Ks)
            Δθ_Treshold = 10^-4.
            Sorptivity_Time =Array{Float64}(iT_TransStead)
            # Sorptivity_Cumul = 0.
				# Sorptivity = 0.
				
				Sorptivity_Time[iT_TransStead] = 0.
				Sorptivity_Time[1] = 0. # for the very very rare case
				for iT in 1:iT_TransStead-1
                    if θs - θ_Time[iT]  >  Δθ_Treshold
						Sorptivity_Time[iT] = SORPTIVITY(θ_Time[iT], θs, θs, θr, Hkg, σ, Ks)
					else
						Sorptivity_Time[iT] = Sorptivity_Time[max(iT-1,1)]
					end
                end
                
            # for iT in iT_TransStead:-1:2
            #     if θ_Time[iT] - θ_Time[iT-1] >  Δθ_Treshold
				# 		  Sorptivity_Time[iT] = SORPTIVITY(θ_Time[iT-1], θ_Time[iT], θ_Ini, θs, θr, Hkg, σ, Ks)
				# 		  + Sorptivity_Time[iT+1]
						  
							
				# 		  Sorptivity = SORPTIVITY(θ_Time[iT-1], θs, θ_Ini, θs, θr, Hkg, σ, Ks)

				# 		  println( "$iT $(θ_Time[iT-1]) $(θ_Time[iT])  $Sorptivity   $(Sorptivity_Time[iT])")

            #         Err = Sorptivity_Time[iT] -Sorptivity

            #         # println("Sorptivity Err = $Err")
            #     else
            #         Sorptivity_Time[iT] = Sorptivity_Cumul
				# 	 end
					 
            # end
            # Sorptivity_Time[1] = 0.

            return Sorptivity_Time
        end



        function SORPTIVITY_DUAL(θ_Ini, θ_Ini2, σ_Synt, θs, θr, Hkg, σ, Ks)
            println("$θ_Ini, $θ_Ini2, $σ_Synt, $θs, $θr, $Hkg, $σ, $Ks")
            S_1 = SORPTIVITY(θ_Ini, θs, θr, Hkg, σ, Ks)
            S_2 = SORPTIVITY(θ_Ini2, θs, θr, Hkg, σ, Ks)

           
            Hkg_Synt = SORPTIVITY_2_Hkg(S_1, θ_Ini, θs, θr, σ_Synt, Ks)

            S_1_Synt = SORPTIVITY(θ_Ini, θs, θr, Hkg_Synt, σ_Synt, Ks)
            S_2_Synt = SORPTIVITY(θ_Ini2, θs, θr, Hkg_Synt, σ_Synt, Ks)

            println("S_1=  $S_1  S_1_Synt= $S_1_Synt  S_2= $S_2   S_2_Synt= $S_2_Synt     ")
            # S_1=  0.5606979358854783  S_1_Synt= 0.5606979368630025  S_2= 1.2020525284191446   S_2_Synt= 1.1426394002426912
        end
    end


    #= =============== VAN GENUCHTEN =============== =#
    module vg
        using QuadGK, wrc, diffusivity, kunsat, cst, wrc, param, sorptivity, BlackBoxOptim, Optim
        export  DθDh_CUMUL, SORPTIVITY_2_Hvg
    
        function SORPTIVITY(θ_Start, θ_End, θs, θr, Hvg, N, Ks, Km)
            Se_Max = 0.9999999

            θ_Start = θ_Start + 10.^-15.
            # Se_Start = wrc.se.θ_2_Se(θ_Start, θs, θr)

            θ_End = wrc.se.Se_2_θ(Se_Max, θs, θr)
            # Se_End = wrc.se.θ_2_Se(θ_End, θs, θr)
         

            function SORPTIVITY_FUNC(θ, θ_Start, θ_End, θs, θr, Hvg, N, Ks, Km)
                Sorptivity = (θ_End + θ - 2.0 * θ_Start) * diffusivity.vg.DIFFUSIVITY(θ, θs, θr, Hvg, N, Ks, Km)
                return Sorptivity
            end

            Sorptivity =  quadgk(θ -> SORPTIVITY_FUNC(θ, θ_Start, θ_End, θs, θr, Hvg, N, Ks, Km), θ_Start, θ_End, abstol=10^-20.)[1]
            return Sorptivity = Sorptivity^0.5 
        end


        function SORPTIVITY2(θ_Start, θ_End, θ_Ini, θs, θr, Hvg, N, Ks, Km)
            θ_Start = θ_Start + 10.^-15.
            
            Se_Max = wrc.vg.H_2_Se(param.H_Min_Vg, Hvg, N, Km)
            # SorptSeMax = min(max(1.0 -10.^(-5.3223273667*N + 5.4377861098),0.75),0.999999999)
            
            # θ_End = wrc.se.Se_2_θ(Se_Max, θs, θr)
            # θ_End = min(θ_End, θs-10.^-10.)

            θ_End = min( wrc.se.Se_2_θ(Se_Max, θs, θr), θs-10.^-10., θ_End)
        
            function SORPTIVITY_FUNC(θ, θ_Ini, θs, θr, Hvg, N, Ks, Km)
                Sorptivity = (θs + θ - 2.0 * θ_Ini) * diffusivity.vg.DIFFUSIVITY(θ, θs, θr, Hvg, N, Ks, Km)
                return Sorptivity
            end

            Sorptivity =  quadgk(θ -> SORPTIVITY_FUNC(θ, θ_Ini, θs, θr, Hvg, N, Ks, Km), θ_Start, θ_End, abstol=10^-15.)[1]
            return Sorptivity = Sorptivity^0.5 
        end


    # SORPTIVITY 2 Hkg VAN GENUCHTEN HYDRAULIC param
     function SORPTIVITY_2_Hvg(Sorptivity, θ_Ini, θs, θr, N, Ks, Km)
        
            function SORPTIVITY_ROOT(Sorptivity, θ_Ini, θs, θr, Hvg, N, Ks, Km)
                 Sorptivity_Root = (log10(1. + Sorptivity) - log10(1. + sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg, N, Ks, Km)))^2.
            end

            Optimization =  optimize(Hvg ->  SORPTIVITY_ROOT(Sorptivity, θ_Ini, θs, θr, 10.^Hvg, N, Ks, Km), log10(param.Hvg_Min), log10(param.Hvg_Max))
            Hvg = 10.^Optim.minimizer(Optimization)[1]
            return Hvg
        end


        function DθDh_CUMUL(θs, θr, Hvg, N, Km)
            H_Min = 0.
            H_Max = 1E+10
            m = 1. - Km / N
            DθDh_Cumul =  quadgk(H -> wrc.vg.DθDh(H, θs, θr, Hvg, N, Km), H_Min, H_Max, abstol=10^-8.)[1]
            return DθDh_Cumul
        end
    end
end
