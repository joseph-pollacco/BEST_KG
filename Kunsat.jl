   
 DIR_Working = pwd() # Saving the current working directory
 push!(LOAD_PATH, DIR_Working) # All the subroutines are in the current working directory
 include("Cst.jl")  

module kunsat

    #= =============== BROOKS AND COOREY =============== =#
    module bc
        export SE_2_KR
        function SE_2_KR(Se, N, Km=2, p=1) # Brooks and Corey UHC
            ɳ = 2/(N-Km)+2+p
            return BC_SE2Kr = Se^ɳ
        end
    end
    
    
    #= =============== KOSUGI =============== =#
    module kg
        export  KUNSAT, KRθini_2_σ, KUNSAT_INTEGRAL
        using SpecialFunctions,  QuadGK

        function KUNSAT(Se, θs, θr, σ, Ks, θs_Mac, σ_Mac)
            Ks_Mat = Ks * (θs_Mac - θr) / (θs - θr)
            Kunsat_Mat = Ks_Mat * sqrt(Se) * (0.5*erfc( erfcinv(2.0*Se) + σ/sqrt(2.0) ))^2.0

            Ks_Mac = Ks * (θs - θs_Mac) / (θs - θr)
            Kunsat_Mac = Ks_Mac * sqrt(Se) * (0.5*erfc( erfcinv(2.0*Se) + σ_Mac/sqrt(2.0) ))^2.0

             Kunsat = Kunsat_Mat + Kunsat_Mac
            return Kunsat
         end

        #  function Kr_2_Se(Ks, θs, θr, σ, Ks)
            
        #      =abs (sqrt(Se) * (0.5*erfc( erfcinv(2.0*Se) + σ/sqrt(2.0) ))^2.0 - Kr)

        #     return Se
        #  end

         
        function KUNSAT_INTEGRAL(θ_Ini, θ_Final, θs, θr, σ, Ks)
            Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)
            Se_Final = wrc.se.θ_2_Se(θ_Final, θs, θr)
            Kunsat_Integral = quadgk( Se -> KUNSAT(Se, θs, θr, σ, Ks, θs, σ), Se_Ini, Se_Final, abstol=10^-8.)[1]
            return Kunsat_Integral   
        end

         
		function KRθini_2_σ(Se_Ini, Kr_θini)
			σ = (erfcinv(2. * (Kr_θini / sqrt(Se_Ini)) ^0.5) - erfcinv(2.0*Se_Ini)) * sqrt(2.0)
		return σ
		end
    end


    #= =============== VAN GENUCHTEN =============== =#
    module vg
    using Optim, SpecialFunctions
    export KUNSAT, KRθini_2_N


        function KUNSAT(Se, N, Ks, Km=1., L=0.5)
            # L = Km/2+Km-1
            M = 1 - Km/N
            # kunsat = Ks * (Se.^L) * ( 1. - (1. - Se.^(1/M) ).^M ).^(2/Km)
            kunsat = Ks * (Se.^L) * ( 1.0 - (1.0 - Se.^(1.0/M) ).^M ).^2.0
            return kunsat
        end


        function KRθini_2_N(Se, Kr, Ks, Km=1., L=0.5)
        
            function OF(Se, Kr, N, Km)
                OF_N = (KUNSAT(Se, N, Ks, Km) / Ks - Kr) ^ 2.
                return OF_N
            end

			# N_Min and N_Max depends on Km
            if Km == 1
                N_Min = param.N_Km1_Min
                N_Max = param.N_Km1_Max
            elseif  Km == 2
                N_Min = param.N_Km2_Min
                N_Max = param.N_Km2_Max
            end
            Optimization =  optimize(N ->  OF(Se, Kr, N, Km), N_Min, N_Max, GoldenSection())
            N = Optim.minimizer(Optimization)[1]
            return N
        end
    end
end
