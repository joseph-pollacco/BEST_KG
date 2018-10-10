module relationship
	using cst, param, BlackBoxOptim, sorptivity, Optim, stats
    export σ_2_Hkg, Hkg_2_σ, OPTIMIZATION_σ_2_Hkg
    


    function Hkg_2_σ(Hkg)
        σ = cst.Pσ_1*( log(Hkg) -1. )^cst.Pσ_2 # Hkg[mm]
        return σ = max(min(σ, param.σ_Max), param.σ_Min)
    end
    


    function σ_2_Hkg(σ)
        Hkg = exp(exp(inv(cst.Pσ_2) * log(σ / cst.Pσ_1)) + 1.) #[mm]
        Hkg = max(min(Hkg , param.Hkg_Max), param.Hkg_Min)
        return  Hkg
    end


    function OPTIMIZATION_σ_2_Hkg(SampleTrue, Hkg_Obs, σ_Obs)

        function OF_σ_2_Hkg(SampleTrue, Hkg_Obs, σ_Obs, Pσ_1, Pσ_2)
            Of = 0
            i = 1
            for σ in σ_Obs
                # if SampleTrue[i] == 1
                    Hkg_Sim = exp( exp(inv(Pσ_2) * log(σ / Pσ_1)) + 1. ) #[mm]
                    Hkg_Sim = max(min(Hkg_Sim , param.Hkg_Max), param.Hkg_Min)
                    Of += (log(Hkg_Sim) - log(Hkg_Obs[i]))^2.  
                # end
                i+=1
            end
            return Of 
        end

        Optimization = BlackBoxOptim.bboptimize(Pσ ->  OF_σ_2_Hkg(SampleTrue, Hkg_Obs, σ_Obs, Pσ[1], Pσ[2]); SearchRange = [(0., 100.), (0., 100.)], NumDimensions=2, TraceMode=:silent)
			
        Pσ_1 = BlackBoxOptim.best_candidate(Optimization)[1]
        Pσ_2 = BlackBoxOptim.best_candidate(Optimization)[2]

        println("Pσ_1 = $Pσ_1, Pσ_2 = $Pσ_2")

        return Pσ_1, Pσ_2
    end
    
    
	# Derive Hm_mac which delineates macropore with matrx flow
    function Func_Hm_mac(Hkg, σ)
        Hm_mac = exp( log(Hkg) - σ * 8.^0.5 )
        return Hm_mac
    end



    # Derive the paremater hm_mac automatically
    function Func_Hkg_mac(Hm_mac)
        Hkg_mac = exp(log(Hm_mac) / 2.0)
        return Hm_mac
    end	
end