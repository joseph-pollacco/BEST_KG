
module diffusivity
	#= =============== KOSUGI =============== =#
	module kg
		using wrc, kunsat, param
		export DIFFUSIVITY
		
		function DIFFUSIVITY(θ, θs, θr, Hkg, σ, Ks)
			Se = wrc.se.θ_2_Se(θ, θs, θr)
			
			Kunsat = kunsat.kg.KUNSAT(Se, θs, θr, σ, Ks, θs, σ)
		
			H = wrc.kg.Se_2_H(Se, Hkg, σ)
		
			DθDh = wrc.kg.DθDh(H, θs, θr, Hkg, σ)
				
			return Diffusivity = Kunsat /  DθDh
		# Diffusivity2 = sqrt(Se*π/8.)*Hkg*σ*Ks*((erfc( erfcinv((2.0 *Se)) + σ/sqrt(2.) ))^2.)/(exp(-(2.+σ*sqrt(2.))*erfcinv(2.0 *Se)))end
		end
	end
	 
	   
    #= =============== VAN GENUCHTEN =============== =#
	module vg
		using wrc, kunsat, param
		export DIFFUSIVITY
		
		# This gave the same results as D. Moret-Fernández, B. Latorre / Journal of Hydrology 544 (2017) 352–362
        function DIFFUSIVITY(θ, θs, θr, Hvg, N, Ks, Km) # van genuchten diffusivity        
			Se = wrc.se.θ_2_Se(θ, θs, θr)

            Kunsat = kunsat.vg.KUNSAT(Se, N, Ks, Km)

            H = wrc.vg.Se_2_H(Se, Hvg, N, Km)
            
            DθDh = wrc.vg.DθDh(H, θs, θr, Hvg, N, Km)

            return Diffusivity = Kunsat / (DθDh + 1.e-20)
		end

		
		
		function DIFFUSIVITY_2(θ, θs, θr, Hvg, N, Ks, Km)
			M = 1 - Km/N
			
			Se = wrc.se.θ_2_Se(θ, θs, θr)

			# Jesus diffusivity something is not right
			# Diffusivity = ((1.0 -M)*Ks*Hvg)/(Km*M*(θs-θr)) * (Se.^ (((Km^3)*M+Km-3)/(2*M))) * ( (1.0 -Se.^(1/M) ).^((1-M-Km)/(Km)) - ((-1)^Km)*(1.0 -Se.^(1/M)).^((M+1-Km)/(Km)) + 2.0 *(2.0 -Km) )

			# D. Moret-Fernández, B. Latorre / Journal of Hydrology 544 (2017) 352–362
			Diffusivity=(((1.0 -M)*Ks*Hvg)/(M*(θs-θr))) * (Se.^ (0.5-(1.0 /M))) * ( ((1.0 -Se.^(1.0 /M)).^(-M))  +  ((1.0 -Se.^(1/M)).^M)  -2. )
			
            return Diffusivity
        end
	end
end