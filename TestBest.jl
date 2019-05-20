module testbest
	using sorptivity
	function DθDH_ACCURACY(iS, θs, θr, Hkg, σ, Hvg, N, Km)
		# Testing van Genuchten DθDh
		DθDh_Accur_Vg = sorptivity.vg.DθDh_CUMUL(θs, θr, Hvg, N, Km)
		DθDh_Accur_Vg = floor(Int, 100. * DθDh_Accur_Vg / (θs - θr))
		
		# Testing Kosugi DθDh
		DθDh_Accur_Kg = sorptivity.kg.DθDh_CUMUL(θs, θr, Hkg, σ)
		DθDh_Accur_Kg = floor(Int, 100 .0 * DθDh_Accur_Kg / (θs - θr))
		
		println("iS : $iS   DθDh_Accur_Vg : $DθDh_Accur_Vg %   DθDh_Accur_Kg : $DθDh_Accur_Kg %")
	end
end