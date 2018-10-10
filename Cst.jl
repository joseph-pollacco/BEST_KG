#= CONSTANTS =#
module cst
		bd = 0.916 # bulk density [g cm-3] 0.95  Sam's dataset 0.916 # in
		pd = 2.65 # particle density [g cm-3] 2.59 # Sam's dataset 2.65  # in (g cm-3)
		γ = 0.75 # Shape param
		β = 0.6  # 0.60 Ratio defined using the ratio between 2 estimators of the Sorptivity
		ϵ = 0.00000001
		Y = 0.149 * 10.0 ^2. #[mm^2]
		Kconst = (10. / (24. * 60. * 60.)) *( 1.03663 * 10. ^9.) #  convert from cm/day to mm s−1.
		Pσ_1 = 0.5920
		Pσ_2 = 0.7679
end
