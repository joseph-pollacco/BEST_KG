module stats

using param
export NASH_SUTCLIFFE, NASH_SUTCLIFFE_OF, INFILTRATION_2_iSTEADYTRANSIT, RELATIVEerr

	function NASH_SUTCLIFFE(Obs, Sim)
		N = length(Obs)
      Obs_Mean = mean(Obs[1:N])
      Obs_Mean_Err =  sum((Obs_Mean - Obs[1:N]).^2.)
      Err = sum((Sim[1:N] - Obs[1:N]).^2.)
      return 1. - Err / Obs_Mean_Err 
   end



	function NASH_SUTCLIFFE_OF(Obs, Sim, Power=2.)
		N = length(Obs)
      Obs_Mean = mean(Obs[1:N])
      Obs_Mean_Err =  sum(abs.(Obs_Mean - Obs[1:N]).^Power)
      Err = sum(abs.(Sim[1:N] - Obs[1:N]).^Power)
      return Err / Obs_Mean_Err 
   end



	function RELATIVEerr(Obs, Sim)
		Err = 1. - abs(Obs - Sim ) / Obs
		return Err
	end


	
   function REGRESSION_STEADY_STATE_INF(Inf_data, T_data, good_fit=3) # linear regression analysis of the last data pairs (cumulative infiltration-time)
		R_squared = zeros(length(Inf_data)-1)
		for i in 1:(length(Inf_data)-1)
			R_squared[i] = (cor(T_data[i:length(Inf_data)],Inf_data[i:length(Inf_data)]))^2
			println(R_squared[i])
		end
		if length(R_squared[R_squared.==1.0]) == good_fit # for when perfect fit at the last infiltration data
			j = length(Inf_data)-good_fit
		else
			j = find(x -> x == maximum(R_squared[R_squared.<1.0]), R_squared) # position of the maximun (< than 1)
			j = j[1] # convert Array{Int64,1} into Int64 to be used as a position index
		end
		intercept = linreg( T_data[j:length(Inf_data)], Inf_data[j:length(Inf_data)])[1] # intercept of the fit, units (mm)
		slope = linreg( T_data[j:length(Inf_data)],Inf_data[j:length(Inf_data)] )[2] # slope of the fit, units (mm s-1)
		R² = (cor(T_data[j:length(Inf_data)],Inf_data[j:length(Inf_data)]))^2 # R^2 of the fit or just R² = R_squared[j]

		println( "The j = , $j")
		
		return j, slope, intercept, R² # steady-state infiltration rate units (mm s-1) and intercept units (mm) if Inf_data in (mm) and T_data (s)
	   
	end
	


	# iT for steady transit
	function INFILTRATION_2_iSTEADYTRANSIT(ΔInf_SteadyTransit, Time, Infiltration, Slope_N=6) # linear regression analysis of the last data pairs (cumulative infiltration-time)
		N_Data = length(Time)
		Slope = zeros(N_Data)
		Intercept = zeros(N_Data)
		Inf_Model = zeros(N_Data)
		Err = zeros(N_Data)
	
		iSeadyTransit = 5.
		FlagSteadyState = false
		for i in N_Data-Slope_N:-1:1
		
			if i-1 >= 1

				# Intercept[i], Slope[i] = linreg( Time[i-Slope_N+1:i], Infiltration[i-Slope_N+1:i])
				Intercept[i], Slope[i] = linreg(Time[i: N_Data], Infiltration[i: N_Data])

				Inf_Model[i-1] = Time[i-1]* Slope[i] + Intercept[i]
				
				Err[i-1] = abs( Infiltration[i-1] - Inf_Model[i-1]) / (Time[i] - Time[i-1])
				if Err[i-1]  >= ΔInf_SteadyTransit && !FlagSteadyState # To catch only the very beginning
					FlagSteadyState = true
					iSeadyTransit = i-1
				end

				# println(i-1, " , ", Slope[i], " , " ,Intercept[i], " , ", Inf_Model[i-1], " , " , Infiltration[i-1], " , ", Err[i-1])	
			end
		end
		
		Time_TransStead = Time[iSeadyTransit]

		println("iSeadyTransit= , $iSeadyTransit, $Time_TransStead")

		return iSeadyTransit, Time_TransStead
	end


end