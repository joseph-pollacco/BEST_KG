module convertHydroModel
    using cst, SpecialFunctions, BlackBoxOptim, param, wrc, kunsat, Optim, diffusivity, stats, Sorptivity
    export VANGENUCHTEN_2_KOSUGI, KOSUGI_2_VANGENUCHTEN

       #= Converting vangenuchten params Hvg, N to Sigma, Hkg =#
       function VANGENUCHTEN_2_KOSUGI(Hvg, N, Km=1.)
        
        Optimization = bboptimize(Param ->  OBJECTIVE_FUNCTION_WRC_Kunsat(Hvg, N, Km, 10.^Param[1], Param[2]) ; SearchRange =[ (log10(param.Hkg_Min), log10(param.Hkg_Max)), (param.σ_Min, param.σ_Max)], NumDimensions=2, TraceMode=:silent)
        # Values of the optimal hydraulic params
        Hkg = 10.^(best_candidate(Optimization)[1])
        σ = best_candidate(Optimization)[2]

        #Statistics
        NSE, NSE_Hse, NSE_Kunsat = NASH_SUTCLIFFE_WRC_Kunsat(Hvg, N, Km, Hkg, σ)

        return Hkg, σ, NSE, NSE_Hse, NSE_Kunsat
    end
   

   
    #= Converting Sigma, Hkg to vangenuchten params to Hvg, N =#
    function KOSUGI_2_VANGENUCHTEN(Hkg, σ; Km=1.)
        if Km == 1
            N_Min = param.N_Km1_Min
            N_Max = param.N_Km1_Max
        elseif  Km == 2
            N_Min = param.N_Km2_Min
            N_Max = param.N_Km2_Max
        end   
        
        Optimization = bboptimize(Param ->  OBJECTIVE_FUNCTION_WRC_Kunsat(10.^Param[1], Param[2], Km, Hkg, σ) ; SearchRange =[ (log10(param.Hvg_Min), log10(param.Hvg_Max)), (N_Min, N_Max)], NumDimensions=2, TraceMode=:silent)
        # Values of the optimal hydraulic params
        Hvg = 10.^(best_candidate(Optimization)[1])
        N = best_candidate(Optimization)[2]

        #Statistics
        NSE, NSE_Hse, NSE_Kunsat = NASH_SUTCLIFFE_WRC_Kunsat(Hvg, N, Km, Hkg, σ)

        return Hvg, N, NSE, NSE_Hse, NSE_Kunsat
    end



    #= Converting Sigma, Hkg to vangenuchten params to Hvg, N =#
    function KOSUGI_2_VANGENUCHTEN_2STEPS(Hkg, σ; Km=1)
        # We optimize first N by fitting Kr than we fit hr with θ(h)
        			# N_Min and N_Max depends on Km
        if Km == 1
            N_Min = param.N_Km1_Min
            N_Max = param.N_Km1_Max
        elseif  Km == 2
            N_Min = param.N_Km2_Min
            N_Max = param.N_Km2_Max
        end

        # STEP 1: Optimize N
        Optimization =  optimize(N ->  OBJECTIVE_FUNCTION_Kunsat(N, Km, σ), N_Min, N_Max )
        N = Optim.minimizer(Optimization)[1]
        
        # STEP 2: Optimize Hvg
        # Hvg_Log_Min = max(log10(Hvg / 1000.), log10(param.Hvg_Min))
        # Hvg_Log_Max = min(log10(Hvg * 1000.), log10(param.Hvg_Max))
        Optimization =  optimize(Hvg -> OBJECTIVE_FUNCTION_WRC(10.^ Hvg, N, Km, Hkg, σ), log10(param.Hvg_Min), log10(param.Hvg_Max), GoldenSection() )
        OF_WRC = Optim.minimum(Optimization)
        Hvg = 10.^(Optim.minimizer(Optimization)[1])
        return Hvg, N
    end



    #= Converting Hvg, N vangenuchten params to Hkg, σ =#
    function VANGENUCHTEN_2_KOSUGI_2STEPS(Hvg, N, Km)
        # We optimize first σ by fitting Kr than we fit hr with θ(h)

        # STEP 1: Optimize σ
        Optimization =  optimize(σ ->  OBJECTIVE_FUNCTION_Kunsat(N, Km, σ), param.σ_Min, param.σ_Max )
        σ = Optim.minimizer(Optimization)[1]
        
        # STEP 2: Optimize Hkg
        Hkg_Log_Min = max(log10(Hvg / 1000.), log10(param.Hvg_Min))
        Hkg_Log_Max = min(log10(Hvg * 1000.), log10(param.Hvg_Max))
        Optimization =  optimize(Hkg -> OBJECTIVE_FUNCTION_WRC(Hvg, N, Km, 10.^Hkg, σ), Hkg_Log_Min, Hkg_Log_Max, GoldenSection() )
        OF_WRC = Optim.minimum(Optimization)
        Hkg = 10.^(Optim.minimizer(Optimization)[1])
        return  Hkg, σ
    end

    function OBJECTIVE_FUNCTION_Kunsat(N, Km, σ)
        θs= 1. # Can be any value does not matter
        θr = 0. # Can be any value does not matter
        Ks = 1. # We are computing Ks
        Se = linspace(0.001, 0.99999, 1000) # Range of Se
        OF = 0.
        for i_Se in Se
            Kr_Kg = kunsat.kg.KUNSAT(i_Se, θs, θr, σ, Ks, θs, σ)
            Kr_Vg = kunsat.vg.KUNSAT(i_Se, N, Ks, Km)
            OF_Kr = ((100. * (1. + Kr_Kg)) - (100.*(1. + Kr_Vg)))^2.
            OF = OF_Kr + OF
        end
        return OF
    end


    function OBJECTIVE_FUNCTION_WRC(Hvg, N, Km, Hkg, σ)
            
        H_Linspace = (linspace(-5, 10., 1000)) # Range of H
        OF = 0.
        for iH_Log10 in  H_Linspace
            H =  10.^iH_Log10
            Se_Kg = wrc.kg.H_2_Se(H, Hkg, σ)
            Se_Vg = wrc.vg.H_2_Se(H, Hvg, N, Km)
            
            OF_Se = ((1. + Se_Kg) - (1. + Se_Vg))^2.
            OF = OF + OF_Se
        end
        return OF
    end

    function OBJECTIVE_FUNCTION_WRC_Kunsat(Hvg, N, Km, Hkg, σ)
        θs= 1. # Can be any value does not matter
        θr = 0. # Can be any value does not matter
        Ks = 1. # We are computing Ks
        Se_Kg = Array{Float64}(1001)
        Se_Vg = Array{Float64}(1001)
            
        H_Linspace = (linspace(-5, 10., 1000)) # Range of H
        i=0
        for iH_Log10 in  H_Linspace
            i=i+1
            H =  10.0^iH_Log10
            Se_Kg[i] = wrc.kg.H_2_Se(H, Hkg, σ)
            Se_Vg[i] = wrc.vg.H_2_Se(H, Hvg, N, Km)
        end
        OF_Hse =stats.NASH_SUTCLIFFE_OF((1.+Se_Kg[1: 1000]),(1.+Se_Vg[1: 1000]))

        Kr_Kg= Array{Float64}(1001)
        Kr_Vg= Array{Float64}(1001)
        Se = linspace(0.0, 1., 1000) # Range of Se
        i=0
        for i_Se in Se
            i=i+1
            Kr_Kg[i] = kunsat.kg.KUNSAT(i_Se, θs, θr, σ, Ks, θs, σ)
            Kr_Vg[i] = kunsat.vg.KUNSAT(i_Se, N, Ks, Km)
        end
        OF_K = stats.NASH_SUTCLIFFE_OF(log.(1+Kr_Kg[1:1000]), log.(1+Kr_Vg[1:1000]))

        WOF = 0.5*OF_Hse +  0.5*OF_K

        return WOF
    end

    function NASH_SUTCLIFFE_WRC_Kunsat(Hvg, N, Km, Hkg, σ)
        θs= 1. # Can be any value does not matter
        θr = 0. # Can be any value does not matter
        Ks = 1. # We are computing Ks
        Se_Kg = Array{Float64}(1001)
        Se_Vg = Array{Float64}(1001)
            
        H_Linspace = (linspace(-5, 10., 1000)) # Range of H
        i=0
        for iH_Log10 in  H_Linspace
            i=i+1
            H =  10.0^iH_Log10
            Se_Kg[i] = wrc.kg.H_2_Se(H, Hkg, σ)
            Se_Vg[i] = wrc.vg.H_2_Se(H, Hvg, N, Km)
        end
        NSE_Hse =abs(stats.NASH_SUTCLIFFE( (1.+Se_Kg[1: 1000]),(1.+Se_Vg[1: 1000])))

        Kr_Kg= Array{Float64}(1001)
        Kr_Vg= Array{Float64}(1001)
        Se = linspace(0.001, 0.99999, 1000) # Range of Se
        i=0
        for i_Se in Se
            i=i+1
            Kr_Kg[i] = kunsat.kg.KUNSAT(i_Se, θs, θr, σ, Ks, θs, σ)
            Kr_Vg[i] = kunsat.vg.KUNSAT(i_Se, N, Ks, Km)
        end
        NSE_Kunsat = abs(stats.NASH_SUTCLIFFE( log.(1+Kr_Kg[1:1000]), log.(1+Kr_Vg[1:1000])))

        NSE = (NSE_Hse +  NSE_Kunsat) / 2.

        return NSE, NSE_Hse, NSE_Kunsat
    end

end
