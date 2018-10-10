module option
    DownloadPackage = false
    
    σ_Known = false

    # Data available
    module DataAvailable
        HydraulicParam = true ##  [false] or [true]
        H_Ini = true ##  [false] or [true]
    end

    # θ(h) and k(θ) model used
    module HydroModel
        HydraulicParam = "vangenuchten" ## [kosugi] or [vangenuchten] or [brooks]
    end

    # Inf model used
    module infiltration
        Model = "QuasiExactNormFast" #[QuasiExactNormSlow],  or [QuasiExact] or [QuasiExactNormFast]
    end

    module Plot
        WRC_Kunsat = true
    end

    module output
        TestingDθDh = false
        TestingSorpt = false
        Best = true
        QuasiExact = true
    end
end
