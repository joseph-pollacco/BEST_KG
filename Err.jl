module err
   module kg
      using diffusivity, cst, kunsat, wrc, sorptivity
      export ERR_QE

            function ERR_QE(θs, θr, Hkg, σ, Ks)
                  Err_Bad = 10.0^-5.
                  Se_Ini = 0.25
                  θ_Ini = wrc.se.Se_2_θ(Se_Ini, θs, θr)
                  Se_1 = 0.5 # Soil moisture we want to monitor
                  θ_1 = wrc.se.Se_2_θ(Se_1, θs, θr)
						
						Sorpt = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg, σ, Ks)

                  function DIFF(θ, θs, θr, Hkg, σ, Ks)
							return Diff = diffusivity.kg.DIFFUSIVITY(θ, θs, θr, Hkg, σ, Ks)
                  end

                  Integral = quadgk( θ -> DIFF(θ, θs, θr, Hkg, σ, Ks), θ_1, θs-10.0^-8, abstol=10.0^-8.)[1]

                  LeftTerm = (kunsat.kg.KUNSAT(Se_1, θs, θr, σ, Ks, θs, σ) - kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)) / (Ks - kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ))

                  Err = abs(((θ_1 - θ_Ini) / (θs - θ_Ini)) * (1. - 2.*cst.β * (θs - θ_Ini) * Integral / Sorpt^2.) - LeftTerm)
                  return Err
            end





         function ERR_QE2(θs, θr, Hkg, σ, Ks)
            Err_Bad = 10.0^-5.

            Se_Ini_1 = 0.25
            θ_Ini_1 = wrc.se.Se_2_θ(Se_Ini_1, θs, θr)
            
            Se_1 = 0.5 # Soil moisture we want to monitor
            θ_1 = wrc.se.Se_2_θ(Se_1, θs, θr)
            
            function fθ(θ_Ini, θ, θs)
               return (2.* (θ - θ_Ini)) / (θs + θ - 2.*θ_Ini)
            end

            function INTEGRAL(θ, θ_Ini, θs, θr, Hkg, σ, Ks)
               Integral = (θ - θ_Ini) * diffusivity.kg.DIFFUSIVITY(θ, θs, θr, Hkg, σ, Ks) / fθ(θ_Ini, θ, θs)
            end

            Integral_1 = quadgk( θ -> INTEGRAL(θ, θ_Ini_1, θs, θr, Hkg, σ, Ks), θ_1, θs-10.0^-8, abstol=10.0^-8.)[1]

            Integral_2 = quadgk( θ -> INTEGRAL(θ, θ_Ini_1, θs, θr, Hkg, σ, Ks), θ_Ini_1, θs-10.0^-8, abstol=10.0^-8.)[1]

            LeftTerm = (kunsat.kg.KUNSAT(Se_1, θs, θr, σ, Ks, θs, σ) - kunsat.kg.KUNSAT(Se_Ini_1, θs, θr, σ, Ks, θs, σ)) / (Ks - kunsat.kg.KUNSAT(Se_Ini_1, θs, θr, σ, Ks, θs, σ))

            Err = fθ(θ_Ini_1, θ_1, θs) * (1. - cst.β * Integral_1 / Integral_2) - LeftTerm
            return Err
         end
   end
end