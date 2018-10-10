   
    #= =============== Se =============== =#
module wrc
    module se
        using cst
        export  θ_2_Se, Se_2_θ

        function θ_2_Se(θ, θs, θr=0.)
            Se = (θ-θr) / (θs-θr)
            Se = max(min(Se, 1.), 0.)
            return Se
        end


        function Se_2_θ(Se, θs, θr=0)
            θ = Se * (θs-θr) + θr
            θ = max(min(θ, θs), θr)
            return  θ
        end

    end


	#= =============== KOSUGI =============== =#
    module kg
        using cst, SpecialFunctions, wrc
        export H_2_Se, Se_2_H, H_2_θdual, DθDr_Dual, DθDh

        function H_2_Se(H, Hkg, σ) # Kosugi WRC
            return Se = 0.5*erfc((log(H)-log(Hkg))/(σ*sqrt(2.)))
         end


        function Se_2_H(Se, Hkg, σ) # Kosugi WRC
             H = exp(erfcinv(2.*Se)*σ*sqrt(2.) + log(Hkg))
             return H
        end


        function H_2_θdual(H, θs, θr, Hkg, σ, θs_Mac, Hkg_Mac, σ_Mac) # Kosugi WR
            θ_Mac = (θs - θs_Mac) * 0.5 * erfc((log(H) - log(Hkg_Mac)) / (σ_Mac*sqrt(2.)))
            θ_Mat = (θs_Mac - θr) * 0.5 * erfc((log(H) - log(Hkg)) / (σ*sqrt(2.))) + θr
            θ = θ_Mac + θ_Mat
            return θ
        end


        function DθDr_Dual(r, θs, θr, Hkg, σ, θs_Mac, Hkg_mac, σ_Mac) # Kosugi
            #error with exp()
            rm = cst.Y / Hkg
            rm_mac = cst.Y / Hkg_mac

            PDF_mac = (θs-θs_Mac) * (exp( -((log( r/rm_mac ))^2.) / (2 *σ_Mac^2.))) / (r *σ_Mac * sqrt( 2.*π ))

            PDF_mat =  (θs_Mac-θr) * (exp( -(log( r/rm )^2.) / (2 * σ^2.))) / (r * σ * sqrt( 2.*π ))

            PDFnorm = PDF_mac + PDF_mat
            return PDFnorm
        end


        function DθDh(H, θs, θr, Hkg, σ) # Kosugi
            Dθ_Dh =(θs-θr)* exp( -((log(H/Hkg))^2. ) / (2.*σ^2.) ) / ( H*σ*sqrt(2.*π) )
            return Dθ_Dh
    end
        
#         function DθDh_R(H, θs, θr, Hkg, σ) # Kosugi
# 	         r = cst.Y / (H + cst.ϵ)
# 	         rm = cst.Y / Hkg
# 	         Dθ_Dr =  (θs-θr) * exp( -((log(r/rm))^2.) / (2.*σ^2.) ) / ( r*σ*sqrt(2.*π) )
# 	         Dθ_Dr_h =  (θs-θr) * exp( -((log(Hkg/H))^2.) / (2.*σ^2.) ) / ( (cst.Y/H)*σ*sqrt(2.*π) )
# 	         Dθ_Dh2 = Dθ_Dr_h * cst.Y / H^2.
#              return Dθ_Dh
#         end
    end

  
    #= =============== VAN GENUCHTEN =============== =#
    module vg
        using cst, SpecialFunctions, wrc
        export H_2_Se, Se_2_H, H_2_θ, DθDh

        function H_2_Se(H, Hvg, N, Km=1.) # van Genuchten WRC
            M = 1 - Km/N
            Se = (1. + (H/Hvg)^N )^(-M)
            return Se
        end


        function H_2_θ(H, θs, θr, Hvg, N, Km=1.) # van Genuchten WRC
            M = 1 - Km/N
            Se = (1. + (H/Hvg)^N )^(-M)
            θ = wrc.se.Se_2_θ(Se, θs, θr)
            return θ
        end


        function Se_2_H(Se, Hvg, N, Km=1.) # van Genuchten WRC
            M = 1. - Km / N
            H = Hvg * exp(log(exp(log(Se) / -M) - 1.) / N)
            return H
        end


        function DθDh(H, θs, θr, Hvg, N, Km)
            # The integration should be performed with H and not with Se  quadgk( H->  (,-Inf,0)[1]
            M = 1.- Km/N # van Genuchten

            DθDh = M*(θs-θr)/ (Hvg*(1-M))*(H/Hvg).^(N*M).*(1.+(H/Hvg).^N).^(-M-1.)
            return DθDh
            
            ## Alternative which does not work
            # H = Se_2_H(Se, Hvg, N, Km)
            # DθDha = M * (θs-θr) * (Se^inv(M)) * ((1. - (Se^inv(M)) )^ M) / (Hvg* (1.-M))
            # DθDh =(θs-θr)*M*N*(1./Hvg)*(H/Hvg).^(N-1.).*(1.+(H/Hvg).^N).^(-M-1.)
        end


    end #end module vg

end
