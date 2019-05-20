
module PSD
    using cst
    export  Particle_radius2h, Particle_radius2Pore_radius, Particle_radius2θ

    function Particle_radius2h(θs, Vi, Ri, Ms, ρp, τ) # soil matric potential of the i-th fraction
        h_Ri = Ri # define h_Ri

        for i in 1:length(Ri)
            h_Ri[i] = cst.Y*(( Ms/(θs*ρp)*sum(Vi.0 *Ri.^(-τ)) )^(1/(3-τ)))./Ri[i]
        end
        return h_Ri
    end

    function Particle_radius2Pore_radius(θs, Vi, Ri, Ms, ρp, τ) # pore radius
        return ri = Ri./ ( Ms/(θs*ρp)*sum(Vi.0 *Ri.^(-τ)) )^(1/(3-τ))
    end

    function Particle_radius2θ(θs, Vi, Ri, τ) # soil water content of the i-th fraction
        θ_Ri = Ri # define θ_Ri
        for i in 1:length(Ri)
            θ_Ri[i] = ( θs/(sum(Vi.0 *Ri.^(-τ))) )*sum(Vi[1:i].0 *Ri[1:i].^(-τ))
        end
        return θ_Ri
    end

end
