module ptype
    export  R2Ptype

    function R2Ptype(Vi, Ri, ϕ, θs, Ms, ρp, τ) # deriving Ptype from particle radius and volume fractions
        return Ptype = (6*(1-ϕ)/π)*( Ms/θs/ρp*sum(Vi.0 *Ri^(-τ)) )^(1/(3-τ))
    end

end
