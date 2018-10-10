module shape_param  # computes the shape paramter from Kr_θini
    using SpecialFunctions
    export  Kr_Se_ini2σ, Kr_Se_ini2η

    function Kr_Se_ini2σ(Se_ini, Kr_θini) # Kosugi shape param σ
        return σ = sqrt(2) * ( erfcinv(2*sqrt(Kr_θini/sqrt(Se_ini))) - erfcinv(2*Se_ini) )
    end

    function Kr_Se_ini2η(Se_ini, Kr_θini) # Brooks and Corey shape param η
        return η = log(Kr_θini)/log(Se_ini)
    end
end
