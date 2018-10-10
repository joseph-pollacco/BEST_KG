module porosity
    export  FUNC_Φ, FUNC_Φ2θs

    function FUNC_Φ(ρb, ρp=2.65) # ρb and ρp in (g cm-3)
      return Φ = 1 - ρb/ρp
    end

    function FUNC_Φ2θs(Φ, P_Φ=0.95) # correcting factor P_Φ[-]
      return θs = Φ * P_Φ
    end

end
