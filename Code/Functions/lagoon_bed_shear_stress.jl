function lagoon_bed_shear_stress(bf,Df,ko,g)

    Hs = wave_height(bf,Df,Uref,g)
    Tp = wave_period(bf,Df,Uref,g)
    k = wave_number(1/Tp,Df,g)
    Um = (pi*Hs/Tp/sinh(k*Df))
    aw = Tp*Um/(2*pi)
    fw = 0.4*(aw/ko)^-0.75
    τ = 1/2*1020*fw*Um^2

    return τ

end
