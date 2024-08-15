function marsh_bed_shear_stress(bf,Dm,Uref,Bfrac,ko,g)

    if Bfrac == 0  
    
        Hs = wave_height(bf,Dm,Uref,g)
        Tp = wave_period(bf,Dm,Uref,g)
        k = wave_number(1/Tp,Dm,g)
        Um = (pi*Hs/Tp/sinh(k*Dm))
        aw = Tp*Um/(2*pi)
        fw = 0.4*(aw/ko)^-0.75
        τ = 1/2*1020*fw*Um^2
    
    else
    
        τ = 0;
    
    end

    return τ

end