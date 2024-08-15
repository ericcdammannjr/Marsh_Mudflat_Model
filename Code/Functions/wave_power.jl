function wave_power(bf,Db,Uref,r,g)

    Hs = wave_height(bf,Db,Uref,g)
    Tp = wave_period(bf,Db,Uref,g)
    k = wave_number(1/Tp,Db,g)
    cg = 2*pi/k/Tp*0.5*(1+2*k*Db/(sinh(2*k*Db)))
    W = cg*9800/16*abs(Hs)^2

    return W
    
end