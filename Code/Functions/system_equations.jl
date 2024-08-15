function system_equations(du,u,p,t)
    
    # Model Parameters
    
    bf0 = p[1]
    Bpeak = p[2]
    Co = p[3]
    χref = p[4]
    Dmin = p[5]
    g = p[6]
    ka = p[7]
    ke = p[8]
    ko = p[9]
    L = p[10]
    λ = p[11]
    νGp = p[12]
    ϕ = p[13]
    r = p[14]
    R = p[15]
    ρ = p[16]
    ρo = p[17]
    Tω = p[18]
    τcr = p[19]
    Uref = p[20]
    ws = p[21]
    x = p[22]

    # State Variables

    bf = u[1]
    df = u[2]
    dm = u[3]

    # Width of the Marsh

    bm = L-bf

    # Organic Sediment Production
    
    Dmax = 0.7167*r-0.483
    AA = 0.25*(-Dmin-Dmax)*(Dmax-3*Dmin)
    B = Bpeak*(dm-Dmax)*(dm-Dmin)/AA
    
    if B <= 1e-3
    
       B = 0;
    
    end
    
    Bfrac = (B/Bpeak)
    AMC = (182.5)*B*(νGp)/(365*24*60*60)
    Rref = AMC*χref
    O = (1/ϕ)*(Rref/ρo)

    # Average Depths

    db = dm+(df-dm)*(1-exp(-x*0.1/df))
    Db = (db+(db-min(1,db/r)*r))/2
    Df = (df+(df-min(1,df/r)*r))/2
    Dm = (dm+(dm-min(1,dm/r)*r))/2

    # Sediment Concentration in the Lagoon
    
    τ1 = lagoon_bed_shear_stress(bf,Df,ko,g)
    S1 = max((τ1-τcr)/τcr,0)
    Cr = ρ*λ*S1/(1+λ*S1)

    # Sediment Concentration in the Marsh
    
    if Dm > 1e-4
        
        τ2 =  marsh_bed_shear_stress(bf,Dm,Uref,Bfrac,ko,g)
         
    else 
             
        τ2 = 0
         
    end

    S2 = max((τ2-τcr)/τcr,0)
    Cm = ρ*λ*S2/(1+λ*S2)

    # Sediment Flux Between the Lagoon and Marsh
    
    Fm = (Cr-Cm)*min(r,dm)/Tω
    
    # Sediment Flux Between the Open Ocean and Lagoon

    Fc = (Cr-Co)*min(r,df)/Tω 

    # Marsh Platform Erosion rate
    
    W = wave_power(bf,Db,Uref,r,g)
    Be = ke*W/(db-dm)
    
    # Marsh Platform Progradation Rate
    
    Ba = ka*Cr*ws/ρ

    # Morphodynamic Equations
    
    du[1] = Be-Ba                                       
    du[2] = -(Be-Ba)*(df-dm)/bf+(bm/bf)*Fm/ρ+Fc/ρ+R
    du[3] = -Fm/ρ-O+R

    nothing

end