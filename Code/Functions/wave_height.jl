function wave_height(fetch,h,Uref,g)

    delta = h*g/Uref^2
    chi = fetch*g./Uref^2
    epsilon = 3.64*10^-3*(tanh(0.493*delta^0.75)*tanh(3.13*10^-3*chi^0.57/tanh(0.493*delta^0.75)))^1.74
    Hs = 4*sqrt(Uref^4*epsilon/g^2)
    
    return Hs

end