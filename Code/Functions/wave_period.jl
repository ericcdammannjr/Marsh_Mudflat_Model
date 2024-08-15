function wave_period(fetch,h,Uref,g)

    delta = h*g/Uref^2
    chi = fetch*g./Uref^2
    ni = 0.133*(tanh(0.331*delta^1.01)*tanh(5.215*10^-4*chi^0.73/tanh(0.331*delta^1.01)))^-0.37
    Tp = Uref/ni/g

    return Tp

end