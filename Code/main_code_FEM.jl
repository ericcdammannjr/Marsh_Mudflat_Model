using Plots
pyplot()

include("Functions/lagoon_bed_shear_stress.jl")
include("Functions/marsh_bed_shear_stress.jl")
include("Functions/wave_height.jl")
include("Functions/wave_number.jl")
include("Functions/wave_period.jl")
include("Functions/wave_power.jl")

# Computation Parameters

t0 = 0
tmax = 150
tmax_s = tmax*365*24*60*60
n = 100000
dt = (tmax_s-t0)/n
dtyrs = (tmax-t0)/n
t = t0:dtyrs:tmax

# Model Parameters

Bpeak = 2.500
Co = 0.1
χref = 0.15
Dmin = 0.0
g = 9.8
ka = 2.0
ke = 0.1/(365*24*60*60)
ko = 0.001
L = 10000
λ = 0.0001
νGp = 0.0138
ϕ = 0.377
r = 1.4
R = 5*(10^-3)/(365*24*60*60)
ρ = 1000
ρo = 1000
Tω = 12.5*(60*60)
τcr = 0.1
Uref = 8
ws = 0.5*10^-3
x = 10

# Initial Conditions

bf0 = 9000
df0 = 3.278
dm0 = (0.7167*r-0.483)/2
bm0 = L-bf0

# Vector Preallocation 

bf = zeros(n+1); bf[1] = bf0
df = zeros(n+1); df[1] = df0
dm = zeros(n+1); dm[1] = dm0
bm = zeros(n+1); bm[1] = bm0
dbfdt = zeros(n)
ddfdt = zeros(n)
ddmdt = zeros(n)
Ba = zeros(n)
Be = zeros(n)
O = zeros(n)
Rref = zeros(n)

for i in 1:n

    # Organic Sediment Production

    Dmax = 0.7167*r-0.483
    AA = 0.25*(-Dmin-Dmax)*(Dmax-3*Dmin)
    B = Bpeak*(dm[i]-Dmax)*(dm[i]-Dmin)/AA

    if B <= 1e-3

        B = 0;

    end

    Bfrac = (B/Bpeak)
    AMC = (182.5)*B*(νGp)/(365*24*60*60)
    Rref[i] = AMC*χref
    O[i] = (1/ϕ)*(Rref[i]/ρo)

    # Average Depths

    db = dm[i]+(df[i]-dm[i])*(1-exp(-x*0.1/df[i]))
    Db = (db+(db-min(1,db/r)*r))/2
    Df = (df[i]+(df[i]-min(1,df[i]/r)*r))/2
    Dm = (dm[i]+(dm[i]-min(1,dm[i]/r)*r))/2

    # Sediment Concentration in the Lagoon

    τ1 = lagoon_bed_shear_stress(bf[i],Df,ko,g)
    S1 = max((τ1-τcr)/τcr,0)
    Cr = ρ*λ*S1/(1+λ*S1)

    # Sediment Concentration in the Marsh

    if Dm > 1e-4
    
        τ2 =  marsh_bed_shear_stress(bf[i],Dm,Uref,Bfrac,ko,g)
     
    else 
         
        τ2 = 0
     
    end

    S2 = max((τ2-τcr)/τcr,0)
    Cm = ρ*λ*S2/(1+λ*S2)

    # Sediment Flux Between the Lagoon and Marsh

    Fm = (Cr-Cm)*min(r,dm[i])/Tω

    # Sediment Flux Between the Open Ocean and Lagoon

    Fc = (Cr-Co)*min(r,df[i])/Tω 

    # Marsh Platform Erosion rate

    W = wave_power(bf[i],Db,Uref,r,g)
    Be[i] = ke*W/(db-dm[i])

    # Marsh Platform Progradation Rate

    Ba[i] = ka*Cr*ws/ρ

    # Morphodynamic Equations

    dbfdt[i] = Be[i]-Ba[i]                                      
    ddfdt[i] = -(Be[i]-Ba[i])*(df[i]-dm[i])/bf[i]+(bm[i]/bf[i])*Fm/ρ+Fc/ρ+R
    ddmdt[i] = -Fm/ρ-O[i]+R

    # Foward Euler

    bf[i+1] = bf[i]+dt*dbfdt[i]
    df[i+1] = df[i]+dt*ddfdt[i]
    dm[i+1] = dm[i]+dt*ddmdt[i]
    bm[i+1] = L-bf[i+1]

end

# Plotting Results

plot(t,[bf[:] df[:] bm[:] dm[:]],

    layout=(2,2),
    title = ["Mudflat Width" "Mudflat Depth" "Marsh Width" "Depth of Marsh Below MHW"],
    xlabel = "Years",
    ylabel = "Meters",
    linecolor = :blue,
    linewidth = 0.5,
    titlefontsize = 10,
    xtickfontsize = 6,
    ytickfontsize = 8,
    legend = false,
    grid = false

)