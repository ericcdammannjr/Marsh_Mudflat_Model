using DifferentialEquations, Plots
pyplot()

include("Functions/behavior.jl")
include("Functions/lagoon_bed_shear_stress.jl")
include("Functions/lagoon_depth.jl")
include("Functions/marsh_bed_shear_stress.jl")
include("Functions/system_equations.jl")
include("Functions/wave_height.jl")
include("Functions/wave_number.jl")
include("Functions/wave_period.jl")
include("Functions/wave_power.jl")

# Computational Paramters

years = 150
t0 = 1
tf = years*365*24*60*60
tspan = (t0,tf)
tyears = range(t0,years,years)

# Model Paramters 

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
p = [bf0 Bpeak Co χref Dmin g ka ke ko L λ νGp ϕ r R ρ ρo Tω τcr Uref ws x]
prob_df0 = NonlinearProblem(lagoon_depth,[2],p)
df0 = solve(prob_df0, LevenbergMarquardt(), abstol=10^-28, reltol=10^-28, maxiters=10000)
dm0 = (0.7167*r-0.483)/2
u0 = [bf0 df0 dm0]

# Solution to System of ODEs

prob = ODEProblem(system_equations,u0,tspan,p)
sol = solve(prob,Rosenbrock23(),abstol=10^-6,reltol=10^-6,callback = behavior(),saveat = range(t0,tf,years))

# Solution Plot

plot([sol[1,:] sol[2,:] L.-sol[1,:] sol[3,:]],

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