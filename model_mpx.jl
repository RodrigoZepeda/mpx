#https://github.com/epirecipes/sir-julia/blob/master/markdown/ode_turing/ode_turing.md
using DifferentialEquations, CSV, Turing, Plots, DataFrames, StatsPlots, 
        DataFramesMeta, Parameters, MCMCChains, Pipe

# Open file and read cases 
@info "Loading data from Mexico City"
df = CSV.read("data/MPX_DGE.csv", DataFrame)
df = @pipe df |>
    filter(:Estado  => n -> n  == "Ciudad de México", _) |>
    select(_, :Fecha_reporte, :Casos, :Incidencia) |>
    sort(_, :Fecha_reporte)
 
cases = df[2:end,:Incidencia]

# Establish model (SEIR)
function mpx!(dy, y, params, t)
    
    #Parameters
    #--------------------------------------------------
    @unpack ρ, ν, σ, θ, δ, h, ϕ = params;

    #Model compartments
    #--------------------------------------------------
    S   = y[1];   #Susceptible
    E   = y[2];   #Exposed
    I   = y[3];   #Infected
    R   = y[4];   #Recovered
    V   = y[5];   #Vaccinated
    PSS = y[6];   #Susceptible-Susceptible pairing
    PSE = y[7];   #Susceptible-Exposed pairing
    PSI = y[8];   #Susceptible-Infected pairing
    PSR = y[9];   #Susceptible-Recovered pairing
    PSV = y[10];  #Susceptible-Vaccinated pairing
    PEE = y[11];  #Exposed-Exposed pairing
    PEI = y[12];  #Exposed-Infected pairing
    PER = y[13];  #Exposed-Recovered pairing
    PEV = y[14];  #Exposed-Vaccinated pairing
    PII = y[15];  #Infected-Exposed pairing
    PIR = y[16];  #Infected-Recovered pairing
    PIV = y[17];  #Infected-Vaccinated pairing
    PRR = y[18];  #Recovered-Recovered pairing
    PRV = y[19];  #Recovered-Vaccinated pairing
    PVV = y[20];  #Vaccinated-Vaccinated pairing
    Tot = y[21];  #Total infected (cummulative). 
    N   = sum(y[1:5]) + 2*sum(y[6:(end-1)]); 

    #Model
    #--------------------------------------------------
    @inbounds begin
        dy[1]  = -(ρ + ν)*S + σ*(2*PSS + PSE + PSI + PSR + PSV);          #dS
        dy[2]  = -(ρ + θ)*E + σ*(PSE + 2*PEE + PEI + PER + PEV);          #dE
        dy[3]  = -(ρ + δ)*I + θ*E + σ*(PSI + PEI + 2*PII + PIR + PIV);    #dI
        dy[4]  = -ρ*R + δ*I + σ*(PSR + PER + PIR + 2*PRR + PRV);          #dR
        dy[5]  = -ρ*V + ν*S + σ*(PSV + PEV + PIV + PRV + 2*PVV);          #dV
        dy[6]  = 0.5*ρ*S^2/N - (σ + 2*ν)*PSS;                             #dPSS
        dy[7]  = ρ*S*E/N - (σ + θ + ν)*PSE;                               #dPSE
        dy[8]  = ρ*(1 - h)*S*I/N + θ*PSE - (σ + ϕ*h + δ + ν)*PSI;         #dPSI
        dy[9]  = ρ*S*R/N + δ*PSI - (σ + ν)*PSR;                           #dPSR
        dy[10] = ρ*S*V/N  + ν*PSS - (σ + ν)*PSV;                          #dPSV
        dy[11] = 0.5*ρ*E^2/N - (σ + 2*θ)*PEE;                             #dPEE
        dy[12] = ρ*E*I/N + ρ*h*S*I/N + ϕ*h*PSI + θ*PEE - (σ + θ + δ)*PEI; #dPEI
        dy[13] = ρ*E*R/N + δ*PEI - (σ + θ)*PER;                           #dPER
        dy[14] = ρ*E*V/N + ν*PSE - (σ + θ)*PEV;                           #dPEV
        dy[15] = 0.5*ρ*I^2/N + θ*PEI - (σ + 2*δ)*PII;                     #dPII
        dy[16] = ρ*I*R/N + δ*PII + θ*PER - (σ + δ)*PIR;                   #dPIR
        dy[17] = ρ*I*V/N + θ*PEV + ν*PSI - (σ + δ)*PIV;                   #dPIV
        dy[18] = 0.5*ρ*R^2/N + δ*PIR - σ*PRR;                             #dPRR 
        dy[19] = ρ*R*V/N + δ*PIV + ν*PSR - σ*PRV;                         #dPRV
        dy[20] = 0.5*ρ*V^2/N + ν*PSV - σ*PVV;                             #dPVV
        dy[21] = θ*(E + PSE + PEE + PEI + PER + PEV);                     #dCummulative
    end
    return nothing
end;

@model bayes_sir(cases) = begin
    
    #Total population
    tpop = 311000

    # Calculate number of timepoints
    l     = length(cases)    
    tspan = (0.0, float(l))

    #Parameters
    ρ          ~ truncated(Normal(0.0, 10.0), lower = 0.0, upper = Inf)
    σ          ~ truncated(Normal(0.0, 10.0), lower = 0.0, upper = Inf)
    θ          ~ truncated(Normal(0.0, 1.0), lower = 0.0, upper = 1.0)
    δ          ~ truncated(Normal(0.0, 1.0), lower = 0.0, upper = 1.0)
    h          ~ truncated(Normal(0.0, 1.0), lower = 0.0, upper = 1.0)
    ϕ          ~ truncated(Normal(0.0, 1.0), lower = 0.0, upper = 1.0)
    prop_ss    ~ Beta(1.0, 1.0)
    prop_ii    ~ Beta(1.0, 1.0)
    #p_detect   ~ Beta(2.0, 1.0)
    φ          ~ InverseGamma(2, 3)

    #Parameters for ODE
    params = (ρ = ρ, σ = σ, ν = 0, δ = δ, θ = θ, h = h, ϕ = ϕ)

    #Initial values
    cases_init = cases[1]
    S   = (tpop - cases_init)*prop_ss; 
    E   = 0.0
    I   = cases_init*prop_ii; 
    R   = 0.0; 
    V   = 0.0; 
    PSS = (tpop - cases_init)*(1.0 - prop_ss)/2.0; 
    PSE = 0.0; 
    PSI = cases_init*(1.0 - prop_ii)/2.0; 
    PSR = 0.0; 
    PSV = 0.0; 
    PEE = 0.0; 
    PEI = 0.0; 
    PER = 0.0; 
    PEV = 0.0;
    PII = 0.0; 
    PIR = 0.0; 
    PIV = 0.0; 
    PRR = 0.0; 
    PRV = 0.0; 
    PVV = 0.0; 
    Tot = dI + dPSI + dPEI + 2*dPII + dPIR + dPIV;
    y0  = [S, E, I, R, V, PSS, PSE, PSI, PSR, PSV, PEE, PEI, PER, PEV, PII, PIR, PIV, PRR, PRV, PVV, Tot]
    
    #ODE Solution
    prob = ODEProblem(mpx!, y0, tspan, params)
    sol  = solve(prob, Rodas4(); p=params, saveat=1.0, save_idxs = length(y0))

    #Break inference 
    #https://turing.ml/dev/docs/using-turing/advanced#update-the-accumulated-log-probability-in-the-model-definition
    if !(sol.retcode == :Success)
        @Turing.addlogprob! -Inf
        return
    end

    #Get incidence
    μ   = max.(1.0*(sol.u[2:end] - sol.u[1:(end-1)]), 1.e-10)
    
    #Negative Binomial cases
    #r = Number of successes; p = Probability of success. 
    #Hence mean = (1 - p)r/p and variance = (1 - p)r/p^2
    variance  = μ .+ μ.^2 ./ φ
    p_negbin  = μ ./ variance
    r_negbin  = μ .* p_negbin ./ (1 .- p_negbin)

    cases .~ NegativeBinomial.(r_negbin, p_negbin)
    
end;

ode_nuts = sample(bayes_sir(cases), NUTS(), MCMCThreads(), 200, 4; progress=true)

#Model description
diags = describe(ode_nuts)
print(DataFrame(diags[1]))
print(DataFrame(MCMCChains.gelmandiag(ode_nuts)))

nsim = 1000
tpop = 311000
muestra = DataFrame(sample(ode_nuts, nsim))
function prob_func(prob_result,i,repeat)
    cases_init = cases[1]
    remake(prob_result; 
        p  = (ρ = muestra[i,:ρ], θ = muestra[i,:θ], ν = 0.0, σ = muestra[i,:σ], 
                δ = muestra[i,:δ], h = muestra[i,:h], ϕ = muestra[i,:ϕ]),
        u0 = [(tpop - cases_init)*muestra[i,:prop_ss], 0.0, cases_init*muestra[i,:prop_ii], 
                0.0, 0.0, (tpop - cases_init)*(1.0 - muestra[i,:prop_ss])/2.0, 0.0, 
                cases_init*(1.0 - muestra[i,:prop_ii])/2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, cases_init])
end


@info "Simulating posterior"
tspan       = (0.0, length(cases))
y0          = [0.5*(tpop - cases[1]), tpop*cases[1],(1 - 0.5)*(tpop - cases[1])/2, 0, 0.5*cases[1]/2, 0, 0, 0, 0, cases[1]]    
prob_result = ODEProblem(mpx!, y0, tspan, params)
sol         = solve(EnsembleProblem(prob_result, prob_func = prob_func), Rodas4(), EnsembleSerial(), trajectories = nsim, save_idxs = 21)

summ = EnsembleSummary(sol,quantiles=[0.025,0.5,0.975])

@info "Plotting"
Plots.plot(summ, title = "Viruela símica en México", xlabel = "t", ylabel = "Casos acumulados")
StatsPlots.scatter!(2:length(cases), df[2:end,:Casos])


sol         = solve(EnsembleProblem(prob_result, prob_func = prob_func), Rodas4(), EnsembleSerial(), trajectories = nsim, save_idxs = 21)

summ = EnsembleSummary(sol,quantiles=[0.025,0.5,0.975])

@info "Plotting"
Plots.plot(summ, title = "Viruela símica en México", xlabel = "t", ylabel = "Casos acumulados")

