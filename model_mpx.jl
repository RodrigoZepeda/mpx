#https://github.com/epirecipes/sir-julia/blob/master/markdown/ode_turing/ode_turing.md
using DifferentialEquations, CSV, Turing, Plots, DataFrames, StatsPlots, 
        DataFramesMeta, Parameters, MCMCChains, Pipe, Random

Random.seed!(29890225)

# Open file and read cases 
@info "Loading data from Mexico City"
df = CSV.read("data/MPX_DGE.csv", DataFrame)
df = @pipe df |>
    filter(:Estado  => n -> n  == "Ciudad de México", _) |>
    select(_, :Fecha_reporte, :Casos, :Incidencia) |>
    sort(_, :Fecha_reporte)
 
cases = df[2:end,:Incidencia]

include("additional_functions/model.jl")

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
    p_detect   ~ Beta(2.0, 1.0)
    φ          ~ InverseGamma(2, 3)

    #Parameters for ODE
    params = (ρ = ρ, σ = σ, ν1 = 0, ν2 = 0, tvac = Inf, δ = δ, θ = θ, h = h, ϕ = ϕ)

    #Initial values
    cases_init = cases[1]/p_detect;
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
    Tot = I + PSI + PEI + 2*PII + PIR + PIV;
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
    μ   = max.(p_detect*(sol.u[2:end] - sol.u[1:(end-1)]), 1.e-10)
    
    #Negative Binomial cases
    #r = Number of successes; p = Probability of success. 
    #Hence mean = (1 - p)r/p and variance = (1 - p)r/p^2
    variance  = μ .+ μ.^2 ./ φ
    p_negbin  = μ ./ variance
    r_negbin  = μ .* p_negbin ./ (1 .- p_negbin)

    cases .~ NegativeBinomial.(r_negbin, p_negbin)
    
end;

nsim = 1000 #change to 100 for experiments

ode_nuts = sample(bayes_sir(cases), NUTS(), MCMCThreads(), nsim, 4; progress=true)

#Model description
diags = describe(ode_nuts)
CSV.write("diagnostics/rhat.csv",DataFrame(diags[1]))
CSV.write("diagnostics/GelmanRubin.csv",DataFrame(MCMCChains.gelmandiag(ode_nuts)))


tpop = 311000
muestra = DataFrame(sample(ode_nuts, nsim))
CSV.write("data/model_fitted_params.csv", muestra)

function prob_func(prob_result,i,repeat)
    cases_init = cases[1]/muestra[i,:p_detect]
    remake(prob_result; 
        p  = (ρ = muestra[i,:ρ], θ = muestra[i,:θ], ν1 = 0.0, ν2 = 0.0, tvac = Inf, σ = muestra[i,:σ], 
                δ = muestra[i,:δ], h = muestra[i,:h], ϕ = muestra[i,:ϕ]),
        u0 = [(tpop - cases_init)*muestra[i,:prop_ss], 0.0, cases_init*muestra[i,:prop_ii], 
                0.0, 0.0, (tpop - cases_init)*(1.0 - muestra[i,:prop_ss])/2.0, 0.0, 
                cases_init*(1.0 - muestra[i,:prop_ii])/2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, cases_init])
end


@info "Simulating posterior"
tspan       = (0.0, 2*length(cases))
y0          = [0.5*(tpop - cases[1]), tpop*cases[1],(1 - 0.5)*(tpop - cases[1])/2, 0, 0.5*cases[1]/2, 0, 0, 0, 0, cases[1]]    
prob_result = ODEProblem(mpx!, y0, tspan, params)
sol         = solve(EnsembleProblem(prob_result, prob_func = prob_func), Rodas4(), EnsembleSerial(), trajectories = nsim)

summ = EnsembleSummary(sol, 0:(1/7):2*length(cases),  quantiles=[0.025,0.5,0.975])

@info "Plotting"
Plots.plot(summ, title = "Viruela símica en México", 
            xlabel = "t", ylabel = "Casos acumulados")
StatsPlots.scatter!(2:length(cases), df[2:end,:Casos])


#Save ensemble EnsembleSummary
CSV.write("data/model_summary.csv", DataFrame(summ.med))
CSV.write("data/model_summary_qlow.csv", DataFrame(summ.qlow))
CSV.write("data/model_summary_qhigh.csv", DataFrame(summ.qhigh))


@info "Simulating vaccine scenarios"
nu = vcat(exp.(LinRange(log(0.01),log(1), 15)), 0.0)
for j in 1:length(nu)
    function prob_func(prob_result,i,repeat)
        cases_init = cases[1]/muestra[i,:p_detect]
        remake(prob_result; 
            p  = (ρ = muestra[i,:ρ], θ = muestra[i,:θ], ν1 = 0.0, ν2 = nu[j], tvac = 5, σ = muestra[i,:σ], 
                    δ = muestra[i,:δ], h = muestra[i,:h], ϕ = muestra[i,:ϕ]),
            u0 = [(tpop - cases_init)*muestra[i,:prop_ss], 0.0, cases_init*muestra[i,:prop_ii], 
                    0.0, 0.0, (tpop - cases_init)*(1.0 - muestra[i,:prop_ss])/2.0, 0.0, 
                    cases_init*(1.0 - muestra[i,:prop_ii])/2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, cases_init])
    end

    tspan       = (0.0, 2*length(cases))

    y0          = [0.5*(tpop - cases[1]), tpop*cases[1],(1 - 0.5)*(tpop - cases[1])/2, 0, 0.5*cases[1]/2, 0, 0, 0, 0, cases[1]]    

    prob_result = ODEProblem(mpx!, y0, tspan, params)

    sol         = solve(EnsembleProblem(prob_result, prob_func = prob_func), Rodas4(), EnsembleSerial(), trajectories = nsim)

    summ = EnsembleSummary(sol,0:1:length(cases),  quantiles=[0.025,0.5,0.975])

    #Save ensemble EnsembleSummary
    CSV.write("simulations/model_sims_point" * string(nu[j]) * ".csv", DataFrame(summ.med))
    CSV.write("simulations/model_sims_low" * string(nu[j]) * ".csv", DataFrame(summ.qlow))
    CSV.write("simulations/model_sims_high" * string(nu[j]) * ".csv", DataFrame(summ.qhigh))
end
