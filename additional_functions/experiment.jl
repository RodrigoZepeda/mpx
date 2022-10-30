#https://github.com/epirecipes/sir-julia/blob/master/markdown/ode_turing/ode_turing.md
using DifferentialEquations, CSV, Turing, Plots, DataFrames, StatsPlots, 
        DataFramesMeta, Parameters, MCMCChains, Pipe

#Model
include("model.jl")

#Initial conditions
dS   = 186582; dE = 0.0; dI = 18; dR = 0.0; dV = 0.0; dPSS = 124388/2; dPSE = 0.0; 
dPSI = 12; dPSR = 0.0; dPSV = 0.0; dPEE = 0.0; dPEI = 0.0; dPER = 0.0; dPEV = 0.0;
dPII = 2; dPIR = 0.0; dPIV = 0.0; dPRR = 0.0; dPRV = 0.0; dPVV = 0.0; 

y0     = Float64[dS, dE, dI, dR, dV, dPSS, dPSE, dPSI, dPSR, dPSV, dPEE, dPEI, dPER, dPEV, dPII, dPIR, dPIV, dPRR, dPRV, dPVV, dI + dPSI + dPEI + 2*dPII + dPIR + dPIV];
N      = sum(y0) + sum(y0[6:(end-1)]);

#Example of solving only one problem
@info "Solving"
tspan = (0, 20)
model_parameters = (ρ = 3.9456, ν1 = 0.0,  ν2 = 0.1, tvac = Inf, σ = 19.26, δ = 0.67, θ = 0.79, h = 0.85, ϕ = 0.46) #rho, sigma, delta, h, phi
prob_result = ODEProblem(mpx!, y0, tspan, model_parameters)
sol         = solve(prob_result, save_idxs = [21]) #Chose idx 3,8,15,16 
Plots.plot(sol, title = "Modelo SEQIR", label = "Total I(t)", xlabel = "t")

model_parameters = (ρ = 3.9456, ν1 = 0.0,  ν2 = 0.1, tvac = 5, σ = 19.26, δ = 0.67, θ = 0.79, h = 0.85, ϕ = 0.46) #rho, sigma, delta, h, phi
prob_result = ODEProblem(mpx!, y0, tspan, model_parameters)
sol         = solve(prob_result, save_idxs = [21]) #Chose idx 3,8,15,16 
Plots.plot!(sol, title = "Modelo SEQIR", label = "Total I(t)", xlabel = "t")


# Example of solving multiple simulations
ρ   = rand(truncated(Normal(30.0, 10.0), lower = 1.0, upper = Inf), nsim)
σ   = rand(truncated(Normal(15.0, 10.0), lower = 1.0, upper = Inf), nsim)
θ   = rand(truncated(Normal(0.0, 1.0), lower = 0.0, upper = 1.0), nsim)
δ   = rand(truncated(Normal(0.0, 1.0), lower = 0.0, upper = 1.0), nsim)
ϕ   = rand(truncated(Normal(0.0, 1.0), lower = 0.0, upper = 1.0), nsim)
h   = rand(Beta(1,1), nsim)
function prob_func(prob_result,i,repeat)
    remake(prob_result, p = (ρ = ρ[i], σ = σ[i], ν1 = 0.0, δ = δ[i], θ = θ[i], h = h[i], ϕ = ϕ[i]))
end


@info "Simulating posterior"
tspan       = (0.0, 10)
prob_result = ODEProblem(mpx!, y0, tspan)
sol         = solve(EnsembleProblem(prob_result, prob_func = prob_func), Rodas4(), EnsembleSerial(), trajectories = nsim, save_idxs = 21)

@info "Plotting"
Plots.plot(sol, title = "Modelo SEQIR", label = "Total I(t)", xlabel = "t")


