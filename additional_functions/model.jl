# Establish model (SEIR)
function mpx!(dy, y, params, t)
    
    #Parameters
    #--------------------------------------------------
    @unpack ρ, ν1, ν2, tvac, σ, θ, δ, h, ϕ = params;

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

    if (t < tvac)
        ν = ν1;
    else
        ν = ν2;
    end

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