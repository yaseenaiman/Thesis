function [W, phase] = combinerGen(SP)

a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(1j*pi*(0:N_r-1).'*cos(aoa));

maxPhase = 2*pi;
if SP.RFonly
    if SP.beamforming
        phase = (randi(2^SP.N_Q, SP.N_s, 1))/2^SP.N_Q * maxPhase;
        W_temp = a_r(SP.N_r, phase);
    else
        phase = (randi(2^SP.N_Q, SP.N_r, SP.N_s)-1)/2^SP.N_Q * maxPhase;
        W_temp = exp(1j*phase);
    end
    W = sqrt(1/SP.N_r)*W_temp;
else
    if SP.beamforming
        phase = (randi(2^SP.N_Q, SP.L_r, 1))/2^SP.N_Q * maxPhase;
        W_RF = a_r(SP.N_r, phase);
    else
        phase = (randi(2^SP.N_Q, SP.N_r, SP.L_r)-1)/2^SP.N_Q * maxPhase;
        W_RF = exp(1j*phase);
    end
    W_RF = sqrt(SP.N_r)*W_RF;
    W_BB = randn(SP.L_r) + 1j*randn(SP.L_r);
    W    = W_RF*W_BB;
end