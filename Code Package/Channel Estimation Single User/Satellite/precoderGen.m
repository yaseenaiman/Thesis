function [F, phase] = precoderGen(SP)

a_t = @(N_t, aod) 1/sqrt(N_t)*exp(1j*pi*(0:N_t-1).'*cos(aod));

maxPhase = 2*pi;
if SP.RFonly
    if SP.beamforming
        phase = (randi(2^SP.N_Q, SP.N_s, 1))/2^SP.N_Q * maxPhase;
        F_temp = a_t(SP.N_t, phase);
    else
        phase = (randi(2^SP.N_Q, SP.N_t, SP.N_s)-1)/2^SP.N_Q * maxPhase;
        F_temp = exp(1j*phase);
    end
    F = sqrt(SP.N_s)*F_temp/norm(F_temp, 'fro');
else
    if SP.beamforming
        phase = (randi(2^SP.N_Q, SP.L_t, 1))/2^SP.N_Q * maxPhase;
        F_RF = a_t(SP.N_t, phase);
    else
        phase = (randi(2^SP.N_Q, SP.N_t, SP.L_t)-1)/2^SP.N_Q * maxPhase;
        F_RF = exp(1j*phase);
    end
    F_BB_temp = randn(SP.L_t) + 1j*randn(SP.L_t);
    F_BB      = sqrt(SP.N_s)*F_BB_temp/norm(F_RF*F_BB_temp, 'fro');
    F         = F_RF*F_BB;
end