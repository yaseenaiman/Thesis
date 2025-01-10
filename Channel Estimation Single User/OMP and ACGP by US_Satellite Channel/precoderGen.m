function [F, phase] = precoderGen(SP)

a_t = @(N_t, aod) 1/sqrt(N_t)*exp(1j*pi*(0:N_t-1).'*cos(aod));

maxPhase = 2*pi;

        phase = (randi(2^SP.N_Q, SP.N_t, SP.N_s)-1)/2^SP.N_Q * maxPhase;
        F_temp = exp(1j*phase);
    
    F = sqrt(SP.N_s)*F_temp/norm(F_temp, 'fro');

end