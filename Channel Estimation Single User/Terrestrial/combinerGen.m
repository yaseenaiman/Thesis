function [W, phase] = combinerGen(SP)

a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(1j*pi*(0:N_r-1).'*cos(aoa));

maxPhase = 2*pi;

        phase = (randi(2^SP.N_Q, SP.N_r, SP.N_s)-1)/2^SP.N_Q * maxPhase;
        W_temp = exp(1j*phase);
    W = sqrt(1/SP.N_r)*W_temp;

end