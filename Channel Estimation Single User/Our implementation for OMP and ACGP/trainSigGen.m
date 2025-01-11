function S = trainSigGen(N, N_s)

% S = zeros(N_s, N+N_c-1);
    S_rp = exp(1j*(2*pi*rand(N_s, N)));
    S = 1/sqrt(N_s)*S_rp;
end