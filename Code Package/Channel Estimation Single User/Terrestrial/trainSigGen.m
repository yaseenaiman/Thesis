function S = trainSigGen(N, N_c, N_s, type)

% S = zeros(N_s, N+N_c-1);

if strcmp(type, 'Hadamard')
    S_had = hadamard(N);
    S_had = 1/sqrt(N_s)*S_had(2:N_s+1, :);
    S = S_had;
elseif strcmp(type, 'ZC')
    S_ZC = zeros(N_s, N);
    zcSeq = 1/sqrt(N_s)*ZCgen(N);
    for i = 1:N_s
        S_ZC(i, :) = circshift(zcSeq, i);
    end
    S = S_ZC;
elseif strcmp(type, 'random')
    S_rand = sign(randn(N_s, N)) + 1j*sign(randn(N_s,N));
    S = 1/sqrt(2*N_s)*S_rand;
elseif strcmp(type, 'random phase')
    S_rp = exp(1j*(2*pi*rand(N_s, N)));
    S = 1/sqrt(N_s)*S_rp;
end