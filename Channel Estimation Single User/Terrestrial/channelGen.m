function [H, A_t, A_r, AoD, AoA] = channelGen(SP)


    a_t = @(N_t, aod) 1/sqrt(N_t)*exp(-1j*(0:N_t-1).'*aod);
    a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(-1j*(0:N_r-1).'*aoa);

[~, ~, DictAoD, DictAoA] = dictSteeringMatGen(SP.N_t, SP.N_r, SP.G_t, SP.G_r);

    tau = (SP.N_c-1)*rand(SP.N_p, 1); % zeros(N_p, 1);
alpha = sqrt(1/2)*(randn(SP.N_p,1) + 1j*randn(SP.N_p, 1)); % channel gain


    AoD = 2*pi*rand(1, SP.N_p);
    AoA = 2*pi*rand(1, SP.N_p);

A_t = a_t(SP.N_t, AoD);
A_r = a_r(SP.N_r, AoA);
    
H = zeros(SP.N_r, SP.N_t, SP.N_c);
for d = 1:SP.N_c
    p     = RC_shaping(1, SP.beta, (d-1)-tau);
    Delta = diag(p.*alpha);
    H_temp = A_r*Delta*A_t';
    H(:,:,d) = sqrt(SP.N_r*SP.N_t/SP.N_p)*H_temp; %/norm(H_temp, 'fro');
    disp('')
end
