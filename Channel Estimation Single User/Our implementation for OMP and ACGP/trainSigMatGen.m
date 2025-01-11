function S = trainSigMatGen(N, N_c, N_s)
s_mat = trainSigGen(N, N_s);

% N-symbol long received signal
S = zeros(N, N_c*N_s);
s_mat_temp = s_mat.';
for d = 1:N_c
    S(:, (d-1)*N_s+1:d*N_s) = s_mat_temp;
    s_mat_temp = [zeros(1, N_s) ; s_mat_temp];
    s_mat_temp = s_mat_temp(1:N,:);
end
disp('')
