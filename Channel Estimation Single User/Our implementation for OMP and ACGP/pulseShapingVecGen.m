function p_mat = pulseShapingVecGen(G_c, N_c, RC_beta)

p_mat = zeros(G_c, N_c);

for n=1:G_c
    for d=1:N_c
        p_mat(n,d) = RC_shaping(1, RC_beta, (d-1)-(n-1)*N_c/G_c);
    end
end