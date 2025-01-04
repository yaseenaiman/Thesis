function [Psi, A_tx, A_rx, DictAoD, DictAoA] = PsiMatGen(SP)

[A_tx, A_rx, DictAoD, DictAoA] = dictSteeringMatGen(SP.N_t, SP.N_r, SP.G_t, SP.G_r);
DictMat = kron(conj(A_tx), A_rx);
p_mat = pulseShapingVecGen(SP.G_c, SP.N_c, SP.beta);
Psi = zeros(SP.N_c*SP.N_r*SP.N_t, SP.G_c*SP.G_r*SP.G_t);

Psi = kron(eye(SP.N_c), DictMat);

disp('')


