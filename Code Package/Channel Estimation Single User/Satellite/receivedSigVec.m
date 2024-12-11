function [y, y_withoutNoise] = receivedSigVec(SP, H, F, W, S)

v_mat = SP.sigma*sqrt(1/2)*(randn(SP.N_r, SP.N) + 1j*randn(SP.N_r, SP.N));
H_2d_mat = reshape(H, SP.N_r, []);
r_mat = sqrt(SP.rho)*H_2d_mat*kron(eye(SP.N_c), F)*S.' + v_mat;
y_mat = W'*r_mat;
y = y_mat(:); % measurement vector

y_withoutNoise = W'*sqrt(SP.rho)*H_2d_mat*kron(eye(SP.N_c), F)*S.';

