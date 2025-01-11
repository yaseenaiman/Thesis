function [x_OMP]= omp_yas(x,PhiPsi,m)

% Initialization
r_OMP = x;
y_quant=x;
x_OMP = zeros(size(PhiPsi,2), 1);
iterations =m;
    Gamma = []; % Index set 

for n = 1:iterations
      % OMP step
    [~, idx] = max(abs(PhiPsi' * r_OMP)); % Select column
     Gamma = union(Gamma, idx); % Update the active set

        % Step 3: Extract the submatrix of Psi corresponding to Gamma
        Psi_Gamma = PhiPsi(:, Gamma);
    
     % Step 4: Solve the least squares problem
        x_OMP(Gamma) = (Psi_Gamma' * Psi_Gamma) \ (Psi_Gamma' * y_quant);
            r_OMP = y_quant - Psi_Gamma * x_OMP(Gamma); % Update residual
    

end