function [x_ACGP]= acgp_yas(x,PhiPsi,m)

r_ACGP = x;
x_ACGP = zeros(size(PhiPsi,2), 1);
iterations =m;
Gamma_ACDP= []; % Index set
            d_acgp = zeros(size(PhiPsi,2), 1);

for nn=1:iterations
    % ACGP step
    g = PhiPsi' * r_ACGP; % Gradient
    [~, idx_a] = max(abs(g)); % Select column
    Gamma_ACDP = union(Gamma_ACDP, idx_a); % Update the active set
 if nn == 1
            d_acgp = zeros(size(PhiPsi,2), 1);
            d_acgp(Gamma_ACDP) = g(Gamma_ACDP); % Unit vector for the selected index

 else
            % Step 6: Compute b
                c1 = PhiPsi(:, Gamma_ACDP) * d_acgp(Gamma_ACDP); % c vector

            b = -(c1' * PhiPsi(:, Gamma_ACDP)*g(Gamma_ACDP)) / (c1' * c1);

            % Step 7: Update the direction
            d_acgp(Gamma_ACDP) = g(Gamma_ACDP) + b * d_acgp(Gamma_ACDP);
        end
    c = PhiPsi(:, Gamma_ACDP) * d_acgp(Gamma_ACDP); % c vector
    alpha = (r_ACGP' * c) / (c' * c); % Step size
    x_ACGP(Gamma_ACDP) = x_ACGP(Gamma_ACDP) + alpha * d_acgp(Gamma_ACDP); % Update solution
    r_ACGP = r_ACGP - alpha * c; % Update residual
end
end