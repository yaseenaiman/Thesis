function SP = top_sim_function_residual(SP)

[Psi, A_tx, A_rx, DictAoD, DictAoA] = PsiMatGen(SP); 
for snrInd = 1:length(SP.SNR_db_array)
    SP.rho = SP.SNR_lin_array(snrInd)*SP.sigma^2;
    % Simulation iteration
    for simInd = 1:SP.iter 
        fprintf('Iteration: %d \n', simInd)
        y_mat = zeros(SP.L_r*SP.N, SP.M);
        y_mat = mat2cell(y_mat, SP.L_r*SP.N, ones(1,SP.M));
        Phi = zeros(SP.N*SP.M*SP.L_r, SP.N_c*SP.N_r*SP.N_t);
        Phi = mat2cell(Phi, SP.N*SP.L_r*ones(1, SP.M), SP.N_c*SP.N_r*SP.N_t);
        receivedSigLen = SP.L_r*SP.N;
        PhiSize = SP.N*SP.L_r;
        [H, A_t, A_r, AoD, AoA] = channelGen(SP);
        % Frame iteration
        for frameInd = 1:SP.M 
            %% Receive vector generation
            [F, precoderPhase] = precoderGen(SP);
            [W, combinerPhase] = combinerGen(SP);
            S = trainSigMatGen(SP.N, SP.N_c, SP.N_s, SP.trainType);
            [y, y_withoutNoise] = receivedSigVec(SP, H, F, W, S);
            y_mat{frameInd} = y;
           
            %% Measurement matrix generation
            Phi_frame = kron(S*kron(eye(SP.N_c), F.'), W');
            Phi{frameInd} = ...
                sqrt(SP.rho)*Phi_frame;
        end
        y_mat = cell2mat(y_mat);
        y_vec = y_mat(:);
  y_quant = y_vec;
        for iif=1:SP.M
        if iif ==1
            y_quant(1:((SP.N_r*SP.L_r)-round(0.2*(SP.N_r*SP.L_r)))*iif) = quantization(y_vec(1:((SP.N_r*SP.L_r)-round(0.2*(SP.N_r*SP.L_r)))*iif), SP.ADCbits, max(abs([real(y_vec);imag(y_vec)])));
        else 
                    y_quant((iif-1)*SP.N_r*SP.L_r:((iif-1)*SP.N_r*SP.L_r+((SP.N_r*SP.L_r)-round(0.2*(SP.N_r*SP.L_r))))) = quantization(y_vec((iif-1)*SP.N_r*SP.L_r:((iif-1)*SP.N_r*SP.L_r+((SP.N_r*SP.L_r)-round(0.2*(SP.N_r*SP.L_r))))), SP.ADCbits, max(abs([real(y_vec);imag(y_vec)])));
    
            
        end
        end
%         y_quant = quantization(y_vec, SP.ADCbits, max(abs([real(y_vec);imag(y_vec)])));
     
        Phi = cell2mat(Phi);
        h = H(:);
        PhiPsi = Phi*Psi;
        
       % Author ACGP
 m1 = size(PhiPsi);
    m=m1(2); 
        [x_hat_LS, ~,~] = greed_nomp(y_quant, PhiPsi, m,'maxIter',10); 

        h_hat_LS = Psi*x_hat_LS;
        NMSE_LS(simInd, snrInd) = 10*log10(norm(h-h_hat_LS)^2/norm(h)^2);
        SP.NMSE_LS = NMSE_LS;
        mean(SP.NMSE_LS)
        


         disp('Starting OMP by Author')
        [x_hat_OMP_author] = OMPa(PhiPsi, y_quant, SP.N_p*9); 
        h_hat_OMP_author = Psi*x_hat_OMP_author;
        NMSE_OMP_author(simInd, snrInd) = 10*log10(norm(h-h_hat_OMP_author)^2/norm(h)^2);
        SP.NMSE_author= NMSE_OMP_author;
        mean(SP.NMSE_author)
         
        
 disp('Starting ACGP by Yaseen')
 [x_ACGP] = acgp_yas(y_quant, PhiPsi, 10); 
x_hat_NOMP = x_ACGP;
        h_hat_NOMP = Psi*x_hat_NOMP;
        NMSE_NOMP(simInd, snrInd) = 10*log10(norm(h-h_hat_NOMP)^2/norm(h)^2);
        SP.NMSE_NOMP= NMSE_NOMP;
        mean(SP.NMSE_NOMP)

        % OMP
 disp('Starting OMP by Yaseen')
  [x_OMP] = omp_yas(y_quant, PhiPsi, SP.N_p*9); 
 x_hat_OMP = x_OMP;
h_hat_OMP = Psi*x_hat_OMP;
        NMSE_OMP(simInd, snrInd) = 10*log10(norm(h-h_hat_OMP)^2/norm(h)^2);
 SP.NMSE_OMP = NMSE_OMP;
        
        
        % MMSE
disp('Starting MMSE')
% Assuming noise variance is sigma^2
noise_var = SP.sigma^2;
% Identity matrix of size (number of columns in PhiPsi)
I = eye(size(PhiPsi, 2));
% MMSE estimation
x_hat_MMSE = (PhiPsi' * PhiPsi + noise_var * I) \ (PhiPsi' * y_quant);
h_hat_MMSE = Psi * x_hat_MMSE;
NMSE_MMSE(simInd, snrInd) = 10 * log10(norm(h - h_hat_MMSE)^2 / norm(h)^2);
SP.NMSE_MMSE = NMSE_MMSE;
mean(SP.NMSE_MMSE)

    end
  

end
SP.LS = mean(SP.NMSE_LS);
SP.OMP = mean(SP.NMSE_OMP);
SP.NOMP = mean(SP.NMSE_NOMP);
SP.Omp_author = mean(SP.NMSE_author);
SP.MMSE = mean(SP.NMSE_MMSE);
% Plotting the rates
plot(SP.SNR_db_array ,SP.LS,'b','LineWidth',1.5);
 hold on, plot(SP.SNR_db_array ,SP.OMP,'r -o','LineWidth',1.5);
hold on; plot(SP.SNR_db_array ,SP.NOMP,'g -o','LineWidth',1.5);
plot(SP.SNR_db_array, SP.MMSE, 'm', 'LineWidth', 1.5); % Add MMSE 
plot(SP.SNR_db_array, SP.Omp_author, 'black -o', 'LineWidth', 1.5); % Add OMP by Author

set(gcf,'color','w');
set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)
    
legend('ACGP by Author','OMP by Yaseen', 'ACGP by Yaseen', 'MMSE','OMP by Author')
xlabel('SNR dB','FontSize',12);grid;
ylabel('MMSE dB','FontSize',12);



