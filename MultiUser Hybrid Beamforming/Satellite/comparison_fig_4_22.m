
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear all;clc;warning('off');
% ----------------------------- System Parameters -------------------------
Num_users=8; % Number of users
TX_ant=64; %Number of UPA TX antennas
TX_ant_w=sqrt(TX_ant); % width of the planar array
TX_ant_h=sqrt(TX_ant); % hieght of the planar array 

NtRf = 16;%number of transmitter RF chains
RX_ant=16; %Number of UPA RX antennas
RX_ant_w=sqrt(RX_ant); % width 
RX_ant_h=sqrt(RX_ant); % hieght

% ----------------------------- Channel Parameters ------------------------
Num_paths=5; %Number of channel paths
fc=28e9; % Carrier Frequency 
Re=6378*1000; %Earth Radius
Hight = 560*1000; % Satellite Altitude
Rs = Re+Hight; %Total distance between the center of the earth and the satellite

% ----------------------------- Simulation Parameters ---------------------
SNR_dB_range=-20:5:20;  % SNR in dB
Rate_SU=zeros(1,length(SNR_dB_range)); % single user spectral efficiency (no-interference)
Rate_BS=zeros(1,length(SNR_dB_range));% analog-only beamsteering spectral efficiency
Rate_HP=zeros(1,length(SNR_dB_range)); % hybrid ZF precoding spectral efficiency
Rate_HP_MSE=zeros(1,length(SNR_dB_range)); %  hybrid MMSE precoder spectral efficiency 
Rate_HP_Kalman=zeros(1,length(SNR_dB_range)); % hybrid Kalman precoder spectral efficiency
Rate_HP_FD_ZF=zeros(1,length(SNR_dB_range)); % fully digital hybrid ZF precoding spectral efficiency
Rate_HP_FD_MSE=zeros(1,length(SNR_dB_range)); % fully digital hybrid MMSE precoding spectral efficiency
Rate_HP_ayash = zeros(1,length(SNR_dB_range)); %  OMP_Ayash Spectral Efficiency
Rate_HP_Yu = zeros(1,length(SNR_dB_range)); %  Novel Algorithm Spectral Efficiency 
Rate_HP_Original = zeros(1,length(SNR_dB_range)); % PE-AltMin Spectral Efficiency
ITER=500; % Number of iterations
    


%% RF codebook (Design Space is Steering_Vec)
% design space parameters 
Kd=pi;  % Assuming: K=2pi/lambda, D=lambda/2
Num_Qbits = 6; %Number of quantization bits
Num_Directions=2^Num_Qbits; %Possible directions
Step=2*pi/Num_Directions;
Antennas_index=0:1:TX_ant-1;
Theta_Quantized=0:Step:2*pi-Step;

for i=1:1:length(Theta_Quantized)
    Steering_Vec(:,i)=sqrt(1/TX_ant)*exp(1j*Antennas_index*Theta_Quantized(i));
end


%% --------------- Simulation starts ---------------------------------------
for iter=1:1:ITER
    
    % Generate user channels 
    [H,a_TX,a_RX]=generate_channels(Num_users,TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths, fc, Hight, Re, Rs); 
    % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions

        
    % Intialization of the hybrid precoder and combiner
    Frf=zeros(TX_ant,Num_users); % BS RF precoders 
    Wrf=zeros(RX_ant,Num_users); % MS RF precoders 
    
    for u=1:1:Num_users
        Frf(:,u)=a_TX(:,u); %Ideal FRF, assuming AODs are known
        Wrf(:,u)=a_RX(:,u); %Ideal WRF, assuming AOAs are known
    end      
    
  
    
    %% Effective Channel Assuming Wrf (Combiners of users)is known or  AOAs are perfectly known 
    % when performing beamforming at the transmitter end. 
    
 for u=1:1:Num_users
        Channel=zeros(RX_ant,TX_ant);
        Channel(:,:)= H(u,:,:);
        He_fd(u,:)=Wrf(:,u)'*Channel;     
 end
    
    
 %% Fully-digital zero-forcing precoding to minimize interference amongst users
    
 Ffd=He_fd'*pinv(He_fd*He_fd'); %Fully digital beamforming
    for u=1:1:Num_users % Normalization of the fully digital precoders
        Ffd(:,u)=Ffd(:,u)/sqrt((Ffd(:,u))'*(Ffd(:,u)));
    end 
    %% Spectral efficiency calculations
    count=0;
    for SNR_dB=SNR_dB_range
        count=count+1;
        SNR=10^(.1*SNR_dB)/Num_users; % SNR value
        
        sigma2=1/SNR;
        
        
        
 %% Fully Digital MMSE, it is another technique to minimize 
 % the interference amongst users
 FfdMSE=inv(He_fd'*He_fd+Num_users*sigma2*eye(TX_ant))*He_fd';
 for u=1:1:Num_users % Normalization of the precoder
 FfdMSE(:,u)=FfdMSE(:,u)/sqrt((FfdMSE(:,u))'*(FfdMSE(:,u)));
 end
    
    
    %%  Maximum-Ratio Transmission (MRT) for interference cancellation
    
 FfdMRT=He_fd';
 for u=1:1:Num_users % Normalization of the  precoder
 FfdMRT(:,u)=FfdMRT(:,u)/sqrt((FfdMRT(:,u))'*(FfdMRT(:,u)));
 end
 %% 
% Fopt=  Ffd; % Zero Forcing
Fopt=  FfdMRT; %  We are relying on MRT approach to minimize 
 % the interference amongs users, then, we are implementing 
 % various Hybrid Beamforming algorithms using it. Since the user does not
 % have the ability to cancel the interference, this is carried out at the
 % Satellite station side
% Fopt=  FfdMSE; % MMSE
%% Novel Algorithm Combining OMP by Ayash and PE AltMin by Xianghao Yu

[FRF_Novel, FBB_Novel] = Novel_Algorithm(Fopt, NtRf, Steering_Vec);
FBB_Novel = sqrt(Num_users) * FBB_Novel / norm(FRF_Novel * FBB_Novel,'fro'); %Normalization
F_Novel =FRF_Novel*FBB_Novel;


%% MMSE baseband precoder, Hybrid Beamforming
FRFMSE = DiscreteRFGeneration(Fopt,Steering_Vec,NtRf); %Generating FRF Using Fully digital MMSE Fopt
FbbMSE=inv((He_fd*FRFMSE)'*(He_fd*FRFMSE)+Num_users*sigma2*FRFMSE'*FRFMSE)*(He_fd*FRFMSE)';
       
        for u=1:1:Num_users % Normalization of the hybrid precoders
            FbbMSE(:,u)=FbbMSE(:,u)/sqrt((FRFMSE*FbbMSE(:,u))'*(FRFMSE*FbbMSE(:,u)));
        end
        
        
%% Baseband zero-forcing precoding, Hybrid Beamforming 
 FRF_ZF = DiscreteRFGeneration(Ffd,Steering_Vec,NtRf); %Generation of FRF depending on ZF fully digital Fopt
 Fbb=(He_fd*FRF_ZF)'*((He_fd*FRF_ZF)*(He_fd*FRF_ZF)')^(-1); 
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb(:,u)=Fbb(:,u)/sqrt((FRF_ZF*Fbb(:,u))'*(FRF_ZF*Fbb(:,u)));
    end
    
    
    %% Kalman Beamforming
       
 FRFKalman = DiscreteRFGeneration(Fopt,Steering_Vec,NtRf); %Depnding on Fully digital MMSE Fopt

  Fbbk=eye(NtRf,Num_users);
  RN1=Fbbk*Fbbk';
  Qm=eye(Num_users)*sigma2;
  ITERK=10; % Number of Kalman iterations
        
  for ii=1:ITERK
    Hk=He_fd*FRFKalman;
    K=RN1*Hk'*pinv(Hk*RN1*Hk'+Qm);
    errk=(eye(Num_users)-Hk*Fbbk);
    errk=errk/norm(errk);
    Fbbk=Fbbk+K*errk; 
    RN=RN1-K*Hk*RN1;
    RN1=RN; 
  end
        
        for u=1:1:Num_users % Normalization of the hybrid precoders
            Fbbk(:,u)=Fbbk(:,u)/sqrt((FRFKalman*Fbbk(:,u))'*(FRFKalman*Fbbk(:,u)));
        end
       Fkalman= FRFKalman*Fbbk;

    %%  OMP Hybrid Beamforming using Fopt
    Frf_ayash=[];
    Fres_ayash = Fopt;
    At = Steering_Vec;
    tic
for m=1:1:NtRf
    % Selecting the best RF beamforming vector
    Epsi=At'*Fres_ayash;
    [val,Ind_Direction]=max(diag(Epsi*Epsi'));
    Frf_ayash=[Frf_ayash At(:,Ind_Direction)];
    % Digital precoding
    Fbb_ayash=pinv(Frf_ayash)*Fopt;
    Fres_ayash=(Fopt-Frf_ayash*Fbb_ayash)/norm(Fopt - Frf_ayash * Fbb_ayash,'fro');
end
    Fbb_ayash = sqrt(Num_users)*Fbb_ayash/(norm(Frf_ayash*Fbb_ayash,'fro'));
     F_ayash = Frf_ayash*Fbb_ayash;



%% AltMin Original Approach Assuming Continous Phase shifters (No Codebook)
 [FRF_YuOriginal, FBB_YuOriginal] = PE_AltMin( Fopt, NtRf);
 FBB_YuOriginal = sqrt(Num_users) * FBB_YuOriginal / norm(FRF_YuOriginal * FBB_YuOriginal,'fro');
F_YuOriginal =FRF_YuOriginal*FBB_YuOriginal;


%% Single User (No Interference) and Analog beamforming calculations
        for u=1:1:Num_users
            Int_set=[]; % interference index
            for i=1:1:Num_users
                if(i~=u)
                    Int_set=[Int_set i]; 
                end
            end
            Channel=zeros(RX_ant,TX_ant);
            Channel(:,:)= H(u,:,:);
            [U_channel S_channel V_channel]=svd(Channel);
            
            % Single-user rate
            Rate_SU(count)=Rate_SU(count)+log2(1+SNR*S_channel(1,1)^2)/(Num_users*ITER);
            
            % Analog-only beamforming
            SINR_BS=(SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf(:,Int_set)).^2))+1);
            Rate_BS(count)=Rate_BS(count)+log2(1+SINR_BS)/(Num_users*ITER);            
        end
    
        % ZF Hybrid Precoding
        Rate_HP(count)=Rate_HP(count)+log2(det(eye(Num_users)+SNR*((He_fd*FRF_ZF)*(Fbb*Fbb')*(He_fd*FRF_ZF)')))/(Num_users*ITER);
            
        % Hybrid Precoding MMSE
         Rate_HP_MSE(count)=Rate_HP_MSE(count)+log2(det(eye(Num_users)+SNR*((He_fd*FRFMSE)*(FbbMSE*FbbMSE')*(He_fd*FRFMSE)')))/(Num_users*ITER);

        % New Algorithm
         Rate_HP_Yu(count)=Rate_HP_Yu(count)+log2(det(eye(Num_users)+SNR*(He_fd*(F_Novel*F_Novel')*He_fd')))/(Num_users*ITER);

          % AltMin Algorithm
         Rate_HP_Original(count)=Rate_HP_Original(count)+log2(det(eye(Num_users)+SNR*(He_fd*(F_YuOriginal*F_YuOriginal')*He_fd')))/(Num_users*ITER);

         
        % Hybrid Precoding Kalman
        Rate_HP_Kalman(count)=Rate_HP_Kalman(count)+log2(det(eye(Num_users)+SNR*(He_fd*(Fkalman*Fkalman')*He_fd')))/(Num_users*ITER);

        % ZF fully digital precoding
        Rate_HP_FD_ZF(count)=Rate_HP_FD_ZF(count)+log2(det(eye(Num_users)+SNR*(He_fd*(Ffd*Ffd')*He_fd')))/(Num_users*ITER);
        
        % MSE fully digital precoding
        Rate_HP_FD_MSE(count)=Rate_HP_FD_MSE(count)+log2(det(eye(Num_users)+SNR*(He_fd*(FfdMSE*FfdMSE')*He_fd')))/(Num_users*ITER);
        
        % OMP Algorithm 
        Rate_HP_ayash(count)=Rate_HP_ayash(count)+log2(det(eye(Num_users)+SNR*(He_fd*(F_ayash*F_ayash')*He_fd')))/(Num_users*ITER);
        
 
    end % End of SNR loop
end % End of ITER loop

%Plotting the spectral efficiency
    figure,plot(SNR_dB_range,Rate_SU,'-mv','LineWidth',1.5); hold on;
        hold on; plot(SNR_dB_range,Rate_HP_FD_MSE,'--k','LineWidth',2);
        hold on; plot(SNR_dB_range,Rate_HP_Original,'->','Color',[0.5 0 0.8],'LineWidth',2);

                hold on; plot(SNR_dB_range,Rate_HP_Yu,'-^g','LineWidth',2);
    hold on; plot(SNR_dB_range,Rate_HP_ayash,'-+c','LineWidth',2);
    hold on; plot(SNR_dB_range,Rate_HP_Kalman,'-yo','LineWidth',2);
    hold on; plot(SNR_dB_range,Rate_HP_MSE,'-b*','LineWidth',1.5);
    plot(SNR_dB_range,Rate_HP,'-rs','LineWidth',1.5);


    set(gcf,'color','w');
    set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
    set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)


    
  
    hold on; plot(SNR_dB_range,Rate_BS,'-ko','LineWidth',1.5);
    legend('Single-user (No Interference)','MMSE Fully-Digital Precoding','PE-AltMin','Proposed Algorithm','OMP','Kalman Hybrid Precoding','MMSE Hybrid Precoding','ZF Hybrid Precoding','Analog-only Beamsteering');
xlabel('SNR (dB)','FontSize',12);
ylabel('Spectral Efficiency (bps/ Hz)','FontSize',12);
grid;


function [H,a_TX,a_RX]=generate_channels(Num_users,TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths,fc, Hight, Re, Rs )

H=zeros(Num_users,RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);  % One user channel
a_TX=zeros(TX_ant_w*TX_ant_h,Num_users); % TX steering vector
a_RX=zeros(RX_ant_w*RX_ant_h,Num_users); % RX steering vector
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);

% Define the maximum angle ?_max
theta_max = 30*pi/180; % Example value, adjust as needed

% Number of points to generate
num_points = Num_users; % Example number of points

% Initialize arrays to store generated points
xi = zeros(num_points, 1);
yi = zeros(num_points, 1);

% Generate points within the circle region
for i = 1:num_points
    % Generate positive random angles within the range [0, ?_max]
theta_x = theta_max * rand();
theta_y = theta_max * rand();
    
    % Convert angles to Cartesian coordinates
    x = sin(theta_y)*cos(theta_x);
    y = sin(theta_y)*sin(theta_x);
    
    % Check if the point falls within the circle region
    if x^2 + y^2 <= sin(theta_max)^2
        xi(i) = x;
        yi(i) = y;
    end
end
% Convert Cartesian coordinates to angles
eta_x = acos(xi ./sqrt((xi).^2 + (yi).^2)); % Angle in the x-direction
eta_y = asin(yi ./ sqrt((xi).^2 + (yi).^2)); % Angle in the y-direction

% Constructing the channels
for u=1:1:Num_users

   
    AoA_el(u,:)=(pi*rand(1,1)-pi/2)*ones(1,Num_paths);
    AoA_az(u,:)=(2*pi*rand(1,1))*ones(1,Num_paths);
     AoD_el(u,:)=acos(eta_y(u))*ones(1,Num_paths);
AoD_az(u,:)=acos(eta_x(u)/sin(AoD_el(u)))*ones(1,Num_paths);
   

Nadir(u) = acos(sin(AoD_el(u))*sin(AoD_az(u)));% angle between the satellite's position directly above the Earth's surface (at the sub-satellite point) and the line connecting the satellite to the user. In simpler terms, it represents how far the user is from the satellite's "straight-down" view
alpha_new(u) = acos((Rs/Re)*sin(Nadir(u)));% elevation angle

%   % Convert alpha_new(u) from radians to degrees
% angle_deg = rad2deg(alpha_new(u));
% angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];
% 
% % Find the index of the angle range that angle_deg falls into
% index = find(angle_deg >= angles, 1, 'last');
% 
% % Assign the corresponding sigma_beta value based on the index
% switch index
%     case 1
%         sigma_beta(u) = 3.5;
%     case 2
%         sigma_beta(u) = 3.4;
%     case 3
%         sigma_beta(u) = 2.9;
%     case 4
%         sigma_beta(u) = 3;
%     case 5
%         sigma_beta(u) = 3.1;
%     case 6
%         sigma_beta(u) = 2.7;
%     case 7
%         sigma_beta(u) = 2.5;
%     case 8
%         sigma_beta(u) = 2.3;
%     case 9
%         sigma_beta(u) = 1.2;
%     otherwise
%         error('Invalid angle range');
% end
% % Generate random PDP values
% Beta(u) = abs(sigma_beta (u)* randn(1, 1));
D(u)=sqrt(Re^2*(sin(alpha_new(u)))^2+Hight^2+2*Hight*Re)-Re*sin(alpha_new(u));
% PL(u)= 20 * log10(fc) + 20 * log10(D(u)) +32.45 +Beta(u) + 2; % 2dB is the ionspheric loss

 % Parameters for the Gaussian distribution
% mean_eta = (10^(-PL(u)/10));  
% Generate random numbers from the exponential distribution
% eta_pdp = (exprnd(mean_eta, Num_paths, 1));
% eta_norm(u,:) = (eta_pdp/ sum(eta_pdp))';%total power is 1.
% gamma(u,:) = sqrt(1/2)*(((eta_norm(u,:))').*(randn(Num_paths, 1) + 1i * randn(Num_paths, 1)));

    
    
    
    % Calculate the propagation delay based on altitude and speed of light
propagation_delay = (1/(3e8))*D(u);
std_deviation = 0.1 * propagation_delay;  % Standard Deviation of the propagation delay

     tau =  propagation_delay+(std_deviation * randn(Num_paths, 1)) ;%generated in seconds
  A = 10; % Rician Factor 
 sigmaR = 0.1; % Adjust the spread (standard deviation) as needed

% Generate the real and imaginary parts from Gaussian distributions
real_part(1) =   A+ sigmaR * randn(1, 1);
real_part(2:Num_paths) = sigmaR * randn(Num_paths-1, 1);
imaginary_part = sigmaR * randn(Num_paths, 1);
% Real and imaginary parts to form Rician-distributed alpha
alpha(u,:) = (real_part' + 1j * imaginary_part);
%    alpha(u,:)=sqrt(Beta(u))*gamma(u,:);

 Temp_Channel=zeros(RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);

  for l=1:1:Num_paths
  a_TX(:,u)=transpose(sqrt(1/(TX_ant_w*TX_ant_h))*exp(1j*pi*(ind_TX_w*sin(AoD_az(u,l))*sin(AoD_el(u,l))+ind_TX_h*cos(AoD_el(u,l))) ));
  a_RX(:,u)=transpose(sqrt(1/(RX_ant_w*RX_ant_h))*exp(1j*pi*(ind_RX_w*sin(AoA_az(u,l))*sin(AoA_el(u,l))+ind_RX_h*cos(AoA_el(u,l))) ));
  Temp_Channel=Temp_Channel+ sqrt((TX_ant_w*TX_ant_h/Num_paths)*(RX_ant_w*RX_ant_h))*alpha(u,l)*a_RX(:,u)*a_TX(:,u)';
   

    end
    H(u,:,:)=Temp_Channel;
end

end
