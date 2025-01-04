

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear all;clc;warning('off');
% ----------------------------- System Parameters -------------------------
Num_users=2:2:16; % Number of users
TX_ant=64; %Number of UPA TX antennas
TX_ant_w=sqrt(TX_ant); % width of the planar array
TX_ant_h=sqrt(TX_ant); % hieght of the planar array 
NtRf = 16;%number of transmitter RF chains
RX_ant=16; %Number of UPA RX antennas
RX_ant_w=sqrt(RX_ant); % width 
RX_ant_h=sqrt(RX_ant); % hieght

% ----------------------------- Channel Parameters ------------------------
Num_paths=10; %Number of channel paths

% ----------------------------- Simulation Parameters ---------------------
SNR_dB_range=20;  % SNR in dB
time_SU=zeros(1,length(Num_users)); % single user spectral efficiency (no-interference)
time_BS=zeros(1,length(Num_users));% analog-only beamsteering spectral efficiency
time_HP=zeros(1,length(Num_users)); % hybrid ZF precoding spectral efficiency
time_HP_MSE=zeros(1,length(Num_users)); %  hybrid MMSE precoder spectral efficiency 
time_HP_Kalman=zeros(1,length(Num_users)); % hybrid Kalman precoder spectral efficiency
time_HP_FD_ZF=zeros(1,length(Num_users)); % fully digital  ZF precoding spectral efficiency
time_HP_FD_MSE=zeros(1,length(Num_users)); % fully digital  MMSE precoding spectral efficiency
time_HP_ayash = zeros(1,length(Num_users)); %  OMP_Ayash Spectral Efficiency
time_HP_Yu = zeros(1,length(Num_users)); %  Proposed Algorithm Spectral Efficiency 
time_HP_Original = zeros(1,length(Num_users)); % PE-AltMin Spectral Efficiency
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
for count=1:length(Num_users)
for iter=1:1:ITER
    
    %% Spectral efficiency calculations
  

    % Generate user channels 
    [H,a_TX,a_RX]=generate_channels(Num_users(count),TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths); 
    % H is a 3-dimensional matrix, with Num_users(count),RX_ant,TX_ant dimensions

        
    % Intialization of the hybrid precoder and combiner
    Frf=zeros(TX_ant,Num_users(count)); % BS RF precoders 
    Wrf=zeros(RX_ant,Num_users(count)); % MS RF precoders 
    
    for u=1:1:Num_users(count)
        Frf(:,u)=a_TX(:,u); %Ideal FRF, assuming AODs are known
        Wrf(:,u)=a_RX(:,u); %Ideal WRF, assuming AOAs are known
    end      
    
  
    
    %% Effective Channel Assuming Wrf (Combiners of users)is known or  AOAs are perfectly known 
    % when performing beamforming at the transmitter end. 
    
 for u=1:1:Num_users(count)
        Channel=zeros(RX_ant,TX_ant);
        Channel(:,:)= H(u,:,:);
        He_fd(u,:)=Wrf(:,u)'*Channel;     
 end
    
    

    
        SNR=10^(.1*SNR_dB_range)/Num_users(count); % SNR value
        
        sigma2=1/SNR;
        
        
        
    %%  Maximum-Ratio Transmission (MRT) for interference cancellation
    tic
 FfdMRT=He_fd';
 for u=1:1:Num_users(count) % Normalization of the  precoder
 FfdMRT(:,u)=FfdMRT(:,u)/sqrt((FfdMRT(:,u))'*(FfdMRT(:,u)));
 end
 mrtTime=toc;
 %% 
Fopt=  FfdMRT; %  We are relying on MRT approach to minimize 
 % the interference amongs users, then, we are implementing 
 % various Hybrid Beamforming algorithms using it. Since the user does not
 % have the ability to cancel the interference, this is carried out at the
 % Satellite station side

%% Novel Algorithm Combining OMP by Ayash and PE AltMin by Xianghao Yu
tic
[FRF_Novel, FBB_Novel] = Novel_Algorithm(Fopt, NtRf, Steering_Vec);
FBB_Novel = sqrt(Num_users(count)) * FBB_Novel / norm(FRF_Novel * FBB_Novel,'fro'); %Normalization
F_Novel =FRF_Novel*FBB_Novel;
Novel=toc;

%% MMSE baseband precoder, Hybrid Beamforming
tic
FRFMSE = DiscreteRFGeneration(Fopt,Steering_Vec,NtRf); %Generating FRF Using Fully digital MMSE Fopt
FbbMSE=inv((He_fd*FRFMSE)'*(He_fd*FRFMSE)+Num_users(count)*sigma2*FRFMSE'*FRFMSE)*(He_fd*FRFMSE)';
       
        for u=1:1:Num_users(count) % Normalization of the hybrid precoders
            FbbMSE(:,u)=FbbMSE(:,u)/sqrt((FRFMSE*FbbMSE(:,u))'*(FRFMSE*FbbMSE(:,u)));
        end
        MMSEtime = toc;
        
%% Baseband zero-forcing precoding, Hybrid Beamforming 
tic
 FRF_ZF = DiscreteRFGeneration(Fopt,Steering_Vec,NtRf); %Generation of FRF depending on ZF fully digital Fopt
 Fbb=(He_fd*FRF_ZF)'*((He_fd*FRF_ZF)*(He_fd*FRF_ZF)')^(-1); 
    for u=1:1:Num_users(count) % Normalization of the hybrid precoders
        Fbb(:,u)=Fbb(:,u)/sqrt((FRF_ZF*Fbb(:,u))'*(FRF_ZF*Fbb(:,u)));
    end
    
    ZFhy=toc;
    %% Kalman Beamforming
       tic
 FRFKalman = DiscreteRFGeneration(Fopt,Steering_Vec,NtRf); %Depnding on Fully digital MMSE Fopt

  Fbbk=eye(NtRf,Num_users(count));
  RN1=Fbbk*Fbbk';
  Qm=eye(Num_users(count))*sigma2;
  ITERK=10; % Number of Kalman iterations
        
  for ii=1:ITERK
    Hk=He_fd*FRFKalman;
    K=RN1*Hk'*pinv(Hk*RN1*Hk'+Qm);
    errk=(eye(Num_users(count))-Hk*Fbbk);
    errk=errk/norm(errk);
    Fbbk=Fbbk+K*errk; 
    RN=RN1-K*Hk*RN1;
    RN1=RN; 
  end
        
        for u=1:1:Num_users(count) % Normalization of the hybrid precoders
            Fbbk(:,u)=Fbbk(:,u)/sqrt((FRFKalman*Fbbk(:,u))'*(FRFKalman*Fbbk(:,u)));
        end
       Fkalman= FRFKalman*Fbbk;
KalTime=toc;
    %%  OMP Hybrid Beamforming 
    tic
    Frf_ayash=[];
    Fres_ayash = Fopt;
    At = Steering_Vec;
for m=1:1:NtRf
    % Selecting the best RF beamforming vector
    Epsi=At'*Fres_ayash;
    [val,Ind_Direction]=max(diag(Epsi*Epsi'));
    Frf_ayash=[Frf_ayash At(:,Ind_Direction)];
    % Digital precoding
    Fbb_ayash=pinv(Frf_ayash)*Fopt;
    Fres_ayash=(Fopt-Frf_ayash*Fbb_ayash)/norm(Fopt - Frf_ayash * Fbb_ayash,'fro');
end
    Fbb_ayash = sqrt(Num_users(count))*Fbb_ayash/(norm(Frf_ayash*Fbb_ayash,'fro'));
     F_ayash = Frf_ayash*Fbb_ayash;
AyashTime=toc;


%% AltMin Original Approach Assuming Continous Phase shifters (No Codebook)
tic
 [FRF_YuOriginal, FBB_YuOriginal] = PE_AltMin( Fopt, NtRf);
 FBB_YuOriginal = sqrt(Num_users(count)) * FBB_YuOriginal / norm(FRF_YuOriginal * FBB_YuOriginal,'fro');
F_YuOriginal =FRF_YuOriginal*FBB_YuOriginal;

AltminTime=toc;

    
        % ZF Hybrid Precoding
        time_HP(count)=time_HP(count)+((ZFhy+mrtTime)/ITER);
            
        % Hybrid Precoding MMSE
         time_HP_MSE(count)=time_HP_MSE(count)+((MMSEtime+mrtTime)/ITER);

        % New Algorithm
         time_HP_Yu(count)=time_HP_Yu(count)+((Novel+mrtTime)/ITER);

          % AltMin Algorithm
         time_HP_Original(count)=time_HP_Original(count)+((AltminTime+mrtTime)/ITER);

         
        % Hybrid Precoding Kalman
        time_HP_Kalman(count)=time_HP_Kalman(count)+((KalTime+mrtTime)/ITER);

        
        % OMP Algorithm 
        time_HP_ayash(count)=time_HP_ayash(count)+((AyashTime+mrtTime)/ITER);
        
 
    end % End of ITER loop
end % End of count loop

%Plotting the spectral efficiency
      figure,  plot(Num_users,time_HP_Original*1e3,'->','Color',[0.5 0 0.8],'LineWidth',2);

                hold on; plot(Num_users,time_HP_Yu*1e3,'-^g','LineWidth',2);
    hold on; plot(Num_users,time_HP_ayash*1e3,'-+c','LineWidth',2);
    hold on; plot(Num_users,time_HP_Kalman*1e3,'-yo','LineWidth',2);
    hold on; plot(Num_users,time_HP_MSE*1e3,'-b*','LineWidth',1.5);
    plot(Num_users,time_HP*1e3,'-rs','LineWidth',1.5);


    set(gcf,'color','w');
    set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
    set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)


    
  
    legend('PE-AltMin','Proposed Algorithm','OMP','Kalman Hybrid Precoding','MMSE Hybrid Precoding','ZF Hybrid Precoding');
xlabel('Number of Users','FontSize',12);
ylabel('Time (ms)','FontSize',12);
grid;



function [H,a_TX,a_RX]=generate_channels(Num_users,TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths)

H=zeros(Num_users,RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);  % One user channel
a_TX=zeros(TX_ant_w*TX_ant_h,Num_users); % TX steering vector
a_RX=zeros(RX_ant_w*RX_ant_h,Num_users); % RX steering vector
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);
% Constructing the channels
for u=1:1:Num_users
    AoD_el(u,:)=pi*rand(1,Num_paths)-pi/2;
    AoD_az(u,:)=2*pi*rand(1,Num_paths);
    AoA_el(u,:)=pi*rand(1,Num_paths)-pi/2;
    AoA_az(u,:)=2*pi*rand(1,Num_paths);
    alpha(u,:)=sqrt(1/Num_paths)*sqrt(1/2)*(randn(1,Num_paths)+1j*randn(1,Num_paths));

    Temp_Channel=zeros(RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);
    for l=1:1:Num_paths
        a_TX(:,u)=transpose(sqrt(1/(TX_ant_w*TX_ant_h))*exp(1j*pi*(ind_TX_w*sin(AoD_az(u,l))*sin(AoD_el(u,l))+ind_TX_h*cos(AoD_el(u,l))) ));
        a_RX(:,u)=transpose(sqrt(1/(RX_ant_w*RX_ant_h))*exp(1j*pi*(ind_RX_w*sin(AoA_az(u,l))*sin(AoA_el(u,l))+ind_RX_h*cos(AoA_el(u,l))) ));
        Temp_Channel=Temp_Channel+sqrt((TX_ant_w*TX_ant_h)*(RX_ant_w*RX_ant_h))*alpha(u,l)*a_RX(:,u)*a_TX(:,u)';
    end
    H(u,:,:)=Temp_Channel;
end

end
