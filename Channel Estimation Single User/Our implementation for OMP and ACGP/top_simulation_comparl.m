clear
clc

%% System Parameters
SP = struct();
SP.N_t = 32;
SP.N_r = 16;
SP.L_r = 4;
SP.L_t = SP.L_r;
SP.N_s = SP.L_r;
SP.fc = 28e9;
SP.Re = 6378*1000;
SP.Hight = 560*1000;
SP.Rs = SP.Re+SP.Hight;
SP.N_c = 1;
SP.N_p = 5;
SP.N_Q = 6;
SP.M   = 80;
SP.N   = 16; % symbols per frame
SP.G_t = 2*SP.N_t;
SP.G_r = 2*SP.N_r;
SP.G_c = 4;
SP.ADCbits = 4;
SP.virtualChannel   = 'spatial'; % 'spatial', 'physical'
SP.RFonly           = 1; % RF precodier/combiner only? or with baseband precoder/combiner as well?
SP.trainType        = 'random phase'; % 'Hadamard', 'ZC', 'random', 'random phase'
SP.integerTapDelay  = 0;
SP.BeamAlign        = 0; % AoA and AoD aligned with dictionary grid
SP.beamforming      = 0; % precoder and combiner beamforming
SP.bandwidth        = 500e6;
SP.sigma            = sqrt(10^(.1*(-174+10*log10(SP.bandwidth))));
SP.beta             = 0.8;
SP.iter             = 30;
SP.SNR_db_array     = -20:5:15;
SP.SNR_lin_array    = 10.^(.1*SP.SNR_db_array);
SP.NMSE_LS          = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_OMP         = zeros(SP.iter, length(SP.SNR_db_array));
 SP.NMSE_NOMP = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_author = zeros(SP.iter, length(SP.SNR_db_array));

SP.acgp_res = zeros(SP.iter, SP.N_p*9);
SP.omp_res = zeros(SP.iter, SP.N_p*9);



SP.LS =zeros(1,length(SP.SNR_db_array)); % Achievable rate with estimated channel
SP.OMP = zeros(1,length(SP.SNR_db_array));
SP.NOMP = zeros(1,length(SP.SNR_db_array));
SP.Omp_author = zeros(1,length(SP.SNR_db_array));


%% Simulation
SP_result = top_sim_function_com(SP);


