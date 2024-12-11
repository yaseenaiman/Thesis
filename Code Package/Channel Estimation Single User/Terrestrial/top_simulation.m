clear
clc

%% System Parameters
SP = struct();
SP.N_t = 32;
SP.N_r = 16;
SP.L_r = 1;
SP.L_t = SP.L_r;
SP.N_s = SP.L_r;
SP.N_c = 4;
SP.N_p = 2;
SP.N_Q = 6;
SP.M   = 100;
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
SP.SNR_db_array     = -10:5:10;
SP.SNR_lin_array    = 10.^(.1*SP.SNR_db_array);
SP.NMSE_LS          = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_OMP         = zeros(SP.iter, length(SP.SNR_db_array));
 SP.NMSE_NOMP = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_pcgp = zeros(SP.iter, length(SP.SNR_db_array));

SP.LS =zeros(1,length(SP.SNR_db_array)); % Achievable rate with estimated channel
SP.OMP = zeros(1,length(SP.SNR_db_array));
SP.NOMP = zeros(1,length(SP.SNR_db_array));
SP.pcgp = zeros(1,length(SP.SNR_db_array));

%% Simulation
SP_result = top_sim_function(SP);


