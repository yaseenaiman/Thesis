clear
clc

%% System Parameters
SP = struct();
SP.N_t = 32; % number of transmit antennas
SP.N_r = 16; % number of receive antennas 
SP.L_r = 4;   %nunber of RF chains for the receiver
SP.L_t = SP.L_r; % number of RF chains for the transmitter
SP.N_s = SP.L_r; % number of data streams
SP.N_c = 4; % number of clusters
SP.N_p = 2;  % number of paths per a cluster
SP.N_Q = 6;  % number of phase shifters quantization bits
SP.M   = 80; %number of frames
SP.N   = 16; % symbols per frame, in the theortical analysis, 1 symbol per frame is assumed
SP.G_t = 2*SP.N_t; % Determines the resolution of the quantized angles for the angle of departure (AoD) dictionary
SP.G_r = 2*SP.N_r;  % Determines the resolution of the quantized angles for the angle of arrival (AoA) dictionary
SP.G_c = 4;  % Number of time samples per symbol period
SP.ADCbits = 1; % ADC number of bits
SP.bandwidth        = 500e6; % System BW
SP.sigma            = sqrt(10^(.1*(-174+10*log10(SP.bandwidth)))); % Thermal Noise Standard Deviatiom
SP.beta             = 0.8; %refers to the roll-off factor (?) of the raised cosine (RC) pulse shaping filter
SP.iter             = 30;  % Montecarlo iterations
SP.SNR_db_array     = -20:5:15; % SNR Range
SP.SNR_lin_array    = 10.^(.1*SP.SNR_db_array); % SNR range as a ratio
SP.NMSE_LS          = zeros(SP.iter, length(SP.SNR_db_array)); % NMSE for LS
SP.NMSE_OMP         = zeros(SP.iter, length(SP.SNR_db_array));  % NMSE for OMP
SP.NMSE_NOMP = zeros(SP.iter, length(SP.SNR_db_array)); % NMSE for ACGP

SP.LS =zeros(1,length(SP.SNR_db_array)); 
SP.OMP = zeros(1,length(SP.SNR_db_array));
SP.NOMP = zeros(1,length(SP.SNR_db_array));


%% Simulation
SP_result = top_sim_function(SP);

%% If  G_c =1, you only take one sample per symbol, which is insufficient to capture a smoothly shaped pulse.
%% If  G_c =4, you divide each symbol period into 4 time steps, capturing more detail about the pulse shape.