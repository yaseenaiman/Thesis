function [H, A_t, A_r, AoD, AoA] = channelGen(SP)


    a_t = @(N_t, aod) 1/sqrt(N_t)*exp(-1j*(0:N_t-1).'*aod);
    a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(-1j*(0:N_r-1).'*aoa);




% Define the maximum angle ?_max
theta_max = 30*pi/180; % Example value, adjust as needed

% Number of points to generate
num_points = 1; % 1 user

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

   
     AoA_el=(pi*rand(1,1)-pi/2)*ones(1,SP.N_p);
    AoA=(2*pi*rand(1,1))*ones(1,SP.N_p);
      AoD_el=acos(eta_y)'.*ones(1,SP.N_p);
AoD=acos(eta_x/sin(AoD_el(1)))*ones(1,SP.N_p);
   
Nadir = acos(sin(AoD_el(1))*sin(AoD(1))'); % angle between the satellite's position directly above the Earth's surface (at the sub-satellite point) and the line connecting the satellite to the user. In simpler terms, it represents how far the user is from the satellite's "straight-down" view
alpha_new = acos((SP.Rs /SP.Re)*sin(Nadir)); % elevation angle
%   % Convert alpha_new(u) from radians to degrees
% angle_deg = rad2deg(alpha_new);
% angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];
% 
% % Find the index of the angle range that angle_deg falls into
% index = find(angle_deg >= angles, 1, 'last');
% 
% % Assign the corresponding sigma_beta value based on the index
% switch index
%     case 1
%         sigma_beta = 3.5;
%     case 2
%         sigma_beta= 3.4;
%     case 3
%         sigma_beta = 2.9;
%     case 4
%         sigma_beta = 3;
%     case 5
%         sigma_beta = 3.1;
%     case 6
%         sigma_beta = 2.7;
%     case 7
%         sigma_beta = 2.5;
%     case 8
%         sigma_beta = 2.3;
%     case 9
%         sigma_beta= 1.2;
%     otherwise
%         error('Invalid angle range');
% end
% % Generate random PDP values
% Beta = abs(sigma_beta * randn(1, 1));
D=sqrt(SP.Re^2*(sin(alpha_new))^2+SP.Hight^2+2*SP.Hight*SP.Re)-SP.Re*sin(alpha_new);
% PL= 20 * log10(SP.fc) + 20 * log10(D) +32.45 +Beta + 2; % 2dB is the ionspheric loss
% 
%  % Parameters for the Gaussian distribution
% mean_eta = (10^(-PL/10));  
% % Generate random numbers from the exponential distribution
% eta_pdp = (exprnd(mean_eta, SP.N_p, 1));
% eta_norm = (eta_pdp/ sum(eta_pdp))';%total power is 1.
% gamma = sqrt(1/2)*(((eta_norm)').*(randn(SP.N_p, 1) + 1i * randn(SP.N_p, 1)));
% 
%     
%     
    
    % Calculate the propagation delay based on altitude and speed of light
propagation_delay = (1/(3e8))*D;
std_deviation = 0.1 * propagation_delay;  % Standard Deviation of the propagation delay

     tau =  propagation_delay+(std_deviation * randn(SP.N_p, 1)) ;%generated in seconds
 

A = 10; % Rician Factor 
sigmaR = 0.1; % Adjust the spread (standard deviation) as needed

 % Generate the real and imaginary parts from Gaussian distributions
real_part(1) =   A+ sigmaR * randn(1, 1);
real_part(2:SP.N_p) = sigmaR * randn(SP.N_p-1, 1);
imaginary_part = sigmaR * randn(SP.N_p, 1);
% Real and imaginary parts to form Rician-distributed alpha
alpha = (real_part' + 1j * imaginary_part);




A_t = a_t(SP.N_t, AoD);
A_r = a_r(SP.N_r, AoA);
    
H = zeros(SP.N_r, SP.N_t, SP.N_c);
for d = 1:SP.N_c
    p     = RC_shaping(1, SP.beta, (d-1)-tau);
    Delta = diag(p.*alpha);
    H_temp = A_r*Delta*A_t';
    H(:,:,d) = sqrt(SP.N_r*SP.N_t/SP.N_p)*H_temp; %/norm(H_temp, 'fro');
    disp('')
end
