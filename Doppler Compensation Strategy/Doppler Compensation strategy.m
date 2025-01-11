% Constants
clear;
Re = 6.37e6;          % Earth's radius in meters
h = 355e3;           % Altitude in meters
R = Re+h;           % Satellite LEO circular orbit radius in meters
i = deg2rad(60);      % Inclination in radians
we = 7.292115e-5;     % Earth's rotation rate in rad/second
ws = 9.3077e-4;       % Satellite's angular velocity in rad/sec
Thetao = deg2rad(20); % Initial Theta in radians
Theta = deg2rad(20:1:90); % Range of Theta values in radians
c = 3e8;              % Speed of light in meters/second
fc=28*1e9;
% Calculate visibility (tao) for each Theta value
tao = zeros(size(Theta));
for idx = 1:length(Theta)
    xx= cos(acos((Re / R) * cos(Thetao)) - Thetao);
    yy =cos(acos((Re / R) * cos(Theta(idx))) - Theta(idx));
    arg = acos( xx/ yy);
    tao(idx) = 2 / (ws - (we * cos(i))) * arg;
end

% Plot visibility (tao) against Theta
figure;
plot(rad2deg(Theta), tao);
xlabel('Maximum Elevation Angle (deg)');
ylabel('Visibility Duration (s)');
legend('Satellite Altitude = 355 km');
grid on;
set(gcf,'color','w');
    set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
    set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)

  
  Thetat1 = (Theta(39));
  t1= tao(39)/2;
  t=0:0.1:round(tao(39));
  Doppler = zeros(1,length(t));
  alphat1 = acos(Re/R*cos(Thetat1))-Thetat1;
for idx =1: length(t)
      psi  = (ws-we*cos(i))*(t(idx)-t1);
      op = -Re*R *sin(psi)*cos(acos((Re / R) * cos(Thetat1)) - Thetat1);
      on= c*sqrt(Re^2+R^2 - 2*Re*R*cos(psi)*cos(acos((Re / R) * cos(Thetat1)) - Thetat1));
      wf= ws-we*cos(i);  
      DeltaF(idx)= fc*op/on*wf;
      alpha(idx) = acos(cos(psi)*cos(alphat1));
      s(idx)=sqrt(Re^2+R^2-2*Re*R*cos(alpha(idx)));
    theta(idx)=  acos(R*sin(alpha(idx))/s(idx));
end
  
 % Plot visibility (tao) against Theta
figure;
plot(t, DeltaF/1000);
xlabel('Time (s)');
ylabel('Doppler (kHz)');
legend('Carrier Frequency (fc) = 28 GHz ' );
grid on;
set(gcf,'color','w');
    set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
    set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)


     % Plot visibility (tao) against Theta
figure;
plot(t, s/1000);
xlabel('Time (s)');
ylabel('Slant Range (km)');
grid on;
set(gcf,'color','w');
    set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
    set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)

       % Plot visibility (tao) against Theta
figure;
plot(t, rad2deg(theta));
xlabel('Time (s)');
ylabel('Elevation Angle (Degree)');
grid on;
set(gcf,'color','w');
    set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
    set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)
