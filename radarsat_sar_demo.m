% RADARSAT SAR Processor Demo
% Dependencies Image processing toolbox, SATData.mat
% Imaged Area Coords:  49.02153897564771, -123.15890344471491

clear;                                                                                   
close all;                                                                               
load SATData.mat;
c = 3e8;
% SAR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = 993.4627e3; %meters
y0 = -27.466e3;
r0 = sqrt(x0^2 + y0^2); % (m)                                                                         
lambda = 0.05657; % (m)  
PRF = 1256.98; % Hertz                                              
V = 7062; %(m/s) velocity of satellite   
BW = 300e6; %LFM Bandwidth Hz
tau_p = (1/32.317)*1e-6; %uncompressed pulse width
tau_res = tau_p; % 1/BW???????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_lfm = BW / tau_p;
                                                                                         
del_y = 1; % (m) desired resolution                                                      
theta_B = del_y / r0; %beam width of SAR array
L = 8624; %parameter from satellite SAR data set
T = 1/PRF; % Need T < ( lambda * r0) / (2 * V * w) to prevent aliasing 
w = (lambda * r0) / (2 * V * T);
l = w;

T_L = L/V; %How long the SAR collects data                                             
K_L = floor(T_L/(2*T))+1;                                                                                                   
                                                                                                                         
a0 = (2 * V^2) / (lambda * r0);                                                          
k = -K_L:K_L-1;
r_res = tau_res * c / 2;

%RADARSAT sar data has had gross doppler removed in preprocessing
%fdg = (2*V*y0)/(lambda*r0); %gross doppler removal

rmax = sqrt((x0+l/2)^2 + (abs(y0)+w/2+L/2)^2);
tau_min = 2*(x0 - l/2)/c;
tau_max = 2*rmax/c;
tau_r0 = 2*r0/c;

m_min = (tau_min-tau_r0) / tau_res;
m_max = (tau_max-tau_r0) / tau_res;
m = m_min:m_max-2;

tau = repmat((tau_r0 + m * tau_res)', [1 length(k)]);

V_BB = v;

% Perform RCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau_p_cell = tau_p;
f = k' * 1/(K_L*tau_p_cell); %fft tap frequencies    % 1536 x 1
rmin = sqrt(x0^2 + (abs(y0) - L/2 - w/2)^2);
delta_tau = 2*(sqrt(x0^2 + (y0-V*k*T).^2) - rmin)/c; % 1 x 1536
V_BB_fft = fftshift(fft(V_BB, 2*(K_L)),1);               % 768 x 1536
RCMC = exp(1j*2*pi*f*delta_tau);
V_BB_corrected = ifft(fftshift(RCMC.*V_BB_fft,1));
V_BB_corrected = V_BB_corrected(1:length(m),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image Sharpening Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v_h = exp(1j*2*pi*(V^2/(lambda*r0))*(k*T).^2);
% v_h = repmat(v_h, [length(m) 1]); %attempted fix for v_h, does not work
v_deltaf = exp(1j*2*pi*( (2*y0*V)/(lambda*r0) * (m'*(c*tau_res/2)/r0)) * k*T);
v_q = exp(-1j*pi*(2*V^2 / (lambda*r0) * m'*(c*tau_res/2)/r0) * (k*T).^2);
V_BB_corrected = V_BB_corrected .* v_q .* v_deltaf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sar_image = zeros(length(m), length(k));
for i=1:length(m)
    quadratic_phase_correction = (exp(1j * pi * a0 * T^2 * k.^2));%calculate quadratic phase correction
    V_o = fftshift(fft(V_BB_corrected(i,:) .* quadratic_phase_correction)); %take fft across m for each k
    sar_image(i,:) = V_o;
end

sar_image = abs(sar_image);
sar_image = abs(1 - sar_image);

del_f = 1 / (T * K_L); % (Hz) frequency resolution 
f = (-K_L/2 : (K_L/2 - 1) ) * del_f;
x = m * c * tau_res / 2;
y = ((lambda*r0*f)/(2*V));

delta_y_corr = (y0 * m * (c*tau_res/2)) / r0;

imagesc(y,x,sar_image');
axis off;

colormap gray;
set(gca, 'Ydir', 'normal', 'Xdir', 'normal');
set(gca, 'Clim', [3000, 15000]) %adjust contrast to get decent image
title("RADARSAT data");
grid on;


