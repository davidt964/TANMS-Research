%=============================Clear Cache==================================
clc; close all; clear all;

%==============================Constants===================================
mu_0 = 4*pi * 1E-7;             %Permeability of free space in [H/m]
Z_0 = 50;                       %Impedance of input transmission line in [ohms]
eta = 377;                      %Free space radiation impedance in [ohms]
f = 20E3;                       %Operating frequency of antenna in [Hz]
f_2 = 175E3;                    %Operating frequency of Emily's antenna in [Hz]
P_in = 100E-3;                  %Input power of antenna in [W]
c = 3E8;                        %Speed of light in [m/s]
a = 0.02256;                    %Radius of loop antenna in [m]
S = pi*a^2;                     %Area of loop antenna in [m^2]
lambda = c/f;                   %Operating wavelength of antenna in [m]
lambda_2 = c/f_2;               %Operating wavelength of Emily's antenna in [m]
k = 2 * pi / lambda;            %Wavenumber in [m^-1]
k_2 = 2 * pi / lambda_2;        %Wavenumber  of Emily's antenna in [m^-1]
r = 0.40;                       %Distance of data measurement from center of antenna in [m]
R_r = 31171 * (S^2) / (lambda^4);
                                %Radiation resistance of the antenna
R_r2 = 31171 * (S^2) / (lambda_2^4);
                                %Radiation resistance of Emily's antenna
N_pts = 100;                    %Number of points


%================================Arrays====================================
theta = linspace(0, 2*pi, N_pts);      %Array to store angular values

%%For 30 kHz
B_r = zeros(1, N_pts);                 %Array to store Br
B_r_neg = zeros(1, N_pts);             %Array to store negative Br
B_r_norm = zeros(1, N_pts);            %Array to store normalized Br
B_r_norm_neg = zeros(1, N_pts);        %Array to store normalized negative Br

%%For 175 kHz
B_r_175 = zeros(1, N_pts);             %Array to store Br for 175 kHz antenna
B_r_neg_175 = zeros(1, N_pts);         %Array to store negative Br for 175 kHz antenna
B_r_norm_175 = zeros(1, N_pts);        %Array to store normalized Br for 175 kHz antenna
B_r_norm_neg_175 = zeros(1, N_pts);    %Array to store normalized negative Br for 175 kHz antenna

test = mu_0/sqrt(2) * sqrt((3)/(pi * eta * (k_2)^4) * (Z_0 * R_r2 * P_in)/(R_r2 + Z_0)^2) * 2/r^3;
%===========================Analytical Solution============================
%Note that this analytical solution is only an approximation
%Absolute data
for i=1:N_pts
    B_r(i) = mu_0/sqrt(2) * sqrt( ((3)/(pi * eta * k^4)) * ((Z_0 * R_r * P_in)/(R_r + Z_0)^2)    ) * (2 * cos(theta(i)))/(r^3);
    B_r_neg(i) = -mu_0/sqrt(2) * sqrt((3)/(pi * eta * k^4) * (Z_0 * R_r * P_in)/(R_r + Z_0)^2) * ((2 * cos(theta(i)))/r^3);
    B_r_175(i) = mu_0/sqrt(2) * sqrt((3)/(pi * eta * (k_2)^4) * (Z_0 * R_r2 * P_in)/(R_r2 + Z_0)^2) * ((2 * cos(theta(i)))/r^3);
    B_r_neg_175(i) = -mu_0/sqrt(2) * sqrt((3)/(pi * eta * (k_2)^4) * (Z_0 * R_r2 * P_in)/(R_r2 + Z_0)^2) * (2 * cos(theta(i))/r^3);
%     fprintf('%.10f  %.10f %.5f\n', B_r(i), B_r_175(i), k); %demonstrates
%     that the data for the two frequencies are identical to one another
end

%Normalized data
max_Br = max(B_r);
max_neg_Br = max(B_r_neg);
max_Br_175 = max(B_r_175);
max_neg_Br_175 = max(B_r_neg_175);

for i = 1:N_pts
    B_r_norm(i) = B_r(i) / max_Br;
    B_r_norm_neg(i) = B_r_neg(i) / max_neg_Br;
    B_r_norm_175(i) = B_r_175(i) / max_Br_175;
    B_r_norm_neg_175(i) = B_r_neg_175(i) / max_neg_Br_175;
end

%==================================Plots===================================
%%Polar Plots
set(gcf, 'Position',  [400, 100, 1000, 800]);
figure(1);
polaraxes('ThetaZeroLocation', 'top');
hold on;
polarplot(theta, B_r_norm, 'm', 'LineWidth', 2);
polarplot(theta, B_r_norm_neg, 'm', 'LineWidth', 2);
title(['Normalized $B_r$ vs. $\theta$ for $f =$ ', num2str(f), ' Hz at $r = $ ', num2str(r), 'm'], 'Interpreter', 'latex');
pax = gca;
pax.FontSize = 14;
pax.ThetaDir = 'clockwise';
pax.RLim = [0 1];
textbox_size = [.375 .1 .1 .2];
annotation('textbox', textbox_size, 'String', ['Max $B_r = $ ', num2str(max_Br), ' T'], 'Interpreter', 'latex', 'FitBoxToText', 'on');
pax.RAxisLocation = 270;
hold off;

figure(2);
polaraxes('ThetaZeroLocation', 'top');
hold on;
set(gcf, 'Position',  [400, 100, 1000, 800]);
polarplot(theta, B_r_norm_175, 'm', 'LineWidth', 2);
polarplot(theta, B_r_norm_neg_175, 'm', 'LineWidth', 2);
title(['Normalized $B_r$ vs. $\theta$ for $f =$ ', num2str(f_2), ' Hz at $r = $ ', num2str(r), 'm'], 'Interpreter', 'latex');
pax = gca;
pax.FontSize = 14;
pax.ThetaDir = 'clockwise';
pax.RLim = [0 1];
textbox_size = [.375 .1 .1 .2];
annotation('textbox', textbox_size, 'String', ['Max $B_r = $ ', num2str(max_Br_175), ' T'], 'Interpreter', 'latex', 'FitBoxToText', 'on');
pax.RAxisLocation = 270;

%%3D Radiation Pattern (uncertain if correct)
d = dipole('Length', 0.01, 'Width', .00034);
pattern(d, 175E3);