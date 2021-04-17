%==========================================================================
%Magnetic Damping Simulation
%David Tran
%{
Description: Simulates the magnetic damping of a multiferroic antennna,
generating Mx vs Hx plots.
%}

%=============================Clear Cache==================================
clc; close all; clear all;


%=============================Constants====================================
mu_0 = 4 * pi * 1E-7;                  %Free Space Permeability [N/A^2]
M_s = 4.85E5;                          %Saturation Magnetization [A/m]
gamma = 1.759 / (2 * pi) * 1E11;       %Gyromagnetic Ratio [rad * A * m / N /s]
beta_xx = 1 / 2.1152 / mu_0;           %Inverse Permeability Tensor (xx-direction) [A^2/N]
alpha = 0.045;                         %Gilbert Damping Parameter [Dimensionless]
alpha_2  = 1;
omega = 2 * pi * 400E6;                %Angular Frequency [rad/s]
B = mu_0 * M_s / 50;                   %Magnetic Field Strength [N/A/m]
a = (gamma * mu_0 * M_s) / alpha;      %Inverse Relaxation Time [1/s]
a_2 = (gamma * mu_0 * M_s) / alpha_2;  %Set damping term to a higher value to provide visualizations of results

%==============================Arrays======================================
t = linspace(0, 2*pi/omega, 1000);       %Array to store time values
Hx = B * beta_xx * sin(omega * t) + B * ((1 / mu_0) - beta_xx) * (omega / a)/(1 + (omega / a)^2) * (cos(omega * t) - exp(-a * t) + (omega / a)* sin(omega * t));
Mx = (1 / mu_0) * (B * sin(omega * t)) - Hx;
Bx = mu_0 * (Mx + Hx);

Hx_2 = B * beta_xx * sin(omega * t) + B * ((1 / mu_0) - beta_xx) * (omega / a_2)/(1 + (omega / a_2)^2) * (cos(omega * t) - exp(-a_2 * t) + (omega / a_2)* sin(omega * t));
Mx_2 = (1 / mu_0) * (B * sin(omega * t)) - Hx_2;
Bx_2 = mu_0 * (Mx_2 + Hx_2);

%==============================Plot data===================================
figure(1);
hold on;
subplot(2,1,1);
xlim([-5000 5000]);
ylim([-6000 6000]);
plot(Hx, Mx, 'LineWidth', 1.25, 'Color', 'c');

hold on;
plot(Hx_2, Mx_2, 'LineWidth', 1.25, 'Color', 'k');
xlabel('$H_x$ (A/m)', 'Interpreter', 'latex');
ylabel('$M_x$ (A/m)', 'Interpreter', 'latex');
legend({'$\alpha = 0.045$', '$\alpha = 1.0$'} , 'Interpreter', 'latex', 'Location', 'southeast');
title('$M_x$ vs. $H_x$', 'Interpreter', 'latex');
grid on;


hold on;
subplot(2,1,2);
plot(Bx, Mx, 'LineWidth', 1.25, 'Color', 'c');


hold on;
plot(Bx_2, Mx_2, 'LineWidth', 1.25, 'Color', 'k');
xlabel('$H_x$ (A/m)', 'Interpreter', 'latex');
ylabel('$B_x$ (N/A/m)', 'Interpreter', 'latex');
legend({'$\alpha = 0.045$', '$\alpha = 1.0$'} , 'Interpreter', 'latex', 'Location', 'southeast');
title('$B_x$ vs. $H_x$', 'Interpreter', 'latex');
grid on;
