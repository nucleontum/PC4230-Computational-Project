%% QM3 Project
% Numerical simulation of the evolution of a wavepacket in a 1D harmonic trap 
% using fast fourier transport (fft) method
% 
% Unit of energy: hbar*omega, where h_bar is the Planck constant and omega is 
% the frequency of the trap
% 
% Unit of length: l=sqrt(h_bar/(m*omega)), where sqrt(...) is the square root 
% function and m is the mass of the particle
% 
% Unit of momentum: hbar/l
% 
% energy unit: hbar\omega,  Hamiltonian --> dimensionless
% 
% time dimensionless: omega*t    i d/dxt | >= dimension H |>
% 
% dimensionless time = 2pi. one classical period

clc
clear
close all
%% Parameters of Simulation

a = -20;                        % Left end point 
b = +20;                        % Right end point 
L = b-a;                        % Width of the space
N = 512;                        % No. of cells
X = a+L*(0:N-1)/N;              % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum

T= 50*pi;                       % Duration of the evolution
M = 10^3;                       % Total number of steps in the evolution
dt = T/M;                       % Time step

X0 = 2.0;                       % Location of centre of wavepacket
sigma = 1.0;                    % Width of the initial wavepacket

A = 0.01;                       % Perturbation Amplitude
w = 0.99;                         % Perturbation Frequency
%% Defining the Hamiltionian
% We define 2 vectors that store the split-step propagators in position and 
% momentum space respectively.
% 
% Note that UV is only placed here for completeness. Since it depends on time, 
% it must be defined in the loop over time later in the code.

%UV = exp(-1i*(X.^2/2 + A*sin(X)*cos(w*time))*dt/2);   % One-step propagator in position space, only taking diagonal form
UT = exp(-1i*(P.^2/2)*dt);                             % One-step propagator in momentum space
%% Defining & Initialising the HO Eigenstates

poly_0 = hermiteH(0,X);
ho_0_prep = poly_0.*exp(-(X(1:N)-X0).^2/(2*sigma^2));   % Ground state
ho_0=ho_0_prep/sqrt(sum(abs(ho_0_prep).^2));            % State normalisation

poly_1 = hermiteH(1,X);
ho_1_prep = poly_1.*exp(-(X(1:N)-X0).^2/(2*sigma^2));   % First excited state
ho_1=ho_1_prep/sqrt(sum(abs(ho_1_prep).^2));            % State normalisation

poly_2 = hermiteH(2,X);
ho_2_prep = poly_2.*exp(-(X(1:N)-X0).^2/(2*sigma^2));   % Second excited state
ho_2=ho_2_prep/sqrt(sum(abs(ho_2_prep).^2));            % State normalisation

poly_3 = hermiteH(3,X);
ho_3_prep = poly_3.*exp(-(X(1:N)-X0).^2/(2*sigma^2));   % Third excited state
ho_3=ho_3_prep/sqrt(sum(abs(ho_3_prep).^2));            % State normalisation

poly_4 = hermiteH(4,X);
ho_4_prep = poly_4.*exp(-(X(1:N)-X0).^2/(2*sigma^2));   % Fourth excited state
ho_4=ho_4_prep/sqrt(sum(abs(ho_4_prep).^2));            % State normalisation
%% Animation of the Probability Density of the Wavefunction in Unperturbed Hamiltonian

UV = exp(-1i*(X.^2/2)*dt/2);

figure;
wave = plot(X(1:N),abs(psi(1:N)).^2);   % plotting initial state
ylim([0 0.15])
hold on 
plot (P(1:N),abs(psi(1:N)).^2) 
psi_0=ho_2;

drawnow
for m = 1:M
    psi_1 = UV.*psi_0;
    phi_2 = fft(psi_1); %wavefunction in momentum space
    phi_3 = UT.*phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV.*psi_3;
    psi_0 = psi_4;      %prepare a new cycle
    set(wave, 'YData', abs(psi_0(1:N)).^2)
    pause(0.05)
end
psi=psi_0;              %final state updated 
%% Transition Probability (1)

psi_0 = ho_0;           % this is the state which evolves
phi_0 = ho_1;           % this is the state which we wish to compare psi_0 against

t = [0];
y = [0];
P = [0];

for m = 1:M

    time = t(end) + dt;
    num_prob = abs(dot(psi_0, phi_0))^2;
    %theor_prob = (A^2/8)*exp(-1/2)*abs((exp(1i*(1-w)*time)-1)/(1-w))^2;       % hardcoded the transition probability derived from theory
    theor_prob = (A^2/8)*exp(-1/2)*abs( (exp(1i*(1-w)*time)-1)/(1-w) + (exp(1i*(1+w)*time)-1)/(1+w) )^2;       % Non-RWA version

    t(end+1) = time;
    y(end+1) = num_prob;
    P(end+1) = theor_prob;
    
    UV = exp(-1i*(X.^2/2 + A*sin(X)*cos(w*time))*dt/2);

    psi_1 = UV.*psi_0;
    psi_2 = fft(psi_1);
    psi_3 = UT.*psi_2;
    psi_4 = ifft(psi_3);
    psi_5 = UV.*psi_4;

    psi_0 = psi_5;
end


str1 = sprintf('A = %.2f, w = %.2f (sim)', A, w);
plot(t, y, ...
    'Color', '#E74C3C', ...
    'LineWidth', 2, ...
    'DisplayName', str1 ...
    );                % plot of simulated transition pobability
hold on

str2 = sprintf('A = %.2f, w = %.2f (theory)', A, w);
plot(t, P, ...
    'Color', '#9B59B6', ...
    'LineWidth', 2, ...
    'DisplayName', str2 ...
    );              % plot of theoretically derived transition probability
hold on


legend('FontSize',12);
xlabel('Time, T', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability, P', 'FontSize', 14, 'FontWeight', 'bold');

ax = gca;
ax.FontSize = 12;

grid on
%% Transition Probability (2)

psi_0 = ho_0;           % this is the state which evolves
phi_0 = ho_1;           % this is the state which we wish to compare psi_0 against

t = [0];
y = [0];
P = [0];

for m = 1:M

    time = t(end) + dt;
    num_prob = abs(dot(psi_0, phi_0))^2;
    %theor_prob = (A^2/8)*exp(-1/2)*abs((exp(1i*(1-w)*time)-1)/(1-w))^2;       % hardcoded the transition probability derived from theory
    theor_prob = (A^2/8)*exp(-1/2)*abs( (exp(1i*(1-w)*time)-1)/(1-w) + (exp(1i*(1+w)*time)-1)/(1+w) )^2;       % Non-RWA version

    t(end+1) = time;
    y(end+1) = num_prob;
    P(end+1) = theor_prob;
    
    UV = exp(-1i*(X.^2/2 + A*sin(X)*cos(w*time))*dt/2);

    psi_1 = UV.*psi_0;
    psi_2 = fft(psi_1);
    psi_3 = UT.*psi_2;
    psi_4 = ifft(psi_3);
    psi_5 = UV.*psi_4;

    psi_0 = psi_5;
end


str1 = sprintf('A = %.2f, w = %.2f (sim)', A, w);
plot(t, y, ...
    'Color', '#2980B9', ...
    'LineWidth', 2, ...
    'DisplayName', str1 ...
    );              % plot of simulated transition pobability
hold on

str2 = sprintf('A = %.2f, w = %.2f (theory)', A, w);
plot(t, P, ...
    'Color', '#1ABC9C', ...
    'LineWidth', 2, ...
    'DisplayName', str2 ...
    );              % plot of theoretically derived transition probability
hold on


legend('FontSize',12);
xlabel('Time, T', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability, P', 'FontSize', 14, 'FontWeight', 'bold');

ax = gca;
ax.FontSize = 12;

grid on

%ylim([0,0.1])
%% Transition Probability (3)

psi_0 = ho_0;           % this is the state which evolves
phi_0 = ho_1;           % this is the state which we wish to compare psi_0 against

t = [0];
y = [0];
P = [0];

for m = 1:M

    time = t(end) + dt;
    num_prob = abs(dot(psi_0, phi_0))^2;
    theor_prob = (A^2/8)*exp(-sigma/2)*sigma^(3)*abs((exp(1i*(1-w)*time)-1)/(1-w))^2;       % hardcoded the transition probability derived from theory

    t(end+1) = time;
    y(end+1) = num_prob;
    P(end+1) = theor_prob;
    
    UV = exp(-1i*(X.^2/2 + A*sin(X)*cos(w*time))*dt/2);

    psi_1 = UV.*psi_0;
    psi_2 = fft(psi_1);
    psi_3 = UT.*psi_2;
    psi_4 = ifft(psi_3);
    psi_5 = UV.*psi_4;


    psi_0 = psi_5;
end


str1 = sprintf('A = %.2f, w = %.2f (sim)', A, w);
plot(t, y, ...
    'Color', '#F1C40F', ...
    'LineWidth', 2, ...
    'DisplayName', str1 ...
    );              % plot of simulated transition pobability
hold on

str2 = sprintf('A = %.2f, w = %.2f (theory)', A, w);
plot(t, P, ...
    'Color', '#27AE60', ...
    'LineWidth', 2, ...
    'DisplayName', str2 ...
    );              % plot of theoretically derived transition probability
hold on


legend('FontSize',12);
xlabel('Time, T', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability, P', 'FontSize', 14, 'FontWeight', 'bold');

ax = gca;
ax.FontSize = 12;

grid on
%% Fourier Transform of Probability 

hold off

K = fft(P);

fs = 1/dt;
f = (0:length(K)-1)*fs/length(K);

plot(f,K, ...
    'Color', '#27AE60', ...
    'LineWidth', 2 ...
    );

xlabel('Frequency, f', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Wave Amplitude, K', 'FontSize', 14, 'FontWeight', 'bold');

xlim([0,1])
ylim([-0.06,0.06])

grid on