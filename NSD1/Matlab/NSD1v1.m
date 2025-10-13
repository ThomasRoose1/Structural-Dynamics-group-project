%% NSD1: Vibration analysis of a beam with midplane stretching
clear; clc; close all; tic;

%% Load workspace if possible
workspace = "NSD1_workspace.mat";
if exist(workspace, 'file')
    load(workspace);
    fprintf('Loaded ' + workspace + '\n')
end

%% parameters
rho = 2700;                 %[kg/m^3] mass density
E = 70*10^9;                %[GPa] Elastic modulus
b = 0.02;                   %[m] width
h = 0.02;                   %[m] height
L = 1;                      %[m] length
c = 2.58;                   %[Ns/m] linear viscous damping
F = 939;                    %[N] Transversal force amplitude

A = b*h;                    %[m^2] area
m = 3*rho*A*L / 8;          %[kg] mass
I = (b*h^3)/12;             %[m^4] second moment of area
k1 = 2*pi^4*E*I/L^3;        %[N/m^3] equivalent linear stiffness
k3 = pi^4*E*A/(8*L^3);      %[N/m^3] equivalent cubic stiffness

%% (a) compute and plot the FRF of the linearized Duffing equation
f = 1:1:700;
w = f *(2*pi);
FRF = zeros(1,length(w));
for i = 1:length(w)
    den = -w(i)^2*m + 1i*w(i)*c + k1;
    FRF(:,i) = inv(den);
end

% Compute the resonance frequecny of the undamped Duffing equation
f_resonance = sqrt(k1/m) / (2*pi);

% plot both magnitude and phase
fig = figure;
semilogy(f, abs(FRF));
xlabel('Frequency [Hz]');
ylabel('Magnitude [N/m]');
title('Linear FRF $\frac{\hat{V}(\frac{L}{2},f)}{\hat{F}(f)}$', 'Interpreter', 'latex');
grid on;

% Export the plot as pdf
exportgraphics(fig, "Export_graphics\NSD1p2_linearFRF" +'.pdf','Resolution',1200, Padding=5);
clear fig

%% (b) Perform a slow sweep of the linearized Duffing equation with f = [1 700] Hz and compare with FRF
% Multiply the FRF with the given force to obtain the linear displacement
FRF_m = F * FRF;

% Define a frequency range from 1 to 50 Hz with frequency resolution of
% 0.25 Hz
f_start = 1; f_end = 700;
f_Delta = 1; 

% Define the number of samples of the transient and stable
N_t = 200; N_s = 200;

% Define initial conditions as 0 for both displacement and velocity
x0 = [0;0];

% Define ode options
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);

% % Define the system
alpha = k1/m; beta = 0; gamma = F/m; Delta = c/m;
sys = @(t,x,f) duffing(t,x,alpha,beta,gamma,Delta,f);

% Perform the sweep
[freq, A] = sweep(sys, x0, f_start, f_end, f_Delta, N_t, N_s, odeopts);

%% Plot results linear FRF + sweep up analysis
% plot both magnitude and phase
fig = figure; 
semilogy(f, abs(FRF_m), 'b');    
hold on;
semilogy(freq, A, 'r--');           
xlabel('Frequency [Hz]');
ylabel('Magnitude [m]');
title('FRF and sweep up analysis of linearized Duffing equation', 'Interpreter','latex');
legend('linear FRF', 'sweep up');
grid on;

% Export the plot as pdf
exportgraphics(fig, "Export_graphics\NSD1p2_linearFRF_displacement_sweepup" +'.pdf','Resolution',1200, Padding=5);
clear fig

%% (c) sweep up and down on f = [1 50]
f_start = 1; f_end = 50;
f_Delta = 0.25;

% Define the number of samples of the transient and stable
N_t = 200; N_s = 200;

% Define initial conditions as 0 for both displacement and velocity
x0 = [0;0];

% Define the system
alpha = k1/m; beta = k3/m; gamma = F/m; Delta = c/m;
sys_nl = @(t,x,f) duffing(t,x,alpha,beta,gamma,Delta,f);

% Perform the sweep
[freq_up, A_up] = sweep(sys_nl, x0, f_start, f_end, f_Delta, N_t, N_s, odeopts);

[freq_down, A_down] = sweep(sys_nl, x0, f_end, f_start, f_Delta, N_t, N_s, odeopts);

% plot both magnitude and phase
fig = figure;
semilogy(freq_up, A_up, 'b');
hold on;
semilogy(freq_down, A_down, 'r--');
xlabel('Frequency [Hz]');
ylabel('Magnitude [m]');
title('Sweep up and down of Duffing equation with $f \in [1, \, 50]$ Hz', 'Interpreter', 'latex');
legend('sweep up', 'sweep down');
grid on;

% Export the plot as pdf
exportgraphics(fig, "Export_graphics\NSD1p2_sweep_1_to_50" +'.pdf','Resolution',1200, Padding=5);
clear fig

%% (d) sweep up and down on f = [50 600]
f_start = 50; f_end = 600;
f_Delta = 1;

% Define the number of samples of the transient and stable
N_t = 200; N_s = 200;

% Define initial conditions as 0 for both displacement and velocity
x0 = [0;0];

% Define the system
alpha = k1/m; beta = k3/m; gamma = F/m; Delta = c/m;
sys_nl = @(t,x,f) duffing(t,x,alpha,beta,gamma,Delta,f);

% Perform the sweep
[freq_up, A_up] = sweep(sys_nl, x0, f_start, f_end, f_Delta, N_t, N_s, odeopts);

[freq_down, A_down] = sweep(sys_nl, x0, f_end, f_start, f_Delta, N_t, N_s, odeopts);

% plot both magnitude and phase
fig = figure; 
semilogy(freq_up, A_up, 'b'); 
hold on;
semilogy(freq_down, A_down, 'r--');
xlabel('Frequency [Hz]');
ylabel('Magnitude [m]');
title('Sweep up and down of Duffing equation with $f \in [50, \, 600]$ Hz', 'Interpreter', 'latex');
legend('sweep up', 'sweep down');
grid on;

% Export the plot as pdf
exportgraphics(fig, "Export_graphics\NSD1p2_sweep_50_to_600" +'.pdf','Resolution',1200, Padding=5);
clear fig

%% (e) sweep up with f = [400 450]and sweep down with f = [350 400] 
f_start = 400; f_end = 450;
f_Delta = 0.25;

% Define the number of samples of the transient and stable
N_t = 200; N_s = 200;

% Define initial condition for both displacement and velocity
x0 = [0.0176;8];

% Define the system
alpha = k1/m; beta = k3/m; gamma = F/m; Delta = c/m;
sys_nl = @(t,x,f) duffing(t,x,alpha,beta,gamma,Delta,f);

% Perform the sweep
[freq_up, A_up] = sweep(sys_nl, x0, f_start, f_end, f_Delta, N_t, N_s, odeopts);

% Define frequencies for sweep down
f_start = 400; f_end = 350;

[freq_down, A_down] = sweep(sys_nl, x0, f_start, f_end, f_Delta, N_t, N_s, odeopts);

% plot both magnitude and phase
fig = figure; 
semilogy(freq_up, A_up, 'b'); 
hold on;
semilogy(freq_down, A_down, 'r--');
xlabel('Frequency [Hz]');
ylabel('Magnitude [m]');
yscale([])
title('Sweep up and down of Duffing equation with $f \in [350, \, 450]$ Hz', 'Interpreter', 'latex');
legend('sweep up', 'sweep down');
grid on;

% Export the plot as pdf
exportgraphics(fig, "Export_graphics\NSD1p2_sweep_e" +'.pdf','Resolution',1200, Padding=5);
clear fig

%% (f) not finished

%% (g) - (k)
option = 'j';

% General parameters
N_s = 200; N_o = 200;
alpha = k1/m; beta = k3/m; gamma = F/m; Delta = c/m;
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);

switch option
    case 'g'
        % Specific parameters
        f_e = 34; N_t = 200; x0 = [0;0]; f_range = [0 500];
        sys_nl = @(t,x) duffing(t,x,alpha,beta,gamma,Delta,f_e);
    case 'h'
        % Specific parameters
        f_e = 37; N_t = 200; x0 = [0;0]; f_range = [0 500];
        sys_nl = @(t,x) duffing(t,x,alpha,beta,gamma,Delta,f_e); 
    case 'i'
        % Specific parameters
        f_e = 400; N_t = 800; x0 = [0.1;0]; f_range = [0 1000];
        sys_nl = @(t,x) duffing(t,x,alpha,beta,gamma,Delta,f_e);
    case 'j'
        % Specific parameters
        f_e = 400; N_t = 400; x0 = [0.0176;8]; f_range = [0 1000];
        sys_nl = @(t,x) duffing(t,x,alpha,beta,gamma,Delta,f_e);    
    case 'k'
        % Specific parameters
        f_e = 400; N_t = 1400; x0 = [0.018;8]; f_range = [0 1000];
        sys_nl = @(t,x) duffing(t,x,alpha,beta,gamma,Delta,f_e);    
end

% Call the simulate and plot, poincare_section and frequency_spectrum
% functions
T_e = 1 / f_e;
[t_ss, x_ss] = simulate_and_plot(sys_nl, x0, T_e, N_o, N_t, N_s, odeopts);
poincare_section(t_ss,x_ss,T_e,N_s);
[f1, Mag1s] = frequency_spectrum(t_ss, x_ss, T_e, f_range);


%% Save workspace for later use
save(workspace);
fprintf('Saved workspace as ' + workspace + '\n')    

%% Print elapsed time
elapsedTime = toc;
fprintf('Script completed in %.2f minutes.\n', elapsedTime/60);

