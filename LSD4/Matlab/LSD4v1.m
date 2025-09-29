%% 4DM90 Structural dynamics LSD 4: Impedence coupling
clear; clc; close all;

%% Parameters
m1 = 2; m2 = 2; m3 = 1; m4 = 1;
k1 = 100; k2 = 50; k3 = 50; k4 = 10;
b1 = 0.8; b2 = 0.8; b3 = 0.4; b4 = 0.4;

% Frequency range
f = 0.1:0.001:1.6;
w = f * (2*pi);

%% subsystem 1 matrices
M1 = [m2 0;
     0  m1];

B1 = [b2     -b2;
     -b2    b1+b2];

K1 = [k2     -k2;
     -k2    k1+k2];

%% FRF calculation for subsystem 1
n = length(w);
H1 = zeros(2,2,n);

for i = 1:n
    A = -w(i)^2*M1 + 1i*w(i)*B1 + K1;
    H1(:,:,i) = inv(A);
end

%% Plot all FRFs of subsystem 1
index2 = [1 1; 1 2; 2 1; 2 2];
plotFullFRF(f,H1,index2,"FRF of subsystem 1");

%% Subsystem 2 matrices
M2 = [0 0 0;
    0 m3 0;
    0 0 m4];
B2 = [b3 -b3 0;
    -b3 b3+b4 -b4;
    0 -b4 b4];
K2 = [k3 -k3 0;
      -k3 k3+k4 -k4;
      0 -k4 k4];

%% FRF calculation for subsystem 2
H2 = zeros(3,3,n);

for i = 1:n
    A = -w(i)^2*M2 + 1i*w(i)*B2 + K2;
    H2(:,:,i) = A\eye(3,3);
end

%% Plot all FRFs of subsystem 2
index3 = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3];
plotFullFRF(f,H2,index3,"FRF of subsystem 2");

%% Complete system matrices
M = diag([m1,m2,m3,m4]);
B = [b1+b2      -b2     0       0;
      -b2       b2+b3   -b3     0;
      0         -b3     b3+b4   -b4;
      0         0       -b4     b4;];

K = [k1+k2      -k2     0       0;
     -k2        k2+k3   -k3     0;
     0          -k3     k3+k4   -k4;
     0          0       -k4     k4];

%% FRF calculation for complete system
H = zeros(4,4,n);

for i = 1:n
    A = -w(i)^2*M + 1i*w(i)*B + K;
    H(:,:,i) = A\eye(4,4);
end

%% Plot all FRFs of complete system
index4 = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];
plotFullFRF(f,H,index4,"FRF of complete system");

%% Determine the complete FRF using impedence coupling
% Reorganise FRF matrices
H1_BB = H1(1,1,:);
H1_BI = H1(1,2,:);
H1_IB = H1(2,1,:);
H1_II = H1(2,2,:);

H1_ic = [H1_BB, H1_BI;
         H1_IB, H1_II];

H2_BB = H2(1,1,:);
H2_BI = H2(1,2:3,:);    
H2_IB = H2(2:3,1,:);   
H2_II = H2(2:3,2:3,:); 

H2_ic = [H2_BB, H2_BI;
         H2_IB, H2_II];

H_ic = impedanceCoupling(H1_ic, H2_ic);
%% Plot the FRF of the impendence coupling system
index_ic = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];
% index_ic = ["B" "B"; "B" "I"; "B" "I"; "B" "I"; "I" "B"; "I" "I"; "I" "I"; "I" "I"; "I" "B"; "I" "I"; "I" "I"; "I" "I"; "I" "B"; "I" "I"; "I" "I"; "I" "I"];
plotFullFRF(f,H_ic,index_ic,"FRFs of impedence coupled system");

%% Add noise to the subsystem and complete system
% Define Gaussian noise with standard deviation of 0.001
sd = 0.001;
noise1 = sd * randn(size(H1));
noise2 = sd * randn(size(H2));
noise_complete = sd * randn(size(H));

% Add the noise to each subsystem
H1_noise = H1 + noise1; 
H2_noise = H2 + noise2;
H_noise = H + noise_complete;

%% Plot all FRFs of the subsystems and complete system with noise
index2 = [1 1; 1 2; 2 1; 2 2];
plotFullFRF(f, H1_noise, index2, "FRF of subsystem 1 with noise");

index3 = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3];
plotFullFRF(f,H2_noise,index3,"FRF of subsystem 2 with noise");

index4 = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];
plotFullFRF(f,H_noise,index4,"FRF of complete system with noise");

%% Determine the complete FRF using impedence coupling based on subsystems with noise
% Reorganise FRF matrices
H1_BB_noise = H1_noise(1,1,:);
H1_BI_noise = H1_noise(1,2,:);
H1_IB_noise = H1_noise(2,1,:);
H1_II_noise = H1_noise(2,2,:);

H1_ic_noise = [H1_BB_noise, H1_BI_noise;
               H1_IB_noise, H1_II_noise];

H2_BB_noise = H2_noise(1,1,:);
H2_BI_noise = H2_noise(1,2:3,:);    
H2_IB_noise = H2_noise(2:3,1,:);   
H2_II_noise = H2_noise(2:3,2:3,:); 

H2_ic_noise = [H2_BB_noise, H2_BI_noise;
               H2_IB_noise, H2_II_noise];

H_ic_noise = impedanceCoupling(H1_ic_noise, H2_ic_noise);

%% Plot the impendence coupled FRF with noise
index_ic = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];
name_index = [];
plotFullFRF(f,H_ic_noise,index_ic,"FRFs of impedence coupled system with noise");


%% Function definitions
function plotFullFRF(f, H, index, plotTitle)
    % plot all FRFs of the given FRF matrix in a square grid on log log
    % scale. 
    fig = figure;
    for i = 1:size(index,1)
        subplot(size(H,1), size(H,1), i);
        semilogx(f, squeeze(mag2db(abs(H(index(i,1), index(i,2), :)))));
        xlabel('Frequency [Hz]');
        ylabel('|H_{' + string(index(i,1)) + string(index(i,2)) + '}| [dB]');
        grid on;
    end
    sgtitle(plotTitle);
    cd Export_graphics\;
    exportgraphics(fig, string(plotTitle) +'.pdf','Resolution',1200, Padding=5);
    cd ../;
end

function H_ic = impedanceCoupling(H1, H2)
%   ImpedanceCoupling of two subsystems using the 1 inverse metho
%
%   Inputs:
%       H1 : (n1 x n1 x N) FRF of subsystem 1
%       H2 : (n2 x n2 x N) FRF of subsystem 2
%            Both must share 1 common boundary DOF (row/col 1).
%            Internal DOFs are in rows/cols 2:end.
%
%   Output:
%       H_ic : ( (n1-1 + n2-1 + 1) x (n1-1 + n2-1 + 1) x N ) coupled FRF
%
%   The result is the impedance-coupled FRF for all frequencies.

    N = size(H1,3);      % number of frequency points
    n1i = size(H1,1)-1;  % # of internal DOFs in H1
    n2i = size(H2,1)-1;  % # of internal DOFs in H2

    % Preallocate result
    H_ic = zeros(1+n1i+n2i, 1+n1i+n2i, N);

    for k = 1:N
        % Partition H1
        H1BB = H1(1,1,k);
        H1BI = H1(1,2:end,k);
        H1IB = H1(2:end,1,k);
        H1II = H1(2:end,2:end,k);

        % Partition H2
        H2BB = H2(1,1,k);
        H2BI = H2(1,2:end,k);
        H2IB = H2(2:end,1,k);
        H2II = H2(2:end,2:end,k);

        % First block matrix
        first_term = [ H1BB      H1BI        zeros(1,n2i);
                       H1IB      H1II        zeros(n1i,n2i);
                       zeros(n2i,1) zeros(n2i,n1i) H2II ];

        % Coupling correction (outer product)
        v = [ H1BB;
              H1IB;
             -H2IB ];

        second_term = v * ( (H1BB + H2BB) \ transpose(v) );

        % Store result
        H_ic(:,:,k) = first_term - second_term;
    end
end
