%% LSD 4 4DM90 Structural dynamics
clear; clc; close all;

%% Parameters
m1 = 2; m2 = 2; m3 = 1; m4 = 1;
k1 = 100; k2 = 50; k3 = 50; k4 = 10;
b1 = 0.8; b2 = 0.8; b3 = 0.4; b4 = 0.4;

% Frequency range
f = 0.1:0.01:1.6;
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
figure; hold on;
index1 = [1 1; 1 2; 2 1; 2 2];
for i = 1:4
    subplot(2,2,i);
    loglog(f,squeeze(abs(H1(index1(i,1),index1(i,2),:))));
    xlabel('Frequency [Hz]');
    ylabel('|H_{' + string(index1(i,1)) + string(index1(i,2)) + '}|');
    grid on;
end
sgtitle('FRFs of subsystem 1');
hold off;

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
figure; hold on;
index2 = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3];
for i = 1:9
    subplot(3,3,i);
    loglog(f,squeeze(abs(H2(index2(i,1),index2(i,2),:))));
    xlabel('Frequency [Hz]');
    ylabel('|H_{' + string(index2(i,1)) + string(index2(i,2)) + '}|');
    grid on;
end
sgtitle('FRFs of subsystem 2');
hold off;

%% Complete system matrices
M = diag([m1,m2,m3,m4]);
B = [b1+b2      -b2     0       0;
      -b2       b2+b3   -b3     0;
      0         -b3     b3+b4   -b4;
      0         0       -b4     b4;]

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
figure; hold on;
index = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];
for i = 1:16
    subplot(4,4,i);
    loglog(f,squeeze(abs(H(index(i,1),index(i,2),:))));
    xlabel('Frequency [Hz]');
    ylabel('|H_{' + string(index(i,1)) + string(index(i,2)) + '}|');
    grid on;
end
sgtitle('FRFs of complete system');
hold off;
