%% 4DM90 Structural dynamics LSD 5: DFT/FFT
clear; clc; close all;

%% Define signals
delta_t = 0.01;
T = 5;
delta_f = 1 / T;
N = T / delta_t;
% N = 500;
t = (0:N-1)*delta_t;
fs = 1 / delta_t;

signals = zeros(6,length(t));

% A harmonic signal with f = 1 Hz
signals(1,:) = sin((2*pi) * 1 * t);

% A signal being the sum of two harmonic functions with f1 = 1 Hz and f2 =
% 2 Hz and 90° phase shift
signals(2,:) = signals(1,:) + sin(((2*pi) * 2 * t) + deg2rad(90)); 

% A square wave signal, where each square has a length of 1 sec.
count = 0;
high = true;
for i = 1:length(t)
    if count == 100
        if high == true
            high = false;
        elseif high == false
            high = true;
        end
        count = 0;
    end
    if high == true
        signals(3,i) = 1;
    end
    count = count + 1;
end

% A sawtooth signal
for i = 1:length(t)
    signals(4,i) = mod(i, 100) / 100; % Normalizing to range [0, 1]
end

% A Dirac delta function
signals(5,1) = 1;

% A normally distributed (or Gaussian-distributed white noise signal with
% standard devitation of 1
signals(6,:) = randn(1, length(t));

%% Plot all signals in subplots
figure; 
for i = 1:size(signals,1)
    subplot(2,3,i);
    plot(t,signals(i,:));
    xlabel("time [s]");
    ylabel("amplitude [m]");
    title('Signal '+ string(i));
    ylim([1.5*min(signals(i,:)) 1.5*max(signals(i,:))]);
    grid on;
end

%% Perfrom the fft of all signals 
fftSignals = 1/N*fft(signals, [], 2);
f = (0:N-1) * delta_f; % Frequency vector

% Plot the magnitude of the FFT results
figure;
for i = 1:size(fftSignals, 1)
    subplot(2, 3, i);
    h = stem(f, abs(fftSignals(i,:)));
    h.Marker = '.';        % use a dot instead of a circle
    h.MarkerSize = 6;      % smaller dot
    h.LineWidth = 0.1;       % thinner stem lines
    xlabel("Frequency [Hz]");
    ylabel("Magnitude [m]");
    title('FFT of Signal ' + string(i));
    xlim([0 max(f)]);
    ylim([0 1.2*max(abs(fftSignals(i,:)))]);
    grid on;
end


