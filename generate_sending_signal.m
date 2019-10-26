function [signal, signal_complex] = generate_sending_signal(data)
%SENDER 此处显示有关此函数的摘要
%   data should be a binary array  size:n*1
%   return signal that will be sended
%   此处显示详细说明

%% Parameters

start_flag = [0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1]';
end_flag   = [1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0]';

R = 1e3; % [bits/sec]
duration = 12.8; % [sec]
DataL = R*duration/2;  % Data length in symbols
NFFT = 128;

M = 4; % 4 QAM
snr = 20; % [dB]

Nsym = 4;           % Filter order in symbol durations
beta = 0.5;         % Roll-off factor
sampsPerSym = 10;    % Upsampling factor
L = sampsPerSym*Nsym + 1; % Raised cosine filter order
Fs = R * sampsPerSym;   % Sampling frequency
f1 = Fs/5;
cyclic_prefix_signal = [];

%% Raised cosine filter design

% shape = 'Raised Cosine';
% Specifications of the raised cosine filter with given order in symbols
% b = rcosdesign(beta, Nsym, sampsPerSym);
% rcosSpec = fdesign.pulseshaping(sampsPerSym, shape, 'Nsym,beta', Nsym, beta);
% rcosFlt = design(rcosSpec);
% rcosFlt.Numerator = rcosFlt.Numerator / max(rcosFlt.Numerator);

%% modulate

data = [start_flag; data; end_flag];
reshaped_data = zeros(length(data)/2, 2);
for i = 1:length(data)/2
    reshaped_data(i, 1) = data(i*2-1);
    reshaped_data(i, 2) = data(i*2);
end
dec_data = bi2de(reshaped_data, 'left-msb');
modulated_data = qammod(dec_data, 4);

%% Serial to Parallel
for i = 1:length(modulated_data)/NFFT
    ifft_signal = ifft(modulated_data((i-1)*NFFT+1:i*NFFT));
    % cyclic prefix
    prefix = [ifft_signal(end-L+1:end); ifft_signal];
    cyclic_prefix_signal = [cyclic_prefix_signal; prefix];
end

%% D/A
% signal_complex = filter(rcosFlt, upsample([cyclic_prefix_signal; zeros(Nsym/2,1)], sampsPerSym));
% signal_complex = upfirdn(cyclic_prefix_signal, b, sampsPerSym);
signal_complex = rcupflt([cyclic_prefix_signal; zeros(Nsym/2,1)]);
fltDelay = Nsym / (2*R); % Filter group delay, since raised cosine filter is linear phase and symmetric.
signal_complex = signal_complex(fltDelay*Fs+1:end); % Correct for propagation delay by removing filter transients

I = real(signal_complex);
Q = imag(signal_complex);
t = (0 : length(signal_complex)-1) / Fs;
signal = I'.*cos(2*pi*t*f1) - Q'.*sin(2*pi*t*f1);
end
