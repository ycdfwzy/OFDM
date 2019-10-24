function [data] = decode_received_signal(received_signal)
%DECODE_RECEIVED_SIGNAL 此处显示有关此函数的摘要
%   return decoded data from received data
%   此处显示详细说明

%% Parameters

R = 1e5; % [bits/sec]
duration = 0.128; % [sec]
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

fltDelay = Nsym / (2*R); % Filter group delay, since raised cosine filter is linear phase and symmetric.

%% Raised cosine filter design

shape = 'Raised Cosine';
% Specifications of the raised cosine filter with given order in symbols
rcosSpec = fdesign.pulseshaping(sampsPerSym, shape, 'Nsym,beta', Nsym, beta);
rcosFlt = design(rcosSpec);
rcosFlt.Numerator = rcosFlt.Numerator / max(rcosFlt.Numerator);

%% A/D
t = (0 : length(received_signal)-1) / Fs;
signal_complex = received_signal.*cos(2*pi*f1*t) - 1i * received_signal.*sin(2*pi*f1*t);

filtered_signal = filter(rcosFlt, [signal_complex.'; zeros(fltDelay*Fs,1)]); 
reconstructed_signal = downsample(filtered_signal, sampsPerSym);
reconstructed_signal = reconstructed_signal(Nsym/2+1:end); % Correct for propagation delay by removing filter transients

%% Serial to Parallel
fft_signal = [];
for i = 1:length(reconstructed_signal)/(NFFT+L)
    tmp = reconstructed_signal((i-1)*(NFFT+L)+1 : i*(NFFT+L));
    removing_cyclic_prefix = tmp(end-NFFT+1:end); % removing cyclic prefix
    fft_signal_cur = fft(removing_cyclic_prefix); % fft
    fft_signal = [fft_signal; fft_signal_cur];
end

demodulated_data = qamdemod(fft_signal, M);
data = de2bi(demodulated_data, 'left-msb');
data = data(:);
end

