clear all;
close all;
clc;

[~, standard_flag] = generate_sending_signal([]);
bernoulli_binary_generator = randi([0,1], 12800, 1);
[sending_signal, ~] = generate_sending_signal(bernoulli_binary_generator);
maxx = max(abs(sending_signal));
minx = min(sending_signal);
% calcFreq(sending_signal, 7000);
noised_signal = awgn([sending_signal(end-1312:end), sending_signal, sending_signal(1: 1000)], 20, 'measured');
% calcFreq(noised_signal, 100000);
% noised_signal = noised_signal./sum(noised_signal.^2);
received_data = decode_received_signal(noised_signal, standard_flag);
BER = sum(xor(bernoulli_binary_generator, received_data))/length(bernoulli_binary_generator);
disp(BER);

% test_data1 = randi([0,1], 1280, 1);
% test_data2 = randi([0,1], 1280, 1);
% test_data3 = [];
% [sending_signal1, signal_complex1] = generate_sending_signal(test_data1);
% [sending_signal2, signal_complex2] = generate_sending_signal(test_data2);
% [sending_signal3, signal_complex3] = generate_sending_signal(test_data3);
% noised_signal2 = [rand(1, 1000), sending_signal2]; % awgn(sending_signal2, 20, 'measured');
% decode_received_signal(noised_signal2, signal_complex3);

