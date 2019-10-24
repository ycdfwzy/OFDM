clear all;
close all;
clc;

bernoulli_binary_generator = randi([0,1], 12800, 1);
sending_signal = generate_sending_signal(bernoulli_binary_generator);
audiowrite('test.wav', sending_signal, 48000);
% figure(1);
% plot(t(1:1000),y1(1:1000));
noised_signal = awgn(sending_signal, 20, 'measured');
received_data = decode_received_signal(noised_signal);
BER = sum(xor(bernoulli_binary_generator, received_data));
[y,Fs] = audioread('test.wav');
sound(y,Fs);