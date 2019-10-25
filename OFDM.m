clear all;
close all;
clc;

bernoulli_binary_generator = randi([0,1], 12800, 1);
sending_signal = generate_sending_signal(bernoulli_binary_generator);
sending_signal = sending_signal./sum(sending_signal.^2);
% calcFreq(sending_signal, 7000);
% [y,Fs] = audioread('C:\Users\ycdfwzy\Desktop\20191024_161158.m4a');
% figure(1);
% plot(t(1:1000),y1(1:1000));
noised_signal = awgn(sending_signal, 20, 'measured');
% calcFreq(noised_signal, 100000);
noised_signal = noised_signal./sum(noised_signal.^2);
received_data = decode_received_signal(noised_signal);
BER = sum(xor(bernoulli_binary_generator, received_data))/length(received_data);
audiowrite('test.wav', sending_signal, 7000);