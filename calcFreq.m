function [] = calcFreq(data, Fs)
%CALCFREQ �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    N = length(data);
    n = 0:N-1;
    y = fft(data);
    mag = abs(y);
    frequencies = n * Fs/N;
    plot(frequencies, mag);    %�����Ƶ�ʱ仯�����
    xlabel('Feq/Hz');
    ylabel('Am');
    grid on;
end

