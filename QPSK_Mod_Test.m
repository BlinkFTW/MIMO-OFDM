%%
% Demo QPSK modulation
clear all; close all; clc

inStr = 'Hello';
decStr = double(inStr);

binDat = Binary_cnvrt(decStr);

% Map 0 bits to -1
for j = 1:length(binDat)
    if binDat(j) == 0
        binDat(j) = -1;
    end
end

sym_data = QPSK_bin2symb(binDat);

fs = 1000;
fc = 10;
N = 128;    % symbol length
t_end = 2;  % symbol period in seconds
t = 0:1/fs:t_end-1/fs;
k = 0:N-1;

% QPSK modulated data stream
Q_data = zeros(1,N*length(sym_data));
for j = 0:length(sym_data)-1
    Q_data(N*j+1:N*j+N) = exp(1i*(2*pi*fc.*k/fs+angle(sym_data(j+1))));
end
plot(t,real(Q_data))
% I = cos(2*pi*fc.*k/fs);
% Q = sin(2*pi*fc.*k/fs);
% s10 = 1-1i;
% c10 = exp(1i*(2*pi*fc.*t+angle(s10)));%*exp(angle(s10));

% figure(3)
% subplot(211)
% plot(k,I+Q)
% title('phase shift of -\pi/4')
% subplot(212)
% plot(k,c10)

%%
% f = -fs/2:1/t(end):fs/2;
% Xf = fftshift(fft(c10));
% Xm = real(Xf);
% Xa = angle(Xf);
% figure(5)
% stem(f,Xm)
% title('Magnitude')

% %%
% cc10 = ifft(Xf);
% figure(6)
% plot(t,real(cc10))


%%
% figure(1)
% subplot(211)
% plot(t,I)
% title('I component')
% subplot(212)
% plot(t,Q)
% title('Q component')

% symDat = zeros(1,length(binDat)*length(t)/2);
% for i = 1:2:length(binDat)
%     
% end

% figure(2)
% subplot(221)
% plot(t,Q-I)
% title('01')
% subplot(222)
% plot(t,I+Q)
% title('11')
% subplot(223)
% plot(t,-I-Q)
% title('00')
% subplot(224)
% plot(t,I-Q)
% title('10')
