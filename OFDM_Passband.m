%% OFDM Passband code
% Process:
% - input text data
% - convert text to binary
% - convert binary to symbol
% - OFDM idft
clear all; close all; clc

% Load Function libraries
QPSK = QPSK_lib;
OFDM_lib;
%% Input text data

% Load text file
fid = fopen('usdeclar.txt');
cAr = fread(fid,inf);   % Read as ASCII decimal characters
fclose(fid);
% 
% load('txtDocs.mat','moby');
% cAr = moby;
% % Load 32 char string
% load('Test_Strings.mat','QPSKstr'); % This will allow for a single OFDM char
% cAr = double(QPSKstr(1:end-3));
% %cAr = double('A');
%% Convert to binary

binDat = Binary_cnvrt(cAr);

%% Map to QPSK symbol

symDat = QPSK.bin2symb(binDat,1); % symbols are [1 -1 j -j]
                                % correspond to [0 180 90 270] degrees
%% OFDM symbols

N = 128;
O_data = idft(symDat,N);

%% No CP

xp = O_data;

%% Parallel to Serial

xDim = size(xp);        % Matrix dimensions
                        % xDim(1): OFDM symbol length
                        % xDim(2): # of symbols
xn = reshape(xp,1,xDim(1)*xDim(2));
n = 1:length(xn);
% pad zeros between each sample, then use Nyquist root cosine filter?

%% D/A conversion - can take a while for large data

os = 1;  % Over Sample = (Fs/Rs)
xt = interp(xn,os);
figure(1)
stem(0:length(xn)-1,xn)
figure(2)
plot(0:length(xt)-1,xt)

% Check signal's spectrum s.t. it is band limited
Xf = fftshift(fft(xt));
figure(3)
plot((0:length(Xf)-1)-length(Xf)/2,sqrt(Xf.*conj(Xf)))
%% Power Amplification

%%
% %Bandpass Modulation - from SAM proj. Fall 2012
% carrier = sqrt(Rb/Fs)*cos(2*pi*Fc*t);
% bp_op = bb_op*carrier;


%% %%%%%%%%%%%%%%%%%% Channel %%%%%%%%%%%%%%%%%%%%%
% Average TX power
Pt = mean(xt.*conj(xt));    % E{|x[i]|^2}

% Rayleigh (block) fading with AWGN
h = 1;  % Channel Gain
SNR_dB = 12.598;
SNR = 10^(SNR_dB/10);
Es = 2;   % Symbol Energy for QPSK
Pn = (h*Pt)/SNR;    % Average Noise Power
noise = sqrt(Pn/Es).*(randn(1,length(xt))+1i*randn(1,length(xt)));
NoB = mean(noise.*conj(noise)); % Measured Noise Power
%% %%%%%%%%%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%%

y = h.*xt + noise;

%% A/D conversion

yn = y(1:os:end);   % Decimate by factor of 10

% Assume receiver knows xDim(2)
yp = reshape(yn,xDim(1),xDim(2));

%% FFT then Parallel to Serial

ySym = dft(yp,N);
% Assume Receiver knows the length of unencoded message
ySym = ySym(1:length(symDat));  % Remove extra data added from idft from TX

%% Show Symbol Constellation
figure(4)
QPSK.Constellation(ySym, SNR_dB);

%% De-map symbols and convert to binary

[bin2, symDat2] = QPSK.sym2bin(ySym);

%% Calculate SER & BER

SER = sum(symDat2 ~= symDat)/length(symDat);
BER = sum(binDat ~= bin2)/length(binDat);

% %% Convert to ASCII character data - can take a while for large data
% 
% dOut = bin2ascii(bin2)';
% 
% %% Number of error characters
% 
% ChER = sum(dOut ~= cAr)/length(dOut);
