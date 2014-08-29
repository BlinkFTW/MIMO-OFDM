%% Linear Minimum Mean Square Error
% System Parameters
% FFT N = 128
% M = 2
clear all; clc

% Load Function libraries
QPSK = QPSK_lib;
OFDM_lib;
%% Data
% Load text file
fid = fopen('usdeclar.txt');
cAr = fread(fid,inf);   % Read as ASCII decimal characters
fclose(fid);

% % % Load 32 char string
% load('Test_Strings.mat','QPSKstr'); % This will allow for a single OFDM char
% cAr = double(QPSKstr);
% 
%% Convert to binary
% 
binDat = Binary_cnvrt(cAr);

%binDat = ones(1,69480);

%% Map to QPSK symbols

symDat = QPSK.bin2symb(binDat,1); 

%% LMMSE Encoding
% |s1 s3|
% |s2 s4|
a = reshape(symDat,2,length(symDat)/2);

%% Split into Mt channels and perform IDFT

N = 128;
ch1 = idft(a(1,:),N);
ch2 = idft(a(2,:),N);

%% Parallel to Serial for each channel

chDim = size(ch1);
ch1 = reshape(ch1,1,chDim(1)*chDim(2));
ch2 = reshape(ch2,1,chDim(1)*chDim(2));
TX_2 = [ch1;ch2];

%% Channel
Mt = 2;
Mr = 2;
Pt = mean(mean(TX_2.'.*conj(TX_2.')));  % Es/Mt

%%
SNR_dB = 20;
SNR = 10^(SNR_dB/10);
H = sqrt(1/2).*(randn(Mr,Mt)+1i*randn(Mr,Mt));  % Channel Coeff
%H = eye(Mr,Mt);
Pn = Pt/SNR;
% load('ZF_n.mat');
n = sqrt(Pn/2).*(randn(Mr,length(ch1))+1i.*randn(Mr,length(ch1)));
N0 = mean(mean(n.'.*conj(n.')));

%% Receiver

RX_2 = H*TX_2+n;

%% Serial to Parallel

y1 = reshape(RX_2(1,:),chDim(1),chDim(2));
y2 = reshape(RX_2(2,:),chDim(1),chDim(2));

%% DFT

ySym1 = dft(y1,chDim(1));
ySym2 = dft(y2,chDim(1));

ySym1 = ySym1(1:length(symDat)/2);  % Remove extra data added from idft from TX
ySym2 = ySym2(1:length(symDat)/2);

%% LMMSE Decoding

v2 = Mt/SNR;    % Mt*N0/Es
g = (H*H'+v2*eye(Mt))\H;
z = g'*[ySym1; ySym2];

xh = reshape(z,1,length(ySym1)*2);  % Estimation of x

%% Demapper
% De-map symbols and convert to binary

[bin2, symDat2] = QPSK.sym2bin(xh);

% Show Received symbol constellation
close all
figure(2)
QPSK.Constellation(xh, SNR_dB);
%% Calculate SER & BER

SER = mean(symDat2 ~= symDat);
BER = mean(binDat ~= bin2);
