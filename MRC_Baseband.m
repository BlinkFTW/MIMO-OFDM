%% SIMO Maximal Ratio Combining Simulation (Baseband)
% Eric Kang
% April 26, 2014
%
% Parameters:
% N = 128;      % Number of Sub-carriers (Powers of 2)
% SNR_dB = 10;  % Eb/N0 (dB)
% Mrx = 10;     % Number of RX antenna
%
clear all; clc

% Load Function libraries
QPSK = QPSK_lib;
MIMO = MIMO_lib;
OFDM_lib;

%% %%%%%%%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%%
%% Input data

% Load text file
fid = fopen('usdeclar.txt');    % http://www.constitution.org/usdeclar.txt
cAr = fread(fid,inf);   % Read as ASCII decimal characters
fclose(fid);            % Each character is 1 byte

%% Convert to binary

binDat = Binary_cnvrt(cAr);

%% Map to QPSK symbols

symDat = QPSK.bin2symb(binDat,1); 

%% OFDM symbols

N = 128;
xp = idft(symDat,N);

%% Parallel to Serial

xDim = size(xp);        % Matrix dimensions
% xDim(1): OFDM symbol length
% xDim(2): Total # of symbols
xn = reshape(xp,1,xDim(1)*xDim(2));

%% %%%%%%%%%% Channel + Receiver %%%%%%%%%%%%%

SNR_dB = 12.598;
Mrx = 2;

% Maximum Ratio Combining
% Simulates channel propagation and MRC combiner
% x:    Transmitted signal
% M:    Number of RX antenna

Pt = mean(xn.*conj(xn));
%% Channel
snr = 10^(SNR_dB/10);
h = sqrt(1/2).*(randn(1,Mrx)+1i*randn(1,Mrx)).';  % Channel Coeff

% h = ones(Mrx,1);
Pn = Pt/snr;
noise = sqrt(Pn/2).*(randn(Mrx,length(xn))+1i.*randn(Mrx,length(xn)));

%% Rx
y = h*xn+noise;
a = conj(h).';    % Combining Coeff
z = a*y;    % Combiner
z_eq = 1/sum(a.'.*conj(a.')).*z; % Equalization

% y = MIMO.MRC(xn,Mrx,SNR_dB);    % y is Serial

%% FFT then Parallel to Serial

yp = reshape(z_eq,xDim(1),xDim(2));    % Serial to Parallel
ySym = dft(yp,N);
% Assume Receiver knows the length of unencoded message
ySym = ySym(1:length(symDat));  % Remove extra data added from idft from TX

%% Show Symbol Constellation

close all
figure(1)
QPSK.Constellation(ySym, SNR_dB);

%% De-map symbols and convert to binary

[bin2, symDat2] = QPSK.sym2bin(ySym);

%% Calculate SER & BER

SER = mean(symDat2 ~= symDat);
BER = mean(binDat ~= bin2);

%% Convert to ASCII character data

dOut = bin2ascii(bin2)';

%% Number of error characters

ChER = mean(dOut ~= cAr);
