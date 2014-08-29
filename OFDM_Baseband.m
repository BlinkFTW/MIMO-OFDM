%% OFDM Baseband code
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

% % Load 32 char string
% load('Test_Strings.mat','QPSKstr'); % This will allow for a single OFDM char
% cAr = double(QPSKstr);
% %cAr = double('A');
%% Convert to binary
binDat = Binary_cnvrt(cAr);

%% Map to QPSK symbol
symDat = QPSK.bin2symb(binDat,1); % symbols are [1 -1 j -j]
% correspond to [0 180 90 270] degrees
%% OFDM symbols
N = 128;
O_data = idft(symDat,N);

%% Cyclic Prefix
% cpS = floor(N*0.25);    % CP size
% xp = [O_data(end-cpS+1:end,:); O_data]; % Parallel data

% No CP
xp = O_data;
%% Parallel to Serial
xDim = size(xp);        % Matrix dimensions
% xDim(1): OFDM symbol length
% xDim(2): # of symbols
xn = reshape(xp,1,xDim(1)*xDim(2));
n = 1:length(xn);
% pad zeros between each sample, then use Nyquist root cosine filter?

% Average TX power
Pt = mean(xn.*conj(xn));    % E{|x[i]|^2}

%%
SNR_dB = 12.598;
SNR = 10^(SNR_dB/10);

%% %%%%%%%%%%%%%%%%%% Channel %%%%%%%%%%%%%%%%%%%%%
% Rayleigh (block) fading with AWGN
h = 1;  % Channel Gain
% h = sqrt(1/2).*(randn(1)+1i*randn(1));

% Pn = Pt*(h.*conj(h))/SNR;    % Average Noise Power
Pn = Pt/SNR;    % Average Noise Power
noise = sqrt(Pn/2).*(randn(1,length(xn))+1i*randn(1,length(xn)));

%% %%%%%%%%%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%%

y = h.*xn + noise;

%% FFT then Parallel to Serial

yp = reshape(y,xDim(1),xDim(2));
ySym = dft(yp,N);
% Assume Receiver knows the length of unencoded message
ySym = ySym(1:length(symDat));  % Remove extra data added from idft from TX
% yn = reshape(ySym,1,xDim(1)*xDim(2));

%% Show Symbol Constellation
figure(2)
QPSK.Constellation(ySym, SNR_dB);

%% De-map symbols and convert to binary

[bin2, symDat2] = QPSK.sym2bin(ySym);

%% Calculate SER & BER

SER = mean(symDat2 ~= symDat);
BER = mean(binDat ~= bin2);
% 
% %% Convert to ASCII character data
% 
% dOut = bin2ascii(bin2)';
% 
% %% Number of error characters
% 
% ChER = mean(dOut ~= cAr);

%% Compare STDev to Noise Variance
