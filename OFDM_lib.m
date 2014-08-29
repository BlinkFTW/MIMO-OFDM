function OFDM_lib
% Eric Kang
% April 18, 2014
% Contains all OFDM related functions
assignin('base','idft',@idft);
assignin('base','idftMIMO',@idftMIMO);
assignin('base','dft',@dft);
assignin('base','dftMIMO',@dftMIMO);
end

function sk = idft(c_ni,N_pt)
% This function computs the IFFT of the input serial symbol stream made up of
% symbols c_ni
%
% sk - output 
% c_ni - input serial complex symbol stream
% N_pt - specifies the number of IFFT points to compute and MUST BE A POWER
% OF TWO

%% Setup idft matrix

% w/o the use of guardbands
% symNum = ceil(length(c_ni)/N_pt);
% x = N_pt*symNum-length(c_ni);
% c_ni = [c_ni zeros(1,x)];
% iMat = reshape(c_ni,N_pt,symNum);

us = N_pt-floor(N_pt*0.1);  % Usable sub-carriers to assure a guardband that's 10% of BW
                            % from [LTE in a Nutshell]
symNum = ceil(length(c_ni)/us); % Number of OFDM symbols

% Append x zeros to c_ni
x = us*symNum-length(c_ni);
c_ni = [c_ni zeros(1,x)];

temp = reshape(c_ni,us,symNum);

iMat = zeros(N_pt,symNum);
% fG = ceil((N_pt*0.1-1)/2);  % Front guardband size
% bG = floor((N_pt*0.1-1)/2); % Back guardband size
% % Insert data around designated gaurdbands
% iMat(fG+1:fG+us/2,:) = temp(1:us/2,:);
% iMat(fG+us/2+2:end-bG,:) = temp(us/2+1:end,:);

iMat(1:us/2,:) = temp(1:us/2,:);
iMat(us/2+1+floor(N_pt*0.1):end,:) = temp(us/2+1:end,:);
%% perform ifft
sk = ifft(iMat);
end

function TX_M = idftMIMO(c_ni,N_pt)
% This function computs the IFFT of the input MxN symbol stream made up of
% complex valued symbols
%
% sk - output 
% c_ni - input MxN complex symbol stream
% N_pt - specifies the number of IFFT points to compute and MUST BE A POWER
% OF TWO

%% Setup idft matrix
cSize = size(c_ni);
us = N_pt-floor(N_pt*0.1);  % Usable sub-carriers to assure a guardband that's 10% of BW
                            % from [LTE in a Nutshell]
symNum = ceil(cSize(2)/us); % Number of OFDM symbols
TX_M = zeros(cSize(1),cSize(2));
for i= 1:cSize(1)
    % Append x zeros to c_ni
    x = us*symNum-cSize(2);
    c_ni(i,:) = [c_ni(i,:) zeros(1,x)];
    
    temp = reshape(c_ni(i,:),us,symNum);
    
    iMat = zeros(N_pt,symNum);
    % fG = ceil((N_pt*0.1-1)/2);  % Front guardband size
    % bG = floor((N_pt*0.1-1)/2); % Back guardband size
    % % Insert data around designated gaurdbands
    % iMat(fG+1:fG+us/2,:) = temp(1:us/2,:);
    % iMat(fG+us/2+2:end-bG,:) = temp(us/2+1:end,:);
    
    iMat(1:us/2,:) = temp(1:us/2,:);
    iMat(us/2+1+floor(N_pt*0.1):end,:) = temp(us/2+1:end,:);
    %% perform ifft
    sk = ifft(iMat);
    %% Parallel to Serial
    TX_M(i,:) = reshape(sk,1,N_pt,symNum);
end
end

function c_ni = dft(sk,N_pt)
% Computes the FFT of the input ODFM data

iMat = fft(sk);
% fG = ceil((N_pt*0.1-1)/2);  % Front guardband size
% bG = floor((N_pt*0.1-1)/2); % Back guardband size

us = N_pt-floor(N_pt*0.1);  % Usable sub-carriers to assure a guardband that's 10% of BW
                            % from [LTE in a Nutshell]
symNum = length(iMat(1,:)); % Number of OFDM symbols

temp = zeros(us,symNum);
% temp(1:us/2,:) = iMat(fG+1:fG+us/2,:);
% temp(us/2+1:end,:) = iMat(fG+us/2+2:end-bG,:);

temp(1:us/2,:) = iMat(1:us/2,:);
temp(us/2+1:end,:) = iMat(us/2+1+floor(N_pt*0.1):end,:);

c_ni = reshape(temp,1,us*symNum);

end

function c_ni = dftMIMO(RX_M,N_pt)
% Computes the N_pt FFT of the input MxN ODFM data 

ySize = size(RX_M);
c_ni = zeros(ySize(1),ySize(2));
us = N_pt-floor(N_pt*0.1);  % Usable sub-carriers to assure a guardband that's 10% of BW
symNum = ceil(ySize(2)/us); % Number of OFDM symbols
for i = 1:ySize(1)
    %% Serial to Parallel for each row of RX_M
    sk = reshape(RX_M(i,:),N_pt,symNum);
    iMat = fft(sk);
    % fG = ceil((N_pt*0.1-1)/2);  % Front guardband size
    % bG = floor((N_pt*0.1-1)/2); % Back guardband size
    
    temp = zeros(us,symNum);
    % temp(1:us/2,:) = iMat(fG+1:fG+us/2,:);
    % temp(us/2+1:end,:) = iMat(fG+us/2+2:end-bG,:);
    
    temp(1:us/2,:) = iMat(1:us/2,:);
    temp(us/2+1:end,:) = iMat(us/2+1+floor(N_pt*0.1):end,:);
    
    c_ni(i,:) = reshape(temp,1,us*symNum);
end
end
