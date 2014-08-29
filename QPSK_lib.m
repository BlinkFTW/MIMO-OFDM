function QPSK = QPSK_lib
% Eric Kang
% April 18, 2014
% Contains all functions regarding QPSK
QPSK.bin2symb = @bin2symb;
QPSK.sym2bin = @sym2bin;
QPSK.Constellation = @Constellation;
QPSK.stdev = @stdev;
end

function symb = bin2symb(binDat,m)
%%
% Convert binary array to QPSK symbols
% binDat -> binary array
% m      -> Constellation mapping mode
%           m = 1 gives pi/4 constellation
%           m = 2 gives pi/2 constellation
% symb -> QPSK symbol array (complex)

% inStr = 'Hello';
% decStr = double(inStr);
% 
% binDat = Binary_cnvrt(decStr);

if mod(length(binDat),2) ~= 0
    error('Binary string length must be EVEN for QPSK')
end

% Convert to bipolar binary
for i = 1:length(binDat)
    if binDat(i) == 0
        binDat(i) = -1;
    end
end

symb = zeros(1,length(binDat)/2);

if m == 1
    % pi/4 constellation
    for j = 1:length(symb)
        bi = j*2-1; % binary array index
        symb(j) = binDat(bi) + 1i.*binDat(bi+1); 
    end
else
    % pi/2 constellation table
    for j = 1:length(symb)
        bi = j*2-1; % binary array index
        if binDat(bi) == -1 && binDat(bi+1) == -1
            symb(j) = 1;
        elseif binDat(bi) == 1 && binDat(bi+1) == 1
            symb(j) = -1;
        elseif binDat(bi) == -1 && binDat(bi+1) == 1
            symb(j) = 1i;
        else
            symb(j) = -1i;
        end
    end
end
end

function [binDat, mSym] = sym2bin(symb)
% Map demodulated baseband symbols to QPSK symbols
% Output Binary 

% Histogram of symbols show boundaries range from -pi to pi
% Decision boundaries will be [0 90 -90]

an = angle(symb);
mSym = zeros(1,length(an));
for i = 1:length(symb)
    if an(i) >= 0
        if an(i) >= pi/2  % Map to 01
            mSym(i) = -1+1i;
        else                % Map to 11
            mSym(i) = 1+1i;
        end
    else
        if an(i) <= -pi/2     % Map to 00
            mSym(i) = -1-1i;
        else                % Map to 10
            mSym(i) = 1-1i;
        end
    end
end

binDat = zeros(1,length(an)*2);
binDat(1:2:end) = real(mSym);
binDat(2:2:end) = imag(mSym);

for i = 1:length(binDat)
    if binDat(i) == -1
        binDat(i) = 0;
    end
end
end

function Constellation(c_ni,SNR_dB)
% Plots the QPSK symbol Constellation map
% c_ni:     Complex symbol stream
% SNR_dB:   SNR at receiver in dB

M = 31;
k = 0:M;
s = sqrt(2).*exp(1i*2*pi*k./M); % PSK base circle
p = [1+1i -1+1i 1-1i -1-1i];    % QPSK symbols
% Boundaries
b = ceil(max(abs(max(real(c_ni))),abs(max(imag(c_ni)))));
v = [b*1i -b*1i];   % Vertical Line
hold on
plot(c_ni,'r.') % Received Symbols
xlabel('In-phase Component')
ylabel('Quatrature Component')
title(['Received Symbols for SNR = ' num2str(SNR_dB) 'dB'])
plot(p,'ko')    % Symbol Constellations
plot(v,'k')
plot(-b:b,zeros(1,2*b+1),'k')   % Horizontal Boundary
plot(s,'--');
legend('RX symbols','QPSK Symbols','Decision Boundary')
grid on
hold off

end

function std = stdev(symb)
% Calculate the average symbol variance from QPSK symbols
an = angle(symb);
x = zeros(1,length(an));
for i = 1:length(symb)
    if an(i) >= 0
        if an(i) >= pi/2  % Map to 01
            x(i) = sqrt((-1-real(symb(i)))^2+(1-imag(symb(i)))^2);
        else                % Map to 11
            x(i) = sqrt((1-real(symb(i)))^2+(1-imag(symb(i)))^2);
        end
    else
        if an(i) <= -pi/2     % Map to 00
            x(i) = sqrt((-1-real(symb(i)))^2+(-1-imag(symb(i)))^2);
        else                % Map to 10
            x(i) = sqrt((1-real(symb(i)))^2+(-1-imag(symb(i)))^2);
        end
    end
end

std = mean(x);
% vr = std^2;
end
