function MIMO = MIMO_lib
% Eric Kang
% April 26, 2014
% Contains all functions regarding MIMO
MIMO.MRC = @MRC;

end

function z_eq = MRC(x,M,snr_db)
% Maximum Ratio Combining
% Simulates channel propagation and MRC combiner
% x:    Transmitted signal
% M:    Number of RX antenna

Pt = mean(x.*conj(x));
%% Channel
snr = 10^(snr_db/10);
h = sqrt(1/2).*(randn(M,1)+1i*randn(M,1));  % Rayleigh Channel Coeff
%h = ones(M,1); % AWGN Channel
Pn = Pt/snr;
noise = sqrt(Pn/2).*(randn(M,length(x))+1i.*randn(M,length(x)));

%% Rx
y = h*x+noise;
a = conj(h).';    % Combining Coeff
z = a*y;    % Combiner
z_eq = 1/sum(a.'.*conj(a.')).*z; % Equalization
end
