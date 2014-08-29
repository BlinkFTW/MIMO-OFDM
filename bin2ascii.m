function decStr = bin2ascii(binStr)
% Converts binary array into decimal string
% Every 8 'bits' of array are converted to decimal
% use bin2dec
decStr = zeros(1,length(binStr)/8);
% Concatinate binary to string in groups of 8

for i = 0:length(decStr)-1
    decStr(i+1) = bin2dec(num2str(binStr(i*8+1:i*8+8)));
end

end
