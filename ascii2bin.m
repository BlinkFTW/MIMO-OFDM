function y = ascii2bin(x_int)
% This function takes a single ASCII character as input and outputs its
% corresponding ASCII value in a binary array
%
% x = Single ASCII char in DEC form
% y = 8 bit binary number

% x_int = double(x);  % Integer equivalent

x_binStr = dec2bin(x_int); %Binary # in string form

y = zeros(1,length(x_binStr));
for i = 1:length(y)
    if x_binStr(i) == '1'
        y(i+1) = 1;
    else
        y(i+1) = 0;
    end
end

% Append any necessary zeros s.t. y is 8-bits
y = [zeros(1,8-length(y)) y];

end
