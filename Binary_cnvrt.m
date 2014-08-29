function bStrm = Binary_cnvrt(txtStr)
%%
% Eric Kang
% March 29, 2014
%
% Convert the ascii characters to binary values
% data will contain the binary information for cAr
%
% txtStr = String of ASCII characters

txtStr = double(txtStr);    % Converts to decimal
bStrm = zeros(1,length(txtStr)*8);
for i = 0:length(txtStr)-1
    bStrm(i*8+1:i*8+8) = ascii2bin(txtStr(i+1));
end

end
