function [binvec] = dec2binarray(dec, numBits)

    % Convert the decimal number to a binary string.
    switch nargin
        case 1
           binstr = dec2bin(dec);
        case 2
           binstr = dec2bin(dec,numBits);
    end
    
    binvec = logical(str2num([fliplr(binstr);blanks(length(binstr))]')');
    
end