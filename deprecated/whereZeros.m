function y = whereZeros(x)

fprintf(2, ['Warning: ''' mfilename ''' is a horrible function and no code like this should exist.\n']);

% derivative of non-linear line filtering response

y = -(x+1/2)*exp(-1/2*(x+1/2)^2/(x^2)) + ...
   +(x-1/2)*exp(-1/2*(x-1/2)^2/(x^2));