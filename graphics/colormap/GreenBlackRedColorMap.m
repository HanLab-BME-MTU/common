function x = GreenBlackRedColorMap(n, r)
if ~exist('n', 'var')
    n = 64;
    r = 1;
end
x = zeros(n, 3);

x((129-n):128, 1) =(((1:n)/n).^r)' ;
x(1:n, 2) = (((n:-1:1)/n).^r)';
colormap(x)


% if ~exist('n', 'var')
%     n = 32;
%     r = 1;
% end
% x = zeros(64, 3);
% 
% x((65-n):64, 1) =(((1:n)/n).^r)' ;
% x(1:n, 2) = (((n:-1:1)/n).^r)';
% colormap(x)