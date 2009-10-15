function W = WTATrou(I, varargin)
% W = WTATROU(I) computes the A Trou Wavelet Transform. A description of
% the algorithm can be found in:
% J.-L. Starck, F. Murtagh, A. Bijaoui, "Image Processing and Data
% Analysis: The Multiscale Approach", Cambridge Press, Cambridge, 2000. 
%
% Input:
% I: the input image of size N x M.
% W: the wavelet coefficients which consists of a array of size N x M x K+1
% where K = ceil(max(log2(N), log2(M))). The coefficients are organized as
% follows:
% W(:, :, 1:K) corresponds to the wavelet coefficients (also called detail
% images) at scale k = 1...K
% W(:, :, K+1) corresponds to the last approximation image A_K.
%
% You can use plotWTATrou(W) to display the wavelet coefficients.
%
% Sylvain Berlemont, 2009

[N, M] = size(I);

J = ceil(max(log2(N), log2(M)));

W = zeros(N, M, J + 1);

I = double(I);
lastA = I;

for k = 1:J
    newA = convolve(I, k);
    W(:, :, k) = lastA - newA;
    lastA = newA;
end

W(:, :, J + 1) = lastA;

    function F = convolve(I, k)
        [N, M] = size(I);
        k1 = 2^(k - 1);
        k2 = 2^k;
        
        tmp = padarray(I, [k2 0], 'replicate');
        
        % Convolve the columns
        for i = k2+1:k2+N
            I(i - k2, :) = 6 * tmp(i, :) + 4 * tmp(i + k1, :) + ...
                4 * tmp(i - k1, :) + tmp(i + k2, :) + tmp(i - k2, :);
        end

        tmp = padarray(I * .0625, [0, k2], 'replicate');
        
        % Convolve the rows
        for i = k2+1:k2+M
            I(:, i - k2) = 6 * tmp(:, i) + 4 * tmp(:, i + k1) + ...
                4 * tmp(:, i - k1) + tmp(:, i + k2) + tmp(:, i - k2);
        end

        F = I * .0625;