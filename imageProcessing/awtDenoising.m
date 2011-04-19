function Irec = awtDenoising(I, varargin)
% Irec = AWTDENOISING(I) reconstructs the original image I from a
% soft thresholding of its A Trou wavelet coefficients.
%
% A description of the algorithm can be found in:
% "Olivo-Marin J.C. 2002. Extraction of spots in biological images using
% multiscale products. Pattern Recognit. 35: 1989�1996."
%
% Irec = AWTDENOISING(I, nBands) uses up to nBands (inclusive) of the
% A Trou Wavelet transform to reconstruct image I. The default value is
% ceil(max(log2(N), log2(M))), where [N, M] = size(I).
%
% Irec = AWTDENOISING(..., includeLoBand) allows to add the
% approximation A_K (lowest band) to the reconstructed image. The default
% value is 1 (true).
%
% Irec = AWTDENOISING(..., nSigma) allows to specify the number of
% standard deviations being used in the soft threshold. The default value
% is 3.
%
% Output:
% Irec is the reconstructed image.
%
% Sylvain Berlemont, 2009

[N, M] = size(I);

K = ceil(max(log2(N), log2(M)));

nBands = K;

if nargin > 1 && ~isempty(varargin{1})
    nBands = varargin{1};
    
    if nBands < 1 || nBands > K
        error('invalid range for nBands parameter.');
    end
end

includeLoBand = 1;

if nargin > 2 && ~isempty(varargin{2})
    includeLoBand = varargin{2};
end

nSigma = 3;

if nargin > 3 && ~isempty(varargin{3})
    nSigma = varargin{3};
end

W = awt(I, nBands);

Irec = zeros(size(I));

if includeLoBand
    Irec = W(:, :, nBands + 1);
end

madFactor = nSigma / norminv(0.75, 0, 1);
for k = 1:nBands
    S = W(:, :, k);
    S(abs(S) < madFactor * mad(S(:), 1)) = 0;
    Irec = Irec + S;
end
