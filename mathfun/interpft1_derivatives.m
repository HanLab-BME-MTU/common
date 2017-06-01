function [ vq ] = interpft1_derivatives( v, xq, derivs , freq)
%interpft1_derivatives Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 4)
    freq = false;
end

    K = floor(size(v,1)/2);
    freqM = ifftshift(-K:K).'*1i;
    if(~freq)
        v_hat = fft(v);
    else
        v_hat = v;
    end
    derivDim = ndims(v)+1;
    freqMs = arrayfun(@(x) freqM.^x,derivs,'UniformOutput',false);
    freqMs = cat(derivDim,freqMs{:});
    
    v_hat = bsxfun(@times,v_hat,freqMs);
    
    xqrep = ones(1,derivDim);
    xqrep(derivDim) = length(derivs);
    xq = repmat(xq,xqrep);
    
    vq = interpft1([0 2*pi],v_hat, xq, 'horner_freq');
    

end

