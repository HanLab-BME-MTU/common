function [norm_values,interval_int,x_roots] = normft1(x,dim,TOL,dofft)
% normft calculates the 1-norm of a fourier series with regular samples
%
% INPUT
% x - regularly sampled values to be interpolated.
%     Values are considered to be sampled at (0:length(x)-1)*2*pi/length(x)
% dim - dimension along which to find maxima
% TOL - tolerance to determine if log(abs(root)) is zero to determine if
%       the root of the derivative is real.
%       If tolerance is negative, then use 10*abs(log(abs(root))) as the
%       tolerance only if no roots are found with tolerance at -TOL.
% dofft - logical. If true, transforms dim with fft. Default: true
%
% OUTPUT
% norm_values - value of the 1-norm
% interval_int - value of the integral between roots, NaN if no roots
% x_roots - roots between 0 and 2*pi, NaN if no roots
%
% See also norm, chebfun.norm
%
% Author: Mark Kittisopikul, April 2017
    
%     original_size = size(x);
    if(nargin > 1 && ~isempty(dim))
        x = shiftdim(x,dim-1);
        unshift = ndims(x) - dim + 1;
    else
        if(isrow(x))
            % If the input is a row vector, transpose it without conjugation
            dim = 2;
            unshift = 1;
            x = x.';
        else
            dim = 1;
            unshift = 0;
        end
    end
    
    if(nargin < 3 || isempty(TOL))
    	% Tolerance for log(abs(root)) to be near zero, in which case the root is real
        % Set negative so that tolerance adapts if no roots are found
        TOL = -eps(class(x))*1e2;
    end
    if(nargin < 4 || isempty(dofft))
        dofft = true;
    end



    output_size = size(x);
%     output_size(1) = output_size(1) - 1;
    

    s = size(x);
    scale_factor = s(1);
    
    if(s(1) == 1)
        norm_values = shiftdim(abs(x),unshift)*2*pi;
        interval_int = norm_values;
        x_roots = NaN(s);
        return;
    end

    % Calculate fft and nyquist frequency
    if(dofft)
        x_h = fft(x);
    else
        x_h = x;
    end
    nyquist = ceil((s(1)+1)/2);

    % If there is an even number of fourier coefficients, split the nyquist frequency
    if(~rem(s(1),2))
        % even number of coefficients
        % split nyquist frequency
        x_h(nyquist,:) = x_h(nyquist,:)/2;
        x_h = x_h([1:nyquist nyquist nyquist+1:end],:);
        x_h = reshape(x_h,[s(1)+1 s(2:end)]);
        output_size(1) = output_size(1) + 1;
    end
    % Wave number, unnormalized by number of points
    freq = [0:nyquist-1 -nyquist+1:1:-1]';

    % calculate integral, not including 2*pi factor here
    qx_h = bsxfun(@rdivide,x_h,freq * 1i);
    qx_h(1,:,:) = 0;
    
    % use companion matrix approach
    x_mean = x_h(1,:,:)./s(1);
    x_h = -fftshift(x_h,1);
    x_h = x_h(:,:);
%     r = zeros(output_size,'like',dx_h);
    output_size1 = output_size(1);
    nProblems = prod(output_size(2:end)); 
    batchSize = min(1024,nProblems);
    nBatches = ceil(nProblems/batchSize);
    % Only use parallel workers if a pool already exists
    nWorkers = ~isempty(gcp('nocreate'))*nBatches;
    in = ones(1,nBatches)*batchSize;
    in(end) = in(end) + nProblems - sum(in);
    x_h = mat2cell(x_h,size(x_h,1),in);
    r = cell(1,nBatches);
    parfor (i=1:nBatches, nWorkers)
        r{i} = zeros(output_size1,size(x_h{i},2),'like',x_h{i});
        for j = 1:size(x_h{i},2);
            try
                % roots outputs only column vectors which may be shorter than
                % expected
                x_h_roots = roots(x_h{i}(:,j));
                x_h_roots(end+1:output_size1) = 0;
                r{i}(:,j) = x_h_roots;
            catch err
                switch(err.identifier)
                    case 'MATLAB:ROOTS:NonFiniteInput'
                        r{i}(:,j) = NaN(output_size1,1);
                end
            end
        end
    end
    r = [r{:}];
    r = reshape(r,output_size);
    % magnitude
    magnitude = abs(log(abs(r)));
    % keep only the real answers
    real_map = (magnitude <= abs(TOL));
    % If tolerance is negative and no roots are found, then use the
    % root that is closest to being real
    if(TOL < 0)
        no_roots = ~any(real_map);
        real_map(:,no_roots) = bsxfun(@le,magnitude(:,no_roots),min(magnitude(:,no_roots))*10);
    end
%     r(imaginary_map) = NaN;
%     real_map = ~imaginary_map;
%     clear imaginary_map
    
   
    % In the call to roots the coefficients were entered in reverse order (negative to positive)
    % rather than positive to negative. Therefore, take the negative of the angle..
    % angle will return angle between -pi and pi
    
    r = -angle(r(real_map));
    
    % Map angles to between 0 and 2 pi, moving negative values up
    % a period
    neg_roots = r < 0;
    r(neg_roots) = r(neg_roots) + 2*pi;
    
    x_roots = NaN(output_size);
    x_roots(real_map) = r;
    
    % Sort in order within the periodic domain
    x_roots = sort(x_roots);
    % Prepend the last root and calculate integral at roots
    x_roots_wrapped = [nanmax(x_roots); x_roots];
    
    % Use Horner's method 
    interval_int = interpft1([0 2*pi],qx_h,x_roots_wrapped,'horner_freq');
    % Take the difference of adjacent roots, wrapping around periodic
    % boundary
    x_roots_wrapped(1,:,:) = x_roots_wrapped(1,:,:)-2*pi;
    interval_int = diff(interval_int) + bsxfun(@times,diff(x_roots_wrapped),x_mean);
    
    % Norm is the sum of the absolute value of the integral over each
    % interval
    norm_values = nansum(abs(interval_int));
    norm_is_zero = norm_values == 0;
    norm_values(norm_is_zero) = abs(x_mean(norm_is_zero))*2*pi;
    
    % Realign dimensions
    norm_values = shiftdim(norm_values,unshift);
    if(nargout > 1)
        interval_int = shiftdim(interval_int,unshift);
        if(nargout > 2)
            x_roots = shiftdim(x_roots,unshift);
        end
    end   
    
end
