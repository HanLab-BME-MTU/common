function [maxima,minima,maxima_value,minima_value,other,other_value] = interpft_extrema(x,dim,sorted)
% interpft_extrema finds the extrema of the function when interpolated by fourier transform
% This is the equivalent of doing sinc interpolation.
%
% INPUT
% x - regularly sampled values to be interpolated.
%     Values are considered to be sampled at (0:length(x)-1)*2*pi/length(x)
% dim - dimension along which to find maxima
% sorted - logical value about whether to sort the maxima by value
%
% OUTPUT
% maxima - angle between 0 and 2*pi indicating location of local maxima
% minima - angle between 0 and 2*pi indicating location of local minima
% maxima_value - value of interpolated function at maxima
% minima_value - value of interpolated function at minima
% other - angle between 0 and 2*pi where second derivative is zero
% other_value - value of interpolated function at other
%
% If no output is requested, then the maxima and minima will be plotted.
% Maxima are indicated by a red vertical line and red circle.
% Minima are indicated by a green vertical line and a green circle
%
% 1D Example
%
% r = rand(7,1);
% figure;
% plot((0:length(r)-1)/length(r)*2*pi,r,'ko');
% hold on;
% plot((0:199)/200*2*pi,interpft(r,200),'k');
% interpft_extrema(r);
% hold off;
%
% 2D Example
% r = rand(11,3);
% figure;
% plot((0:size(r,1)-1)/size(r,1)*2*pi,r,'ko');
% hold on;
% plot((0:199)/200*2*pi,interpft(r,200),'k');
% interpft_extrema(r);
% hold off;
% 
% This function works by calculating a fourier series of degree length(x)
% that fits through the input points. Then the fourier series is then considered
% as a trigonometric polynomial and the roots are solved by evaluating the eigenvalues
% of a companion matrix via the builtin function roots.
%
% See also interpft, roots
%
% Author: Mark Kittisopikul, May 2016

    % Tolerance for log(abs(root)) to be near zero, in which case the root is real
    TOL = eps*1e2;
    
%     original_size = size(x);
    if(nargin > 1)
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
    
    if(nargin < 2)
        sorted = false;
    end

    output_size = size(x);
    output_size(1) = output_size(1) - 1;

    s = size(x);
    scale_factor = s(1);

    % Calculate fft and nyquist frequency
    x_h = fft(x);
    nyquist = ceil((s(1)+1)/2);

    % If there is an even number of fourier coefficients, split the nyquist frequency
    if(~rem(s(1),2))
        % even number of coefficients
        % split nyquist frequency
        x_h(nyquist,:) = x_h(nyquist,:)/2;
        x_h = x_h([1:nyquist nyquist nyquist+1:end],:);
        output_size(1) = output_size(1) + 1;
    end
    % Wave number, unnormalized by number of points
    freq = [0:nyquist-1 -nyquist+1:1:-1]';

    % calculate derivatives, unscaled by 2*pi since we are only concerned with the signs of the derivatives
    dx_h = bsxfun(@times,x_h,freq * 1i);
    dx2_h = bsxfun(@times,x_h,-freq.^2);
    
    % use companion matrix approach
    if(dx_h(end) ~= 0)
%         r = roots(fftshift(-dx_h./dx_h(end)));
        dx_h = -fftshift(dx_h,1);
        r = zeros(output_size,'like',dx_h);
        % FIXME
        r = r(:,:).';
        dx_h = dx_h(:,:).';
        if(~isempty(gcp('nocreate')))
            parfor i=1:size(r,1)
                r(i,:) = roots(dx_h(i,:));
            end
        else
            for i=1:size(r,1);
                r(i,:) = roots(dx_h(i,:));
            end
        end
%         dx_h = dx_h.';
        r = r.';
        r = reshape(r,output_size);
%       dx_h = reshape(dx_h.',output_size);
        % keep only the real answers
%         r = r(abs(log(abs(r))) < TOL);
        imaginary_map = abs(log(abs(r))) > TOL;
%         real_map = ~imaginary_map;
        r(imaginary_map) = NaN;
    else
        % no roots
        maxima = [];
        minima = [];
        maxima_value = [];
        minima_value = [];
        other = [];
        other_value = [];
        return;
    end
   
    % In the call to roots the coefficients were entered in reverse order (negative to positive)
    % rather than positive to negative. Therefore, take the negative of the angle..
    % angle will return angle between -pi and pi
    extrema = -angle(r);
    
    % Map angles to between 0 and 2 pi, moving negative values up
    % a period
    neg_extrema = extrema < 0;
    extrema(neg_extrema) = extrema(neg_extrema) + 2*pi;
    
    % calculate angles multiplied by wave number
    theta = bsxfun(@times,extrema,shiftdim(freq,-ndims(extrema)));
    % waves
    waves = exp(1i*theta);
    % evaluate each wave by fourier coeffient
%     dx2 = waves(:,:)*dx2_h;
    ndims_waves = ndims(waves);
    dim_permute = [ndims_waves 2:ndims_waves-1 1];
    dx2 = real(sum(bsxfun(@times,waves,permute(dx2_h,dim_permute)),ndims_waves));
    
    % dx2 could be equal to 0, meaning inconclusive type
    % local maxima is when dx == 0, dx2 < 0
%     maxima = extrema(dx2 < 0);
    output_template = NaN(output_size);
    
    maxima = output_template;
    minima = output_template;
    
    maxima_map = dx2 < 0;
    minima_map = dx2 > 0;
    
    maxima(maxima_map) = extrema(maxima_map);
    % local minima is when dx == 0, dx2 > 0
%     minima = extrema(dx2 > 0);
    minima(minima_map) = extrema(minima_map);
    
    if(nargout > 2 || nargout == 0 || sorted)
	% calculate the value of the extrema if needed
%         maxima_value = waves*x_h/scale_factor;
        extrema_value = real(sum(bsxfun(@times,waves,permute(x_h,dim_permute)),ndims_waves))/scale_factor;
        maxima_value = output_template;
        minima_value = output_template;
        minima_value(minima_map) = extrema_value(minima_map);
        maxima_value(maxima_map) = extrema_value(maxima_map);
        
        if(sorted)

            
            maxima_value_inf = maxima_value;
            maxima_value_inf(isnan(maxima_value_inf)) = -Inf;
            [~,maxima,maxima_value] = sortMatrices(maxima_value_inf,maxima,maxima_value,'descend');
            clear maxima_value_inf;
            numMax = sum(maxima_map);
            numMax = max(numMax(:));
            maxima = maxima(1:numMax,:,:);
            maxima_value = maxima_value(1:numMax,:,:);
            
            minima_value_inf = minima_value;
            minima_value_inf(isnan(minima_value_inf)) = Inf;
            [~,minima,minima_value] = sortMatrices(minima_value_inf,minima,minima_value);
            clear minima_value_inf;
            numMin = sum(minima_map);
            numMin = max(numMin(:));
            minima = minima(1:numMin,:,:);
            minima_value = minima_value(1:numMin,:,:);
                       
        end
        
%         maxima_value = reshape(maxima_value,output_size);
%         minima_value = reshape(minima_value,output_size);
        maxima_value = shiftdim(maxima_value,unshift);
        minima_value = shiftdim(minima_value,unshift);
        
        if(nargout > 4)
            % calculate roots which are not maxima or minima
            other = output_template;
            other_map = dx2 == 0;
            other(other_map) = extrema(other_map);
            
            other_value = output_template;
            other_value(other_map) = extrema_value(other_map);
            
            if(sorted)
                other_value_inf = other_value;
                other_value_inf(isnan(maxima_value_inf)) = -Inf;
                [~,other,other_value] = sortMatrices(other_value_inf,other,other_value,'descend');
                clear other_value_inf;
            end
            
%             other = reshape(other,output_size);
%             other_value = reshape(other_value,output_size);
            other = shiftdim(other,unshift);
            other_value = shiftdim(other_value,unshift);
        end
    end
    

    
    if(nargout == 0)
        % plot if no outputs requested
        if(~isempty(extrema))
	    % Maxima will be green
            real_maxima = maxima(maxima_map);
            plot([real_maxima real_maxima]',ylim,'g');
            plot(real_maxima,real(maxima_value(maxima_map)),'go');
	    % Minima will be red
            real_minima = minima(minima_map);
            plot([real_minima real_minima]',ylim,'r');
            plot(real_minima,real(minima_value(minima_map)),'ro');
        else
            warning('No extrema');
        end
    end
    
%     maxima = reshape(maxima,output_size);
%     minima = reshape(minima,output_size);
    maxima = shiftdim(maxima,unshift);
    minima = shiftdim(minima,unshift);
    
end
