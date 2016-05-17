function [maxima,minima,maxima_value,minima_value,other,other_value] = interpft_extrema(x)
% interpft_extrema finds the extrema of the function when interpolated by fourier transform
% This is the equivalent of doing sinc interpolation.
%
% INPUT
% x - regularly sampled values to be interpolated.
%     Values are considered to be sampled at (0:length(x)-1)*2*pi/length(x)
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
% Example
%
% r = rand(7,1);
% figure;
% plot((0:length(r)-1)/length(r)*2*pi,r,'ko');
% hold on;
% plot((0:199)/200*2*pi,interpft(r,200),'k');
% interpft_extrema(r);
% hold off;
% 
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

    s = size(x);
    scale_factor = s(1);
    % If the input is a row vector, transpose it without conjugation
    if(s(1) == 1 && length(s) == 2)
        x = x.';
        s = size(x);
        scale_factor = s(1);
    end

    % Calculate fft and nyquist frequency
    x_h = fft(x);
    nyquist = ceil((s(1)+1)/2);

    % If there is an even number of fourier coefficients, split the nyquist frequency
    if(~rem(s(1),2))
        % even number of coefficients
        % split nyquist frequency
        x_h(nyquist,:) = x_h(nyquist,:)/2;
        x_h = x_h([1:nyquist nyquist nyquist+1:end],:);
        s(1) = s(1) + 1;
    end
    % Wave number, unnormalized by number of points
    freq = [0:nyquist-1 -nyquist+1:1:-1]';

    % calculate derivatives, unscaled by 2*pi since we are only concerned with the signs of the derivatives
    dx_h = x_h .* freq * 1i;
    dx2_h = -x_h .* freq.^2;
    
    % use companion matrix approach
    if(dx_h(end) ~= 0)
        r = roots(fftshift(-dx_h./dx_h(end)));
        % keep only the real answers
        r = r(abs(log(abs(r))) < TOL);
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
    theta = bsxfun(@times,extrema,freq');
    % waves
    waves = exp(1i*theta);
    % evaluate each wave by fourier coeffient
    dx2 = waves*(dx2_h);
    
    % dx2 could be equal to 0, meaning inconclusive type
    % local maxima is when dx == 0, dx2 < 0
    maxima = extrema(dx2 < 0);
    % local minima is when dx == 0, dx2 > 0
    minima = extrema(dx2 > 0);
    
    if(nargout > 2 || nargout == 0)
	% calculate the value of the extrema if needed
        maxima_value = waves*x_h/scale_factor;
        minima_value = maxima_value(dx2 > 0);
	if(nargout > 4)
		% calculate roots which are not maxima or minima
		other = extrema(dx2 == 0);
		other_value = maxima_value(dx2 == 0);
	end
        maxima_value = maxima_value(dx2 < 0);
    end
    
    if(nargout == 0)
        % plot if no outputs requested
        if(~isempty(extrema))
	    % Maxima will be green
            plot([maxima maxima]',ylim,'g');
            plot(maxima,real(maxima_value),'go');
	    % Minima will be red
            plot([minima minima]',ylim,'r');
            plot(minima,real(minima_value),'ro');
        else
            warning('No extrema');
        end
    end

end
