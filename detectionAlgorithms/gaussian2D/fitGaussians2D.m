
% Inputs:   frame : ny * nx * nc array, where ny and nx are the height and width,
%                   respectively, and nc is the number of channels
%            mask :      
%

% Output: pmat:
%               rows: x, y, A, s, c, x_pstd, y_pstd, A_pstd, 


% Usage for a single-channel frame with mask and fixed sigma:
% fitGaussians2D(frame, x_v, y_v, 'sigma', sigma_v, 'mask', mask);


% Francois Aguet, March 28 2011

function output = fitGaussians2D(frame, x, y, A, sigma, c, mode, varargin)

np = length(x);

output = struct('x', [], 'y', [], 'A', [], 's', [], 'c', [],...
                'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 's_pstd', [], 'c_pstd', [],...
                'c_maskStd', [], 'c_resStd', [], 'pval_KS', [], 'pval_A', []);

xi = round(x);
yi = round(y);
[ny,nx] = size(frame);

if mod(length(varargin),2)~=0
    error('Optional arguments need to be entered as pairs.');
end

idx = find(strcmpi(varargin, 'mask'));
if ~isempty(idx)
    mask = varargin{idx+1};
    labels = bwlabel(mask);
else
    %mask = [];
    labels = zeros(size(frame));
end

iRange = [squeeze(min(min(frame))) squeeze(max(max(frame)))];

estIdx = regexpi('xyAsc', ['[' mode ']']);
fixIdx = setdiff(1:5, estIdx);

             
% initialize output arrays
output.x = NaN(1,np);
output.y = NaN(1,np);
output.A = NaN(1,np);


output.x_pstd = NaN(1,np);
output.y_pstd = NaN(1,np);
output.A_pstd = NaN(1,np);
% if fitSigma
    output.s = NaN(1,np);
    output.s_pstd = NaN(1,np);
% end
output.c = NaN(1,np);
output.c_pstd = NaN(1,np);
output.c_maskStd = NaN(1,np);
output.c_resStd = NaN(1,np);

output.pval_A = NaN(1,np);
output.pval_KS = NaN(1,np);

sigma_max = max(sigma);
w2 = ceil(2*sigma_max);
w3 = ceil(3*sigma_max);
w4 = ceil(4*sigma_max);
dw = w4-w3;

% mask template: ring with inner radius w3, outer radius w4
[x,y] = meshgrid(-w4:w4);
r = sqrt(x.^2+y.^2);
bandMask = zeros(size(r));
bandMask(r<=w4 & r>=w3) = 1;

% mask template: disk with radius w3
[x,y] = meshgrid(-w3:w3);
r = sqrt(x.^2+y.^2);
diskMask = zeros(size(r));
diskMask(r<w3) = 1;

% figure; imagesc(diskMask); axis image;

% indexes for cross-correlation coefficients
n = length(mode);
idx = pcombs(1:n);
i = idx(:,1);
j = idx(:,2);
ij = i+n*(j-1);
ii = i+n*(i-1);
jj = j+n*(j-1);


for p = 1:np
    
    if (xi(p)>w4 && xi(p)<=nx-w4 && yi(p)>w4 && yi(p)<=ny-w4)
          window = frame(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        
        % label mask
        maskWindow = labels(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
        
        % estimate background
        cmask = bandMask;
        cmask(maskWindow~=0) = 0;
        window = frame(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        if isempty(c)
            c_init = mean(window(cmask==1));
        else
            c_init = c(p);
        end
        %cs = std(window(cmask==1));
        % set any other components to NaN
        window(maskWindow~=0) = NaN;
        
        % standard deviation of the background within mask
        bgStd = nanstd(window(bandMask==1));
        
        % reduce to w = 3*sigma from w = 4*sigma
        window = window(dw+1:end-dw, dw+1:end-dw);
        
        % mask w/ 3*sigma disk
        window(diskMask==0) = NaN;
        npx = sum(isfinite(window(:)));
        
        % fit
        if isempty(A)
            A_init = max(window(:))-c_init;
        else
            A_init = A(p);
        end
        [prm, prmStd, C, res] = fitGaussian2D(window, [0 0 A_init sigma(p) c_init], mode);
        
        K = corrFromC(C,ij,ii,jj);

        dx = prm(1);
        dy = prm(2);
        
        % eliminate points where localization failed or which are close to image border
        if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2 && prm(3)<2*diff(iRange))
            
            output.x(p) = xi(p) + dx;
            output.y(p) = yi(p) + dy;
            output.A(p) = prm(3);
            output.s(p) = prm(4);
            output.c(p) = prm(5);
            
            stdVect = zeros(1,5);
            stdVect(estIdx) = prmStd;
            
            output.x_pstd = stdVect(1);
            output.y_pstd = stdVect(2);
            output.A_pstd = stdVect(3);
            output.s_pstd = stdVect(4);
            output.c_pstd = stdVect(5);
            
            output.c_maskStd = bgStd;
            output.c_resStd = res.std;
            output.pval_KS(p) = res.pval;
            output.pval_A = NaN; % to be added
            
        end
    end
end



function K = corrFromC(C,ij,ii,jj)
n = size(C,1);
K = zeros(n,n);

K(ij) = C(ij) ./ sqrt(C(ii).*C(jj));
% remaining components are redundant
% K = K + K';
% K(sub2ind([n n], 1:n, 1:n)) = 1;
