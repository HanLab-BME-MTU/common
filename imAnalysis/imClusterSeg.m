function [seg_img, obj_val, varargout]=imClusterSeg(img_in, CONTR, varargin)
% IMKMEANSSEG image segmentation with clustering (kmeans or EM)
% 
%             This function does a segmentation of the image using the 
%             kmeans or the EM algorithm. For the kmeans the number of
%             cluster (objects) has to be specified (k_cluster), for the 
%             EM algorithm the maximumal and minimal number of clusters
%             has to be specified (k_min, k_max) the function finds then
%             the optimal number of clusters. 
%             The EM algorithm can be slow, if so use image binning (2 to
%             4), and if running the function on image stacks use previous
%             results (pp, mu) as innitial results for the next image!
%             Watch out, if initial (pp, mu) are given their dimension must
%             match k_max!
%
%             If there is an error seg_img = -99 is returned
%  
%             K-means can be used with different definitions of the
%             distance: 
%             'cityblock'
%             'sqEuclidean'
%             'cosine'
%
% SYNOPSIS    [seg_img, obj_val]=imClusterSeg(img_in, CONTR, varargin)
%
% INPUT       img       : the image
%             DEPTH     : depth of the image
%             CONTR     : flag for control image display
% 
% OUTPUT      seg_img    : segmented image, intensities equal cluster #
%             obj_val    : mean values of the segmented objects
%
%                           
% DEPENDENCES   imFindThresh uses { kmeans,
%                                   mixtures4
%                                 } 
%               imFindThresh is used by {
%                                 } 
%
% Matthias Machacek 02/10/04

%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if DEPTH == 8 | DEPTH == 10 | DEPTH == 12 | DEPTH == 14 | DEPTH == 16 | DEPTH == 255 ...
%       | DEPTH == 1023 | DEPTH == 4095 | DEPTH == 16383 | DEPTH == 65535 
% else 
%    error('unsuported image depth, only images of 8,10,12,14,16 bit are accepted.');
% end
% 
% %convert
% switch DEPTH
%    case 8
%       DEPTH=2^8-1;
%    case 10
%       DEPTH=2^10-1;
%    case 12
%       DEPTH=2^12-1;
%     case 14
%       DEPTH=2^14-1;  
%     case 16
%       DEPTH=2^16-1;
% end


l=length(varargin);
for i=1:l
    if strcmp(varargin{i},'method')
        METHOD=varargin{i+1};
    elseif strcmp(varargin{i},'k_cluster')
        K_CLUSTER=varargin{i+1};
    elseif strcmp(varargin(i),'k_min')
        K_MIN=varargin{i+1};
    elseif strcmp(varargin(i),'k_max')
        K_MAX=varargin{i+1};  
    elseif strcmp(varargin(i),'p0')
        P0=varargin{i+1};  
    elseif strcmp(varargin(i),'mu0')
        MU0=varargin{i+1};          
    elseif strcmp(varargin(i),'distance')
        DISTANCE=varargin{i+1};  
    elseif strcmp(varargin(i),'binning')
        BINNING=varargin{i+1};    
    end
end

%%%%%%%%%  general parameter %%%%%%%%%%%%%%
if ~exist('METHOD','var')
    METHOD = 'em';
end
if ~exist('BINING','var')
    BINNING = 0;
end

%%%%%%%%%%%  kmeans parameter %%%%%%%%%%%%%
if ~exist('K_CLUSTER','var')
   K_CLUSTER=2;
end
if ~exist('DISTANCE','var')
   DISTANCE='cityblock';
   %'cityblock'
   %'sqEuclidean'
   %'cosine'
end

%%%%%%%%%%   EM parameter %%%%%%%%%%%%%%%%
if ~exist('K_MIN','var')
   K_MIN=2;
end
if ~exist('K_MAX','var')
   K_MAX=4;
end
%the probability of the individual normal distributions
if ~exist('P0','var')
   P0=[];
end
%the mean values of the individual normal distributions
if ~exist('MU0','var')
   MU0=[];
end
%%%%%%%%%%%%%%%%%%% End parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_img_org m_img_org] = size(img_in);

%do image binning
if BINNING > 0
    img_bin = img_in(1:BINNING:n_img_org, 1:BINNING:m_img_org);
else 
   img_bin =  img_in;
end

%size of the image 
[n_img m_img] = size(img_bin);


%calculate the histogram
%[img_hist, hist_val] = imhist(img_in,DEPTH);


%reshape image into vector
img_in_vec = reshape(img_bin, n_img*m_img, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  clustering with kmeans  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(METHOD,'kmeans')
    [cluster_index, cluster_c, cluster_d_sum]  = kmeans(img_in_vec, K_CLUSTER, 'Distance', DISTANCE);
    
    obj_val = cluster_c;
    seg_img = reshape(cluster_index, n_img, m_img);
    
    %resize image to original size
    if BINNING > 0
        seg_img = imresize(seg_img, [n_img_org, m_img_org], 'bicubic')   
    end
    
elseif strcmp(METHOD,'em')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  clustering with EM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %check first of the initial values have the right dimension
    if ~isempty(P0)
        if size(P0) ~= size(MU0)
            seg_img = -99;
            obj_val = -99;
            return 
        end
    end
    if ~isempty(MU0)
        if size(MU0) ~= size(MU0)
            seg_img = -99;
            obj_val = -99;
            return 
        end
    end    
    
    if CONTR
        verbose = 1;
    else
        verbose = 0;         
    end
    
    %accuracy
    th = 10^(-4);
    regularize = 0;
    covoption = 0;
    [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(img_in_vec', K_MIN, K_MAX,...
                                                         regularize, th, covoption, P0, MU0, verbose);
   
    
                                                     
    %reshape image into vector
    img_in_vec = reshape(img_in, n_img_org*m_img_org, 1);                                                  
    %calculate the probability of each data point
    for k = 1:bestk
        var = bestcov(:,:,k);
        m = bestmu(:,k);
        ff = ((2*pi*(var+realmin))^(-1/2));
        y(:,k) = bestpp(k) * ff * exp((-1/(2*var))*(img_in_vec-m).^2);   
    end
    
    %find the highest probability for each data point
    [val, index] = max(y,[],2);     
    
    %put vector into matrix shape
    seg_img = reshape(index, n_img, m_img);

    
    varargout(:,1) = bestpp;
    varargout(:,2) = bestmu;
end

if CONTR
    %         fig_hist_h = figure;
    %         plot(hist_val, img_hist);
    %         hold on
    %         plot(cluster_c(1), 100, '.r');
    %         plot(cluster_c(2), 100, '.r');
    figure
    imshow(seg_img,[]);
end
