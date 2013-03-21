function [diffusedVectorField] = diffuseVectorField(referenceVectorField, varargin)
% Diffuses an ND vector field
% 
%     [diffusedVectorField] = diffuseVectorField(referenceVectorField, varargin)
% 
% Computes the diffused version of an N-Dimensional vector field by solving
% a variational optimization problem described in the following paper:
% 
%   Chenyang, X. and J. L. Prince (1998). "Snakes, shapes, and gradient vector 
%   flow." IEEE Transactions on Image Processing, 7(3): 359-369.
% 
% Required Input Arguments:
% 
%       referenceVectorField: Input N-D reference vector field.
%     
%                             Must be a cell array of size equal to the 
%                             dimensionality of the vector field. Each 
%                             element of the cell array must be an ND 
%                             matrix representing the corresponding 
%                             component of the vector field.
%                           
%   
% Optional Input Arguments:
% 
%       regularizationWeight: A postive real number specifying the weight of
%                             the regularization term in the optimization 
%                             function. 
%                                 
%                             This parameter allows you to control the 
%                             amount of smoothing.
%                             
%                             Default: 0.1
%                             Prescribed Range: (0, 0.2]                          
% 
%             maxIterations: A positive integer specifying the maximum number 
%                            of iterations of the gradient descent procedure.
%                            
%                            If you want the diffusion to reach steady state
%                            set this to a very high value. If not, set it
%                            based on the amount of diffusion you want.
% 
%                tolerance:  A small positive real number determines a 
%                            tolerance level for the convergence of the            
%                            gradient descent procedure.                               
% 
%                            Default: 10^-4
% 
%                   spacing: pixel spacing of the input image.
%                            Must be an array of size equal to the number 
%                            of dimensions of the vector field.
%                            
%                            Default: 1 (isotropic spacing)                      
%                            
%             
%                 debugMode: true/false
%                            A bunch of stuff is printed in debug mode. If 
%                            you are just a user of this function and do
%                            not care about what is going on inside, then
%                            dont bother setting it to true. 
%                            
%                            Default: false
%                            
%           flagParallelize: true/false
%           
%                            If set to true, then parts of the algorithm are
%                            run in parallel. This mode will only work if
%                            you have a multi-code CPU and the parallel
%                            computing toolbox of matlab.
%                            
%                            In particular, the euler step of the different
%                            components of the vector field are computed in
%                            parallel.
%                            
%                            Default: false
%                            
% Output Arguments:                           
% 
%       diffusedVectorField: a diffused version of the input vector field. 
% 
% Example:  
% 
%     im = zeros(300,400);
%       
%     [X,Y] = meshgrid(1:size(im,2), 1:size(im,1));
%     
%     xc = round(0.5 * size(im,2));
%     yc = round(0.5 * size(im,1));
%     rad = round(0.25 * min(size(im)));
%     
%     im = ((X - xc).^2/rad^2) + ((Y - yc).^2/rad^2) - 1 <= 0;
%     im = double(bwperim(im));
% 
%     referenceVectorField = cell(1, 2);    
%     [referenceVectorField{:}] = gradient(im);
%     
%     [diffusedVectorField] = diffuseVectorField( referenceVectorField, 'debugMode', true);
% 
%
% Author: Deepak Roy Chittajallu (Created Mar 20, 2013)
%                            
    p = inputParser;
    p.addRequired( 'referenceVectorField', @iscell);
    p.parse(referenceVectorField);
    
    % validate each element of the reference vector field      
    for i = 1:numel(referenceVectorField)              
        assert( isnumeric(referenceVectorField{i}) && ...
                sum(size(referenceVectorField{i}) > 1) >= 2, ...
                'reference vector field should be a matrix of dimension >= 2' );
            
        if i == 1
            imsize = size(referenceVectorField{1});
        else
            assert( all(size(referenceVectorField{i}) == imsize), ...
                    'all components of vector field must be of same size');
        end
    end
    
    imdims = ndims(referenceVectorField{1});
    
    p.addParamValue( 'regularizationWeight', 0.1, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'maxIterations', 10000, @(x) (isnumeric(x) && isscalar(x)) );    
    p.addParamValue( 'tolerance', 10^-4, @(x) (isnumeric(x) && isscalar(x)) );    
    
    p.addParamValue('spacing', ones(1, imdims), ...
                    @(x) (isnumeric(x) && ~isscalar(x) && numel(x) == imdims) );
    p.addParamValue('debugMode', false, @(x) ( islogical(x) & isscalar(x) ) );
    p.addParamValue('flagParallelize', false, @(x) ( islogical(x) & isscalar(x) ) );
    p.parse(referenceVectorField, varargin{:} );
    
    parameters = p.Results;
    
    % get parameters
    regularizationWeight = parameters.regularizationWeight;
    maxIterations = parameters.maxIterations;
    
    spacing = num2cell(parameters.spacing);
    flagDebugMode = parameters.debugMode;
    flagParallelize = parameters.flagParallelize;
    tolerance = parameters.tolerance;
    
    % if parallelization is requested and matlabpool is not open, then open it
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end    
    
    % pad the vector field matrix so that pixels on the boundary are handled
    % properly
    padMethod = 'symmetric';
    padSize = ones(1,imdims); % pad just enough to compute the laplacian at the boundary
    referenceVectorFieldPadded = cell(size(referenceVectorField));
    for i = 1:imdims
        referenceVectorFieldPadded{i} = padarray(referenceVectorField{i}, padSize,padMethod);
    end
    unpadmask = padarray(ones(imsize), padSize, 0 );
    borderCorrectionInd = getMatrixPaddingIndices(imsize, padSize, padMethod, 'both' );
        
    for j = 1:imdims
        borderCorrectionInd{j} = borderCorrectionInd{j} + 1;
    end
    
    % compute the squared magnitude of the vector field
    vmag2 = zeros(size(referenceVectorFieldPadded{1}));
    for j = 1:imdims
        vmag2 = vmag2 + referenceVectorFieldPadded{j}.^2;
    end
    
    % Iteratively diffuse the vector field
    fprintf(1, '\n' );
    iterChange = zeros(1,imdims);   
    
    initialVectorField = referenceVectorFieldPadded;
    vectorField = initialVectorField;
    
    for i = 1:maxIterations        
        if flagParallelize
            parfor j = 1:imdims

                % boundary correction/ensure-mirroring
                vold = vectorField{j}(borderCorrectionInd{:}); 

                % compute euler step
                vectorField{j} = vold + regularizationWeight * (2 * imdims) * del2(vectorField{j}, spacing{:}) - vmag2 .* (vectorField{j} - initialVectorField{j});

                % record maximum change from previous iteration
                vchange = vectorField{j} - vold;            
                iterChange(j) = max(abs(vchange(:)));            

            end                                
        else            
            for j = 1:imdims

                % boundary correction/ensure-mirroring
                vold = vectorField{j}(borderCorrectionInd{:}); 

                % compute euler step
                vectorField{j} = vold + regularizationWeight * (2 * imdims) * del2(vectorField{j}, spacing{:}) - vmag2 .* (vectorField{j} - initialVectorField{j});

                % record maximum change from previous iteration
                vchange = vectorField{j} - vold;            
                iterChange(j) = max(abs(vchange(:)));            

            end        
        end
        
        % print iteration progress
        fprintf(1, '%3d ', i );
        if( rem(i,20) == 0 )
            fprintf(1, '\n' );
        end                        
        
        changeTrend(i) = max(iterChange);
        
        if ~any(iterChange > tolerance)
            fprintf(1, '\nsteady-state reached in %d iterations\n', i );
            break;
        end
    end    
    fprintf(1, '\n' );
    
    if flagDebugMode
        figure, plot(changeTrend);        
    end
    
    % unpad the components of the vector field
    for j = 1:imdims
        vectorField{j} = reshape( vectorField{j}(unpadmask > 0), imsize );
    end    
    
    diffusedVectorField = vectorField;
    
    % display the vector field if debugMode is on and the vector field 2D
    if flagDebugMode && imdims == 2
        
        [X,Y] = meshgrid(1:imsize(2), 1:imsize(1));
        vmag = sqrt(vmag2);
        vmag = reshape( vmag(unpadmask > 0), imsize );
        
        figure;
        
        subplot(1,2,1)
        
            imshow(vmag, []);
            hold on;
                quiver(X, Y, referenceVectorField{:});
            hold off;
            title( 'Reference Vector Field' );
            
        subplot(1,2,2)
        
            imshow(vmag, []);
            hold on;
                quiver(X, Y, diffusedVectorField{:});
            hold off;
            title( 'Diffused Vector Field' );
    end
    
    % close matlapool if it was not open before calling the function
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end
    
end