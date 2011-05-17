function closedContours = closeContours(contoursIn,matIn)
%CLOSECONTOURS takes the input  open contours and closes them by filling in the gaps where they meet the image border
%
% closedContours = closeContours(contoursIn,matIn)
%
%
% Description:
% 
% The contours as returned by the matlab contouring functions (contour.m,
% contourc.m etc) will contain gaps where a given contour meets the edge of
% the matrix is being contoured. This function closes these contours by
% filling in the appropriate areas at the image boundary, so that the
% closed contour encloses areas of the matrix where the value is higher
% than the isovalue for that contour. The contours should be first
% separated using separateContours.m Any input contours which are already
% closed will be un-affected.
% 
% Input: 
% 
%   contoursIn - A 1xM cell array of the M input contours, as returned by
%                separateContours.m
% 
%   matIn -    The matrix the contours were derived from. 
% 
% 
% Output:
% 
%   closedContours - A 1xM array of the closed contours. 
% 
% 
% Hunter Elliott
% 4/2010
% 

%% -------- Parameters ---------- %%

nChk = 5;


%% ----- Input ----- %%

if nargin < 2 || isempty(matIn) || isempty(contoursIn)
    error('Must input both a cell-array of contours and the matrix which the contours came from!')
end

if ~iscell(contoursIn)
    error('The input contours must be separated into a cell array! Try using separateContours.m!')
end


%% -----  Init ----- %%

[M,N] = size(matIn);

nContours = length(contoursIn);

%Get the coordinates of the image border in a  clock-wise fashion, starting
%at origin.
borderCoord = vertcat([1:(N-1)      ones(1,M-1)*N  N:-1:2         ones(1,M-1)],...
                      [ones(1,N-1)  1:(M-1)        ones(1,N-1)*M  M:-1:2]);
    
nB = length(borderCoord);                  
                  

closedContours = cell(nContours,1);

%% ----- Closure ---- %%

for j = 1:nContours
    
    %Verify that this contour needs closure by checking that the first
    %point touches the image border.
    iTouchStart = find(arrayfun(@(x)(all(borderCoord(:,x) == round(contoursIn{j}(:,1)))),1:nB),1,'first');    
    
    if ~isempty(iTouchStart)
        %Find where the last point touches the border
        iTouchEnd = find(arrayfun(@(x)(all(borderCoord(:,x) == round(contoursIn{j}(:,end)))),1:nB),1,'last');
        
        %Check which direction the values increase in
        iBefore = mod(iTouchEnd-nChk:iTouchEnd-1,nB);%Use modulus in case it starts near the origin
        iAfter = mod(iTouchEnd+1:iTouchEnd+nChk,nB);
        iBefore(iBefore==0) = nB;
        iAfter(iAfter==0) = nB;
        iBefore = sub2ind(size(matIn),borderCoord(2,iBefore),borderCoord(1,iBefore));%Convert these to matrix indices
        iAfter = sub2ind(size(matIn),borderCoord(2,iAfter),borderCoord(1,iAfter));
        
        if max(matIn(iBefore)) < max(matIn(iAfter))
            incClockwise = true;
        else
            incClockwise = false;                                    
        end
        
        if incClockwise %If they increase clockwise...
            
            %Check if the border used for closing crosses the origin
            if iTouchEnd < iTouchStart

                %Add these border points to the contour in normal order
                closedContours{j} = [contoursIn{j} borderCoord(:,iTouchEnd+1:iTouchStart-1)];

            elseif iTouchEnd > iTouchStart
                %Add the border points in normal order, looping past the
                %origin
                closedContours{j} = [contoursIn{j} borderCoord(:,iTouchEnd+1:end) borderCoord(:,1:iTouchStart-1)];
                
            end%If they are equal, we don't want to do anything - this curve just touches the border at one point. 
            
        else
            %Check if the border used crosses the origin
            if iTouchEnd > iTouchStart

                %Add these border points to the contour in reverse order
                closedContours{j} = [contoursIn{j} borderCoord(:,iTouchEnd-1:-1:iTouchStart+1)];

            elseif iTouchEnd < iTouchStart
                %Add the border points in reverse order, looping past the
                %origin
                closedContours{j} = [contoursIn{j} borderCoord(:,iTouchEnd-1:-1:1) borderCoord(:,end:-1:iTouchStart+1)];
                
            end%If they are equal, we don't want to do anything - this curve just touches the border at one point. 
            
        end
    else
       closedContours{j} = contoursIn{j}; 
    end
end

