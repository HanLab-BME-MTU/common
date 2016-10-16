function matrix = cell2matInconsistent(cellArray)
%function features = cell2matInconsistent(matrix) makes matrix out of cell
%array according to the largest cell array component and fill empty spaces
%with NaNs. cellArray should be in 1xm format.
% Sangyoon Han 2016 Oct.

numRowsCellArray = size(cellArray,1);
% Get the max nSamples
nSampleRowsAll = max(sum(cellfun(@length,cellArray)));
nSampleRowsMax = max(max(cellfun(@length,cellArray)));
numCols=size(cellArray,2);
  
matrix = NaN(nSampleRowsAll,numCols);
for ii=1:numRowsCellArray
    workingCellArray = cellArray(ii,:);
    curMatrix = NaN(nSampleRowsMax,numCols);
    % Get the corresponding vector
    for jj=1:numCols
        curVec = workingCellArray{1,jj};
        numRowsCurVec=length(curVec);
        curMatrix(1:numRowsCurVec,jj)=curVec;
    end
    % Assign the values
    matrix((ii-1)*nSampleRowsMax+1:(ii)*nSampleRowsMax,:)=curMatrix;
end

end

 
