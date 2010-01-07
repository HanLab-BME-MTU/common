function sampleArray = loadMovieArraySamples(movieArray)

nMovies = length(movieArray);




for j = 1:nMovies
    
    movieArray{j} = setupMovieData(movieArray{j});
    
    
    %Load the activity samples    
    
    load([movieArray{j}.activity.directory filesep movieArray{j}.activity.fileName])
    
    sampleArray(j).activity = allWindowSamples;
    
    clear allWindowSamples

    %Load the whole-cell samples
    if isfield(movieArray{j}.activity,'wholeCell')
        
        load([movieArray{j}.activity.directory filesep movieArray{j}.activity.wholeCell.fileName])    
        sampleArray(j).activity.wholeCell = wholeCellSamples;
        
    end
    %Load the whole-cell samples
%     if isfield(movieArray{j}.activity,'cellCenter')
        
        load([movieArray{j}.activity.directory filesep movieArray{j}.activity.cellCenter.fileName])    
        sampleArray(j).activity.wholeCell.center = cellCenterSamples;
        
%     end


    %Load the protrusion samples
   
    if checkMovieProtrusionSamples(movieArray{j})
        load([movieArray{j}.protrusion.directory filesep movieArray{j}.protrusion.samples.fileName])

        sampleArray(j).protrusion = protrusionSamples;

        clear protrusionSamples
    end
    if isfield(movieArray{j},'stimulation')
        sampleArray(j).stimulation.iFrame = movieArray{j}.stimulation.iFrame;        
    end
    
    
end