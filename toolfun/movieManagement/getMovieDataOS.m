function osString = getMovieDataOS(movieData)

%Determines what OS the directory specifiers in movieData are in and
%returns the computer string for that OS


%Hunter Elliott, 2/2009

osString =[];


%Check the image and analysis directories for OS


%Check for WIN filesep
if any(regexp(movieData.imageDirectory,'\'))
    
    %Make sure the analysis directory matches
    if any(regexp(movieData.analysisDirectory,'\'))
        osString = 'PCWIN';
    else
        errordlg('!!Warning!! OS conflict in movieData!!',mfilename)
        return        
    end
end

%Check for unix-based filesep (Linus Torvalds himself admitted Linux was
%based on UNIX
if any(regexp(movieData.imageDirectory,'/'))
    
    if ~isempty(osString)
        errordlg('!!Warning!! OS conflict in movieData!!',mfilename)
        return        
    end
    
    %Make sure the analysis directory matches
    if any(regexp(movieData.analysisDirectory,'/'))
        osString = 'GLNX';
    else
        errordlg('!!Warning!! OS conflict in movieData!!',mfilename)
        return        
    end
end