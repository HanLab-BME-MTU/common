function movieData = convertMovieDataOS(movieData,toOS)

%Converts all the directory specifiers in the input movieData to a
%specified operating system. If no OS is specified, the current OS is used

%Hunter Elliott, 2/2009







%If the convert-to OS is not specified, use current OS
if nargin < 2 || isempty(toOS)
    toOS = computer;        
end
if nargin < 1
    errordlg('Must input the movieData!!',mfilename)
    return
end


%Check the movieData manually because setupMovieData uses this function
if ~isstruct(movieData) || ~isfield(movieData,'analysisDirectory') || ~isfield(movieData,'imageDirectory')
    errordlg('Invalid movieData!!',mfilename)
    return
end



%Also allow input of other OS strings
if any(strcmpi(toOS,{'Windows_NT', 'WIN','PC','Windows','PCWIN64','win32','win64'}))
    toOS = 'PCWIN';
end
%Include mac in the unix catergory
if any(strcmpi(toOS,{'UNIX','LINUX','LNX','GLNX86','GLNXA64','MACI','LNX86','LNXA64'}))
    toOS = 'GLNX';
end

%Determine what OS the movieData is in currently
fromOS = getMovieDataOS(movieData);

%Check if the from and to OSs are the same, do nothing
if strcmpi(fromOS,toOS)
    return
else %otherwise, convert all specified directories
    
    movieData = dirFieldConvert(movieData);       
    
end


