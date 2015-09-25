function ML=indexLSFMData(moviePaths,moviesRoot,varargin)
% - Batch sort time points in ch0/, ch1/ (and optionnaly deskew, dezip on demand) 
% - Create movieData and analysis folder for each movie
% - Print MIP in analysis folder
% - Create a movieList saved in an anlysis folder below <moviesRoot>
%
% EXAMPLE:  indexLSFMData('/path/to/your/cells/Cell_*/whatever/CAM_0{ch}*.tif','/path/to/your/cells','lateralPixelSize',160,'axialPixelSize',400,'copyFile',true)
%
% INPUT:  - moviePaths is either: 
%              ** A regular expression describing all the file of intereste using the following format: 
%                '/path/to/cell/Cell_*/path/to/tiff/xp_whatever_chanel_{ch}_*.tif'
%                - each images related to a cell must be in a separate folder
%                - the token {ch} represents the location of the channel number. It is mandatory.  
%                
%              ** a  matlab cell of path  that represent cell  directories containing *.tif files  ( or optionnaly compressed  *.tif.bz2) describing
%                 individual timepoint.
%                 - The option 'filePattern' must thus be used to describe channel name for example: 
%                   'Iter_ sample_scan_3p35ms_zp4um_ch{ch}_stack*'
%
%         - moviesRoot: original root for the list of movies:
%               ** contain the analysis folder with the movieList.mat and
%               associated outputdir
%               **  
%         
%
% OPTION: - newDataPath {[]}: if specificied, all the new
%           cell folders are writen under this folder ( analysis file with MD with or
%           without data) WARNING: copyFALSE default witll movie files. 
%         - copyFile {false}: if true copy the file instead of moving them. 
%         - 'writeData' {true}: If false, do not write data, only create MovieData under movieRoots or  NewDataPath



ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('moviePaths', @(x)(iscell(x)||ischar(x)));
ip.addRequired('moviesRoot', @ischar);
ip.addParamValue('movieListName','movieList.mat', @ischar);
ip.addParamValue('deskew',false, @islogical);
ip.addParamValue('writeData',true, @islogical);
ip.addParamValue('copyFile',false, @islogical);
ip.addParamValue('createMIP',true, @islogical);
ip.addParamValue('lateralPixelSize',1, @isfloat);
ip.addParamValue('axialPixelSize',1, @isfloat);
ip.addParamValue('timeInterval',1, @isfloat);
ip.addParamValue('newDataPath',[], @ischar);
ip.addParamValue('chStartIdx',1, @isnumeric);
ip.addParamValue('filePattern','Iter_ sample_scan_3p35ms_zp4um_ch{ch}_stack*', @ischar); 
ip.parse(moviePaths,moviesRoot,varargin{:});

p=ip.Results;


filePattern=p.filePattern;
if(ischar(moviePaths))
    [fileDirRegexp,fileRegexp,ext]=fileparts(moviePaths);
    filePattern=[fileRegexp ext];
    dirs=rdir([fileDirRegexp filesep]);    % filesep is important due to a bug in rdir ...
    moviePaths=unique(cellfun(@(x) fileparts(x),{dirs.name},'unif',0));    
end

MDs=cell(1,length(moviePaths));
parfor cellIdx=1:length(moviePaths)
    cPath=moviePaths{cellIdx};
    channelList=[];
    writeData=p.writeData;
    chIdx=ip.Results.chStartIdx;
    while(true)      
        if(p.deskew)
            outputDirCH=[cPath filesep 'deskew' filesep 'ch' num2str(chIdx)];
        else
            outputDirCH=[cPath filesep 'ch' num2str(chIdx)];
        end        
        if(~isempty(p.newDataPath))
            outputDirCH=strrep(outputDirCH,p.moviesRoot,p.newDataPath);
        end
        
        filelistCH=dir([cPath filesep strrep(filePattern,'{ch}',num2str(chIdx))]);
        
        % If it is empty, the reason might be that this as already been
        % carried out. In which case check out the channel folder. 
        if(isempty(filelistCH))
            writeData=false;
            filelistCH=dir([outputDirCH filesep strrep(filePattern,'{ch}',num2str(chIdx)) ] );
        end
        
        % If this channel does not exist, stop building channels. 
        if(isempty(filelistCH))
           break;
        end
        
        if(~exist(outputDirCH,'dir')) mkdir(outputDirCH); end;
        if(writeData)
            for fileIdx=1:length(filelistCH)
                file=filelistCH(fileIdx).name;
                
                if(p.deskew)
                    writeDeskewedFile([cPath filesep file],outputDirCH);
                else
                    if(p.copyFile)
                        copyfile([cPath filesep file],outputDirCH);
                    else
                        movefile([cPath filesep file],outputDirCH);
                    end
                end
                
            end
        end
        channelList=[channelList Channel(outputDirCH)];
        chIdx=chIdx+1;
    end
    

    if(~exist([cPath filesep 'analysis'],'dir')) mkdir([cPath filesep 'analysis']); end
    %%
    tiffReader=TiffSeriesReader({channelList.channelPath_},'force3D',true);
    %%
    MD=MovieData(channelList,[cPath filesep 'analysis'],'movieDataFileName_','movieData.mat','movieDataPath_',[cPath filesep 'analysis'], ...
                            'pixelSize_',p.lateralPixelSize,'pixelSizeZ_',p.axialPixelSize,'timeInterval_',p.timeInterval);
    MD.setReader(tiffReader);                    
    MD.sanityCheck();
    MD.save();
    if(p.createMIP)
        printMIP(MD);
    end
    MDs{cellIdx}=MD;
end


mkdir([moviesRoot filesep 'analysis']);
ML=MovieList(MDs,[moviesRoot filesep 'analysis'],'movieListFileName_',p.movieListName,'movieListPath_',[moviesRoot filesep 'analysis']);
ML.save();


function writeDeskewedFile(filePath,outputDir)
[~,~,ef]=fileparts(filePath);
written=false;
% This while loop handles server instability
while ~written
    try
        if(strcmp(ef,'.bz2'))
            unzipAndDeskewLatticeTimePoint(filePath,outputDir);
        elseif((strcmp(ef,'.tif')))
            deskewLatticeTimePoint(filePath,outputDir);
        else
            error('unsupported format')
        end
        written=true;
    catch
        written=false;
    end;
end

%% SNIPPETS
%     outputDirFinal=[cPath filesep 'deskew' filesep 'final'];
%     if(~isempty(p.rootPath))
%         outputDirFinal=strrep(outputDirFinal,p.rootPath,p.movieListOutputPath);
%     end
%     mkdir(outputDirFinal);
%     filelistFinal=dir([cPath filesep 'Iter_ sample_scan_3p35ms_zp4um_final*.*']);
%     for fileIdx=1:length(filelistFinal)
%         file=filelistFinal(fileIdx).name;
%         writeDeskewedFile([cPath filesep file],outputDirFinal);
%     end    
