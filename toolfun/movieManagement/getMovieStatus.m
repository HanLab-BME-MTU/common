function varargout = getMovieStatus(movieData)

%statusBinary = getMovieStatus(movieData)
%This displays a graphical representation of the status of the various
%processing steps that may/may not have been performed on the input
%movieData(s). The input movieData may be a cell array of movieData or a
%single movieData.
%If return arguments are requested, the status is also returned as an
%array. NOT IMPLEMENTED YET

%Hunter Elliott, 3/2009

% Input Processing

if ~iscell(movieData) %If a single movieData was input.
    movieData = {movieData};
end

nMovies = length(movieData);


% list of potential processing steps
allProcessing = { ...
    'analysisDirectory',...
    'contours',...
    'windows',...
    'activity',...
    'protrusion',...
    {'protrusion','samples'},...    
     'masks',...
     'FRETprocessing',...
     {'FRETprocessing','fretShadeDirectory'},{'FRETprocessing','donorShadeDirectory'}};

 
 %allProcessing = allProcessing([1 7:end]);
allProcessing = allProcessing([1 2]);
 
nProcesses = length(allProcessing);

warning('off','MATLAB:tex')
%figPos = get(0,'ScreenSize') .* .8;
%figPos(1:2) = figPos(1:2) + 100;
figPos = [ 1172         524         719         626];
statFig = figure('Position',figPos);
hold on
figLoc = get(statFig,'Position');

wtBar = waitbar(0,'Please wait, getting movie stati....');
for j = 1:nProcesses
        
    
    
    for k = 1:nMovies                        
        
        movieData{k} = refreshMovieData(movieData{k});               
        
        subplot(nMovies,nProcesses,(k-1)*nProcesses + j )
        hold on
        axis off
        
        if k == 1
            %Draw the process name
             text(mean(xlim),min(ylim)+1.5,allProcessing{j},'HorizontalAlignment','Center','FontSize',8)                         
        end
        
        %Check if the field exists and fill the axes accordingly       
        if j > 1 
            if ~iscell(allProcessing{j})
                if isfield(movieData{k},allProcessing{j})
                    xs = xlim;
                    ys = ylim;
                    isField = true;
                    if isfield(movieData{k}.(allProcessing{j}),'status') && movieData{k}.(allProcessing{j}).status == 1
                        isGood  = true;

                       
                    else
                        isGood = false;
                        
                    end


                else
                    isField = false;
                end
            else %If its a sub-field
                catStr = {'movieData{k}'};
                isField = true;
                for m = 1:length(allProcessing{j})                    
                                        
                    evString = strcat( 'isfield(' , catStr{1}, ',''', allProcessing{j}{m}, ''')' );
                    if isField && eval(evString) 
                        isField = true;
                        catStr = strcat(catStr,'.',allProcessing{j}(m));
                        evString = strcat('isfield(' , catStr{1}, ',''status'')');
                        if eval(evString) && eval(strcat(catStr{1}, '.status == 1'))
                            isGood = true;
                        else
                            isGood = false;
                        end
                    else
                        isField = false;
                    end    
                    
                end
                
                
            end
            
            if isField && isGood
                 fill([xs(1) xs(1) xs(2) xs(2)],[ys(1) ys(2) ys(2) ys(1)],'b')                                                        
            elseif isField
                plot([xs(1) xs(2)],[ys(1) ys(2)],'r')
                plot([xs(2) xs(1)],[ys(1) ys(2)],'r')                
            end
            
            
            
        end
        
        
        
        %Get axes location in figure coordinates
        currAxesLoc = get(gca,'OuterPosition');        
        
        %Convert these to screen coordinates
        axPos(1) = figLoc(1) + currAxesLoc(1) * figLoc(3);
        axPos(2) = figLoc(2) + currAxesLoc(2) * figLoc(4);
        axPos(3) = figLoc(3) * currAxesLoc(3);
        axPos(4) = figLoc(4) * currAxesLoc(4);
        
        if j == 1
            %Draw the analysis directory next to the subplot
            text(max(xlim),mean(ylim),movieData{k}.analysisDirectory,'HorizontalAlignment','Right','FontSize',8,'Interpreter','none')
        end
        
        
        waitbar( (nMovies*(j-1) + k ) / (nProcesses*nMovies),wtBar);
    
        if j == nProcesses
            text(max(xlim)+.2,mean(ylim),num2str(k),'HorizontalAlignment','Left','FontSize',8)    
        
        end
        
    end
    
    
    
end


    
varargout = {};
close(wtBar);