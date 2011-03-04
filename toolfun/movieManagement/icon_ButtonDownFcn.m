function icon_ButtonDownFcn(hObject, eventdata)
% This function call up a help dialog box when user click any of the icons
% in all GUIs.
%
handles = guidata(hObject);

ud = get(hObject, 'UserData');

if ~isempty(ud) && isfield(ud, 'class')
%% Open original txt help file. 
% Help files are specified by variable "helpFile"
    
%     switch ud.class
%         
%         % ----------- Package -------------------
%         
%         case 'UTrackPackage'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'SegmentationPackage'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'BiosensorPackage'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         % ----------- Process -------------------
%         
%         % Segmentation Package
%         case 'SegmentationProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%         
%         case 'ThresholdProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         % Biosensors Package
%         case 'BackgroundMasksProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'MaskRefinementProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'DarkCurrentCorrectionProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'ShadeCorrectionProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'BackgroundSubtractionProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'TransformationProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'BleedthroughCorrectionProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'RatioProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'PhotobleachCorrectionProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'OutputRatioProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         % UTrack Package            
%         case 'DetectionProcess'
%             helpFile = 'CharlesHelpFile.pdf'; 
%             
%         case 'SubResolutionProcess'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'TrackingProcess'
%             helpFile = 'CharlesHelpFile.pdf';       
%         
%         % ----------- Other Tools (not object) -------------------
%         
%         case 'MovieData'
%             helpFile = 'CharlesHelpFile.pdf';
%         
%         case 'costMatLinearMotionLink2GUI'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'costMatLinearMotionCloseGaps2GUI'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'kalmanInitializationGUI'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'trackingVisualGUI'
%             helpFile = 'CharlesHelpFile.pdf';  
%             
%         case 'detectionVisualGUI'
%             helpFile = 'CharlesHelpFile.pdf';     
%             
%         case 'overlayFeaturesMovie'
%             helpFile = 'CharlesHelpFile.pdf';
%             
%         case 'plotTracks2D'
%             helpFile = 'CharlesHelpFile.pdf';  
%             
%         case 'plotCompTrack'
%             helpFile = 'CharlesHelpFile.pdf';         
%             
%         case 'overlayTracksMovieNew'
%             helpFile = 'CharlesHelpFile.pdf';              
%             
%         otherwise
%             warndlg('The help file has not been created yet.', 'modal');
%             return
%     end
    helpFile=[ud.class '.pdf'];
    if exist(helpFile, 'file')
        open(helpFile)
    else
        warndlg(['Cannot find Help file:', helpFile], 'modal');
        return
    end

else


%% Open GUI based help window
    
    userData = get(handles.figure1, 'UserData');
    % Help dialog from MovieData panel
    splitTag = regexp(get(get(hObject,'parent'), 'tag'), '_','split');

    % Pass handle to userData
    % If from package GUI, call pre-defined help dialog
    % if called from setting GUI, call user-defined help dialog 'msgboxGUI'

    if isfield(userData, 'crtProc')
        % Help dialog from setting panel
        if ~isempty(userData.crtProc)
            userData.helpFig = msgboxGUI('Text', sprintf([get(hObject,'UserData'), ...
                '\n', userfcn_copyright ]),'Title',['Help - ' userData.crtProc.name_] );
        else
            userData.helpFig = msgboxGUI('Text', sprintf([get(hObject,'UserData'), ...
                '\n', userfcn_copyright ]),'Title','Help');
        end


    elseif strcmp(splitTag{1}, 'axes') && length(splitTag) >1

            if strcmpi(splitTag{2}, 'help') % Help icon
                if length(splitTag) < 3
                    % Package help
                    userData.packageHelpFig = msgbox(sprintf(get(hObject,'UserData')), ...
                        ['Help - ' userData.crtPackage.name_], 'custom', get(hObject,'CData'), userData.colormap, 'replace');
                else
                    % Process help
                    procID = str2double(splitTag{3});
                    if ~isnan(procID)

                        procName = regexp(userData.crtPackage.processClassNames_{procID}, 'Process','split');
                        userData.processHelpFig(procID) = msgbox(sprintf(get(hObject,'UserData')), ...
                         ['Help - ' procName{1}], 'custom', get(hObject,'CData'), userData.colormap, 'replace');
                    end
                end

            else % Process status icon

                userData.iconHelpFig = msgbox(get(hObject,'UserData'), ...
                    'Help', 'custom', get(hObject,'CData'), userData.colormap, 'replace');            
            end

    else
        userData.iconHelpFig = msgbox(get(hObject,'UserData'), ...
            'Help', 'custom', get(hObject,'CData'), userData.colormap, 'replace'); 
    end

    set(handles.figure1, 'UserData', userData);
end
