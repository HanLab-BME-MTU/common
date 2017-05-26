function path=sprintfPath(XYProjTemplate,fIdx)
% making sprintf windows-path-proof
% sprintf applied on the filename only
[folder,file,ext]=fileparts(XYProjTemplate);
filename=sprintf(file,fIdx);
path=fullfile(folder,[filename ext]);
 