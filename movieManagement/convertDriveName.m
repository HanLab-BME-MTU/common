function outDrive = convertDriveName(inDrive)

%Converts the input drive name form unix to windows or vice versa using the
%known LCCB drives

%Hunter Elliott


%Matrix with drive matches
driveMat = {'H:' , '~' ;
                  'L:','/mnt/lccb';
                  'M:','/mnt/projects';
                  'P:','/mnt/public';
                  'S:','/mtn/fsm'};
    
              
              
winMatch = strcmpi(inDrive,driveMat);

%If PC
if any(winMatch(:,1))
    
    outDrive = driveMat(winMatch(:,1),2);
    
elseif  any(winMatch(:,2))
    outDrive = driveMat(winMatch(:,2),1);
    
else
    outDrive = [];
    errordlg('Unrecognized drive name!',mfilename)
end
