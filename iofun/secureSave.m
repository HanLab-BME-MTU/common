function secureSave(varargin)
%wrapper for save command
%this version never overwrites an existing file, but adds a number to
%the filename instead

%c: 10/10/02 dT

fname=varargin{1};
[path,body,nr,ext]=getFilenameBody(fname);
if isempty(path)
    path=pwd;
end;
if isempty(ext)
    ext='.mat';
end;
fname=[path filesep body num2str(nr) ext];
%if filename already exists, add a (increasing) number at end
while(exist(fname,'file'))
    if isempty(nr)
        nr=1;
    else
        nr=nr+1;
    end;
    fname=[path  filesep body num2str(nr) ext];
end;
vars=sprintf('''%s'',', varargin{2:end});
vars=vars(1:end-1);
savecmd=['save(''' fname ''',' vars ');'];
evalin('caller',savecmd);
