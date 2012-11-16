function exportDatasetArray(a,varargin)
%exportDatasetArray is R2009a's export fn made backcompatible with R2008a
%
% This function writes dataset arrays into a .txt or .xls file, depending
% on the extension in the given filename. If no filename is given, the file
% will be called dataset.txt.
%
% Kathryn Applegate, 09/2009 - Matlab R2008a
%
%   EXPORT(DS,'File',FILENAME) writes the dataset array DS to a tab-delimited
%   text file, including variable names and observation names (if present).
%   If the observation names exist, the name in the first column of the first
%   line of the file is the first dimension name for the dataset (by default,
%   'Observations').  EXPORT overwrites any existing file named FILENAME.
%
%   EXPORT(DS) writes to a text file whose default name is the name of the
%   dataset array DS, appended by '.txt'.  If EXPORT cannot construct the file
%   name from the dataset array input, it writes to the file 'dataset.txt'.
%   EXPORT overwrites any existing file.
%
%   EXPORT(DS,'File',FILENAME,'Delimiter',DELIM) writes the dataset array DS
%   to a text file, using the delimiter DELIM.  DELIM can be any of ' ', '\t',
%   ',', ';', '|' or their corresponding string names 'space', 'tab', 'comma',
%   'semi', or 'bar'.
%
%   EXPORT(DS,'XLSFile',FILENAME) writes the dataset array DS to an Excel
%   spreadsheet file, including variable names and observation names (if
%   present).  You may also specify the 'Sheet' and 'Range' parameter
%   name/value pairs, with parameter values as accepted by the XLSREAD
%   function.
%
%   EXPORT(DS, ..., 'WriteVarNames',FALSE) does not write the variable names
%   to the text file.  EXPORT(..., 'WriteVarNames',TRUE) writes the names as
%   column headings in the first line of the file, and is the default.
%
%   EXPORT(DS ..., 'WriteObsNames',FALSE) does not write the observation names
%   to the text file.  EXPORT(..., 'WriteObsNames',TRUE) writes the names as
%   the first column of the file, and is the default.
%
%   In some cases, EXPORT creates a text file that does not represent DS
%   exactly, as described below.  If you use DATASET('File',FILENAME) to read
%   that file back in and create a new dataset array, the result may not have
%   exactly the same format or contents as the original dataset array.
%
%   *  EXPORT writes out numeric variables using long g format, and
%      categorical or character variables as unquoted strings.
%   *  For non-character variables that have more than one column, EXPORT
%      writes out multiple delimiter-separated fields on each line, and
%      constructs suitable column headings for the first line of the file.
%   *  EXPORT writes out variables that have more than two dimensions as two
%      dimensional variables, with trailing dimensions collapsed.
%   *  For cell-valued variables, EXPORT writes out the contents of each cell
%      as a single row, in multiple delimiter-separated fields, when the
%      contents are numeric, logical, character, or categorical, and writes
%      out a single empty field otherwise.
%   *  EXPORT writes out both the time and the data fields of timeseries
%      variables as separate columns.
%
%   Save DS as a mat file if you need to import it again as a dataset array.
%
%   Examples:
%      load hospital
%
%      % Write to a comma-delimited file 'hospital.txt'
%      export(hospital,'delimiter',',')
%
%      % Write to a comma-delimited file 'dataset.txt'
%      export(hospital(1:10,1:3),'delimiter',',')
%
%      % Write to a comma-delimited file 'Under40.dat'
%      export(hospital(hospital.Age<40,:),'File','Under40.dat','delimiter',',')
%      
%   See also DATASET.

%   Copyright 2008 The MathWorks, Inc. 
%   $Revision: 1.1.6.3 $  $Date: 2008/12/01 07:41:31 $


fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);


t=cell(length(varargin)/2,2);
c=1;
for iVar=1:length(varargin)/2
    t(iVar,1:2)=varargin(c:c+1);
    c=c+2;
end
clear varargin;
varargin=t;

% check whether 'WriteVarNames' was an input
writevarnames = 1;
idx=find(cellfun(@(x) strcmpi(x,'WriteVarNames'),varargin(:,1)));
if ~isempty(idx)
    if ischar(varargin{idx,2}) && strcmpi(varargin{idx,2},'false')==1
        writevarnames = 0;
    elseif ~ischar(varargin{idx,2}) && varargin{idx,2}~=1
        writevarnames = 0;
    end
end

% check whether 'WriteObsNames' was an input
writeobsnames = 1;
idx=find(cellfun(@(x) strcmpi(x,'WriteObsNames'),varargin(:,1)));
if ~isempty(idx)
    if ischar(varargin{idx,2}) && strcmpi(varargin{idx,2},'false')==1
        writeobsnames = 0;
    elseif ~ischar(varargin{idx,2}) && varargin{idx,2}~=1
        writeobsnames = 0;
    end
end


% check whether 'File' was an input
idx=find(cellfun(@(x) strcmpi(x,'file'),varargin(:,1)));
if isempty(idx)
    fileArg='dataset.txt';
else
    isTxt=~isempty(strfind(varargin{idx,2},'.txt'));
    isXls=~isempty(strfind(varargin{idx,2},'.xls'));
    if isTxt
        fileArg=varargin{idx,2};
        xlsfileArg=[];
    elseif isXls
        fileArg=[];
        xlsfileArg=varargin{idx,2};
    end
end


% put a in the 2009 format
temp=a;
clear a;
a.obsnames=temp.Properties.ObsNames;
a.nobs=length(a.obsnames);
a.varnames=temp.Properties.VarNames;
a.nvars=length(a.varnames);
a.data=cell(1,a.nvars);
for i=1:a.nvars
   a.data(i)={temp.(a.varnames{i})}; 
end
a.props.DimNames=temp.Properties.DimNames;

if ~isempty(fileArg) || isempty(xlsfileArg)
    % Create a default name if needed
    if isempty(fileArg)
        dsname = inputname(1);
        if isempty(dsname)
            dsname = 'dataset';
        end
        fileArg = [dsname '.txt'];
    end
    writeFile(a,fileArg,writevarnames,writeobsnames,varargin{:});
else
    try
        Excel = actxserver('Excel.Application');
    catch me
        error('stats:dataset:export:NoCOMServer', ...
              'Could not start Excel server for export.  Export to a text file instead.');
    end
    Excel.Quit;
    
    writeXLSFile(a,xlsfileArg,writevarnames,writeobsnames,varargin{:});
end

end


%-----------------------------------------------------------------------
function writeFile(a,filename,writevarnames,writeobsnames,varargin)
%WRITEFILE Write a dataset array to a text file

idx=find(cellfun(@(x) strcmpi(x,'Delimiter'),varargin(:,1)));
if isempty(idx)
    delimiter='tab';
else
    delimiter=varargin{idx,2};
end

tab = sprintf('\t');
lf = sprintf('\n');

% Set the delimiter
switch delimiter
case {'tab', '\t'}
  delimiter = tab;
case {'space',' '}
  delimiter = ' ';
case {'comma', ','}
  delimiter = ',';
case {'semi', ';'}
  delimiter = ';';
case {'bar', '|'}
  delimiter = '|';
otherwise
  error('stats:dataset:export:UnrecognizedDelimiter', ...
        sprintf('Unrecognized delimiter character: ''%c''.',delimiter(1)));
end

realDoubleFmt = ['%.15g' delimiter];
complexDoubleFmt = ['%.15g+%.15gi' delimiter];
realSingleFmt = ['%.7g' delimiter];
complexSingleFmt = ['%.7g+%.7gi' delimiter];
realIntegerFmt = ['%d' delimiter];
complexIntegerFmt = ['%d+%di' delimiter];

% Open the file for writing.
fid = fopen(filename,'wt'); % text mode: CRLF -> LF

if fid == -1
   error('stats:dataset:ExportOpenFailed', ...
         'Unable to open file ''%s'' for writing.', filename);
end

adata = a.data;
varType = zeros(1,a.nvars);
ncellCols = cell(1,a.nvars);
for j = 1:a.nvars
    varj = adata{j};
    varType(j) = typeCode(varj);
    if iscell(varj)
        ncellCols{j} = max(cellfun(@ncolsCell,varj),[],1);
    end
end

% Write variable names to the first line of the file as column headers
if writevarnames
    str = '';

    % Start with a column header for observation names
    if writeobsnames
        str = sprintf('%s%s',a.props.DimNames{1},delimiter);
    end

    for j = 1:a.nvars
        varnamej = colHeaders(adata{j},a.varnames(j),ncellCols{j});
        for jj = 1:length(varnamej)
            str = [str sprintf('%s%s',varnamej{jj},delimiter)];
        end
    end

    % Write out the header line
    if ~isempty(str), str(end) = lf; end % replace trailing delimiter with '\n'
    fprintf(fid,'%s',str);
end

% Write each row of the dataset array to the file
buf = ''; bufLen = 0;
for i = 1:a.nobs
    if writeobsnames
        str = sprintf('%s%s',a.obsnames{i},delimiter);
        start = bufLen+1; bufLen = bufLen+length(str);
        buf(start:bufLen) = str;
    end
    for j = 1:a.nvars
        varj = adata{j};
        type = varType(j);

        if type == 7
            % Write out the time field for a timeseries
            str = writeElement(varj.time(i),2);
            start = bufLen+1; bufLen = bufLen+length(str);
            buf(start:bufLen) = str;
            if ~varj.IsTimeFirst, varj = varj'; end
            varj = varj.data;
            type = typeCode(varj);
        end

        if type == 8 % cell variable
            str = writeCell(varj(i,:),ncellCols{j});
        elseif type == 1 % char variable
            str = writeElement(varj(i,:,:),type); % keep each row separate
        else % type < 8, non-cell variable
            str = writeElement(varj(i,:),type); % N-D to 2-D
        end
        start = bufLen+1; bufLen = bufLen+length(str);
        buf(start:bufLen) = str;
    end
    if ~isempty(buf), buf(bufLen) = lf; end % replace trailing delimiter with '\n'
    if mod(i,10) == 0
        fprintf(fid,'%s',buf(1:bufLen));
        bufLen = 0;
    end
end
if mod(a.nobs,10) ~= 0, fprintf(fid,'%s',buf(1:bufLen)); end

fclose(fid);


    function str = writeCell(x,ncellCols)
    % WRITECELL Write a cell-valued dataset array element to the file
        [dum,ncols] = size(x);
        str = '';
        for k = 1:ncols
            cellk = x{k};
            cellType = typeCode(x{k});
            
            % Pad with empty fields if necessary to match other rows.
            pad = repmat(delimiter,1,ncellCols(k)-ncolsCell(cellk));
            if cellType == 1
                % Write the cell as a set of char strings.  
                str = sprintf('%s%s%s',str,writeElement(cellk,cellType),pad);
            else
                % Write the cell as a row of values.  writeElement will write an
                % empty field if the contents are not atomic.
                str = sprintf('%s%s%s',str,writeElement(cellk(:)',cellType),pad);
            end
        end
    end

    function str = writeElement(x,type)
    % WRITEELEMENT Write a non-cell dataset array element to the file
    
    if issparse(x), x = full(x); end % Only one row, no memory issues
    
    switch type
    case 1 % char
        if size(x,1) == 1
            str = [x,delimiter];
        else
            % Write out each row as a separate field in the string, including
            % rows in higher dims.
            x = permute(x,[2 1 3:ndims(x)]);
            x = x(:,:); x(end+1,:) = ',';
            str = x(:)';
        end
    case 2 % double
        if isreal(x)
            str = sprintf(realDoubleFmt,x); % may be multiple values
        else
            str = sprintf(complexDoubleFmt,[real(x); imag(x)]); % may be multiple values
        end
    case 3 % single
        if isreal(x)
            str = sprintf(realSingleFmt,x); % may be multiple values
        else
            str = sprintf(complexSingleFmt,[real(x); imag(x)]); % may be multiple values
        end
    case {4 5} % integer/logical
        if isreal(x)
            str = sprintf(realIntegerFmt,x); % may be multiple values
        else
            str = sprintf(complexIntegerFmt,[real(x); imag(x)]); % may be multiple values
        end
    case 6 % categorical
        if size(x,2) == 1
            str = [char(x) delimiter];
        else
            % Can't use strcat, that drops trailing whitespace, including tabs
            str = cell2mat(cellfun(@(str)[str delimiter],cellstr(x),'UniformOutput',false));  % may be multiple values
        end
    case {0,7,8} % not a known type, or a timeseries, or a cell
        % Don't attempt to write out
        str = sprintf('%s',delimiter);
    end
    end

end % writeFile function


%-----------------------------------------------------------------------
function writeXLSFile(a,filename,writevarnames,writeobsnames,varargin)
%WRITEXLSFILE Write a dataset array to an Excel spreadsheet file

sheet=1;
range='A1';


lr = '';
validRange = false;
if ischar(range) && isvector(range) && size(range,1)==1
    rangeSplit = regexpi(range,':','split');
    if length(rangeSplit)==1 || length(rangeSplit)==2
        ul = rangeSplit{1};
        [row,col] = spec2rowcol(ul);
        validRange = ~(isnan(row) || isnan(col));
        if validRange && length(rangeSplit) > 1
            lr = rangeSplit{2};
            [row2,col2] = spec2rowcol(lr);
            validRange = ~(isnan(row2) || isnan(col2));
        end
    end
end
if ~validRange
   error('stats:dataset:InvalidRange', ...
         'The range must be a string of the form ''A1:B2'' or ''A1''.');
end
if isempty(lr)
    maxCol = Inf;
else
    % convert ur:ll to ul:lr
    tmp = min(col,col2); col2 = max(col,col2); col = tmp;
    tmp = min(row,row2); row2 = max(row,row2); row = tmp;
    lr = [':' rowcol2spec(row2,col2)];
    maxCol = col2;
end

% Write observation names.
vars = cell(a.nobs+writevarnames,0);
ulcol = col;
if writeobsnames
    obsnames = a.obsnames;
    if writevarnames
        obsnames = [a.props.DimNames{1}; obsnames];
    end
    vars = [vars obsnames];
    col = col + 1;
end

for j = 1:a.nvars
    if col > maxCol, break; end
    
    varj = a.data{j};
    varnamej = a.varnames(j);
    ists = false;

    if iscell(varj)
        % xlswrite cannot write out non-scalar-valued cells -- convert cell
        % variables to a cell of the appropriate width containing only
        % scalars.
        [dum,ncols] = size(varj);
        ncellColsj = max(cellfun(@ncolsCell,varj),[],1);
        newNumCols = sum(ncellColsj);
        newVarj = cell(a.nobs,newNumCols);
        
        % Expand out each column of varj into as many columns as needed to
        % have only scalar-valued cells, possibly padded with empty cells.
        cnt = 0;
        for jj = 1:ncols
            varjj = varj(:,jj);
            num = ncellColsj(jj);
            newVarjj = cell(a.nobs,num);
            for i = 1:a.nobs
                varjj_i = varjj{i};
                if ischar(varjj_i)
                    vals = char2cell(varjj{i});
                elseif isnumeric(varjj_i) || islogical(varjj_i)
                    vals = num2cell(varjj_i);
                elseif isa(varjj_i,'categorical')
                    vals = cellstr(varjj_i);
                else
                    vals = cell(0,0); % write out only an empty cell
                end
                newVarjj(i,1:numel(vals)) = vals(:)';
            end
            newVarj(:,cnt+(1:num)) = newVarjj;
            cnt = cnt + num;
        end
        
        varj = newVarj;
        ncols = newNumCols;
        ncellColsj = ones(1,newNumCols);
    else
        ncellColsj = [];
        % xlswrite will convert any input to cell array anyway, may as well do
        % it here in all cases to get correct behavior for character and for
        % cases xlswrite won't handle.
        if ischar(varj)
            varj = char2cell(varj);
        elseif isnumeric(varj) || islogical(varj)
            varj = num2cell(varj(:,:));
        elseif isa(varj,'categorical')
            varj = cellstr(varj(:,:));
        elseif isa(varj,'timeseries')
            ists = true;
            if ~varj.IsTimeFirst, varj = varj'; end
            varj = [num2cell(varj.time) num2cell(varj.data(:,:))];
        else % write out empty cells
            varj = cell(a.nobs,1);
        end
    end

    [dum,ncols] = size(varj);
    if writevarnames
        if ists
            vn = varnamej;
            if ncols > 2
                varnamej = strcat(varnamej,'_',num2str((1:ncols-1)'))';
            end
            varnamej = [strcat(vn,'_time') varnamej];
        elseif ncols > 1
            varnamej = strcat(varnamej,'_',num2str((1:ncols)'))';
        end
        varj = [varnamej; varj];
    end
    vars = [vars varj];
    col = col + ncols;

    if mod(j,10) == 0
        writeXLSVars();
        vars = cell(a.nobs+writevarnames,0);
        ulcol = col;
    end
end
if mod(j,10) ~= 0, writeXLSVars(); end

    function writeXLSVars
        ul = rowcol2spec(row,ulcol);
        [success,msg] = xlswrite(filename,vars,sheet,[ul lr]);
        if ~success
            error('stats:dataset:export:XlswriteFailed', ...
                  'Error writing dataset variable ''%s'' to ''%s'':\n%s.',a.varnames{j},filename,msg.message);
        end
    end

end % writeXLSFile function


%-----------------------------------------------------------------------
function varnamej = colHeaders(varj,varnamej,ncellColsj)
%COLHEADERS Create multiple column headers from a dataset variable name

% Need an extra column header for the time data of a timeseries
ists = isa(varj,'timeseries');
if ists
    timename = strcat(varnamej,'_time');
    if ~varj.IsTimeFirst, varj = varj'; end
    varj = varj.data;
end

% Need multiple column headers if the variable has multiple columns.
if ischar(varj)
    ncols = 1;
else
    [dum,ncols] = size(varj); % Treat N-D as 2-D.
end
if ncols > 1
    varnamej = strcat(varnamej,'_',num2str((1:ncols)'))';
end

% Need multiple column headers if the variable is a cell containing non-scalars.
if iscell(varj) && any(ncellColsj > 1)
    vnj = cell(1,sum(ncellColsj));
    cnt = 0;
    for jj = 1:ncols
        num = ncellColsj(jj);
        vnj(cnt+(1:num)) = strcat(varnamej(jj),'_',num2str((1:num)'))';
        cnt = cnt + num;
    end
    varnamej = vnj;
end

if ists
    varnamej = [timename varnamej];
end

end % colHeaders function


%-----------------------------------------------------------------------
function type = typeCode(x)
% TYPECODE Return the type of a variable, coded as an integer.
if ischar(x)
    type = 1;
elseif isa(x,'double')
    type = 2;
elseif isa(x,'single')
    type = 3;
elseif isinteger(x)
    type = 4;
elseif islogical(x)
    type = 5;
elseif isa(x,'categorical')
    type = 6;
elseif isa(x,'timeseries')
    type = 7;
elseif iscell(x)
    type = 8;
else
    type = 0; % other, in this case not a standard type
end
end


%-----------------------------------------------------------------------
function m = ncolsCell(c)
if ischar(c)
    [n,dum,d] = size(c);
    n = max(n,1); % treat '' as a single row
    m = n*d; % each row is one "column", but keep rows and pages separate
elseif isnumeric(c) || islogical(c) || isa(c,'categorical')
    m = max(numel(c),1); % always write out at least one empty field
else
    m = 1; % other types are written as an empty field
end
end


%-----------------------------------------------------------------------
function cs = char2cell(c)
% Treat each row as a separate string, including rows in higher dims.

[n,dum,d] = size(c); szOut = [n,d];
if ndims(c) > 2
    c = permute(c,[2 1 3:ndims(c)]);
    c = c(:,:)';
end
if isempty(c)
cs=cellstr(c); % added to avoid reshaping error due to size change
else
cs = reshape(cellstr(c),szOut);
end
end

%-----------------------------------------------------------------------
function spec = rowcol2spec(row,col)
mult676 = floor((col-1)/676);
mult26 = floor((col-mult676*676-1)/26);
mod26 = mod(col-1,26);
if col <= 26
    spec = sprintf('%s%d',char('A'+mod26),row);
elseif col <= 702
    spec = sprintf('%s%s%d',char('A'+mult26-1),char('A'+mod26),row);
else
    spec = sprintf('%s%s%s%d',char('A'+mult676-1),char('A'+mult26-1),char('A'+mod26),row);
end
end


%-----------------------------------------------------------------------
function [row,col] = spec2rowcol(corner)
row = NaN; col = NaN;
tok = regexpi(corner,'^([a-z]{1,3})([0-9]+)$','tokens');
if ~isempty(tok)
    tok = tok{1};
    if length(tok) == 2
        c = tok{1}; r = tok{2};
        mults = [676; 26; 1]; % 'a' => 1, 'xfd' => 16384
        col = (lower(c) - 'a' + 1) * mults((end-length(c)+1):end);
        row = str2double(r);
    end
    % Let xlswrite decide if the row or column is too large
end
end
