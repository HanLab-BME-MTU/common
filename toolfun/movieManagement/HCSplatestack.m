function [chplatestack, starti, startsw] = HCSplatestack(varagin)
chplatestack = [];
dirurl = varagin;
file_lists = dir(fullfile(dirurl, '*.TIF')); ct = 0;
lenf = zeros(1, 10);
for i2 = 1:10
lenf(i2) = length(file_lists(i2).name); % see if there is inconsistent double digit naming
end
% if lenf(10) ~= mean(lenf)
%     renameSingledigitfiles(dirurl, file_lists, max(lenf));
% end
file_lists = dir(fullfile(dirurl, '*.TIF')); ct = 0;
%[starti, startsw] = getindexstart(file_lists(1).name);
h = waitbar(0, 'Loading HCS Images');
for i1 = 1:length(file_lists)
        waitbar(i1/length(file_lists));
        [starti, startsw] = getindexstart(file_lists(i1).name);
        wp = file_lists(i1,1).name(starti:max(startsw)+1); %well position
        wpv = double(wp(1))-64; %get numericle order from alphabetic sequence
        wph = str2double(wp(2:3));
        if min(abs(str2double(wp(min(startsw)-starti+3))-(0:9))) == 0 %see if site number is double digit.
        shp = str2double([wp(min(startsw)-starti+2),wp(min(startsw)-starti+3)]);
        chn = str2double(wp(max(startsw)-starti+2));
        else
            shp = str2double(wp(min(startsw)-starti+2));
            chn = str2double(wp(max(startsw)-starti+2));
        end
        ct = ct + 1;  
        if ct == 1
            wpvr = wpv; wphr = wph;
        end
        chplatestack{chn}{wpv-wpvr+1, wph-wphr+1}{shp} = [file_lists(i1,1).name];
end
close(h)

function renameSingledigitfiles(dirurl, file_lists, maxlength)
    [starti, startsw] = getindexstart(file_lists(1).name);
    ss = startsw(1);
    for i3 = 1:length(file_lists)
        if length(file_lists(i3).name) ~= maxlength
        movefile(strcat(dirurl, file_lists(i3).name), strcat(dirurl, file_lists(i3).name(1:ss), '0', file_lists(i3).name(ss+1:end)));
        end
    end


function [starti, startsw] = getindexstart(filename)
% % generate a-z string
for n = 1:26
    strr{n} = char(n -1 + 'a');
end
% generate 0-9 string
for n = 1:10
    stnum{n} = num2str(n-1);
end
index_char = regexp(filename, strr, 'ignorecase');
startsw = []; starti = [];
for i = 1:length(index_char) %it does not check the 'tif'
    if ~isempty(index_char{i})
        while ~isempty(index_char{i})
            if index_char{i}(1)+3<length(filename)
                if sum(strcmp(stnum, filename(index_char{i}(1)+1))) == 0
                    index_char{i}(1) = [];
                elseif sum(strcmp(stnum, filename(index_char{i}(1)+2))) > 0 ...
                        && sum(strcmp(stnum, filename(index_char{i}(1)+3))) == 0 ...
                        && isempty(starti)
                    starti = index_char{i}(1);
                elseif sum(strcmp({'s', 'w', 'S', 'W'}, filename(index_char{i}(1)))) > 0 ...
                        && sum(strcmp(stnum, filename(index_char{i}(1)+1))) > 0
                    startsw = [startsw index_char{i}(1)];
                    break;
                else
                    index_char{i}(1) = [];
                end
            else
                index_char{i}(1) = [];
            end
        end
    end
end
