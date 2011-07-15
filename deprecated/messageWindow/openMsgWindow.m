function openMsgWindow
fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

% Opens/activates display
msgH=findobj('Tag','DISPLAY_MSG');
if isempty(msgH)
   msgH=displayMsg;
else
   figure(msgH)
end

