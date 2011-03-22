
% Values from http://www.olympusfluoview.com/applications/fpcolorpalette.html

% Francois Aguet, October 2010

function lambda = name2wavelength(name)

if iscell(name)
    lambda = cellfun(@(x) convert(x), name);
else
    lambda = convert(name);
end


function lambda = convert(name)

switch lower(name)
    case {'bfp', 'ebfp'}
        lambda = 440e-9;
    case 'cfp'
        lambda = 475e-9;
    case 'egfp'
        lambda = 507e-9;
    case 'gfp'
        lambda = 509e-9;
    case 'alexa488'
        lambda = 519e-9;
    case 'yfp'
        lambda = 527e-9;
    case 'alexa555'
        lambda = 565e-9;
    case {'dtomato', 'tdtomato'}
        lambda = 581e-9;
    case 'dsred'
        lambda = 583e-9;
    case 'tagrfp'
        lambda = 584e-9;
    case 'alexa568'
        lambda = 603e-9;
    case {'rfp', 'mrfp'}
        lambda = 607e-9;
    case 'mcherry'
        lambda = 610e-9;
    case 'texasred'
        lambda = 615e-9;
    case 'alexa647'
        lambda = 665e-9;
    otherwise
        if isnumeric(name)
            lambda = name;
        else
            lambda = NaN;
            %error('Shortcut not valid. Please enter wavelength in [nm].');
        end
end