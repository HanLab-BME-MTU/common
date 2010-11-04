
% Values from http://www.olympusfluoview.com/applications/fpcolorpalette.html

% Francois Aguet, October 2010

function lambda = name2wavelength(name)

switch name
    case 'CFP'
        lambda = 475e-9;
    case 'GFP'
        lambda = 509e-9;
    case 'EGFP'
        lambda = 507e-9;
    case 'Alexa488'
        lambda = 519e-9;
    case 'YFP'
        lambda = 527e-9;
    case 'Alexa555'
        lambda = 565e-9;
    case 'Alexa568'
        lambda = 603e-9;
    case 'RFP'
        lambda = 607e-9;
    case 'mCherry'
        lambda = 610e-9;
    case 'Alexa647'
        lambda = 665e-9;
    otherwise
        if isnumeric(name)
            lambda = name;
        else
            error('Shortcut not valid. Please enter wavelength in [nm].');
        end
end