% Francois Aguet, October 2010

function lambda = name2wavelength(name)

% values from http://www.olympusfluoview.com/applications/fpcolorpalette.html
switch name
    case 'GFP'
        lambda = 509e-9;
    case 'EGFP'
        lambda = 507e-9;
    case 'mCherry'
        lambda = 610e-9;
    otherwise
        if isnumeric(name)
            lambda = name;
        else
            error('Shortcut not valid. Please enter wavelength in [nm].');
        end
end
