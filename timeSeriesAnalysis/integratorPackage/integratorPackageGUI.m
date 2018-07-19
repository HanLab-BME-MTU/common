function varargout = integratorPackageGUI(varargin)
% Launch the GUI for the Windowinf Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%

% Sebastien Besson, Sep 2011

varargout{1} = packageGUI('IntegratorPackage',varargin{:});

end