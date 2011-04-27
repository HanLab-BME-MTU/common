function varargout = segmentationPackageGUI(varargin)
% Launch the GUI for the Segmentation Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Sebastien Besson 4/2011
%

varargout{1} = packageGUI(@SegmentationPackage,varargin{:});

end