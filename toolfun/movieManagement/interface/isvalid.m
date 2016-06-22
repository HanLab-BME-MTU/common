function logic= isvalid(aDouble)
% Dirty trick to handle isvalid absence in <2013b
warning('Using shim isvalid');
logic=true;
