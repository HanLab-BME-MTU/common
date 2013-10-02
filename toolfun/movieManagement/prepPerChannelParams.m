function params = prepPerChannelParams(params,nChan)

if isfield(params,'PerChannelParams') && ~isempty(params.PerChannelParams) && iscell(params.PerChannelParams);
      
    nPCPar = numel(params.PerChannelParams);
    
    for j = 1:nPCPar
        
        nEl = numel(params.(params.PerChannelParams{j}));
        if  nEl == 1
            params.(params.PerChannelParams{j}) = repmat(params.(params.PerChannelParams{j}),[1 nChan]);
        elseif nEl ~= nChan
            error(['The parameter "' params.PerChannelParams{j} '" was designated as a per-channel parameter, but contained ' num2str(nEl) ' elements - this must be specified as either a scalar or have array of size equal to the number of channels!'])
        end                                    
    end        
end
