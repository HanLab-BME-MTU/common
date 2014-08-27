classdef ChannelReader < SubIndexReader
    %ChannelReader Creates a reader confined to a single channel
    %
    % The main advantage of using this is that it simplifies the parameters
    % needed for the getSize* methods as well as other Reader methods

    % Mark Kittisopikul
    % mark.kittisopikul@utsouthwestern.edu
    % Lab of Khuloud Jaqaman
    % UT Southwestern
    
    properties
        channel
    end
    
    methods
        function obj = ChannelReader(varargin)
            channel = [];
            if(nargin > 1)
                reader = varargin{1};
                c = varargin{2};
                assert(isscalar(c),'Channel number must be a scalar value');
                if(isa(reader,'Channel') && nargin == 1)
                   reader = reader.getReader();
                   % use c below
                   c = reader;
                end
                if(isa(c,'Channel'))
                    % save the channel as a property
                    channel = c;
                    % we could use subindex+1 instead
                    c = c.getChannelIndex();
                end
                if(nargin < 3 || varargin{3} == ':')
                    t = 1:reader.getSizeT(c);
                    varargin{3} = t;
                end
                if(nargin < 4 || varargin{4} == ':')
                    z = 1:reader.getSizeZ(c);
                    varargin{4} = z;
                end
            end
            obj = obj@SubIndexReader(varargin{:});
            if(~isempty(channel))
                obj.channel = channel;
            end
        end
        % obj.subIndices{1} 
        function s = getSizeX(obj)
            s = obj.reader.getSizeX(obj.subIndices{1});
        end
        function s = getSizeY(obj)
            s = obj.reader.getSizeY(obj.subIndices{1});
        end
        function b = getBitDepth(obj)
            b = obj.reader.getBitDepth(obj.subIndices{1});
        end
        function n = getImageFileNames(obj)
            n = obj.reader.getImageFileNames(obj.subIndices{1});
        end
        function n = getChannelNames(obj)
            n = obj.reader.getChannelNames(obj.subIndices{1});
        end
    end
    
end

