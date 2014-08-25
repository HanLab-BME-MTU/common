classdef TimeSeriesReader < CellReader
% TimeSeriesReader Reads XYT matrices as if arranged in a CxZ cell array
%

    % Mark Kittisopikul
    % mark.kittisopikul@utsouthwestern.edu
    % Lab of Khuloud Jaqaman
    % UT Southwestern
    methods
        function obj = TimeSeriesReader(reader)
            obj = obj@CellReader(reader);
        end
        function s = getSize(obj)
            s = [obj.reader.getSizeC(obj.sizeParam)
                 obj.reader.getSizeZ(obj.sizeParam)]';
        end
        function subIndices = getLinSub(obj,varargin)
            % linearize c z
            subIndices = varargin;
            subIndices([1 3]) = obj.getLinSub@CellReader(varargin{ 1:2:nargin-1 });
        end
        function matrix = toMatrix(obj)
            % obj.to3D uses toCell, so the dimensions are rearranged
            matrix = reshape(obj.to3D,[obj.getSizeY obj.getSizeX obj.getSizeT obj.size]);
        end
        % need to fix subindexing
        function out = toCell(obj)
            S.type = '{}';
            S.subs = {':',':'};
            S = obj.expandSubs(S);
            out = obj.loadCell(S.subs{:});
        end


    end
    methods ( Access = protected )
        function images = loadCell(obj,varargin)
            varargin(nargin :3) = {NaN};
            [cv,zv tv] = deal(varargin{:});
            if(isnan(tv))
                tv = 1 : obj.getSizeT(obj.sizeParam);
            end
            ndim = nargin - 1;

            subs = {cv tv zv};
            subsNan = cellfun(@(x) any(isnan(x)),subs);
            ctzImages = obj.loadCell@CellReader(subs{~subsNan});
            

            images = cell([length(cv) length(zv) 1]);
            for c = 1:length(cv)
                for z = 1:length(zv)
                        images{c,z} = cat(3, ctzImages{c,:,z} );
                end
            end
        end
        function R = getSubIndexReader(obj,S)
            classfcn = str2func(class(obj));
            S.subs{3} = 1: obj.getSizeT(obj.sizeParam);
            S.subs = S.subs([1 3 2]);
            R = classfcn(SubIndexReader(obj,S(1).subs{:}));
        end
    end
end
