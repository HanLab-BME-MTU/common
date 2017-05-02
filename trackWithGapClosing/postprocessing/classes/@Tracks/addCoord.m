function obj=addCoord(obj,tracks)
   % tracks.addCoord(TracksToAdd)
   % addition coordinate, lifetime must be the same. 
   % TODO: tracks lifetime must overlap entirely with obj
   % Philippe Roudot 2017
   if(length(tracks)~=length(obj))
     error('Added tracks set must have the same size.')
   end

   if(length(obj)==1)
   if(tracks.lifetime~=obj.lifetime)
     error('Added tracks must have the same lifetime')
   end
   obj.x=obj.x+tracks.x;
   obj.y=obj.y+tracks.y;
   obj.z=obj.z+tracks.z;
 else
   arrayfun(@(o,t) o.add(t),obj,tracks )
 end
end
