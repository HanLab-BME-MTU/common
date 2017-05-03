function obj=multCoord(obj,scalar)
   % Philippe Roudot 2017

   if(length(obj)==1)
     obj.x=obj.x*scalar;
     obj.y=obj.y*scalar;
     obj.z=obj.z*scalar;
   else
     arrayfun(@(o,t) o.add(t),obj,tracks )
   end
 end
