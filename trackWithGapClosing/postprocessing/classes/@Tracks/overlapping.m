function obj=overlapping(obj,tr)
 % Does not work with merge and split
 % Philippe Roudot 2017
  if(length(obj)==1)
    [F,idxTr,idxObj] = intersect(tr.f,obj.f);
    obj.startFrame=min(F);
    obj.endFrame=max(F);
    obj.x=obj.x(idxObj);
    obj.y=obj.y(idxObj);
    obj.z=obj.z(idxObj);
  else
    arrayfun(@(o,t) o.overlapping(t),obj,tracks );
  end
end
