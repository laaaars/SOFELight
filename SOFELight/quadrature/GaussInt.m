classdef GaussInt < QuadRule
  properties
  end
  methods % constructor
    function obj = GaussInt(n)
      obj = obj@QuadRule(n);
    end
  end
   methods % initialize
    function initData(obj)
      [obj.points, obj.weights] = obj.getGaussPoints();
      obj.points = (1 + obj.points)/2;
      obj.weights = obj.weights/2;
    end
   end
end