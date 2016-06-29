classdef GaussPoint < QuadRule
  properties
  end
  methods % constructor
    function obj = GaussPoint()
      obj = obj@QuadRule(0);
    end
  end
   methods % initialize
    function initData(obj)
      obj.points = [];
      obj.weights = 1;
    end
   end
end