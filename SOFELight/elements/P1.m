classdef P1 < Element
  methods % constructor
    function obj = P1(dim)
      obj = obj@Element(dim, 2:(dim+1), 2:(dim+1), 1);
      obj.doFTuple = zeros(1,dim+1);
      obj.doFTuple(1) = 1;
    end
  end
  methods % evaluation
    function B = getD0Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1));
      B(1,:) = 1-sum(points,2);
      for i = 1:nP
        B(i+1,:) = points(:,i);
      end
    end
    function B = getD1Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP);
      B(1,:,:) = -1;
      for i = 1:nP
        B(i+1,:,i) = 1;
      end
    end
  end
end