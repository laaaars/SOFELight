classdef Element < SOFEClass
  properties
    dimension,
    nV, nB
    order
    doFTuple
    conformity
  end
  methods % constructor
    function obj = Element(dimension, nV, nB, order)
      obj.dimension = dimension;
      obj.nV = nV;
      obj.nB = nB;
      obj.order = order;
      obj.conformity = 'H1';
    end
  end
  methods % evaluation local
    function B = getBasis(obj, points, order)
      switch order
        case 0
          B = obj.getD0Basis(points); % nBxnPxnC
        case 1
          B = obj.getD1Basis(points); % nBxnPxnCxnD
          if ndims(B)<4, B = permute(B, [1 2 4 3]); end
      end
    end
  end
end

