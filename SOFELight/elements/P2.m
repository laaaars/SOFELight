classdef P2 < Element
  properties
  end
  methods % constructor
    function obj = P2(dim)
      obj = obj@Element(dim, 2:(dim+1), [], 2);
      obj.nB = cumsum(1:dim+1);
      obj.nB = obj.nB(2:dim+1);
      obj.doFTuple = zeros(1, dim+1);
      obj.doFTuple(1:2) = 1;
    end
  end
  methods % evaluation
    function B = getD0Basis(obj, points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1));
      switch nP
        case 1
          B(1,:) = 1-points;
          B(2,:) = points;
          B(3,:) = 4*points.*(1-points);
        case 2
          x = points(:,1); y = points(:,2);
          B(1,:) = 1-x-y;
          B(2,:) = x;
          B(3,:) = y;
          B(4,:) = 4*x.*(1-x-y);
          B(5,:) = 4*x.*y;
          B(6,:) = 4*y.*(1-x-y);
        case 3
          x = points(:,1); y = points(:,2); z = points(:,3);
          B(1,:) = 1-x-y-z;
          B(2,:) = x;
          B(3,:) = y;
          B(4,:) = z;
          B(5,:) = 4*x.*(1-x-y);
          B(6,:) = 4*x.*y;
          B(7,:) = 4*y.*(1-x-y);
          B(8,:) = 4*z.*(1-x-y);
          B(9,:) = 4*x.*z;
          B(10,:) = 4*y.*z;
      end
    end
    function B = getD1Basis(obj,points)
      nP = size(points, 2);
      B = zeros(obj.nB(nP), size(points,1), nP);
      switch nP
        case 1
          B(1,:,1) = -1;
          B(2,:,1) = 1;
          B(3,:,1) = 4-8*points;
        case 2
          x = points(:,1); y = points(:,2);
          B(1,:,:) = -1;
          B(2,:,1) = 1;
          B(3,:,2) = 1;
          B(4,:,1) = 4*(1-2*x-y);
          B(4,:,2) = -4*x;
          B(5,:,1) = 4*y;
          B(5,:,2) = 4*x;
          B(6,:,1) = -4*y;
          B(6,:,2) = 4*(1-x-2*y);
        case 3
          x = points(:,1); y = points(:,2); z = points(:,3);
          B(1,:,:) = -1;
          B(2,:,1) = 1;
          B(3,:,2) = 1;
          B(4,:,3) = 1;
          B(5,:,1) = 4*(1-2*x-y-z);
          B(5,:,2) = -4*x;
          B(5,:,3) = -4*x;
          B(6,:,1) = 4*y;
          B(6,:,2) = 4*x;
          B(7,:,1) = -4*y;
          B(7,:,2) = 4*(1-x-2*y-z);
          B(7,:,3) = -4*y;
          B(8,:,1) = -4*z;
          B(8,:,2) = -4*z;
          B(8,:,3) = 4*(1-x-y-2*z);
          B(9,:,1) = 4*z;
          B(9,:,3) = 4*x;
          B(10,:,2) = 4*z;
          B(10,:,3) = 4*y;
      end
    end
  end
end