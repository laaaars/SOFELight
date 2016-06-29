classdef GaussTri < QuadRule
  properties
  end
  methods
    function obj = GaussTri(n)
      obj = obj@QuadRule(n);
    end
    function initData(obj)
      [p,w] = obj.getGaussPoints(); 
      p = 0.5 + 0.5*p;
      w = w/2;
      wArray = w*w';
      W = wArray(:);
      P = [kron(ones(length(w),1),p) kron(p, ones(length(w),1))];
      % transform to triangle
      obj.weights = W.*(1-P(:,2));
      obj.points(:,1) = P(:,1).*(1-P(:,2));
      obj.points(:,2) = P(:,2);
    end
  end
end