classdef QuadRule < SOFEClass
  properties
    order % exact for polynomials of ORDER
    points
    weights
  end
  methods % constructor
    function obj = QuadRule(n)
      obj.order = n;
      obj.initData();
    end
  end
  methods
    function initData(obj)
      obj.points = [];
    end
  end
  methods
    function [P, W] = getGaussPoints(obj)
      % Generates the points P and weights W for the Gauss-Legendre
      % quadrature of order ORDER
      n = ceil((obj.order + 1)/2); % number of points
      P = zeros(n,1);
      W= P;
      m = (n+1)/2;
      for ii=1:m
        z = cos(pi*(ii-.25)/(n+.5)); % Initial estimate.
        z1 = z+1;
        while abs(z-z1)>eps
          p1 = 1;
          p2 = 0;
          for jj = 1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
          end
          pp = n*(z*p1-p2)/(z^2-1); % The L.P. derivative.
          z1 = z;
          z = z1-p1/pp;
        end
        P(ii) = -z; % Build up the abscissas.
        P(n+1-ii) = z;
        W(ii) = 2/((1-z^2)*(pp^2)); % Build up the weights.
        W(n+1-ii) = W(ii);
      end
    end
  end
end

