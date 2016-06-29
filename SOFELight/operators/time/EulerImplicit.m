classdef EulerImplicit < SOFEClass
  properties
    massOp
    statOp
    initCond
    T, dt
    solution
  end
  methods % constructor
    function obj = EulerImplicit(data, statOp)
      obj.statOp = statOp;
      obj.T = data.T;
      obj.dt = data.dt;
      obj.initCond = data.initCond;
      obj.massOp = Mass(obj.statOp.fes);
    end
  end
  methods
    function compute(obj)
      N = ceil(obj.T/obj.dt);
      obj.solution = cell(N,1);
      % init
      obj.solution{1} = obj.statOp.fes.getInterpolation(obj.initCond, 0);
      obj.massOp.assemble();
      obj.statOp.assemble();
      for k = 1:N
        S = obj.massOp.stiffMat + obj.dt*obj.statOp.stiffMat;
        b = obj.dt*obj.statOp.loadVec + obj.massOp.stiffMat*obj.solution{k} - S*obj.statOp.shift;
        I = obj.statOp.freeDoFs;
        % solve
        obj.solution{k+1} = zeros(size(obj.solution{k}));
        obj.solution{k+1}(~I) = obj.statOp.shift(~I);
        obj.solution{k+1}(I) = obj.statOp.shift(I) + ...
                               Operator.solveLS(S(I, I), b(I), obj.statOp.solver);
        fprintf('timestep: %f / %d\n', k, N);
      end
    end
  end
end