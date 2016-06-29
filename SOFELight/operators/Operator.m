classdef Operator < SOFEClass
  properties
    data, fes
    solution, stiffMat, loadVec, shift, freeDoFs
    solver = 'direct';
  end
  methods % constructor
    function obj = Operator(data, fes)
      obj.data = data;
      obj.fes = fes;
    end
  end
  methods % assemble & solve
    function assemble(obj)
      obj.stiffMat = [];
    end
    function compute(obj)
      t = tic; obj.outputBig('Begin assemble ...');
      obj.assemble();
      fprintf('%d DoFs\n', sum(obj.freeDoFs));
      obj.outputBig(['... assembled (',num2str(toc(t)),' sec)']);      
      t = tic; obj.outputBig('Begin solve ...');
      obj.solve(obj.solver);
      obj.outputBig(['... solved (',num2str(toc(t)),' sec)']);
    end
    function solve(obj, type)
      b = obj.loadVec - obj.stiffMat*obj.shift;
      obj.solution = zeros(size(b));
      obj.solution(~obj.freeDoFs) = obj.shift(~obj.freeDoFs);
      A = obj.stiffMat(obj.freeDoFs,obj.freeDoFs);
      b = b(obj.freeDoFs);
      obj.solution(obj.freeDoFs) = obj.shift(obj.freeDoFs) + Operator.solveLS(A, b, type);
    end    
  end
  methods(Static = true)
    function x = solveLS(A, b, type)
      tol = 1e-12; maxit = 1000;
      switch type
        case 'direct'
          x = A\b;
        case 'cg'
          D = spdiags(diag(A), 0, size(A,1), size(A,1));
          x = pcg(A,b,tol,maxit,D);
        case 'bicgstab'
          setup.type = 'nofill';
          [ML,MU] = ilu(A,setup);
          x = bicgstab(A,b,tol,maxit,ML,MU);
        case 'gmres'
          setup.type = 'nofill';
          [ML,MU] = ilu(A,setup);
          x = gmres(A,b,30,tol,maxit,ML,MU);
        otherwise
          x = zeros(size(b));
      end
    end
  end
end
