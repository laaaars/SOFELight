classdef Mass < Operator
  methods % constructor
    function obj = Mass(fes)
      obj = obj@Operator([], fes);
    end
  end
  methods % assemble & solve
    function assemble(obj)
      % lhs & rhs
      N = obj.fes.getNDoF();
      obj.stiffMat = sparse(N, N);
      nBlock = obj.fes.mesh.nBlock;
      for k = 1:nBlock+1
        idx = obj.fes.mesh.getBlock(0, k);
        if ~isempty(idx)
          IDX = [idx(1):idx(2)]';
          [points, weights] = obj.fes.getQuadData(0);
          trafo = obj.fes.mesh.evalTrafo(points, IDX);
          % LHS
          basis = obj.fes.evalGlobalBasis(points, 0, IDX); % nExnBxnPxnC
          dX = bsxfun(@times, trafo, weights'); % nExnP
          eLHS = sum(bsxfun(@times, sum(bsxfun(@times, permute(basis, [1 2 5 3 4]), ...
                                                       permute(basis, [1 5 2 3 4])), 5), ...
                                                       permute(dX, [1 3 4 2])), 4); % nExnBxnB
          % dofMaps
          rLHS = obj.fes.getDoFMap(0, IDX)';
          rLHS = repmat(rLHS,[1 1 size(rLHS,2)]);
          c = permute(rLHS, [1 3 2]);
        else
          rLHS = []; c = []; eLHS = [];
        end
        %
        obj.stiffMat = obj.stiffMat + sparse(rLHS(:), c(:), eLHS(:), N, N);
        if nBlock > 1
          fprintf('progress assembly: %d / %d\r', k, nBlock+1);
        end
      end
      if nBlock > 1, fprintf('\n'); end
    end
  end
end
