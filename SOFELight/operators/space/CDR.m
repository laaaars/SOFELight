classdef CDR < Operator
  methods % constructor
    function obj = CDR(data, fes)
      obj = obj@Operator(data, fes);
    end
  end
  methods % assemble
    function assemble(obj)
      % lhs & rhs
      N = obj.fes.getNDoF();
      obj.stiffMat = sparse(N, N);
      obj.loadVec = zeros(N, 1);
      nBlock = obj.fes.mesh.nBlock;
      for k = 1:nBlock+1
        idx = obj.fes.mesh.getBlock(0, k);
        if ~isempty(idx)
          IDX = (idx(1):idx(2))';
          [points, weights] = obj.fes.getQuadData(0);
          trafo = obj.fes.mesh.evalTrafo(points, IDX);
          % LHS
          A = obj.fes.mesh.evalFunction(obj.data.a, points, IDX); % nExnP
          B = obj.fes.mesh.evalFunction(obj.data.b, points, IDX); % nExnPxnD
          C = obj.fes.mesh.evalFunction(obj.data.c, points, IDX); % nExnP
          gradBasis = obj.fes.evalGlobalBasis(points, 1, IDX); % nExnBxnPxnCxnD
          bGradBasis = sum(bsxfun(@times, gradBasis, permute(B, [1 4 2 5 3])), 5); % nExnBxnPxnC
          basis = obj.fes.evalGlobalBasis(points, 0, IDX); % nExnBxnPxnC
          dX = bsxfun(@times, trafo, weights'); % nExnP
          dXA = bsxfun(@times, A.*trafo, weights'); % nExnP
          dXC = bsxfun(@times, C.*trafo, weights'); % nExnP
          nE = size(gradBasis, 1); nB = size(gradBasis, 2); nP = size(A,2);
          gradBasis = reshape(gradBasis, nE, nB, nP, []);
          eLHS = sum(bsxfun(@times, sum(bsxfun(@times, permute(gradBasis, [1 2 5 3 4]), ...
                                                       permute(gradBasis, [1 5 2 3 4])), 5), ...
                                                       permute(dXA, [1 3 4 2])), 4); % nExnBxnB
          eLHS = eLHS + sum(bsxfun(@times, sum(bsxfun(@times, permute(basis, [1 2 5 3 4]), ...
                                                              permute(bGradBasis, [1 5 2 3 4])), 5), ...
                                                              permute(dX, [1 3 4 2])), 4); % nExnBxnB
          eLHS = eLHS + sum(bsxfun(@times, sum(bsxfun(@times, permute(basis, [1 2 5 3 4]), ...
                                                              permute(basis, [1 5 2 3 4])), 5), ...
                                                              permute(dXC, [1 3 4 2])), 4); % nExnBxnB
          % RHS
          F = obj.fes.mesh.evalFunction(obj.data.f, points, IDX); % nExnPxnC
          eRHS = sum(bsxfun(@times, sum(bsxfun(@times, permute(F, [1 4 2 3]), ...
                                                       permute(basis, [1 2 3 4])), 4), ...
                                                       permute(dX, [1 3 2])), 3); % nExnB
          % dofMaps
          rRHS = obj.fes.getDoFMap(0, IDX)';
          rLHS = repmat(rRHS,[1 1 nB]);
          cLHS = permute(rLHS, [1 3 2]);
          % assemble
          obj.stiffMat = obj.stiffMat + sparse(rLHS(:), cLHS(:), eLHS(:), N, N);
          obj.loadVec = obj.loadVec + accumarray(rRHS(:), eRHS(:), [N, 1]);
        end
        if nBlock > 1
          fprintf('progress assembly: %d / %d\r', k, nBlock+1);
        end
      end
      if nBlock > 1, fprintf('\n'); end
      % Neumann BC
      IDX = obj.fes.mesh.connectInfo.isBoundary(@(x)~obj.data.fix(x));
      if any(IDX)
        [points, weights] = obj.fes.getQuadData(1);
        eRHSBdy = obj.fes.mesh.evalFunction(obj.data.g, points, IDX); % nExnPxnC
        if ~isempty(points)
          basis = obj.fes.evalGlobalBasis(points, 0); % nExnBxnPxnC
          trafo = obj.fes.mesh.evalTrafo(points, IDX);
          dX = bsxfun(@times, trafo, weights'); % nExnP
          eRHSBdy = sum(bsxfun(@times, sum(bsxfun(@times, permute(eRHSBdy, [1 4 2 3]), ...
                                                          permute(basis, [1 2 3 4])), 4), ...
                                                          permute(dX, [1 3 2])), 3); % nExnB
        end
        % dofMaps
        rRHSBdy = obj.fes.getDoFMap(1, IDX)';
        obj.loadVec = obj.loadVec + accumarray(rRHSBdy(:), eRHSBdy(:), [N, 1]);
      end
      % shift
      obj.shift = zeros(N, 1);
      IDX = obj.fes.mesh.connectInfo.isBoundary(obj.data.fix);
      obj.shift = obj.fes.getInterpolation(obj.data.shift, 1, IDX);
      % freeDofs
      obj.freeDoFs = true(N, 1);
      bFix = obj.fes.mesh.connectInfo.isBoundary(obj.data.fix);
      if any(bFix)
        obj.freeDoFs = ~obj.fes.extractDoFs(1, bFix);
      end 
    end
  end
end
