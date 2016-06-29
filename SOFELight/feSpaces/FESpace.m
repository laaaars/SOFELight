classdef FESpace < SOFEClass
  properties
    mesh
    element
    quadRule
  end
  methods % constructor & get
    function obj = FESpace(mesh, element)
      obj.element = element;
      obj.mesh = mesh;
      obj.quadRule = MeshConnect.getQuadRule(element.dimension, max(2*(element.order),1));
    end
    function [Rp, Rw] = getQuadData(obj, codim)
      Rp = obj.quadRule{codim+1}.points;
      Rw = obj.quadRule{codim+1}.weights;
    end
  end
  methods % evaluation of global basis
    function R = evalGlobalBasis(obj, points, dOrder, varargin) % [I]
      B = obj.element.getBasis(points, dOrder);
      switch dOrder
        case 0
          R = permute(B, [4 1 2 3]); % 1xnBxnPxnC
        case 1
          [~,~,jacInv] = obj.mesh.evalTrafoTriple(points, varargin{:});
          jacInv = permute(jacInv, [1 5 2 3 4]); % nEx1xnPxnDxnD
          dBasis = permute(B, [5 1 2 3 4]); % 1xnBxnPxnCxnD
          R = sum(bsxfun(@times, permute(dBasis, [1 2 3 4 6 5]), ...
                                 permute(jacInv, [1 2 3 6 5 4])), 6); % nExnBxnPxnCxnD
      end
    end
  end
  methods % evaluation of DoFVectors
    function R = evalDoFVector(obj, U, points, order, varargin)
      codim = obj.element.dimension - size(points,2);
      if (codim>0 && order>0)
        error('! Higher order derivatives for traces not supported !');
      end
      basis = obj.evalGlobalBasis(points, order, varargin{:}); % [1/nE]xnBxnPxnCx[nD]
      dMap = obj.getDoFMap(codim, varargin{:});
      value = zeros(size(dMap)); % nBxnE
      I = dMap ~= 0; value(I) = U(dMap(I));
      R = sum(bsxfun(@times, value', basis), 2); % nEx1xnPxnCx[nD]
      R = permute(R, [1 3:(4+order) 2]); % nExnPxnCx[nD]
    end
  end
  methods % dofManager
    function [R, nDoF] = getDoFMap(obj, codim, varargin)
      nD = obj.element.dimension;
      doFs = cell(nD+1,1); R = cell(nD+1,1);
      nDoF = 0;
      for dim = 0:nD
        N = obj.mesh.connectInfo.getNumber(nD-dim);
        currNDoF = obj.element.doFTuple(dim+1)*N;
        doFs{dim+1} = reshape(nDoF+(1:currNDoF), [], N); % nBxnE
        nDoF = nDoF + currNDoF;
      end
      for dim = 0:nD
        R{nD-dim+1} = cell(dim+1,1);
        for d = 0:dim
          cArray = obj.mesh.connectInfo.connectArray{dim+1, d+1};
          R{nD-dim+1}{d+1} = reshape(doFs{d+1}(:, cArray'), [], size(cArray,1));
        end
        R{nD-dim+1} = cell2mat(R{nD-dim+1});
      end
      I = ':'; if nargin > 2, I = varargin{1}; end
      R = R{codim+1}(:,I);
    end
    function nDoF = getNDoF(obj)
      [~, nDoF] = obj.getDoFMap(0);
    end
    function R = extractDoFs(obj, codim, I)
      [dMap, nDoF] = obj.getDoFMap(codim, I);
      R = accumarray(dMap(:), 1, [nDoF, 1])>0;
    end
  end
  methods % interpolation
    function R = getInterpolation(obj, f, codim, varargin) % [I]
      if codim == obj.element.dimension
        P = []; nP = 1; S.lhs = 1;
      else
        points = (0:obj.element.order)'/obj.element.order;
        P = points;
        for i = 2:obj.element.dimension-codim
          P = [kron(ones(length(points),1),P), kron(points,ones(length(points)^(i-1),1))];
        end
        P = P((sum(P,2)<=1),:); nP = size(P,1);
        S.lhs = reshape(obj.element.getBasis(P, 0), [], nP)'; % nPxnB, nB==nP
      end
      S.rhs = reshape(obj.mesh.evalFunction(f, P, varargin{:}), [], nP)'; % nPxnE
      S = S.lhs\S.rhs; % nBxnE
      R = zeros(obj.getNDoF,1);
      R(obj.getDoFMap(codim, varargin{:})) = S;
    end
  end
end