classdef Mesh < SOFEClass
  properties
    element
    connectInfo
    nBlock = 1;
  end
  methods % constructor & more
    function obj = Mesh(nodes, elems)  
      dim =  size(nodes, 2);   
      obj.element = P1(dim);
      obj.connectInfo = MeshConnect.getConnectInfo(nodes, elems);
    end
    function R = getBlock(obj, codim, i)
      nE = obj.connectInfo.getNumber(codim);
      if obj.nBlock>nE, error('!Number of blocks exceeds number of elements!'); end
      nEPart = floor(nE/obj.nBlock);
      nERem = nE-nEPart*obj.nBlock;
      if i<=obj.nBlock
        R = (i-1)*nEPart + [1; nEPart];
      elseif i == obj.nBlock+1
        R = obj.nBlock*nEPart + [1; nERem];
        if R(1)>R(2), R = []; end
      else
        error('!Block index must be <= number of blocks+1!');
      end
    end
  end
  methods % reference map
    function R = evalReferenceMap(obj, points, order, varargin)     
      I = ':'; if nargin > 3, I = varargin{1}; end
      nP = size(points,1); nD = obj.connectInfo.dimension;
      codim = nD - size(points,2);
      nodes = obj.connectInfo.getEntity(nD);
      if isempty(points) % for 1D special case
        R = permute(nodes(I,:), [1 3 2]);
        return
      end
      basis = obj.element.getBasis(points, order);
      nB = size(basis, 1);
      connect = obj.connectInfo.getEntity(codim, I); % linear
      nE = size(connect,1);
      R = reshape(nodes(connect(:),:)', nD*nE, nB)* ...
          reshape(basis,nB,[]); % (nD*nE)x(nP*[...])
      R = permute(reshape(R, nD, nE, nP, []), [2 3 1 4]); % nExnPxnDx[...];
    end
    function [det, R, RInv] = evalTrafoTriple(obj, points, varargin)
      R = obj.evalReferenceMap(points, 1, varargin{:});
      switch size(R,3)
        case 1
          det = R;
        case 2
          det = R(:,:,1,1).*R(:,:,2,2) - R(:,:,1,2).*R(:,:,2,1);
        case 3
          det = R(:,:,1,1).*R(:,:,2,2).*R(:,:,3,3) + ...
                R(:,:,1,2).*R(:,:,2,3).*R(:,:,3,1) + ...
                R(:,:,1,3).*R(:,:,2,1).*R(:,:,3,2) - ...
                R(:,:,1,1).*R(:,:,2,3).*R(:,:,3,2) - ...
                R(:,:,1,2).*R(:,:,2,1).*R(:,:,3,3) - ...
                R(:,:,1,3).*R(:,:,2,2).*R(:,:,3,1);
      end
      if nargout > 2
        switch size(R,3)
          case 1
            RInv = 1;
          case 2
            RInv = -R;
            RInv(:,:,1,1) = R(:,:,2,2);
            RInv(:,:,2,2) = R(:,:,1,1);
          case 3
            RInv(:,:,1,1) = R(:,:,2,2).*R(:,:,3,3) - R(:,:,2,3).*R(:,:,3,2);
            RInv(:,:,2,1) = -(R(:,:,2,1).*R(:,:,3,3) - R(:,:,3,1).*R(:,:,2,3));
            RInv(:,:,3,1) = R(:,:,2,1).*R(:,:,3,2) - R(:,:,2,2).*R(:,:,3,1);
            RInv(:,:,1,2) = -(R(:,:,1,2).*R(:,:,3,3) - R(:,:,1,3).*R(:,:,3,2));
            RInv(:,:,2,2) = R(:,:,1,1).*R(:,:,3,3) - R(:,:,1,3).*R(:,:,3,1);
            RInv(:,:,3,2) = -(R(:,:,1,1).*R(:,:,3,2) - R(:,:,1,2).*R(:,:,3,1));
            RInv(:,:,1,3) = R(:,:,1,2).*R(:,:,2,3) - R(:,:,1,3).*R(:,:,2,2);
            RInv(:,:,2,3) = -(R(:,:,1,1).*R(:,:,2,3) - R(:,:,1,3).*R(:,:,2,1));
            RInv(:,:,3,3) = R(:,:,1,1).*R(:,:,2,2) - R(:,:,1,2).*R(:,:,2,1);
        end
        RInv = bsxfun(@rdivide, RInv, det);
      end
    end
    function R = evalTrafo(obj, points, varargin)
      switch obj.connectInfo.dimension
        case 1
          R = abs(obj.evalTrafoTriple(points, varargin{:}));
        case 2
          switch size(points,2)
            case 1
              R = obj.evalReferenceMap(points, 1, varargin{:});
              R = sqrt(sum(R.^2,3));
            case 2
              R = abs(obj.evalTrafoTriple(points, varargin{:}));
          end
        case 3
          switch size(points,2)
            case 2
              R = obj.evalReferenceMap(points, 1, varargin{:});
              R = sqrt((R(:,:,2,1).*R(:,:,3,2) - R(:,:,3,1).*R(:,:,2,2)).^2 + ...
                       (R(:,:,3,1).*R(:,:,1,2) - R(:,:,1,1).*R(:,:,3,2)).^2 + ...
                       (R(:,:,1,1).*R(:,:,2,2) - R(:,:,2,1).*R(:,:,1,2)).^2);   
            case 3
              R = abs(obj.evalTrafoTriple(points, varargin{:}));
          end
      end
    end
  end
  methods % data eval
    function R = evalFunction(obj, F, points, varargin)
      P = obj.evalReferenceMap(points, 0, varargin{:}); % nExnPxnD
      nE = size(P,1); nP = size(P,2);
      P = reshape(P, nE*nP, []); % (nE*nP)xnD
      R = [];
      if ~isempty(P)
        R = reshape(F(P), nE, nP, []);
      end
    end
  end
  methods % refinement
    function uniformRefine(obj, N)
      for i = 1:N
        obj.connectInfo.uniformRefine();
      end
    end
  end
  methods % display
    function show(obj)
      obj.connectInfo.show();
    end
  end
end