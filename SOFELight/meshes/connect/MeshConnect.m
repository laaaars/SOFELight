classdef MeshConnect < handle
  properties
    dimension
    nodes
    connectArray
  end
  methods % constructor
    function obj = MeshConnect(nodes, elem)
      obj.nodes = nodes;
      obj.dimension = size(nodes,2);
      obj.connectArray = cell(obj.dimension+1);
      obj.connectArray{obj.dimension+1,1} = elem;
    end
  end
  methods % mesh information
    function R = getEntity(obj, codim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      if codim == obj.dimension
        R = obj.nodes(I,:);
      else
        R = obj.connectArray{obj.dimension-codim+1}(I,:);
      end
    end
    function R = getNumber(obj, codim)
      R = size(obj.getEntity(codim), 1);
    end
    function R = getCenter(obj, codim, varargin)
      I = ':'; if nargin > 2, I = varargin{1}; end
      if codim == obj.dimension
        vertices = obj.getEntity(obj.dimension);
        R = vertices(I,:);
        return;
      end
      R = obj.getEntity(codim); R = R(I,:);
      [nE, nV] = size(R);
      vertices = obj.getEntity(obj.dimension);
      R = permute(sum(reshape(vertices(R,:), nE, nV, []),2)/nV,[1 3 2]);
    end
    function R = isBoundary(obj, varargin)
      E2F = obj.getElem2Face();
      R = accumarray(E2F(:),1)==1;
      if nargin > 1
        if ~isempty(varargin{1})
          R = find(R);
          I = varargin{1}(obj.getCenter(1, R));
          R = R(I);
          R = accumarray(R,1,[obj.getNumber(1) 1])>0;
        else
          R = [];
        end
      end
    end
    function R = isSurface(obj, varargin)
      if nargin > 1
        idx = find(varargin{:}(obj.getEntity(obj.dimension)));
      else
        idx = 1:obj.getNumber(obj.dimension);
      end
      E2F = obj.getElem2Face();
      E2F = E2F(any(ismember(obj.getEntity(0), idx),2),:);
      R = accumarray(E2F(:), 1, [obj.getNumber(1) 1])==1;
    end
  end
  methods(Static = true)
    function R = getConnectInfo(nodes, elems)
      switch size(nodes, 2)
        case 1
          R = MeshConnectInt(nodes, elems);
        case 2
          R = MeshConnectTri(nodes, elems);
        case 3
          R = MeshConnectTet(nodes, elems);
      end
    end
    function R = getQuadRule(dim, order)
      switch dim
        case 1
          R{2} = GaussPoint();
          R{1} = GaussInt(order);
        case 2
          R{3} = GaussPoint();
          R{2} = GaussInt(order);
          R{1} = GaussTri(order);
        case 3
          R{4} = GaussPoint();
          R{3} = GaussInt(order);
          R{2} = GaussTri(order);
          R{1} = GaussTet(order);
      end
    end
  end
end