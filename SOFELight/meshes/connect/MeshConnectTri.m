classdef MeshConnectTri < MeshConnect
  methods % constructor
    function obj = MeshConnectTri(nodes, elems)
      obj = obj@MeshConnect(nodes, elems);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectArray{1,1} = (1:size(obj.nodes,1))';
      obj.connectArray{2,1} = [obj.connectArray{3,1}(:,[1,2]); ...
                               obj.connectArray{3,1}(:,[2,3]); ...
                               obj.connectArray{3,1}(:,[1,3])];
      [obj.connectArray{2,1}, ~, e2F] = unique(sort(obj.connectArray{2,1},2),'rows');      
      obj.connectArray{3,2} = reshape(e2F, [], 3);
      for d = 1:2
        obj.connectArray{d+1,d+1} = (1:size(obj.connectArray{d+1,1},1))';
      end
    end
  end
  methods % connectivity information
    function R = getElem2Face(obj)
      R = obj.connectArray{3,2}; 
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(0);
      nF = obj.getNumber(1); nN = obj.getNumber(2);
      obj.nodes = [obj.nodes; obj.getCenter(1)];
      newIndices = nN + (1:nF);
      el = [el newIndices(obj.connectArray{3,2})];
      obj.connectArray{3,1} = [el(:,[1 4 6]);el(:,[4 2 5]);el(:,[6 5 3]);el(:,[5 6 4])];
      obj.updateConnectivity();
    end
  end
  methods % display
    function show(obj, varargin)
      nodes = obj.getEntity(2);
      trisurf(obj.getEntity(0), nodes(:,1), nodes(:,2), zeros(size(nodes,1),1));
      axis equal, view(0,90);
    end
  end
end