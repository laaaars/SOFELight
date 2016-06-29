classdef MeshConnectInt < MeshConnect
  methods % constructor
    function obj = MeshConnectInt(nodes, elem)
      obj = obj@MeshConnect(nodes, elem);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectArray{1,1} = [obj.connectArray{2,1}(:,1); ...
                               obj.connectArray{2,1}(:,2)];
      obj.connectArray{1,1} = unique(obj.connectArray{1,1});
      obj.connectArray{2,2} = (1:size(obj.connectArray{2,1},1))';
    end
  end
  methods % connectivity information
    function R = getElem2Face(obj)
      R = obj.connectArray{2,1};
    end
  end
  methods % display
    function show(obj)
      plot(obj.nodes, zeros(obj.getNumber(1),1), '.');
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(0);
      nE = obj.getNumber(0); nN = obj.getNumber(1);
      obj.nodes = [obj.nodes; obj.getCenter(0)];
      newIndices = nN + (1:nE);
      el = [el newIndices'];
      obj.connectArray{2,1} = [el(:,[1 3]); el(:,[3 2])];
      obj.updateConnectivity();
    end
  end
end