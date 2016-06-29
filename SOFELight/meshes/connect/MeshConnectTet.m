classdef MeshConnectTet < MeshConnect
  methods % constructor
    function obj = MeshConnectTet(nodes, elems)
      obj = obj@MeshConnect(nodes, elems);
      obj.updateConnectivity();
    end
    function updateConnectivity(obj)
      obj.connectArray{1,1} = (1:size(obj.nodes,1))';
      obj.connectArray{3,1} = [obj.connectArray{4,1}(:,[1,2,3]); ...
                               obj.connectArray{4,1}(:,[1,2,4]); ...
                               obj.connectArray{4,1}(:,[2,3,4]); ...
                               obj.connectArray{4,1}(:,[1 3 4])];
      [obj.connectArray{3,1}, ~, e2F] = unique(sort(obj.connectArray{3,1},2),'rows');    
      obj.connectArray{2,1} = [obj.connectArray{4,1}(:,[1,2]); obj.connectArray{4,1}(:,[2,3]); ...
                               obj.connectArray{4,1}(:,[1,3]); obj.connectArray{4,1}(:,[1 4]); ...
                               obj.connectArray{4,1}(:,[2,4]); obj.connectArray{4,1}(:,[3 4])];
      [obj.connectArray{2,1}, ~, e2E] = unique(sort(obj.connectArray{2,1},2),'rows');    
      obj.connectArray{4,3} = reshape(e2F,[], 4);
      obj.connectArray{4,2} = reshape(e2E,[], 6);
      obj.connectArray{3,2} = obj.getFace2Edge();
      for d = 1:obj.dimension
        obj.connectArray{d+1,d+1} = (1:size(obj.connectArray{d+1,1},1))';
      end
    end
  end
  methods % connectivity information
    function R = getElem2Face(obj)
      R = obj.connectArray{4,3};
    end
    function R = getElem2Edge(obj)
      R = obj.connectArray{4,2};
    end
    function R = getFace2Edge(obj)
      R = zeros(obj.getNumber(1), 3);
      e2F = obj.getElem2Face();
      e2E = obj.getElem2Edge();
      orientF = obj.getOrientation(3,2);
      orient1 = permute(orientF(1,:,:),[2 3 1]);
      orient2 = permute(orientF(2,:,:),[2 3 1]);
      [~, ind] = unique(e2F);
      orient1 = orient1(ind); orient2 = orient2(ind);
      [elems, type] = ind2sub([obj.getNumber(0), 4], ind);
      nodeIxAtFace = [1 2 3; 1 5 4; 2 6 5; 3 6 4];
      for t = 1:4
        for k = 0:2
          I = (type == t) & (orient1 == k+1);
          R(I,:) = e2E(elems(I), circshift(nodeIxAtFace(t,:)',-k));
        end
      end
      I = (orient2 == -1);
      R(I,:) = R(I, [3 2 1]);
    end
    function R = getOrientation(obj, dim1, dim2)
      switch dim2
        case 1
          switch dim1
            case 3
              e = obj.getEntity(0);
              R = ones(size(e,1), 6);
              R(e(:,1)>e(:,2),1) = -1;
              R(e(:,2)>e(:,3),2) = -1;
              R(e(:,1)>e(:,3),3) = -1;
              R(e(:,1)>e(:,4),4) = -1;
              R(e(:,2)>e(:,4),5) = -1;
              R(e(:,3)>e(:,4),6) = -1;
            case 2
              R = ones(obj.getNumber(1), 3);
          end
        case 2
          e = obj.getEntity(0);
          R = ones(2, size(e,1), 4); % ! two orient flags
          face = [1 2 3; 1 2 4; 2 3 4; 1 3 4];
          for i = 1:4
              [~,R(1,:,i)] = min(e(:,face(i,:)),[],2);
              [~,P] = sort(e(:,face(i,:)),2);
              even = (P(:,1)<P(:,2) & P(:,2)<P(:,3)) | ...
                     (P(:,2)<P(:,3) & P(:,3)<P(:,1)) | ...
                     (P(:,3)<P(:,1) & P(:,1)<P(:,2));
              R(2,~even,i) = -1;
          end
      end
    end
  end
  methods % refinement
    function uniformRefine(obj)
      el = obj.getEntity(0); nN = obj.getNumber(3);
      obj.nodes = [obj.nodes; obj.getCenter(2)];
      newNodeNumber = (nN+1 : nN+obj.getNumber(2));
      el = [el newNodeNumber(obj.connectArray{4,2})];
      obj.connectArray{4,1} = [el(:,[1 5 7 8]); el(:,[2 6 5 9]); ...
                               el(:,[3 7 6 10]); el(:,[4 9 8 10]); ...
                               el(:,[5 6 7 9]); el(:,[5 7 8 9]); ...
                               el(:,[7 8 9 10]); el(:,[7 9 6 10])];
      obj.updateConnectivity();
    end
  end
  methods % display
    function show(obj, varargin)
      expr = [];
      if nargin > 1, expr = varargin{1}; end
      simpplot(obj.getEntity(3), obj.getEntity(0), expr);
      axis equal
    end
  end
end