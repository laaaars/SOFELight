classdef Visualizer < SOFEClass
  properties
    fes
  end
  methods % constructor
    function obj = Visualizer(fes)
      obj.fes = fes;
    end
  end
  methods % plot
    function showDoFVector(obj, U, N, varargin)
      if numel(U)~=obj.fes.getNDoF
        error('w is no DoFVector!');
      end
      switch obj.fes.element.dimension
        case 1
          points = linspace(0,1,N+1)';
          [P,I] = unique(obj.fes.mesh.evalReferenceMap(points, 0));
          Z = obj.fes.evalDoFVector(U, points, 0);
          plot(P(:), Z(I));
        case 2
          nE = obj.fes.mesh.connectInfo.getNumber(0);
          nV = obj.fes.element.nV(obj.fes.element.dimension);
          [Y,X] = meshgrid(linspace(0,1,N+1), linspace(0,1,N+1));
          if nV == 3
            idx = X(:)+Y(:)<=1;
            X = X(idx); Y = Y(idx);
          end
          value = reshape(obj.fes.evalDoFVector(U,[X(:) Y(:)],0)', [], 1); % (nVxnE)
          vertices = obj.fes.mesh.evalReferenceMap([X(:) Y(:)],0); % nExnVx2
          vertices = reshape(permute(vertices, [2 1 3]), [], 2); % (nV*nE)x2
          if nV == 4
            col1 = reshape(1:N*(N+1),N+1,N);
            col1(N+1,:) = []; col1 = col1(:);
            elem = [col1 col1+1 col1+N+1 col1+N+2]; %nESubxnVx1
            elem = elem(:,[1 2 4 3]);
            offset = permute((0:nE-1)*(N+1)^2, [1 3 2]); % 1x1xnE
          else
            col = reshapeTop(1:(N+2)*(N+1)/2,N+1:-1:1);
            col1 = col; col1(N+1:N:N^2+N+1) = 0; col1(col1 == 0) = [];
            col2 = col; col2(:,1) = []; col2(col2 == 0) = [];
            col3 = col; col3(1,:) = []; col3(:,1) = []; col3(col3==0) = [];
            col4 = col; col4(N+1:N:N^2+N+1) = 0; col4(1,:) = 0; col4(col4==0) = [];
            elem = [col1 col1+1 col2; col3 col3-1 col4]; %nESubxnVx1
            offset = permute((0:nE-1)*(N+1)*(N+2)/2, [1 3 2]); % 1x1xnE
          end
          elem = permute(bsxfun(@plus, elem, offset), [1 3 2]); % nESubxnExnV
          patch('faces', reshape(elem, [], nV), 'vertices', [vertices value], ... 
                'facevertexcdata',value,'facecolor','interp', ...
                'edgecolor','interp');
        case 3
          I = obj.fes.mesh.connectInfo.isSurface(varargin{:});
          nE = sum(I); nV = obj.fes.element.nV(obj.fes.element.dimension-1);
          [Y,X] = meshgrid(linspace(0,1,N+1), linspace(0,1,N+1));
          if nV == 3
            idx = X(:)+Y(:)<=1;
            X = X(idx); Y = Y(idx);
          end
          value = reshape(obj.fes.evalDoFVector(U,[X(:) Y(:)], 0, I)', [], 1); % (nVxnE)
          vertices = obj.fes.mesh.evalReferenceMap([X(:) Y(:)], 0, I); % nExnVx3
          vertices = reshape(permute(vertices, [2 1 3]), [], 3); % (nV*nE)x3
          if nV == 4
            col1 = (1:N*(N+1))'; col1((N+1)*(1:N)) = [];
            faces = [col1 col1+1 col1+N+1 col1+N+2]; %nESubxnVx1
            faces = faces(:, [1 2 4 3]);
            offset = permute((0:nE-1)*(N+1)^2, [1 3 2]); % 1x1xnE
          else
            col = reshapeTop(1:(N+2)*(N+1)/2,N+1:-1:1);
            col1 = col; col1(N+1:N:N^2+N+1) = 0; col1(col1 == 0) = [];
            col2 = col; col2(:,1) = []; col2(col2 == 0) = [];
            col3 = col; col3(1,:) = []; col3(:,1) = []; col3(col3==0) = [];
            col4 = col; col4(N+1:N:N^2+N+1) = 0; col4(1,:) = 0; col4(col4==0) = [];
            col1 = reshape(col1,[],1);
            col2 = reshape(col2,[],1);
            col3 = reshape(col3,[],1);
            col4 = reshape(col4,[],1);
            faces = [col1 col1+1 col2; col3 col3-1 col4]; %nESubxnVx1
            offset = permute((0:nE-1)*(N+1)*(N+2)/2, [1 3 2]); % 1x1xnE
          end
          faces = permute(bsxfun(@plus, faces, offset), [1 3 2]); % nESubxnExnV
          patch('faces', reshape(faces, [], nV), 'vertices', vertices, ... 
                'facevertexcdata',value,'facecolor','interp', ...
                'edgecolor','interp');
          axis tight, axis equal
      end
    end
  end
end