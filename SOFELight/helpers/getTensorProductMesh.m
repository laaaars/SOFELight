function [nodes, elem] = getTensorProductMesh(grid, varargin)
  switch numel(grid)
    case 1
      nodes = grid{1};
      if size(nodes,1) == 1
        nodes = nodes';
      end
      m = numel(nodes);
      elem = [1:(m-1); 2:m]';
    case 2
      m = numel(grid{1});
      n = numel(grid{2});
      nodes = [kron(ones(1,n),grid{1}); ...
               kron(grid{2}, ones(1,m))]';
      elem = [1:m*(n-1)-1; 2:m*(n-1)];
      elem = [elem; elem+m]';
      elem(m:m:end,:) = [];
      if nargin > 1 && varargin{1}
        elem = [elem(:,[1 2 3]); elem(:,[4 3 2])];
      end
    case 3
      m = numel(grid{1});
      n = numel(grid{2});
      p = numel(grid{3});
      nodes = [kron(ones(1,p),kron(ones(1,n),grid{1}))', ...
               kron(ones(1,p),kron(grid{2}, ones(1,m)))', ...
               kron(grid{3},kron(ones(1,n), ones(1,m)))'];
      elem = [1:m*n*(p-1)-m-1; 2:m*n*(p-1)-m]';
      elem = [elem elem+m];
      elem = [elem, elem+m*n];
      delete = bsxfun(@(x,y)x+y,m*(n-1)+1:m*n-1,((0:(p-3))*m*n)');
      elem([m:m:m*n*(p-1)-m-1 delete(:)'],:) = [];
      if nargin > 1 && varargin{1}
        elem = [elem(:,[5 6 1 8]); elem(:,[5 1 7 8]); elem(:,[2 1 6 8]); ...
                elem(:,[3 7 1 8]); elem(:,[4 2 8 1]); elem(:,[3 1 4 8])];
      end
  end
end