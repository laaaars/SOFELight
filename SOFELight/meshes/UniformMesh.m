classdef UniformMesh < Mesh
  methods % constructor
    function obj = UniformMesh(dimension, range, N)
      fprintf('Generate mesh ... ');
      if numel(N)==1
        N = repmat(N,dimension,1);
      end
      grids = cell(1,dimension);
      for d = 1:dimension
        grids{d} = linspace(range(d,1), range(d,2), N(d)+1);
      end
      [nodes, elems] = getTensorProductMesh(grids, 1);
      obj = obj@Mesh(nodes, elems);
      fprintf('DONE\n');
    end
  end
end