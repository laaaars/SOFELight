for dim = 1:3
  figure(dim)
  % MESH
  m = UniformMesh(dim, repmat([0 1],dim,1), 30);
  m.nBlock = 20;
  % FESPACE
  fes = FESpace(m, P1(m.connectInfo.dimension));
  % DATA
  data.a = @(x) 0.05+0*x(:,1);
  data.b = @(x) x;
  data.c = @(x) 1+0*x(:,1);
  data.f = @(x) 1+0*x(:,1);
  data.g = @(x) -0.1+0*x(:,1);
  data.fix = @(x) x(:,1) < 1.0;
  data.shift = @(x) x(:,1).*(1-x(:,1));
  % PDE
  p = CDR(data, fes); p.solver = 'bicgstab';
  p.compute();
  % VISUALIZE
  v = Visualizer(fes); clf
  v.showDoFVector(p.solution, 1, @(x)x(:,3)-x(:,1).^2+x(:,2).^2<0.5);
  drawnow
end
