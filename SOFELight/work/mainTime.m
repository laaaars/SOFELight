% MESH
m = UniformMesh(2, [0 1; 0 1], 30);
m.nBlock = 1;
% FESPACES
fes = FESpace(m, P2(m.connectInfo.dimension));
% DATA
data.a = @(x) 0.1 + 0*x(:,1);
data.b = @(x) x;
data.c = @(x) 1+0*x(:,1);
data.f = @(x) 0+0*x(:,1);
data.fix = @(x) x(:,1) < 0.5;
data.shift = @(x) 0*x(:,1);
data.g = @(x) 0*x(:,1);
data.T = 2; data.dt = 0.05;
data.initCond = @(x)exp(-sum((x-0.25).^2,2)/0.001);
% PROBLEM
p = CDR(data, fes);
p = EulerImplicit(data, p);
p.compute();
% VISUALIZE
v = Visualizer(fes);
N = numel(p.solution);
for k = 1:1:N
  clf
  v.showDoFVector(p.solution{k}, 2);
  fprintf('time: %f / %d\n', (k-1)*data.dt, data.T);
  pause(0.01);
end