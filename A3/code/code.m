
clear all
clc
close all
% mkoctfile("./phigradient.cc")

addpath("../poisson");

scale = -0.5;

L1 = 60 + 84*scale;
L2 = 60;
epsilon = 1;
epsilon0 = 1e-3;
sigma = 1;
r_cut = 6;
m = 23;
noise = 0.1;

dt = 0.0001;
T = 10;
steps = T/dt + 1;

salvaOgni = 1;

H1 = 19;
H2 = 38;
H3 = 19;
W = 180*L1/144;

grd.ncy = L2/r_cut + 2;
grd.ncx = L1/r_cut + 2;
grd.x = linspace (0, L1 + 2*r_cut, grd.ncx+1);
grd.y = linspace (0, L2 + 2*r_cut, grd.ncy+1);

h = 1;
msh.ncx = L1 / h;
msh.ncy = L2 / h;
msh.nx = msh.ncx + 1;
msh.ny = msh.ncy + 1;
msh.x  = linspace (r_cut, L1 + r_cut, msh.nx);
msh.y  = linspace (r_cut, L2 + r_cut, msh.ny);
[X, Y] = meshgrid (msh.x, msh.y);
msh.hx = diff (msh.x);
msh.hy = diff (msh.y);

delta = 0.4;

[X1, Y1] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(delta,L2/4 - delta, H1) + r_cut);
[X2, Y2] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(L2/4 + delta, 3*L2/4 - delta, H2) + r_cut);
[X3, Y3] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(3*L2/4 + delta ,L2 - delta, H3) + r_cut);

N1 = numel(X1);
N2 = numel(X2);
N3 = numel(X3);
Nparticelle = N1 + N2 + N3;

ptcls.x = [[X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)]]';
ptcls.q(1:N1) = -0.5.*ones(1,N1);
ptcls.q((N1+1):(N1+N2)) = 0.5.*ones(1,N2);
ptcls.q((N1+N2+1):(N1+N2+N3)) = -0.5.*ones(1,N3);
randn('seed', 100);
ptcls.v = randn (size (ptcls.x)) * noise;

grd_to_ptcl = init_ptcl_mesh (grd, ptcls);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);

save ("-append",strcat("../posizioni/posizioni_iniziali"),"ptcls")

ptcls_gc.x = [];
ptcls_gc.v = [];
ptcls_gc.q = [];

sx = [2 + grd.ncy : grd.ncy * 2 - 1];
bot = [2 * grd.ncy - 1 : grd.ncy : grd.ncx * grd.ncy - grd.ncy - 1];
dx = [grd.ncx * grd.ncy - 2 * grd.ncy + 2 : grd.ncx * grd.ncy - grd.ncy - 1];
top = [2 + grd.ncy : grd.ncy : grd.ncx * grd.ncy - grd.ncy - grd.ncy + 2];

msh_to_ptcl = init_ptcl_mesh (msh, ptcls);

% COPY OF INTERNAL-BORDER PTCLS IN GC
for ic = sx
  particles = grd_to_ptcl(ic);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = [L1; 0] + ptcls.x(:, k);
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end
end
for ic = dx
  particles = grd_to_ptcl(ic);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = ptcls.x(:, k) - [L1; 0];
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end
end
for ic = top
  particles = grd_to_ptcl(ic);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = [0; L2] + ptcls.x(:, k);
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end
end
for ic = bot
  particles = grd_to_ptcl(ic);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = ptcls.x(:, k) + [0; -L2];
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end
end
%% sx + top
particles = grd_to_ptcl(2 + grd.ncy);
particles = particles{1};
for k = particles
  ptcls_gc.x(:, end + 1) = [L1; L2] + ptcls.x(:, k);
  ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
  ptcls_gc.q(end + 1) = ptcls.q(k);
end
%% sx + bot
particles = grd_to_ptcl(2 * grd.ncy - 1);
particles = particles{1};
for k = particles
  ptcls_gc.x(:, end + 1) = [L1; -L2] + ptcls.x(:, k);
  ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
  ptcls_gc.q(end + 1) = ptcls.q(k);
end
%% dx + top
particles = grd_to_ptcl(grd.ncx * grd.ncy - grd.ncy - grd.ncy + 2);
particles = particles{1};
for k = particles
  ptcls_gc.x(:, end + 1) = [- L1; L2] + ptcls.x(:, k);
  ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
  ptcls_gc.q(end + 1) = ptcls.q(k);
end
%% dx + bot
particles = grd_to_ptcl(grd.ncx * grd.ncy - grd.ncy - 1);
particles = particles{1};
for k = particles
  ptcls_gc.x(:, end + 1) = [- L1; -L2] + ptcls.x(:, k);
  ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
  ptcls_gc.q(end + 1) = ptcls.q(k);
end

ptcls_all.x = [ptcls.x, ptcls_gc.x];
ptcls_all.v = [ptcls.v, ptcls_gc.v];
ptcls_all.q = [ptcls.q, ptcls_gc.q];

grd_to_ptcl = init_ptcl_mesh (grd, ptcls_all);

%espressione di rho per la singola particella centrata in x0,y0 valutata al punto x,y
rho1 = @(x,y,x0,y0) ((1 - sqrt((x - x0).^2 + (y - y0).^2) ./ r_cut) .* 3 ./ (r_cut^2 * pi)).* (((x - x0).^2 + (y - y0).^2) <= r_cut^2);

rho = zeros(msh.ny, msh.nx);
for k = 1:size(ptcls_all.x,2)
  rho = rho + rho1(X, Y, ptcls_all.x(1,k), ptcls_all.x(2,k)) .* ptcls_all.q(k);
end

[phi, L, U, P, Q, R] = poisson (msh, 1, rho);

[dphidpx, dphidpy] = phigradient (ptcls, msh_to_ptcl, msh, phi);

F = force(grd, ptcls_all, grd_to_ptcl, r_cut, sigma, epsilon, epsilon0);
F = F(:, 1:Nparticelle);

F(1,:) = F(1,:) - 4 * pi * (dphidpx .* ptcls.q) / epsilon0;
F(2,:) = F(2,:) - 4 * pi * (dphidpy .* ptcls.q) / epsilon0;

%% Time loop
for ii = 2:1:steps
  if (mod(ii,salvaOgni) == 0)
    save ("-append",strcat("../posizioni/posizioni_",int2str(ii)),"ptcls")
  end
  % VELOCITY - VERLET
  % CALCOLO VELOCITA T+DT/2 CON LA FORZA ALL'ISTANTE T
  v_aux = ptcls.v + 0.5*(F./m)*dt;
  % CALCOLO POSIZIONI T+DT
  ptcls.x = ptcls.x + v_aux*dt;

	% TELEPORT
  for k = 1:1:Nparticelle
    ptcls.x(1,k) = r_cut + mod(ptcls.x(1,k) - r_cut,L1);
    ptcls.x(2,k) = r_cut + mod(ptcls.x(2,k) - r_cut,L2);
  end

  grd_to_ptcl = init_ptcl_mesh (grd, ptcls);

	% COPY OF INTERNAL-BORDER PTCLS IN GC
  for ic = sx
    particles = grd_to_ptcl(ic);
    particles = particles{1};
    for k = particles
      ptcls_gc.x(:, end + 1) = [L1; 0] + ptcls.x(:, k);
      ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
      ptcls_gc.q(end + 1) = ptcls.q(k);
    end
  end
  for ic = dx
    particles = grd_to_ptcl(ic);
    particles = particles{1};
    for k = particles
      ptcls_gc.x(:, end + 1) = ptcls.x(:, k) - [L1; 0];
      ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
      ptcls_gc.q(end + 1) = ptcls.q(k);
    end
  end
  for ic = top
    particles = grd_to_ptcl(ic);
    particles = particles{1};
    for k = particles
      ptcls_gc.x(:, end + 1) = [0; L2] + ptcls.x(:, k);
      ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
      ptcls_gc.q(end + 1) = ptcls.q(k);
    end
  end
  for ic = bot
    particles = grd_to_ptcl(ic);
    particles = particles{1};
    for k = particles
      ptcls_gc.x(:, end + 1) = ptcls.x(:, k) + [0; -L2];
      ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
      ptcls_gc.q(end + 1) = ptcls.q(k);
    end
  end
  %% sx + top
  particles = grd_to_ptcl(2 + grd.ncy);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = [L1; L2] + ptcls.x(:, k);
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end
  %% sx + bot
  particles = grd_to_ptcl(2 * grd.ncy - 1);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = [L1; -L2] + ptcls.x(:, k);
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end
  %% dx + top
  particles = grd_to_ptcl(grd.ncx * grd.ncy - grd.ncy - grd.ncy + 2);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = [- L1; L2] + ptcls.x(:, k);
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end
  %% dx + bot
  particles = grd_to_ptcl(grd.ncx * grd.ncy - grd.ncy - 1);
  particles = particles{1};
  for k = particles
    ptcls_gc.x(:, end + 1) = [- L1; -L2] + ptcls.x(:, k);
    ptcls_gc.v(:, end + 1) = ptcls.v(:, k);
    ptcls_gc.q(end + 1) = ptcls.q(k);
  end

  ptcls_all.x = [ptcls.x, ptcls_gc.x];
  ptcls_all.v = [ptcls.v, ptcls_gc.v];
  ptcls_all.q = [ptcls.q, ptcls_gc.q];

  grd_to_ptcl = init_ptcl_mesh (grd, ptcls_all);

  % CALCOLO FORZE ALL'ISTANTE T+DT
  F = force(grd, ptcls_all, grd_to_ptcl, r_cut, sigma, epsilon, epsilon0);
  F = F(:, 1:Nparticelle);

  ptcls_gc.x = [];
  ptcls_gc.v = [];

  rho = 0 * rho;
  for k = 1:size(ptcls_all.x,2)
    rho = rho + rho1(X, Y, ptcls_all.x(1,k), ptcls_all.x(2,k)) .* ptcls_all.q(k);
  end

  [phi] = poisson (msh, epsilon, rho, L, U, P, Q, R);

  msh_to_ptcl = init_ptcl_mesh (msh, ptcls);

  [dphidpx, dphidpy] = phigradient (ptcls, msh_to_ptcl, msh, phi);
  F(1,:) = F(1,:) - 4 * pi * (dphidpx .* ptcls.q) / epsilon0;
  F(2,:) = F(2,:) - 4 * pi * (dphidpy .* ptcls.q) / epsilon0;

  % CALCOLO VELOCITA T+DT CON LA FORZA ALL'ISTANTE T+DT
  ptcls.v = v_aux + 0.5*(F./m)*dt;

end
