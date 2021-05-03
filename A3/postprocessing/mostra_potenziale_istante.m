
clear
clc

addpath("../codice")
addpath("../poisson")

istante = 2400;


scale = 0;

L1 = 60 + 84*scale;
L2 = 60;
epsilon = 1;
epsilon0 = 1e-3;
sigma = 1;
r_cut = 6;
m = 23;
noise = 0.1;

H1 = 19;
H2 = 38;
H3 = 19;
W = 180*L1/144;


grd.ncy = L2/r_cut + 2;
grd.ncx = L1/r_cut + 2;
grd.x = linspace (0, L1 + 2*r_cut, grd.ncx+1);
grd.y = linspace (0, L2 + 2*r_cut, grd.ncy+1);


delta = 0.4;

[X1, Y1] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(delta,L2/3 - delta, H1) + r_cut);
[X2, Y2] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(L2/3 + delta ,2*L2/3 - delta, H2) + r_cut);
[X3, Y3] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(2*L2/3 + delta ,L2 - delta, H3) + r_cut);



N1 = numel(X1);
N2 = numel(X2);
N3 = numel(X3);
Nparticelle = N1 + N2 + N3;



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

if istante == 1
  aux = load(strcat("../posizioni/posizioni_iniziali"));
else
  aux = load(strcat("../posizioni/posizioni_",int2str(istante)));
end

ptcls.x = aux.ptcls.x;
ptcls.v = aux.ptcls.v;
ptcls.q = aux.ptcls.q;



ptcls_gc.x = [];
ptcls_gc.v = [];
ptcls_gc.q = [];

sx = [2 + grd.ncy : grd.ncy * 2 - 1];
bot = [2 * grd.ncy - 1 : grd.ncy : grd.ncx * grd.ncy - grd.ncy - 1];
dx = [grd.ncx * grd.ncy - 2 * grd.ncy + 2 : grd.ncx * grd.ncy - grd.ncy - 1];
top = [2 + grd.ncy : grd.ncy : grd.ncx * grd.ncy - grd.ncy - grd.ncy + 2];

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






rho1 = @(x,y,x0,y0) ((1 - sqrt((x - x0).^2 + (y - y0).^2) ./ r_cut) .* 3 ./ (r_cut^2 * pi)).* (((x - x0).^2 + (y - y0).^2) <= r_cut^2);


rho = zeros(msh.ny, msh.nx);
for k = 1:size(ptcls_all.x,2)
  rho = rho + rho1(X, Y, ptcls_all.x(1,k), ptcls_all.x(2,k)) .* ptcls_all.q(k);
end
subplot(2,2,1)
mesh (X, Y, reshape (rho, size (X)));
title("rho")
subplot(2,2,2)
[phi, L, U, P, Q, R] = poisson (msh, epsilon, rho);
mesh (X, Y, reshape (phi, size (X)));
title("phi")
subplot(2,2,3)
axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
plot([r_cut, r_cut, (L1 + r_cut), (L1 + r_cut), r_cut], [r_cut, (L2 + r_cut), (L2 + r_cut), r_cut , r_cut])
hold on
axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
plot(ptcls.x(1,1:N1),ptcls.x(2,1:N1), "ob", "markersize", 3, "markerfacecolor", "b")
axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
plot(ptcls.x(1,(N1+1):(N1 + N2)),ptcls.x(2,(N1+1):(N1 + N2)), "or", "markersize", 3, "markerfacecolor", "r")
axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
plot(ptcls.x(1,(N1+N2+1):end),ptcls.x(2,(N1+N2+1):end), "ob", "markersize", 3, "markerfacecolor", "b")
hold off
subplot(2,2,4)
contourf (X, Y, reshape (phi, size (X)));
colorbar
