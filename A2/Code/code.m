
clear all
clc
close all
%% Parameters
scale = 1;

%alfa = 20;
%beta = 16;
%L1 = 250 - alfa*scale;
%L2 = 200 - beta*scale;
L1 = 250;
L2 = 200;
epsilon = 5;
sigma = 1;
r_cut = 2.5*sigma;
m = 1;
noise = 0.1;


%dt = 0.00005;
dt = 0.005;
T = 10;
steps = T/dt + 1;


delta = sigma*2^(1/6);

H1 = (40/scale);
W1 = (40/scale);
H2 = (40/scale);
W2 = (160/scale);

H1_l = (H1-1)*delta;
H2_l = (H2-1)*delta;
W1_l = (W1-1)*delta;
W2_l = (W2-1)*delta;


grd.ncy = L2/r_cut;
grd.ncx = L1/r_cut;
grd.x = linspace (0, L1, grd.ncx+1);
grd.y = linspace (0, L2, grd.ncy+1);

Dist = 5;
yc_1 = 120-20*(1-(1/scale));
yc_2 = 80-Dist+20*(1-(1/scale));

[X1, Y1] = meshgrid (linspace ((L1-W1_l)/2, (L1+W1_l)/2, W1),...
		     linspace ((yc_1)-H1_l/2, (yc_1)+H1_l/2, H1)); 
[X2, Y2] = meshgrid (linspace ((L1-W2_l)/2, (L1+W2_l)/2, W2),...
		     linspace ((yc_2)-H2_l/2, (yc_2)+H2_l/2, H2)); 

N1 = numel(X1);
N2 = numel(X2);
Nparticelle = N1 + N2;

ptcls.x = [[X1(:);X2(:)], [Y1(:);Y2(:)]]';
ptcls.v = randn (size (ptcls.x)) * noise;
ptcls.v(2, 1:N1) = ptcls.v(2, 1:N1) -10;
grd_to_ptcl = init_ptcl_mesh (grd, ptcls);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);


% Initial forces
F = force(grd, ptcls, grd_to_ptcl, r_cut, sigma, epsilon);

%% Time loop
for ii = 2:1:steps
  if (mod(ii,100) == 0)
    save ("-append",strcat("./posizioni/posizioni_",int2str(ii)),"ptcls")
  end
  % VELOCITY - VERLET
  % CALCOLO VELOCITA T+DT/2 CON LA FORZA ALL'ISTANTE T
  v_aux = ptcls.v + 0.5*(F./m)*dt;
  % CALCOLO POSIZIONI T+DT
  ptcls.x = ptcls.x + v_aux*dt;
  % BOUNDARY CONDITIONS
  for k = 1:1:Nparticelle
    qx = ptcls.x(1,k);
    qy = ptcls.x(2,k);
    if qx < 0
      ptcls.x(1,k) = -qx;
      v_aux(1,k) = -v_aux(1,k);
    elseif qx > L1
      ptcls.x(1,k) = L1 - (qx - L1);
      v_aux(1,k) = -v_aux(1,k);
    end
    if qy < 0
      ptcls.x(2,k) = -qy;
      v_aux(2,k) = -v_aux(2,k);
    elseif qy > L2
      ptcls.x(2,k) = L2 - (qy - L2);
      v_aux(2,k) = -v_aux(2,k);
    end
  end
  try
    grd_to_ptcl = init_ptcl_mesh (grd, ptcls);
  catch
    ii
    mostraParticelle(ptcls.x,N1,N2,L1,L2,r_cut);
    break
  end
  % CALCOLO FORZE ALL'ISTANTE T+DT
  F = force(grd, ptcls, grd_to_ptcl, r_cut, sigma, epsilon);
  % CALCOLO VELOCITA T+DT CON LA FORZA ALL'ISTANTE T+DT
  ptcls.v = v_aux + 0.5*(F./m)*dt;
end
