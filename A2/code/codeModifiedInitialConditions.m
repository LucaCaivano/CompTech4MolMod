
clear all
clc
close all



scale = 1;

Simulazione = 1;

L1 = 250;
L2 = 200;
epsilon = 5;
sigma = 1;
r_cut = 2.5*sigma
m = 1;
noise = 0.1;


if Simulazione == 1
dt = 0.005;
elseif Simulazione == 2
dt = 0.00005;
end
NumeroIstante0 = 17700;
T = 20;
T0 = NumeroIstante0*dt;
Tf = T0 + T;
steps = T/dt + 1;

if Simulazione == 1
  salvaOgni = 10;
elseif Simulazione == 2
  salvaOgni = 300;
end


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

if Simulazione == 1
  aux = load(strcat("./posizioni1/posizioni_",int2str(NumeroIstante0)));
elseif Simulazione == 2
  aux = load(strcat("./posizioni2/posizioni_",int2str(NumeroIstante0)));
end  
ptcls = aux.ptcls;


grd_to_ptcl = init_ptcl_mesh (grd, ptcls);
d = cellfun (@numel, grd_to_ptcl, 'UniformOutput', true);



% Initial forces
F = force(grd, ptcls, grd_to_ptcl, r_cut, sigma, epsilon);

%% Time loop
for ii = 2:1:steps
  if (mod(ii,salvaOgni) == 0)
    if Simulazione == 1
      save ("-append",strcat("./posizioni/posizioni_",int2str(NumeroIstante0 + (ii/1))),"ptcls")
    elseif Simulazione == 2
      save ("-append",strcat("./posizioni/posizioni_",int2str(NumeroIstante0 + (ii/1))),"ptcls")
    end
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
