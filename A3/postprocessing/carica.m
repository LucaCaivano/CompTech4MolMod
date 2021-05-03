
scale = -0.5;

dt = 0.00001;
step = 5;
T = 190*dt;
dt = dt*step;
times = dt:dt:T;


L1 = 60 + 84*scale;
L2 = 60;
epsilon = 1;
sigma = 1;
r_cut = 6;


delay = 0.001;

Istanti = numel(times)

H1 = 19;
H2 = 38;
H3 = 19;
W = 180*L1/144;
delta = 0.4;
[X1, Y1] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(delta,L2/4 - delta, H1) + r_cut);
[X2, Y2] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(L2/4 + delta, 3*L2/4 - delta, H2) + r_cut);
[X3, Y3] = meshgrid (linspace(delta,L1-delta,W) + r_cut, linspace(3*L2/4 + delta ,L2 - delta, H3) + r_cut);

N1 = numel(X1);
N2 = numel(X2);
N3 = numel(X3);
Particelle = N1 + N2 + N3;

Posizioni = zeros(2,Particelle,Istanti);

for ii = 2:1:Istanti
  aux = load(strcat("../posizioni/posizioni_",int2str(ii*step)));
  Posizioni(:,:,ii) = aux.ptcls.x;
end

Posizioni_iniziali = zeros(2,Particelle);
aux = load(strcat("../posizioni/posizioni_iniziali"));
Posizioni_iniziali(:,:,ii) = aux.ptcls.x;


figure('position',[1000,0,1000,1250])
pause(1)
plot([r_cut, r_cut, (L1 + r_cut), (L1 + r_cut), r_cut], [r_cut, (L2 + r_cut), (L2 + r_cut), r_cut , r_cut])
hold on
plot(Posizioni_iniziali(1,1:N1,ii),Posizioni_iniziali(2,1:N1,ii), "ob", "markersize", 3, "markerfacecolor", "b")
axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
plot(Posizioni_iniziali(1,(N1+1):(N1 + N2),ii),Posizioni_iniziali(2,(N1+1):(N1+N2),ii), "or", "markersize", 3, "markerfacecolor", "r")
plot(Posizioni_iniziali(1,(N1+N2+1):end,ii),Posizioni_iniziali(2,(N1+N2+1):end,ii), "ob", "markersize", 3, "markerfacecolor", "b")
hold off
pause(delay)

figure('position',[0,0,1000,1250])
pause(1)
for ii = 1:1:Istanti
  axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
  plot([r_cut, r_cut, (L1 + r_cut), (L1 + r_cut), r_cut], [r_cut, (L2 + r_cut), (L2 + r_cut), r_cut , r_cut])
  hold on
  axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
  plot(Posizioni(1,1:N1,ii),Posizioni(2,1:N1,ii), "ob", "markersize", 3, "markerfacecolor", "b")
  axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
  plot(Posizioni(1,(N1+1):(N1 + N2),ii),Posizioni(2,(N1+1):(N1 + N2),ii), "or", "markersize", 3, "markerfacecolor", "r")
  axis([0 (L1 + 2*r_cut) 0 (L2 + 2*r_cut)])
  plot(Posizioni(1,(N1+N2+1):end,ii),Posizioni(2,(N1+N2+1):end,ii), "ob", "markersize", 3, "markerfacecolor", "b")
  hold off
  pause(delay)
end
