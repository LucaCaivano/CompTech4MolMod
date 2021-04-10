
Simulazione = 1;

if Simulazione == 1
  dt = 0.005;
  step = 100;
  T = 17500*dt;
  dt = dt*step;
  times = dt:dt:T;
elseif Simulazione == 2
  dt = 0.00005;
  step = 300;
  T = 318000*dt;
  dt = dt*step;
  times = dt:dt:T;
else
  return
end


delay = 0.1;

Istanti = numel(times)

N1 = 1600;
N2 = 6400;
Particelle = N1 + N2;

Posizioni = zeros(2,Particelle,Istanti);

for ii = 1:1:Istanti
  if Simulazione == 1
    aux = load(strcat("./posizioni1/posizioni_",int2str(ii*step)));
  else
    aux = load(strcat("./posizioni2/posizioni_",int2str(ii*step)));
  end
  Posizioni(:,:,ii) = aux.ptcls.x;
end

figure('position',[1800,0,1000,1250])
pause(1)
for ii = 1:1:Istanti
  plot(Posizioni(1,1:N1,ii),Posizioni(2,1:N1,ii), "or", "markersize", 3, "markerfacecolor", "r")
  axis([0 250 0 200])
  hold on
  plot(Posizioni(1,(N1+1):end,ii),Posizioni(2,(N1+1):end,ii), "ob", "markersize", 3, "markerfacecolor", "b")
  hold off
  pause(delay)
end
