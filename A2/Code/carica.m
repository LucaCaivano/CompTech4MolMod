
dt = 0.005;
step = 100;
dt = dt*step;
T = 11;
times = dt:dt:T;

delay = 1;

Istanti = numel(times);

N1 = 1600;
N2 = 6400;
Particelle = N1 + N2;

Posizioni = zeros(2,Particelle,Istanti);

for ii = 1:1:Istanti

  aux = load(strcat("./posizioni/posizioni_",int2str(ii*step)));

  Posizioni(:,:,ii) = aux.ptcls.x;

end


for ii = 1:1:Istanti  
  plot(Posizioni(1,1:N1,ii),Posizioni(2,1:N1,ii), "or", "markersize", 3, "markerfacecolor", "r")
  title(strcat("T = ",num2str(times(ii))))
  axis([0 250 0 200])
  hold on
  plot(Posizioni(1,(N1+1):end,ii),Posizioni(2,(N1+1):end,ii), "ob", "markersize", 3, "markerfacecolor", "b")
  hold off
  pause(delay)
end
