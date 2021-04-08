
L = 1;
m = 1;
g = 10;

th0 = pi/3;
q0 = L*[sin(th0),-cos(th0)]';
v0 = [0,0]';


T = 10;
dt = 0.01;
err = 10^(-10);


N = T/dt + 1;
q = zeros(2,N);
v = zeros(2,N);
q(:,1) = q0;
v(:,1) = v0;
TT = zeros(1,N);
UU = zeros(1,N);
MM = zeros(1,N);


q(:,2) = q0 + v0*dt;


for i = 3:1:N

  qt = q(:,i-1);
  qtmh = q(:,i-2);
  
  A = 4*((dt^4)*norm(qt)^2)/(m^2);
  B = 4*((dt^2)*( qt(1)*(2*qt(1) - qtmh(1)) + qt(2)*(2*qt(2) - qtmh(2) - g*dt^2) ))/m;
  C = (2*qt(1) - qtmh(1))^2 + (2*qt(2) - qtmh(2) - g*dt^2)^2 - L^2;

  gamma = - C/B;

  do

    aux = (A*gamma^2 + B*gamma + C)/(2*A*gamma+B);
    gamma = gamma - aux;
    
  until (abs(aux) < err)

  q(:,i) = 2*(1+((dt^2)*gamma)/m)*qt - qtmh - [0, g*dt^2]';
  v(:,i-1) = (q(:,i) - qtmh)/(2*dt);
  
endfor

v(:,N) = v(:,N-1); %brutto, da aggiustare

for i = 1:1:N

  TT(i) = 0.5 * m * norm(v(:,i))^2;
  UU(i) = m*g*(q(2,i) + L);
  MM(i) = TT(i) + UU(i);
  
endfor

figure
subplot(2,1,2);
plot(0:dt:T, TT)
grid on
hold on
plot(0:dt:T, UU)
plot(0:dt:T, MM)
legend("Energia cinetica","Energia potenziale","Energia totale")

subplot(2,2,1)
plot(0:dt:T, sqrt((q(1,:)').^2 + (q(2,:)').^2) )
grid on

subplot(2,2,2)
mostraPendolo(q(1,:),q(2,:),[-1,1,-1,1],dt)

