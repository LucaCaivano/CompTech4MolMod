

function [q,v] = shake(NumeroParticelle,n,m,M,V,VP,q0,v0,Force,T,dt,err)

  DimensioneQ = numel(q0);

  NumeroFisse = max(VP(:,1)) - NumeroParticelle;
  M = [M;ones(NumeroFisse,1)];
    
  N = T/dt + 1;
  q = zeros(DimensioneQ,N);
  v = zeros(DimensioneQ,N);
  q(:,1) = q0;
  v(:,1) = v0;

  Masse_aux = zeros(NumeroParticelle*2,1);
  Masse_aux(1:2:(NumeroParticelle*2-1)) = M(1:NumeroParticelle);
  Masse_aux(2:2:NumeroParticelle*2) = M(1:NumeroParticelle);

  %Runge-Kutta 2
  a = Force(q0)./Masse_aux;
  a = [a; zeros(NumeroFisse*2,1)];

  v_aux = v0 + a*dt/2;

  qtilda = q0 + v_aux*dt;

  sol = initial_F(n,m,M,V,VP,dt,q0,qtilda);
  f = F_constraint(n,m,M,V,VP,dt,q0,qtilda);
  [gamma, fval, info] = fsolve(f, sol, optimset("TolFun", err));

  delta_Mat = zeros(DimensioneQ, n+m);

  for h = 1:1:NumeroParticelle
    for k = 1:1:(n+m)

      if (VP(k,1) == h)
	aux = VP(k,2);
	delta_Mat([2*h-1,2*h],k) = 2*(q0([2*h-1,2*h]) - q0([2*aux-1,2*aux]))*(dt^2)/M(h);
      elseif (VP(k,2) == h)
	aux = VP(k,1);
	delta_Mat([2*h-1,2*h],k) = 2*(q0([2*h-1,2*h]) - q0([2*aux-1,2*aux]))*(dt^2)/M(h);	
      endif
      
    endfor
  endfor

  q(:,2) = qtilda + delta_Mat*gamma;

  
  for i = 3:1:N    
    
    qt = q(:,i-1);
    qtmh = q(:,i-2);

    a = Force(qt)./Masse_aux;
    a = [a; zeros(NumeroFisse*2,1)];
    
    qtilda = 2*qt - qtmh + a*(dt^2);

    sol = initial_F(n,m,M,V,VP,dt,qt,qtilda);
    
    f = F_constraint(n,m,M,V,VP,dt,qt,qtilda);

    [gamma, fval, info] = fsolve(f, sol, optimset("TolFun", err));

    
    delta_Mat = zeros(DimensioneQ, n+m);

    for h = 1:1:NumeroParticelle
      for k = 1:1:(n+m)

	if (VP(k,1) == h)
	  aux = VP(k,2);
	  delta_Mat([2*h-1,2*h],k) = 2*(qt([2*h-1,2*h]) - qt([2*aux-1,2*aux]))*(dt^2)/M(h);
	elseif (VP(k,2) == h)
	  aux = VP(k,1);
	  delta_Mat([2*h-1,2*h],k) = 2*(qt([2*h-1,2*h]) - qt([2*aux-1,2*aux]))*(dt^2)/M(h);	
	endif
	
      endfor
    endfor

    q(:,i) = qtilda + delta_Mat*gamma;
    v(:,i-1) = (q(:,i)- q(:,i-2))/(2*dt);
    
  endfor

  v(:,N) = v(:,N-1); %brutto, da aggiustare

endfunction
