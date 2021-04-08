function init = initial_F (n,m,M,V,VP,dt,qt,qtilda)


  A = zeros(n+m,n+m);
  b = zeros(n+m,1);
    
  for k = 1:1:n

    alfa = VP(k,1);
    beta = VP(k,2);

    
    for i = 1:1:(n+m)

      if (V(k,i) == 0)
	continue
      endif

      P = 1;
      
      if i == k
	P = 2*(1/M(alfa) + 1/M(beta))*P*(qt([2*alfa - 1, 2*alfa]) - qt([2*beta - 1, 2*beta]));
      else
	aux = 0;
	if (VP(i,1) == alfa)
	  aux = VP(i,2);
	  P = (2*P*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	elseif (VP(i,2) == alfa)
	  aux = VP(i,1);
	  P = (2*P*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	elseif (VP(i,1) == beta)
	  aux = VP(i,2);
	  P = (2*P*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	else
	  aux = VP(i,1);
	  P = (2*P*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	endif
      endif
            
      A(k,i) = 2*(dt^2)*(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))'*P;
      
    endfor

    b(k) = - norm(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))^2 + V(k,k)^2;
    
  endfor

  for k = (n+1):1:(n+m)

    alfa = VP(k,1);
    beta = VP(k,2);
    
    for i = 1:1:(n+m)

      if (V(k,i) == 0)
	continue
      endif
      
      P = -2*(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))/M(beta);
      if (VP(i,1) == beta)
	aux = VP(i,2);
	P = 2*P'*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
      elseif (VP(i,2) == beta)
	aux = VP(i,1);
	P = 2*P'*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
      else
	continue
      endif

      A(k,i) = (dt^2)*P;      
    endfor

    b(k) = - norm(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))^2 + V(k,k)^2;
    
  endfor

  init = A\b;
  
endfunction
