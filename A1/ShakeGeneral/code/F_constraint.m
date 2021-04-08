

function ret_fun = F_constraint (n,m,M,V,VP,dt,qt,qtilda)

  ret_fun = @(gamma) zeros(n+m,1);

  
  for k = 1:1:n

    alfa = VP(k,1);
    beta = VP(k,2);
    vettore_indicatore = zeros(n+m,1);
    vettore_indicatore(k) = 1;
    
    for i=1:1:(n+m)
      for j=1:1:(n+m)

	if (V(k,i) == 0 || V(k,j) == 0)
	  continue
	endif

	P = 1;	
	if i == k %#1
	  P = 2*(1/M(alfa) + 1/M(beta))*P*(qt([2*alfa - 1, 2*alfa]) - qt([2*beta - 1, 2*beta]));
	else
	  aux = 0;
	  if (VP(i,1) == alfa) %#2
	    aux = VP(i,2);
	    P = (2*P*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	  elseif (VP(i,2) == alfa) %#3
	    aux = VP(i,1);
	    P = (2*P*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	  elseif (VP(i,1) == beta) %#4
	    aux = VP(i,2);
	    P = (2*P*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	  else %#5
	    aux = VP(i,1);
	    P = (2*P*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	  endif
	endif

	if j == k %#6
	  P = 2*(1/M(alfa) + 1/M(beta))*P'*(qt([2*alfa - 1, 2*alfa]) - qt([2*beta - 1, 2*beta]));
	else
	  aux = 0;
	  if (VP(j,1) == alfa) %#7
	    aux = VP(j,2);
	    P = (2*P'*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	  elseif (VP(j,2) == alfa) %#8
	    aux = VP(j,1);
	    P = (2*P'*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	  elseif (VP(j,1) == beta) %#9
	    aux = VP(j,2);
	    P = (2*P'*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	  else %#10
	    aux = VP(j,1);
	    P = (2*P'*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	  endif
	endif

	ret_fun = @(gamma) ret_fun(gamma) + (dt^4)*P*gamma(i)*gamma(j)*vettore_indicatore;
	
      endfor
    endfor    

    for i = 1:1:(n+m)

      if (V(k,i) == 0)
	continue
      endif

      P = 1;
      
      if i == k %#11
	P = 2*(1/M(alfa) + 1/M(beta))*P*(qt([2*alfa - 1, 2*alfa]) - qt([2*beta - 1, 2*beta]));
      else
	aux = 0;
	if (VP(i,1) == alfa) %#12
	  aux = VP(i,2);
	  P = (2*P*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	elseif (VP(i,2) == alfa) %#13
	  aux = VP(i,1);
	  P = (2*P*(qt([2*alfa-1, 2*alfa]) - qt([2*aux-1, 2*aux])))/M(alfa);
	elseif (VP(i,1) == beta) %#14
	  aux = VP(i,2);
	  P = (2*P*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	else %#15
	  aux = VP(i,1);
	  P = (2*P*(qt([2*aux-1, 2*aux]) - qt([2*beta-1, 2*beta])))/M(beta);
	endif
      endif

      
      ret_fun = @(gamma) ret_fun(gamma) + 2*(dt^2)*gamma(i)*(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))'*P*vettore_indicatore;
      
    endfor

    ret_fun = @(gamma) ret_fun(gamma) + (norm(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))^2 - V(k,k)^2)*vettore_indicatore;
    
  endfor


  

  for k = (n+1):1:(n+m)

    alfa = VP(k,1);
    beta = VP(k,2);
    vettore_indicatore = zeros(n+m,1);
    vettore_indicatore(k) = 1;
    
    for i = 1:1:(n+m)
      for j = 1:1:(n+m)

	if(V(k,i) == 0 || V(k,j) == 0)
	  continue
	endif

	P = 1;
	if (VP(i,1) == beta) %#16
	  aux = VP(i,2);
	  P = 2*P*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
	elseif (VP(i,2) == beta) %#17
	  aux = VP(i,1);
	  P = 2*P*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
	else
	  continue
	endif
	if (VP(j,1) == beta) %#18
	  aux = VP(j,2);
	  P = 2*P'*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
	elseif (VP(j,2) == beta) %#19
	  aux = VP(j,1);
	  P = 2*P'*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
	else
	  continue
	endif

	ret_fun = @(gamma) ret_fun(gamma) + gamma(i)*gamma(j)*(dt^4)*P/(M(beta)^2)*vettore_indicatore;
	
      endfor
    endfor

    for i = 1:1:(n+m)

      if (V(k,i) == 0)
	continue
      endif
      
      P = -2*(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))/M(beta);
      if (VP(i,1) == beta) %#20
	aux = VP(i,2);
	P = 2*P'*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
      elseif (VP(i,2) == beta) %#21
	aux = VP(i,1);
	P = 2*P'*(qt([2*beta-1,2*beta]) - qt([2*aux-1,2*aux]));
      else
	continue
      endif

      ret_fun = @(gamma) ret_fun(gamma) + (dt^2)*P*gamma(i)*vettore_indicatore;      
    endfor

    ret_fun = @(gamma) ret_fun(gamma) + (norm(qtilda([2*alfa-1, 2*alfa]) - qtilda([2*beta-1, 2*beta]))^2 - V(k,k)^2)*vettore_indicatore;
    
  endfor

  
endfunction
