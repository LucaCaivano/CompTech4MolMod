
function mostraStruttura (varargin)

  [hax, varargin, nargin] = __plt_get_axis_arg__ ("comet", varargin{:});

%  if (nargin != 6)
%    print_usage ();
%  elseif (nargin == 6)
    VP = varargin{1};
    q = varargin{2};
    N = varargin{3}
    limiti = varargin{4};
    p = varargin{5};
    skip = varargin{6};
%  endif

  prendi = 1:skip:numel(q(1,:));
  q = q(:,prendi);
  p = p*skip;

  i = 1:2:numel(q(:,1));
  x = q(i,:);
  y = q(i+1,:);
  
  oldfig = [];
  if (! isempty (hax))
    oldfig = get (0, "currentfigure");
  endif
  unwind_protect
    hax = newplot (hax);
    limits = limiti;
    num = numel (q(1,:));
    dn = round (num/10);

    %hl = plot (x(:,1), y(:,1), "color", "b", "marker", "o", "linestyle", "none");
    %hold on
    axis (limits);  # set manual limits to speed up plotting

    ## Initialize the timer
    t = p;
    timerid = tic ();

    grid on

    linee = zeros(numel(VP(:,1)),1);
    
    for n = 2:(num)
      m = n - dn;
      m = max ([m, 1]);
      k = min ([n, num]);
      %set (hl(1), "xdata", x(n), "ydata", y(n));

      for iter_linee = 1:1:numel(VP(:,1))

	alfa = VP(iter_linee,1);
	beta = VP(iter_linee,2);
	
	linee(iter_linee) = line([x(alfa,n),x(beta,n)],[y(alfa,n),y(beta,n)],"linestyle","-","color","k");
	
      endfor

      %l1 = line([0,x1(k)],[0,y1(k)],"linestyle","-","color","k");
      %l2 = line([x1(k),x2(k)],[y1(k),y2(k)],"linestyle","-","color","k");
      pause (t - toc (timerid));
      if n < num
	for iter_linee = 1:1:numel(VP(:,1))
	  delete(linee(iter_linee))
	endfor	
      endif
      t += p;
    endfor

  unwind_protect_cleanup
    if (! isempty (oldfig))
      set (0, "currentfigure", oldfig);
    endif
  end_unwind_protect

endfunction

