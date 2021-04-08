
function mostraDoppioPendolo (varargin)

  [hax, varargin, nargin] = __plt_get_axis_arg__ ("comet", varargin{:});

  if (nargin != 7)
    print_usage ();
  elseif (nargin == 7)
    x1 = varargin{1};
    y1 = varargin{2};
    x2 = varargin{3};
    y2 = varargin{4};
    
    limiti = varargin{5};
    p = varargin{6};
    skip = varargin{7};
  endif

  prendi = 1:skip:numel(y1);
  x1 = x1(prendi);
  y1 = y1(prendi);
  x2 = x2(prendi);
  y2 = y2(prendi);
  p = p*skip;
  
  
  oldfig = [];
  if (! isempty (hax))
    oldfig = get (0, "currentfigure");
  endif
  unwind_protect
    hax = newplot (hax);
    limits = limiti;
    num = numel (y1);
    dn = round (num/10);

    hl = plot (x1(1), y1(1), "color", "m", "marker", "none",
               x1(1), y1(1), "color", "c", "marker", "none",
               x1(1), y1(1), "color", "b", "marker", "o");
    hold on
    h2 = plot (x2(1), y2(1), "color", "r", "marker", "none",
               x2(1), y2(1), "color", "g", "marker", "none",
               x2(1), y2(1), "color", "b", "marker", "o");
    axis (limits);  # set manual limits to speed up plotting

    ## Initialize the timer
    t = p;
    timerid = tic ();

    grid on
    
    for n = 2:(num+dn)
      m = n - dn;
      m = max ([m, 1]);
      k = min ([n, num]);
      set (hl(1), "xdata", x1(1:m), "ydata", y1(1:m));
      set (hl(2), "xdata", x1(m:k), "ydata", y1(m:k));
      set (hl(3), "xdata", x1(k),   "ydata", y1(k));
      set (h2(1), "xdata", x2(1:m), "ydata", y2(1:m));
      set (h2(2), "xdata", x2(m:k), "ydata", y2(m:k));
      set (h2(3), "xdata", x2(k),   "ydata", y2(k));
      l1 = line([0,x1(k)],[0,y1(k)],"linestyle","-","color","k");
      l2 = line([x1(k),x2(k)],[y1(k),y2(k)],"linestyle","-","color","k");
      pause (t - toc (timerid));
      if n < num+dn
	delete(l1)
	delete(l2)
      endif
      t += p;
    endfor

  unwind_protect_cleanup
    if (! isempty (oldfig))
      set (0, "currentfigure", oldfig);
    endif
  end_unwind_protect

endfunction

