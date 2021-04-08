
function mostraPendolo (varargin)

  [hax, varargin, nargin] = __plt_get_axis_arg__ ("comet", varargin{:});

  if (nargin != 4)
    print_usage ();
  elseif (nargin == 4)
    x = varargin{1};
    y = varargin{2};
    limiti = varargin{3};
    p = varargin{4};
  endif

  
  oldfig = [];
  if (! isempty (hax))
    oldfig = get (0, "currentfigure");
  endif
  unwind_protect
    hax = newplot (hax);
    limits = limiti;
    num = numel (y);
    dn = round (num/10);

    hl = plot (x(1), y(1), "color", "r", "marker", "none",
               x(1), y(1), "color", "g", "marker", "none",
               x(1), y(1), "color", "b", "marker", "o");
    axis (limits);  # set manual limits to speed up plotting

    ## Initialize the timer
    t = p;
    timerid = tic ();

    grid on
    
    for n = 2:(num+dn)
      m = n - dn;
      m = max ([m, 1]);
      k = min ([n, num]);
      set (hl(1), "xdata", x(1:m), "ydata", y(1:m));
      set (hl(2), "xdata", x(m:k), "ydata", y(m:k));
      set (hl(3), "xdata", x(k),   "ydata", y(k));
      l1 = line([0,x(k)],[0,y(k)],"linestyle","-","color","k");
      pause (t - toc (timerid));
      delete(l1)
      t += p;
    endfor

  unwind_protect_cleanup
    if (! isempty (oldfig))
      set (0, "currentfigure", oldfig);
    endif
  end_unwind_protect

endfunction

