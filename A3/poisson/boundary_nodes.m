function ret = boundary_nodes (nx, ny, idx)

  N  = nx*ny;

  [ix, iy, kk] = ...
  indices (nx, ny);

  ret = [];
  for ii = 1 : numel (idx)
    switch idx(ii)
      case 1
        ret = [ret; kk(iy == 1)];
      case 2
        ret = [ret; kk(iy == ny)];
      case 3
        ret = [ret; kk(ix == 1)];
      case 4
        ret = [ret; kk(ix == nx)];
      case 5
        ret = [ret; kk(iy == 2)];
      case 6
        ret = [ret; kk(iy == (ny-1))];
      case 7
        ret = [ret; kk(ix == 2)];
      case 8
        ret = [ret; kk(ix == (nx-1))];
    endswitch
  endfor

endfunction
