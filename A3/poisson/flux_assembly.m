function  A = flux_assembly (N, ii, nx, ny, hastop, hasbot,
                             haslft, hasrgt, fltop,
                             flbot, fllft, flrgt);

  iv = [ii(hasrgt);      ii(haslft);      ii(hastop);     ii(hasbot)];
  jv = [ii(hasrgt)+ny;   ii(haslft)-ny;   ii(hastop)-1;   ii(hasbot)+1];
  vv = [-flrgt(hasrgt); -fllft(haslft);  -fltop(hastop); -flbot(hasbot)];

  A = sparse (iv, jv, vv, N, N);
  d = sum (A, 1);
  A = A - diag (d);

endfunction
