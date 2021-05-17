
function [E, grad] = func (params, lengths, angles, option = [])

  persistent numeval = 0;

  if ! isempty (option)
    if strcmpi (option, 'resetnumeval')
      numeval = 0;
    elseif strcmpi (option, 'getnumeval')
      E = numeval; grad = 0;
    endif
    return
  endif
  ++numeval;
  
  bohr_in_aa = .5291777721092;
  rad_in_deg = 360 / 2 / pi;
  
  command = sprintf('./set_params.sh %s %3.3d',
                    num2str (params', "%17.17e "), numeval);
  system(command);

  disp("Calling Gaussian ...")
  command = sprintf('g09 c2h6_%3.3d.com', numeval);
  system(command);
  disp("Done!")
  
  [a,b] = system(sprintf('egrep "SCF Done" c2h6_%3.3d.log',numeval));
  es = regexp(b,'=[ ]+([-+.\d]+)','tokens');
  E = str2double(es{1});

  if (nargout > 1)
    disp("Using gradients ...")
    file = fopen(sprintf('c2h6_%3.3d.log',numeval), "r");
    [file_str, count, err_msg] = fscanf(file, '%s');
    es = regexp(file_str,'[a-z][0-9][-0-9.]+-DE\/DX=([-0-9]+.[0-9]+)','tokens');
    es_new = cell2mat(es);
    grad = - str2double(es_new)';
    grad(lengths) /=  bohr_in_aa;
    grad(angles)  /=  rad_in_deg;
  endif

endfunction
