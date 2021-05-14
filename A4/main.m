clear

system("mkdir results");

folder_name = datestr(clock());
folder_name = sprintf('results/"%s"', folder_name);
folder_create = sprintf('mkdir %s', folder_name);
system(folder_create);

global iter
iter = 1;

params = [
  1.8922; %1.5328;
%2.6500   2.5918   2.5274   2.4589   2.3865   2.3102
%2.3102   2.2342   2.1510   2.0649   1.9778   1.8922
%1.8922   1.8094   1.7384   1.6790   1.6330
  ]';

% params = [
%   1.5;
%   109.;
%   120.;
%   1.03;
%   110.;
%   -120.;
%   1.05;
%   1.1;
%   109.;
%   1.03;
%   109.;
%   180.;
%   1.05;
%   109.;
%   120.;
%   1.1;
%   109.;
%   -120.;
%   ]';

noise = 0.5;
rnd_params = params + noise * randn(1, numel(params));

%% option 1 : fminunc
% [x, fval, info, output, grad, hess] = ...
%fminunc (@func, params, optimset ("GradObj", "on", "Display", "iter"));

[x,err,iter_desc, trace_desc] = descent(@func,params,1e-3,5,42,0.1);

folder_move = sprintf('mv *.com *.log %s', folder_name);
system(folder_move);

disp("Done")
