folder_name = fullfile ('results', datestr(clock(), "yyyy-mm-dd_HH:MM:SS"));
mkdir (folder_name);

func ([], [], [], 'resetnumeval');

params = ...
[
 1.0;
 109.;
 120.;
 1.03;
 110.;
 -120.;
 1.05;
 1.1;
 109.;
 1.03;
 109.;
 180.;
 1.05;
 109.;
 120.;
 1.1;
 109.;
 -120.
];


%params per ottenere distanza negativa e verificare i vincoli
%params = ...
%[
% 1.0;
%179.;
% 120.;
 %1.03;
 %110.;
 %-179.;
 %1.05;
 %1.1;
 %109.;
 %1.03;
 %109.;
 %180.;
 %1.05;
 %109.;
 %120.;
 %1.1;
 %179.;
 %-120.
%];

lengths = ...
logical ([
          1 %d1
          0 
          0
          1 %d2
          0
          0
          1 %d3
          1 %d4
          0
          1 %d5
          0
          0
          1 %d6
          0
          0
          1 %d7
          0
          0
          ]);

angles = !lengths;

%M = zeros(numel(params),1);
m = zeros(numel(params),1);
m=m+lengths;
M = [
          Inf %d1
          0 
          0
          Inf %d2
          0
          0
          Inf %d3
          Inf %d4
          0
          Inf %d5
          0
          0
          Inf %d6
          0
          0
          Inf %d7
          0
          0
          ];


angles2D = ...
logical([
          0 
          1 %a1 
          0
          0 
          1 %a2
          0
          0
          0
          1 %a3
          0
          1 %a4
          0
          0
          1 %a5
          0
          0
          1 %a6
          0	 
	 ]);

M = M + 180*angles2D;

angles3D = !(lengths + angles2D);

m = m - 180*angles3D;
M = M + 180*angles3D;


%noise = 0.5;
%rnd_params = params + noise * randn(1, numel(params));

[x,err,iter_desc, trace_desc, H] = ...
descent (@(x) func(x, lengths, angles), params, 1e-3, 5, 3, m, M, eye (numel (params)));

[x,err,iter_desc, trace_desc, H] = ...
descent (@(x) func(x, lengths, angles), x, 1e-4, 20, 2, m, M, H);



folder_move = sprintf ('mv *.com *.log %s', folder_name);
system (folder_move);

nf = func ([], [], [], 'getnumeval');

disp ("Done")
printf ('total number of function evaluations : %d\n', nf)
