folder_name = fullfile ('results', datestr(clock(), "yyyy-mm-dd_HH:MM:SS"));
mkdir (folder_name);

func ([], [], [], 'resetnumeval');

params = ...
[
 2.5;
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

angles = !(lengths);

%noise = 0.5;
%rnd_params = params + noise * randn(1, numel(params));

[x,err,iter_desc, trace_desc, H] = ...
descent (@(x) func(x, lengths, angles), params, 1e-3, 5, 3, eye (numel (params)));

[x,err,iter_desc, trace_desc, H] = ...
descent (@(x) func(x, lengths, angles), x, 1e-4, 50, 2, H);

folder_move = sprintf ('mv *.com *.log %s', folder_name);
system (folder_move);

nf = func ([], [], [], 'getnumeval');

disp ("Done")
printf ('total number of function evaluations : %d\n', nf)
