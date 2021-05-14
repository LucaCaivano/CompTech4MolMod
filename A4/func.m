function [E, grad] = func(params, varargin)

	global iter;
	iter
	params = [params';
			  [  111.3281;
				  120.0138;
				  1.0956;
				  111.3255;
				  -120.009;
				  1.0957;
				  1.0955;
				  111.3244;
				  1.0956;
				  111.3302;
				  180.0395;
				  1.0957;
				  111.3294;
				  120.0185;
				  1.0956;
				  111.3274;
				  -120.0096;
			  ]]'

	command = sprintf('./set_params.sh %s %3.3d', num2str(params, "%17.17f "), iter);
 	system(command);

	disp("Calling Gaussian ...")
	command = sprintf('g09 c2h6_%3.3d.com', iter);
	system(command);
	disp("Done!")

	[a,b] = system(sprintf('egrep "SCF Done" c2h6_%3.3d.log',iter));
	es = regexp(b,'=[ ]+([-+.\d]+)','tokens');
	E = str2double(es{1});

	if (nargout > 1)
	 disp("Using gradients ...")
	 file = fopen(sprintf('c2h6_%3.3d.log',iter), "r");
	 [file_str, count, err_msg] = fscanf(file, '%s');
	 es = regexp(file_str,'[a-z][0-9][-0-9.]+-DE\/DX=([-0-9]+.[0-9]+)','tokens');
	 es_new = cell2mat(es);
	 grad = - str2double(es_new)';
	 grad = grad(1);
	endif

	iter += 1;

endfunction
