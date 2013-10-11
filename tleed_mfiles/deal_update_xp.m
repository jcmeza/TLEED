function [xout,pout,isValid]=deal_update_xp(x,p)
	isValid=0;

        problem_dir=mfilename('fullpath');
        problem_dir=fileparts(problem_dir);

	NMAX=length(p); 
	NDIM=length(x)/NMAX; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	PARM=reshape(x, NMAX,NDIM);
	NTYPE=[];
	for i=1:NMAX
		NTYPE=[NTYPE, p{i}];
	end
	[parm_out,ntype_out]=eval_update_wg(problem_dir, PARM,NTYPE);
	xout=reshape(parm_out, 3*NMAX,1);
	for i=1:NMAX
		pout{i}=ntype_out(i);
	end
	
	if ntype_out(1) ~=0
		isValid=1;
	end	
