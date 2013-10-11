%call:	 [parm_out,ntype_out,rfactor]=tleed_eval(parm_in,ntype_in,rfac_flag)
function [parm_out,ntype_out,rfactor]=tleed_eval(parm_in,ntype_in,rfac_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DELTA=0.4;
	DIR=0;
	RANK=0;
	NMAX=14; 
	NDIM=3; 
	problem_dir=mfilename('fullpath');
	problem_dir=fileparts(problem_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	MINB=parm_in-DELTA*ones(NMAX,NDIM);
	MAXB=parm_in-DELTA*ones(NMAX,NDIM);

	if rfac_flag==0
		[parm_out,ntype_out]=eval_wg(parm_in,ntype_in);
		rfactor=0;
	else
		[parm_out,ntype_out]=eval_wg(parm_in,ntype_in);
		rfactor=GPStleed1_wg(problem_dir,DIR,RANK,parm_in,MINB,MAXB,ntype_in);
	end
