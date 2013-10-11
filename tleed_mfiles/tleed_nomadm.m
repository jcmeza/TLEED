%This is function which returns R-factor of a given struture in the TLEED 
%calucation (specifically for Ni(001)-Li-(5x5), which contains 14 inequivalent atoms). 
%
% call:	[fx, cx, gradfx, gradcx]=tleed_nomadm(x,p)
%	x is a vector variable which indicates the position of the atoms;
%	p is a cell array (contains integers) which indicates the identity
%	  of the atom (Ni:1 or Li:2)

%function [fx, cx, gradfx, gradcx]=tleed_nomadm(x,p)
function fx=tleed_nomadm(x,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some parameters for the tleed optimization problem
% This part could later go to tleed_nomadm__Param.m file
	DELTA=0.4;
	DIR=0;
	RANK=0;
%	problem_dir=getenv('PWD'); %this one doesn't work, it always get the 
%       directory where the matlab invoked.
	problem_dir=mfilename('fullpath');
%	problem_dir=which(mfilename); % either one works!
	problem_dir=fileparts(problem_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	NMAX=length(p); % NMAX is the # of inequivalent atoms in the surface
	NDIM=length(x)/NMAX; %NDIM=3, the dimensionality of the position valuable

	PARM=reshape(x, NMAX,NDIM);
%change the cell arry to an array
	NTYPE=[];
	for i=1:NMAX
		NTYPE=[NTYPE, p{i}];
	end

	MINB=PARM-DELTA*ones(NMAX,NDIM);
	MAXB=PARM-DELTA*ones(NMAX,NDIM);

%	disp('optimization parameters passed to evaluate.f')
%	disp(NTYPE)
%	disp(PARM)

%	y=[PARM, NTYPE'];
%       fid = fopen('exp.txt','r');
%       yp=fscanf(fid, '%f');
%       fclose(fid);
%	yp=reshape(yp, 4,14)';
%	disp(yp)
%	disp('difference from that of previous parameters')
%	disp(y-yp)
%keyboard
	fx=GPStleed1_wg(problem_dir,DIR,RANK,PARM,MINB,MAXB,NTYPE);
	rfactortleed=fx
%pause
%	parmt=read_f('work000/gleedo000')
%	parmg=read_f('work000/gleedo000')
%	disp('difference:parmg(after tleed2)-parmt(before tleed2):')
%	disp(parmg-parmt)
%pause
% write input optimization parameters into a file
%        fid = fopen('exp.txt','w');
%        fprintf(fid,'%12.4f  %12.4f %12.4f %4d\n',y');
%        fclose(fid);
