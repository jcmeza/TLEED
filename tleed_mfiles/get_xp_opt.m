%call:	 [x,p]=get_xp_opt(outputfile)
function [x,p]=get_xp_opt(outputfile)

	ntot=58;
	fid0=fopen('work000/tleed5i000','r');
	%skip 195 lines
	for i=1:194
		tmp=fgetl(fid0);
	end
	Fntype=fscanf(fid0,'%d\n',[1,ntot]); 
%	%note new line '...\n' was needed! (fscanf and fgetl are different 
%	%on the rule of dealing with \n default after the operations!) 
	%skip 1 lines
	tmp=fgetl(fid0);
	%read in 58 lines
	Fparmi=fscanf(fid0,'%f',[3,ntot])';
	fclose(fid0);	

	%fid0=fopen('work000/gleedo000','r');
	fid0=fopen('tmp.txt2','r');
%	%skip 6 lines
	for i=1:6
		tmp=fgetl(fid0);
	end
	%read in 58 lines
	Fparm=fscanf(fid0,'%f',[3,ntot])';
disp(Fparm-Fparmi)
	fclose(fid0);	

	parm=[];
	inc=0;
	for i=1:ntot
        	if Fparm(i,2) >=0.000 && Fparm(i,3) >=0.000 && Fparm(i,3)<=Fparm(i,2)
			inc=inc+1;
                	parm=[parm;Fparm(i,:)];
			ntype{inc}=Fntype(i);
			ntype0(inc)=Fntype(i);
        	end
	end
	%disp(inc)
	if inc~=14
		disp('error in gleedo000')
		p={0 0 0 0 0 0 0 0 0 0 0 0 0 0};
		x=100*ones(42,1);
	else
		p=ntype;
		x=reshape(parm, inc*3,1);
		disp(ntype0)
		disp(parm)
	end
