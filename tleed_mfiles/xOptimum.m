%this function update position parameter with local optimum.
% but it is not able to handle the case when the number of 
%inequivalent atoms in the layer changes from 58 to some other number, e.g., 56 
%call:	x=xOptimum(xfile)
function [parm,intype]=xOptimum
        problem_dir=mfilename('fullpath');
        problem_dir=fileparts(problem_dir);

        xfile=strcat(problem_dir,'/work000/searchs000');
%keyboard 
	ntot=58;
	fid0=fopen(xfile,'r');
	%skip 3 lines
	for i=1:3
		tmp=fgetl(fid0);
	end
	Fparm=fscanf(fid0,'%f',[3,ntot])';

%	%skip 11 lines
	for i=1:11
		tmp=fgetl(fid0);
	end
	%read in 58 lines
	Dparm=fscanf(fid0,'%f',[4,ntot])';
	fclose(fid0);	
	Fparm=Fparm+Dparm(:,1:3);
	parm=[];
	inc=0;
	for i=1:ntot
        	if Fparm(i,2) >=0.000 && Fparm(i,3) >=0.000 && Fparm(i,3)<=Fparm(i,2)
			inc=inc+1;
                	parm=[parm;Fparm(i,:)];
 			intype(inc)=Dparm(i,4);
        	end
	end
	%disp(inc)
	if inc~=14
		disp('error in input xfile' )
		parm=100*ones(14,3);
	else
		[parm,I]=sortrows(parm);
		intype=intype(I);
	end

