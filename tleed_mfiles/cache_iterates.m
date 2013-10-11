%this function adds an iterate on the top of the backup file of iterates 
%in case want to use it as an initial point later. 
%this function is specially written for Ni(001)-(5x5)-Li,shoud be used with caution. 
%
%call :  xp_added=cache_iterates(backupfile)
function xp_added=cache_iterates(backupfile)

	ntot=58;
	fid0=fopen('work000/tleed5i000','r');
	%skip 195 lines
	for i=1:194
		tmp=fgetl(fid0);
	end
	Fntype=fscanf(fid0,'%d\n',[1,ntot]); 
	%note new line '...\n' was needed! (fscanf and fgetl are different 
	%on the rule of dealing with \n default after the operations!) 
	%skip 1 lines
	tmp=fgetl(fid0);
	%read in 58 lines
	Fparm=fscanf(fid0,'%f',[3,ntot])';
	fclose(fid0);	
	%xp_added=[Fntype',Fparm]

	parm=[];
	inc=0;
	for i=1:ntot
        	if Fparm(i,2) >=0.000 && Fparm(i,3) >=0.000 && Fparm(i,3)<=Fparm(i,2)
			inc=inc+1;
                	parm=[parm;Fparm(i,:)];
			ntype{inc}=Fntype(i);
        	end
	end
	%disp(inc)

	iterate.x=reshape(parm, inc*3,1);
	iterate.p=ntype;

%	fid0=fopen('work000/trace000','r');
%	%skip 195 lines
%	for i=1:48
%		tmp=fgetl(fid0);
%	end
%
%	rfactor=fscanf(fid0,'%f\n');
%	fclose(fid0);


	fidt=fopen('ctemp.dat','w');
	fprintf(fidt,'%s','iterate(1).x=[');
	fprintf(fidt,'%f  ',iterate.x);
	fprintf(fidt,'%s\n',']'';');
	fprintf(fidt,'%s','iterate(1).p={'); 
	fprintf(fidt, '%d  ', ntype{1:end});
	%fprintf(fidt,'%s%f\n','}; %rfactor=',rfactor);
	fprintf(fidt,'%s\n','};');
	fclose(fidt);
	Fadd=textread('ctemp.dat','%s','delimiter','\n','whitespace','');

	%if sw==1
	%	%[parmg, ntypeg]=read_f('work000/gleedo000')
	%	[parmg, ntypeg]=read_f
	%	iterate.x=reshape(parmg,ntot,1);
	%end

	%read file into a cell array F {# of line,1} 
	F=textread(file,'%s','delimiter','\n','whitespace','');
	%write to a initial _x0 file
	fid = fopen(file,'w');
	fprintf(fid,'%s\n',F{1:2});
	fprintf(fid, '%s\n',Fadd{1:end});
	%for i=3:length(F)
	%	F{i,1}=strcat('%',F{i,1});
	%end
	fprintf(fid,'%s\n',F{3:end});
	fclose(fid);

	%display the iterate added on top of the cache files.
	p=[];
	for i=1:length(ntype)
		p=[p;ntype{i}];
	end
	xp_added=[parm,p]
	delete('ctemp.dat');
