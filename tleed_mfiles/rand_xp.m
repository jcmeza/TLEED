%call:  iterates=rand_xp(n)
function iterates=rand_xp(n)
pn_bfp;
fid=fopen('x0.m','w')
fprintf(fid,'%s\n','function iterate=tleed_nomadm_x0');
for i=1:n
	%nt=ntype;
	x=parm+randn(14,3)*0.1;
	nt=mod(floor(rand(1,14)*10), 2)+1;
	fprintf(fid,'%s%d%s','iterate(',i,').x=[');
	fprintf(fid,'%16.10f', x);
	fprintf(fid,'%s\n',']'';');
	fprintf(fid,'%s%d%s','iterate(',i,').p={');
	fprintf(fid,'%3d',nt);
	fprintf(fid,'%s\n','};');
end
fclose(fid)
