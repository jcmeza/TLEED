function ani1(xfile)
fid=fopen(xfile,'r');
i=1;
k=1;
ntype_in=fscanf(fid,'%d',[1,14]);
parm_in=fscanf(fid,'%f',[3,14]);
tmp_rfac=fscanf(fid,'%f',[1,1]);
while ~feof(fid)
%k=k+1
%if k==1070
%keyboard
%end
if size(parm_in) ~=[3,14] |size(ntype_in)~=[1,14]
sprintf('%s%f\n','line number=', 16*k)
return
end

if tmp_rfac <1.6000
   [parm_out,ntype_out]=eval_wg(parm_in',ntype_in);
%disp(parm_out)
%disp(ntype_out)
%keyboard
%if ntype_out(1)~=0
   y=[parm_out,ntype_out'];
   outfile=strcat('DATA/run5Step',num2str(i),'.dat');

   fidf=fopen(outfile,'w');
   fprintf(fidf, '%10.4f %10.4f %10.4f  %3d\n', y')  
   fclose(fidf);
   i=i+1;
end

   % read in next 15 lines
   ntype_in=fscanf(fid,'%d',[1,14]);
   parm_in=fscanf(fid,'%f',[3,14]);
   tmp_rfac=fscanf(fid,'%f',[1,1]);
end

fclose(fid)
