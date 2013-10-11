function rfac=ani1(xfile)
n = 14
s = .002

%plot the best known solution
pn_bfp;
z = -parm(:,1);
x = parm(:,2);
y = parm(:,3);
h = plot3(x,y,z,'r.');
axis([-10 10 -10 10 -3 3])
%axis([0 10 0 10 -3 3])
%axis square
grid on
set(h,'EraseMode','xor','MarkerSize',18)

%keyboard
fid=fopen(xfile,'r');

i=1;
parm_in=fscanf(fid,'%f',[3,14]);
tmp_rfac=fscanf(fid,'%f',[1,1]);
%[parm_out,ntype_out]=eval_wg(parm_in',ntype);

while ~feof(fid)
   drawnow
%draw only with the minimum rfactor reached 
if tmp_rfac==0.2151
   rfac(i)=tmp_rfac;
   [parm_out,ntype_out]=eval_wg(parm_in',ntype);
   P{i}=parm_out;
   z = -parm_out(:,1);
   x = parm_out(:,2);
   y = parm_out(:,3);
   set(h,'XData',x,'YData',y,'zData',z)

   i=i+1
end
   % read in next 15 lines
   parm_in=fscanf(fid,'%f',[3,14]);
   tmp_rfac=fscanf(fid,'%f',[1,1]);
end

fclose(fid)
hold on
pn_bfp;
z = -parm(:,1);
x = parm(:,2);
y = parm(:,3);
plot3(x,y,z,'go','MarkerSize',10)

%draw best fit with the best known solution
delta0=(P{1}-parm).*(P{1}-parm);
p0=P{1};
for i=2:length(P)
   delta=(P{i}-parm).*(P{i}-parm);
   if delta <delta0 
       delta0=delta;
       p0=P{i};
   end
end
z = -p0(:,1);
x = p0(:,2);
y = p0(:,3);
plot3(x,y,z,'bx','MarkerSize',10)
%figure
%plot(1:length(rfac),rfac,'bd')
