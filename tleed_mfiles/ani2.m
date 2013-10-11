function rfac=ani1(xfile)
n = 14
s = .002

%plot the best known solution
pn_bfp;
z = parm(:,1);
x = parm(:,2);
y = parm(:,3);
h = plot3(x,y,z,'go','MarkerSize',10);
axis([0 10 0 10 -3 3])
%axis([0 10 0 10 -3 3])
%axis square
grid on
%set(h,'EraseMode','xor','MarkerSize',18)
hold on

%keyboard
%fid=fopen('tw1','r');
fid=fopen(xfile,'r');

i=1;
ntype_in=fscanf(fid,'%d',[1,14]);
parm_in=fscanf(fid,'%f',[3,14]);
tmp_rfac=fscanf(fid,'%f',[1,1]);
%[parm_out,ntype_out]=eval_wg(parm_in',ntype);

while ~feof(fid)
   drawnow
%%draw only with the minimum rfactor reached 
%if tmp_rfac==0.2151
   rfac(i)=tmp_rfac;
   [parm_out,ntype_out]=eval_wg(parm_in',ntype_in);
%   parm_out=parm_in;
%keyboard
   P{i}=parm_out;
   z = parm_out(:,1);
   x = parm_out(:,2);
   y = parm_out(:,3);
   plot3(x,y,z,'r.','MarkerSize',10)
%   set(h,'XData',x,'YData',y,'zData',z)
%hold on

   i=i+1;
%end
   % read in next 15 lines
   ntype_in=fscanf(fid,'%d',[1,14]);
   parm_in=fscanf(fid,'%f',[3,14]);
   tmp_rfac=fscanf(fid,'%f',[1,1]);
end

fclose(fid)
hold on
%pn_bfp;
%z = parm(:,1);
%x = parm(:,2);
%y = parm(:,3);
%plot3(x,y,z,'go','MarkerSize',10)

%draw best fit with the best known solution
rfac0=rfac(1);
ind=1;
p0=P{1};
for i=2:length(P)
   if rfac(i) <= rfac0 
       rfac0=rfac(i);
       ind=i;
       p0=P{i};
   end
end
z = p0(:,1);
x = p0(:,2);
y = p0(:,3);
plot3(x,y,z,'bx','MarkerSize',10)
disp(rfac(ind))

%%draw best fit with the best known solution
%delta0=(P{1}-parm).*(P{1}-parm);
%rfac0=rfac(1);
%ind=1;
%p0=P{1};
%for i=2:length(P)
%   delta=(P{i}-parm).*(P{i}-parm);
%   if (delta <delta0 & rfac(i) <= rfac0) 
%       delta0=delta;
%       rfac0=rfac(i);
%       ind=i;
%       p0=P{i};
%   end
%end
%z = p0(:,1);
%x = p0(:,2);
%y = p0(:,3);
%plot3(x,y,z,'bx','MarkerSize',10)
%disp(rfac(ind))
figure
plot(1:length(rfac),rfac,'bd')
