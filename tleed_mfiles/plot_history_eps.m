%this function plots performance history of NOMADm for tleed problem
%from the logfile of the matlab
%
%call: plot_history(logfile)
function plot_history(logfile)
if exist('tmp.logfile')==2
	!rm tmp.logfile
end
if exist('tmp.rfac')==2
	!rm tmp.rfac
end

%AGL the original file does not have quotes for logfile
%copyfile(logfile, 'tmp.logfile');
copyfile('logfile', 'tmp.logfile');
%AGL the original file does not have ./ in front of get_rfac.awk
%!get_rfac.awk tmp.logfile >tmp.rfac
%!./get_rfac.awk tmp.logfile >tmp.rfac
!./get_TRUTH.awk tmp.logfile >tmp.rfac

fid=fopen('tmp.rfac','r');
rfac=fscanf(fid,'%f', [1,inf]);
fclose(fid);

n=length(rfac)
bfp(1)=rfac(1);
for i=2:n
	if rfac(i) < bfp(i-1)
		bfp(i)=rfac(i);
	else
		bfp(i)=bfp(i-1); 
	end
end

temp=sort(rfac);
for i=1:n
	sorted_rfac(i)=temp(n-i+1);
end

temp=0;
for i=1:n
	if rfac(i)<1.6
		temp=temp+1;
	end
	percent(i)=temp/i;
end

figure
%AGL Set the visible property of this figure to off
%f=figure('vis','off'); 
plot(1:n, bfp,'md-')
title('TLEED-NOMADm: {\bf Real-time performance}')
xlabel('# of function calls')
ylabel('R-factor (best feasible point)')
%text(n,bfp(n),num2str(bfp(n)))
%text(3*n/5,0.6,strcat('R-factor_{min} =  ',num2str(bfp(n))))
%text(3*n/5,0.55,strcat('# of function calls =',num2str(n)))
text(3*n/5,0.5*(bfp(1)+bfp(n)),strcat('R-factor_{min} =  ',num2str(bfp(n))))
text(3*n/5,0.45*(bfp(1)+bfp(n)),strcat('# of function calls =',num2str(n)))
%AGL Save figure as test.fig
%saveas(f,'test.fig')
print -deps2 plot1.eps

figure
plot(1:n, rfac,'bo')
title('TLEED-NOMADm: {\bf Real-time history}')
xlabel('i-th function call')
ylabel('R-factor')
%text(n-100,0.6+rfac(n),num2str(rfac(n)))
text(3*n/5,0.6,strcat('R-factor =  ',num2str(rfac(n))))
text(3*n/5,0.55,strcat(num2str(n),'-th function call'))
print -deps2 plot2.eps

figure
plot(1:n, sorted_rfac,'rx-')
title('TLEED-NOMADm: {\bf Sorted R-factor in descending order for trial points}')
xlabel('index for trial function call')
ylabel('R-factor')
text(n-100,0.1+sorted_rfac(n),num2str(sorted_rfac(n)))
print -deps2 plot3.eps

figure
plot(1:n, percent,'c+-')
title('TLEED-NOMADm: {\bf percentage of the valid function calls}')
xlabel('i-th function call')
ylabel('R-factor')
%text(n-100,0.1+sorted_rfac(n),num2str(sorted_rfac(n)))
print -deps2 plot4.eps

