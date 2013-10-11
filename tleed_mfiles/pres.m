y=parm;
for i=1:100
hold on;
	plot3(y(:,3),y(:,1),y(:,2),'o')
	y=rand(14,3)*0.01+y;
%keyboard
	t = timer('TimerFcn',@mycallback, 'Period', 10000.0);
	wait(t)
end

