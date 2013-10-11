function xx=extend_pos(x)
ith=0;
xx=[];
for II=1:14
	if -x(II,2)~=x(II,2)
		ith=ith+1;
		tmp(ith,:)=x(II,:);
		tmp(ith,2)=-x(II,2);
	end
	if -x(II,3)~=x(II,3)
		ith=ith+1;
		tmp(ith,:)=x(II,:);
		tmp(ith,3)=-x(II,3);
	end
	if -x(II,2)~=x(II,2) & -x(II,3)~=x(II,3)
		ith=ith+1;
		tmp(ith,:)=x(II,:);
		tmp(ith,2)=-x(II,2);
		tmp(ith,3)=-x(II,3);
	end
	xx=[xx;tmp];
end

