#this scrip backup problem files
for file in tleed_nomadm_x0.m tleed_nomadm_Omega.m tleed_nomadm_N.m
do
	cp -p $1/$file $1/$file.$2
done

