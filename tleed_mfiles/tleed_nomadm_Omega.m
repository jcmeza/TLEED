%set feasible region for tleed_nomadm problem. 
% n is the # of the continuous variables, n=42 for tleed_nomadm
function [A,l,u,plist] = tleed_nomadm_Omega(n);
A = eye(n);
%something is wroing in this boundary from the pgaleed1.C

  parmz(0+1)=  -1.8757;
  parmz(1+1)=  -1.7941;
  parmz(2+1)=  -1.8067;
  parmz(3+1)=  -0.3861;
  parmz(4+1)=   0.2472;
  parmz(5+1)=  -0.0461;
  parmz(6+1)=   0.0690;
  parmz(7+1)=   0.1874;
  parmz(8+1)=   1.7112;
  parmz(9+1)=   1.7350;
  parmz(10+1)=  1.7378;
  parmz(11+1)=  1.7467;
  parmz(12+1)=  1.7751;
  parmz(13+1)=  1.7897;
                                                                                              
  parmx(0+1)=   0.0000;
  parmx(1+1)=   3.1141;
  parmx(2+1)=   3.0047;
  parmx(3+1)=   6.2250;
  parmx(4+1)=  -4.0621;
  parmx(5+1)=   1.2552;
  parmx(6+1)=   3.6738;
  parmx(7+1)=  -4.2907;
  parmx(8+1)=   5.0398;
  parmx(9+1)=   0.0000;
  parmx(10+1)=  5.0355;
  parmx(11+1)=  2.4703;
  parmx(12+1)=  2.5445;
  parmx(13+1)=  2.4371;
                                                                                             
  parmy(0+1)=   0.0000;
  parmy(1+1)=   0.0000;
  parmy(2+1)=  -3.0047;
  parmy(3+1)=   1.2913;
  parmy(4+1)=   6.2250;
  parmy(5+1)=   1.2552;
  parmy(6+1)=   1.2125;
  parmy(7+1)=   3.7093;
  parmy(8+1)=   0.0000;
  parmy(9+1)=   0.0000;
  parmy(10+1)=  5.0355;
  parmy(11+1)=  5.0402;
  parmy(12+1)=  0.0000;
  parmy(13+1)=  2.4371;

bfp=[parmz';parmx';parmy'];
l=bfp-0.1;
u=bfp+0.1;

natom=n/3;
for i=1:natom
	plist{i} = {1,2};
end
return



