function [CoordSub] = setupCoordSub()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nsub = 6;
ndim = 3;

CoordSub = zeros(nsub,3);
CoordSub0 = zeros(nsub,3);

data = [0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.3,0.5,0.3,0.5,0.5,0.1,0.1,0.1,0.3,0.3,0.5];
CoordSub0 = reshape(data,nsub,ndim);
AA = 12.45;

for i=1:nsub
    CoordSub(i,1) = 3.5214;
    CoordSub(i,2) = AA*CoordSub0(i,2);
    CoordSub(i,3) = AA*CoordSub0(i,3);
end

end

