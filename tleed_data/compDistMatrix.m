function [DistMatrix, DistMatrix2, zxy, check1, check2] = compDistMatrix(filename,parm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

startrow = 4*(parm-1)+2;
endrow = startrow + 2;
tleedparms = importparms(filename,startrow,endrow);
zxy = tleedparms'; % should be 14 x 3

n = length(zxy);

DistMatrix = zeros(n,n);
for i=1:n
    for j=1:i
        DistMatrix(i,j) = norm(zxy(i,:)-zxy(j,:));
    end
end

check1 = min(min(DistMatrix(DistMatrix>0),[],2));

nsub = 6;
DistMatrix2 = zeros(nsub,n);
CoordSub    = setupCoordSub(); % should be 6 x 3
for i=1:nsub
    for j=1:n
        DistMatrix2(i,j) = norm(CoordSub(i,:)-zxy(j,:));
    end
end

check2 = min(min(DistMatrix2,[],2));

end

