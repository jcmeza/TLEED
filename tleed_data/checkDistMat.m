function [DistMatrix, DistMatrix2, check1, check2] = checkDistMat(zxy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nparms = 14;
n    = length(zxy);

DistMatrix = zeros(n,n);
for i=1:n
    for j=1:i
        DistMatrix(i,j) = norm(zxy(i,:)-zxy(j,:));
    end
end

check1 = min(min(DistMatrix(DistMatrix>0),[],2));

nsub = 6; % number of atoms in substrate
DistMatrix2 = zeros(nsub,n);
CoordSub    = setupCoordSub(); % should be 6 x 3
for i=1:nsub
    for j=1:n
        DistMatrix2(i,j) = norm(CoordSub(i,:)-zxy(j,:));
    end
end

check2 = min(min(DistMatrix2,[],2));

end

