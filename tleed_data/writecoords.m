function [count] = writecoords(kind,zxy)
%WRITECOORDS Summary of this function goes here
%   Detailed explanation goes here
count = 0;

fileID = fopen('tleedx0.dat','w');
count = length(kind)

for i=1:count
    fprintf(fileID, '%d %12.8f %12.8f %12.8f\n', kind(i), zxy(i,1), zxy(i,2), zxy(i,3));
end

end

