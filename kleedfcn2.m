function fx = kleedfcn2(xp)
% kleedfcn matlab/Fortran executable interface
% This matlab function takes care of writing the parameters to 
% the necessary files for KLEED to read and do one run. 
% 
% Assuming KLEEDFCN runs to completion, the output should be stored
% in a file called 'kleedfval.dat', which this function reads
% and passes back to the calling routine

kleed_err = NaN;
invalid_kleed = 1.6;
%invalid_kleed = fbound;

% open tleed input file
fileID = fopen('kleedinputs.dat','w');

n = length(xp);
parms = xp(1:14);
parms = round(parms); % need to make sure parms are integers
x = xp(15:n);

% write x and p to input file
fprintf(fileID,'%4d ', parms);
fprintf(fileID,'\n');
fprintf(fileID,'%14.10f ', x);
fprintf(fileID,'\n');
fclose(fileID);

% call kleed function fortran executable as a black box
% Assumption: all communication through file I/O
% KLEEDFCN will read (x,p) from kleedinputs.dat
% KLEEDFCN will write fval to   kleedfval.dat

system('./kleedfcn.exe');

% read output file and set fx
% return error if file can't be opened
[fileID, message] = fopen('kleedfval.dat','r');


if (message == "") % file exists so presumably KLEED ran to completion
    kleedfval = fscanf(fileID,'%f');
    if (kleedfval < invalid_kleed) % make sure not an invalid structure
        fx = kleedfval; 
    else
        fx = kleed_err; % invalid structure so return NaN
    end
else % trouble reading file - return Nan
    fx = kleed_err; 
end
fclose(fileID);
end

