function fx = tleedfcn(x, parms)
% kleedfcn matlab/Fortran executable interface
% This matlab function takes care of writing the parameters to 
% the necessary files for TLEED to read and do one run. 
% 
% Assuming TLEEDFCN runs to completion, the output should be stored
% in a file called 'tleedfval.dat', which this function reads
% and passes back to the calling routine

tleed_err = 1.6;

% open tleed input file
fileID = fopen('tleedinputs.dat','w');

% write x and p to input file
fprintf(fileID,'%4d ', parms);
fprintf(fileID,'\n');
fprintf(fileID,'%14.10f ', x);
fprintf(fileID,'\n');
fclose(fileID);

% call kleed function fortran executable as a black box
% Assumption: all communication through file I/O
% TLEEDFCN will read (x,p) from tleedinputs.dat
% TLEEDFCN will write fval to   tleedfval.dat

system('./tleedfcn.exe');

% read output file and set fx
% return error if file can't be opened
[fileID, message] = fopen('tleedfval.dat','r');

if (message == "") % file exists so presumably KLEED worked
    tleedfval = fscanf(fileID,'%f');
else % trouble reading file - return error value
    tleedfval = tleed_err;
end
fclose(fileID);

fx = tleedfval;

end
