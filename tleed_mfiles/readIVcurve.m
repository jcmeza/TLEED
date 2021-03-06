function[ivexp,ivth] = readIVcurve(datafile)
%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/meza/MyProjects/TLEED/kleedIV/ivexp1
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/07/14 11:32:09

%% Initialize variables.
file1 = ['/Users/meza/MyProjects/TLEED/kleedIV/ivexp',num2str(datafile)];
file2 = ['/Users/meza/MyProjects/TLEED/kleedIV/ivth',num2str(datafile)];
startRow = 3;
endRow = 617;

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%11f%15f%[^\n\r]';

%% Open the text file.
fileID1 = fopen(file1,'r');
fileID2 = fopen(file2,'r');
%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray1 = textscan(fileID1, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false);
dataArray2 = textscan(fileID2, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID1);
fclose(fileID2);
%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
ivexp = [dataArray1{1:end-1}];
ivth  = [dataArray2{1:end-1}];
%% Clear temporary variables
clearvars filename startRow endRow formatSpec fileID dataArray ans;