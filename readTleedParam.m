%===============================================================================
% tleed_Param:  Parameter file for the tleed_nomadm problem
% ------------------------------------------------------------------------------
% Calls:  setWork
% VARIABLES:
%  Param        = structure of parameters set by the Parameter file
%    .xFileType =   type of file used to test validity (f=fortran, m=matlab)
%    .vData     =   substructure of data used in the validity check
%    
%    .p_dir     =   directory where the problem files are located
%    .nAtoms    =   number of atoms (14)
%    .nDim      =   dimension of the atom position variable (3)
%    .iterate0  =   initial iterate(s)
%      .x       =     continuous  variable values
%      .p       =     categorical variable values
%    .A         =   matrix of linear/bound constraint coefficient
%    .l         =   vector of linear/bound constraint lower bounds 
%    .u         =   vector of linear/bound constraint upper bounds 
%    .plist     =   cell array of lists of allowed values for each .p variable
%    .vData     =   substructure of validity data from KLEED
%  iX0          = index into initial points vector used to establish bounds
%  jX0          = indices of initial points used in optimization
%  M            = matrix of raw initial point data read from file
%  X0           = vector of initial iterates parsed from M
%  point        = temporary storage of an initial iterate
%    .x         =   continuous  variable values
%    .p         =   categorical variable values
%  n            = number of continuous variables in the optimization problem
%===============================================================================
function Param = readTleedParam

 % Set iX0 = scalar index of initial point about which bounds are constructed
 % Set jX0 = scalar/vector of indices of initial points to use in optimization
 iX0 = 1;
 jX0 = 1:10;

 % Data taken from KLEED
 Param.xFileType = 'm';
 CoorSub0 = [0,0,0,0,0,0; 0.1,0.3,0.5,0.3,0.5,0.5; 0.1,0.1,0.1,0.3,0.3,0.5]';
 vData.nSub     = size(CoorSub0,1);
 vData.validN   = 58;
 vData.z0       =  3.5214;
 vData.xUnit    = 12.45;
 vData.dSpace   =  1.90;
 vData.rmin     = [1.0; 1.0; 1.0; 1.0; 1.0];
 vData.rmax     = [1.4; 1.7; 1.5; 1.5; 1.5];
 vData.codes    = {'C','A','D','B','E','F','G'};  % do NOT edit!
 vData.N        = [2,1,4,1,4,4,8];                % do NOT edit!
 vData.CoorSub  = [vData.z0*ones(vData.nSub,1), vData.xUnit*CoorSub0(:,2:3)];
 vData.nTypeSub = 2*ones(vData.nSub,1);
 Param.vData    = vData;

 % Set paths and reset temporary files
 data_dir    = 'tleed_data';
 dataFile    = 'plotinitialx0.dat';
 Param.p_dir = fileparts(mfilename('fullpath'));
 Param.d_dir = [Param.p_dir, filesep, data_dir];
 setWork(Param.d_dir);

 % Read initial points from file
 M  = dlmread(fullfile(Param.d_dir,dataFile))
 X0 = [];
 for k = 1:4:size(M,1)
    point.x = [M(k+1,:),M(k+2,:),M(k+3,:)]';
    point.p = num2cell(M(k,:));
    X0      = [X0; point];  %#ok
 end
 nX0 = length(X0);
 
 % Set NOMADm initial points and feasible region parameters
 if isempty(jX0)
    jX0 = length(X0);
 end
 iX0 = min(iX0,nX0);
 Param.iterate0 = X0(jX0);
 Param.nDim     = 3;
 Param.nAtoms   = length(X0(iX0).p);
 Param.l        = X0(iX0).x - 0.2;
 Param.u        = X0(iX0).x + 0.2;
 Param.A        = eye(Param.nDim*Param.nAtoms);
 for i = 1:Param.nAtoms
 	Param.plist{i} = {1, 2};
 end
return

%===============================================================================
% setWork:  Sets up work directories and folders for running TLEED and KLEED.
% ------------------------------------------------------------------------------
% Called by:  tleed_Param
% VARIABLES:
%   data_dir    = directory where source data files are located
%   file        = structure whose fields are the working directory names
%     .twork000 =   cell array of files to copy into directory of same name
%     .kwork000 =   cell array of files to copy into directory of same name
%     .kleedIV  =   cell array of files to copy into directory of same name
%   work_dir    = cell array of current working directory to be created
%   success     = flag indicating success/failure with directory/file creation
%   message     = error message associated with directory/file creation
%   workfile    = name of current file being copied into current directory
%   source      = full name of current file being copied
%   dest        = full name of current file being created
% ==============================================================================
function setWork(data_dir)

 % Names of working directories and data files to be created/copied
 file.twork000 = {'rfac.d','exp.d','tleed4.i','tleed5.i'};
 file.kwork000 = {'rfac.d','exp.d','kleed4.i','kleed5.i'};
 file.kleedIV  = {};

 % Loop through each working directory
 for work_dir = fieldnames(file)'

    % Create working directory if is it missing
    if ~exist(work_dir{1},'dir')
       [success,message] = mkdir(work_dir{1});
       if ~success, error(message); end
    end

    % Copy the required data files into the working directory
    for workfile = file.(work_dir{1})
       source = fullfile(data_dir,   workfile{1});
       dest   = fullfile(work_dir{1},workfile{1});
       if ~exist(source,'file')
          error(['Missing source file: ',source,' not found.']);
       end
       if ~exist(dest,'file')
          [success,message] = copyfile(source,dest);
          if ~success, error(message); end
       end
    end
 end
return
