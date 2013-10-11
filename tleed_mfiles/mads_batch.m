function [BestF,BestI,RunStats,RunSet] = mads_batch
%MADS_BATCH  Sets up and runs the MADS algorithm without a GUI.
%
%   Syntax:
%      mads_batch
%
%   Description:
%      This function serves as a GUI-free alternative to NOMADm in setting
%      up an optimization problem, setting various algorithm parameters and
%      user options, and calling the MADS optimizer.  It first sets all of the
%      variables to their default values, which are clearly stated in the
%      MADS_DEFAULTS file.  To change a variable from its default value, the
%      user must add a statement to this file to do so.  Some variable
%      statements are included here for convenience, which can be change
%      manually.
%
%   See also MADS_DEFAULTS, MADS

%*******************************************************************************
%   Copyright (c) 2001-2006 by Mark A. Abramson
%
%   This file is part of the NOMADm software package.
%
%   NOMADm is free software; you can redistribute it and/or modify it under the
%   terms of the GNU General Public License as published by the Free Software
%   Foundation; either version 2 of the License, or (at your option) any later
%   version.
%
%   NOMADm is distributed in the hope that it will be useful, but WITHOUT ANY
%   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
%   details.
%
%   You should have received a copy of the GNU General Public License along
%   with NOMADm; if not, write to the Free Software Foundation, Inc., 
%   59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% ------------------------------------------------------------------------------
%   Originally created, 2001.
%   Last modified, 26 August 2006
%
%   Author information:
%   Mark A. Abramson, LtCol, USAF, PhD
%   Air Force Institute of Technology
%   Department of Mathematics and Statistics
%   2950 Hobson Way
%   Wright-Patterson AFB, OH 45433
%   (937) 255-3636 x4524
%   Mark.Abramson@afit.edu
%*******************************************************************************

%*******************************************************************************
% mads_batch:  Runs MADS in batch mode.
% ------------------------------------------------------------------------------
% Calls:     mads_defaults, mads, < user initial points file >
% Variables:
%  Defaults    = structure of MADS default values (see mads_defaults)
%  Options     = structure for options settings (see mads_defaults)
%  problemPath = location of user problem files
%  Problem     = structure of data for optimization problem
%  newPath     = logical indicating if path is not the Matlab path
%  iterate0    = structure of data for the initial iterate (see mads)
%  BestF       = final best feasible solution found
%  BestI       = final least infeasible solution found
%  RunStats    = structure of MADS Run statistics (see mads)
%*******************************************************************************

%*******************************************************************************
% DO NOT MODIFY THIS SECTION
%*******************************************************************************

% Set Options to their default values
clear variables
Defaults = mads_defaults('Truth');
Options  = Defaults.Options;
Problem.nameCache   = 'CACHE';
Problem.typeProblem = 'TRUTH';

%*******************************************************************************
% MODIFY ONLY AFTER THIS POINT AND BEFORE NEXT MESSAGE.
%*******************************************************************************

% Specify Problem Files
problemPath = fullfile(matlabroot,'work','NOMADm','TestProblems','ToyProblems');
Problem.File.F = 'examplePC';             % functions file
Problem.File.O = 'examplePC_Omega';       % linear constraints file
Problem.File.X = 'examplePC_X';           % closed constraints file
Problem.File.I = 'examplePC_x0';          % initial points file
Problem.File.N = 'examplePC_N';           % discrete neighbor file (MVP only)
Problem.File.P = 'examplePC_Param';       % parameter file
Problem.File.C = 'examplePC_Cache.mat';   % previously created Cache file
Problem.File.S = 'examplePC_Session.mat'; % previously created Session file
Problem.File.H = 'examplePC_History.txt'; % iteration history text file
Problem.File.D = 'examplePC_Debug.txt';   % debug log file
Problem.fType  = 'C';                     % type of functions file {M,F,C}
Problem.nc     = 0;                       % number of nonlinear constraints

% Specify Choices for SEARCH
Options.nSearches         = 2;
Options.Search(1).type    = 'LHS';       % For choices, see mads_defaults
Options.Search(1).nIter   = 1;           % Number of iterations for Search #1
Options.Search(1).nPoints = 8;           % Number of poll or sample points
Options.Search(1).sfile   = '';          % filename must include full path
Options.Search(1).file    = '';          % filename must include full path
Options.Search(1).local   = 0;           % flag to turn on trust region
Options.Search(1).merit   = 0;           % flag to penalize clustered data
Options.Search(1).recal = strncmp(Options.Search(1).type(max(1,end-1):end), ...
                                  'NW',2) || ...
                          strncmp(Options.Search(1).type(max(1,end-3):end), ...
                                  'DACE',4);
Options.Search(2).type    = 'DACE';      % For choices, see mads_defaults
Options.Search(2).nIter   = 10;          % Number of iterations for Search #2
Options.Search(2).nPoints = 1;           % Number of poll or sample points
Options.Search(2).sfile   = 'regpoly0';  % filename must include full path
Options.Search(2).file    = '';          % filename must include full path
Options.Search(2).local   = 0;           % flag to turn on trust region
Options.Search(2).merit   = 0;           % flag to penalize clustered data
Options.Search(2).recal = strncmp(Options.Search(1).type(max(1,end-1):end), ...
                                  'NW',2) || ...
                          strncmp(Options.Search(1).type(max(1,end-3):end), ...
                                  'DACE',4);
Options.SurOptimizer      = 'mads';
Options.mvp1Surrogate     = 1;

Options.dace(2).reg       = 'regpoly0';
Options.dace(2).corr      = 'correxp';
Options.dace(2).theta     = 1;
Options.dace(2).lower     = 0.01;
Options.dace(2).upper     = 1000;
Options.dace(2).isotropic = 0;

% Specify Choices for POLL
Options.pollStrategy    = 'Standard_2n'; % For choices, see mads_defaults
Options.pollOrder       = 'Consecutive'; % For choices, see mads_defaults
Options.pollCenter      = 0;             % Poll around n-th filter point
Options.pollComplete    = 0;             % Flag for complete polling
Options.NPollComplete   = 0;             % Flag for complete neighbor polling
Options.EPollComplete   = 0;             % Flag for complete extended polling

% Specify Termination Criteria
Options.Term.delta      = 1e-4;          % minimum mesh size
Options.Term.nIter      = Inf;           % maximum number of iterations
Options.Term.nFunc      = 50000;         % maximum number of function evals
Options.Term.time       = Inf;           % maximum CPU time
Options.Term.nFails     = Inf;           % max number of consecutive Poll fails

% Choices for Mesh Control
Options.delta0          = 0.1;           % initial mesh size
Options.deltaMax        = Inf;           % bound on how coarse the mesh can get
Options.meshRefine      = 0.5;           % mesh refinement factor
Options.meshCoarsen     = 2.0;           % mesh coarsening factor

% Choices for Filter management (for problems with nonlinear constraints)
Options.hmin            = 1e-8;          % minimum infeasible point h-value
Options.hmax            = 1.0;           % maximum h-value of a filter point

% Choices for EXTENDED POLL (for MVP problems)
Options.ePollTriggerF   = 0.01;          % f-value Extended Poll trigger
Options.ePollTriggerH   = 0.01;          % h-value Extended Poll trigger

% MADS flag parameter values
Options.loadCache        = 1;            % load pre-existing Cache file
Options.countCache       = 1;            % count Cache points as function calls
Options.runStochastic    = 0;            % runs problem as a stochastic problem
Options.scale            = 2;            % scale directions using this log base
Options.useFilter        = 1;            % filter (0=none, 1=multi-pt, 2=2-pt)
Options.degeneracyScheme = 'random';     % scheme for degenerate constraints
Options.removeRedundancy = 1;            % discard redundant linear constraints
Options.computeGrad      = 0;            % compute gradient, if available
Options.saveHistory      = 0;            % saves MADS performance to text file
Options.plotHistory      = 1;            % plot MADS performance
Options.plotFilter       = 0;            % plot the filter real-time
Options.plotColor        = 'k';          % color of history plot
Options.debug            = 0;            % turn on status messages for debugging

Options.Sur.Term.delta   = 0.01;         % surrogate minimum mesh size

%*******************************************************************************
% DO NOT MODIFY AFTER THIS POINT
%*******************************************************************************

% Set up figure handles for real-time plots
if (Problem.nc == 0)
   Options.plotFilter = 0;
end
if (Options.plotFilter)
   figure;
   Options.fplothandle = gca;
end
if (Options.plotHistory == 2)
   figure;
   Options.hplothandle = gca;
end

% Set the path, and load any user-provided problem parameters
cwd = pwd;
cd(problemPath);
if (exist(Problem.File.P,'file') == 2)
   Problem.Param = feval(Problem.File.P);
   setappdata(0,'PARAM',Problem.Param);
end

% Get the initial iterates and call the optimizer
if isfield(Problem,'Param') && isfield(Problem.Param,'iterate0')
   iterate0 = Problem.Param.iterate0;
else
   iterate0 = feval(Problem.File.I);
end
[BestF,BestI,RunStats,RunSet] = mads(Problem,iterate0,Options);

% Perform any user-defined post-processing (must have argument)
if (exist(Problem.File.P,'file') == 2) && (nargin(Problem.File.P) < 0)
   Param = feval(Problem.File.P,BestF);
   setappdata(0,'PARAM',Param);
end

cd(cwd);
return
