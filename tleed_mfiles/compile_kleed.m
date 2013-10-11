%this function compiles tleed codes and produces the executable of
%objective function for LEED optimization problem.
!make
%mex -v -fortran -O GPStleed1_wg.f GPStleed1.o tleed1GPS.o tleed2GPS.o tleedlibGPS.o
mex -fortran  GPSkleed_wg.f GPSkleed.o kleedGPS.o  kleedlibGPS.o
%mex -fortran eval_update_wg.f eval_update.o tleed1GPS.o tleed2GPS.o tleedlibGPS.o

             
%in case the output directory of tleed codes is not yet set up
problem_dir=mfilename('fullpath');
problem_dir=fileparts(problem_dir);
work=strcat(problem_dir, '/works000');
if exist(work,'dir')~=7
        mkdir(problem_dir,'works000')
end
