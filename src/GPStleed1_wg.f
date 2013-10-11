#include <fintrf.h>
C This is wrapper file for GPStleed1.f (originally named evaluate.f).
C It converts double precision real numberis (matlab's default) into 
C single precision real numbers (which GPStleed1.f and tleed codes 
C called by it used), otherwise program crashes. 
      SUBROUTINE GPStleed1_wg(problem_dir,DIR,RANK,PARM,MINB,MAXB,NTYPE,
     & FITVAL)
      
      PARAMETER (NMAX=14,NDIM=3)

      REAL*8 PARM(NMAX,NDIM),MINB(NMAX,NDIM),MAXB(NMAX,NDIM),FITVAL
      REAL PARM_f(NMAX,NDIM),MINB_f(NMAX,NDIM),MAXB_f(NMAX,NDIM),
     & FITVAL_f
      INTEGER DIR,RANK,NTYPE(NMAX)
      character*(*) problem_dir

       FITVAL=0.0
       FITVAL_f=FITVAL
       DO 100 I=1,NMAX
              do 200 J=1,3
                     PARM_f(I,J)=PARM(I,J)
                     MINB_f(I,J)=MINB(I,J)
                     MAXB_f(I,J)=MAXB(I,J)
 200          enddo                
 100   enddo          

      write(*,*) 'GPStleed1_wg: before GPStleed1'
       call GPStleed1(problem_dir,DIR,RANK,PARM_f,MINB_f,MAXB_f,NTYPE,
     & FITVAL_f)
       FITVAL=FITVAL_f
      write(*,*) 'GPStleed1_wg: after GPStleed1: fitval = ', fitval_f

       return
       end

C     The gateway subroutine
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
CJCM
      mwPOINTER plhs(*), prhs(*)
      integer nlhs, nrhs
      mwSIZE MXGETM, MXGETN

      integer mxGetPr, mxCreateDoubleMatrix 
      integer strlen
CJCM      integer mxGetM, mxGetN, strlen
CJCM      integer nlhs, nrhs, plhs(*), prhs(*)
C      integer pr_in1, pr_in2, pr_in3, pr_in4, pr_in5, pr_in6,
      integer pr_in2, pr_in3, pr_in4, pr_in5, pr_in6,
     & pr_in7, pl_out1

      character*100 in1
      integer in2, in3, in7(14)
      real*8 in4(14,3),in5(14,3),in6(14,3),rin7(14),out1
       
      integer mnsize

      if(nrhs .ne. 7) then
         call mexErrMsgTxt('Seven inputs required.')
      endif
      if(nlhs .ne. 1) then
         call mexErrMsgTxt('One output required.')
      endif
         
C corresponds to first input(string): problem_dir
C      call mxCopyPtrToCharacter(pr_in1, in1, mnsize) %this doesn't work!
      strlen=mxGetM(prhs(1))*mxGetN(prhs(1))
      status = mxGetString(prhs(1),in1, strlen) 
C corresponds to input DIR 
      pr_in2 = mxGetPr(prhs(2))
      call mxCopyPtrToInteger4(pr_in2, in2, 1)
C corresponds to input RANK 
      pr_in3 = mxGetPr(prhs(3))
      call mxCopyPtrToInteger4(pr_in3, in3, 1)
C corresponds to input PARM
      mnsize =mxGetM(prhs(4))*mxGetN(prhs(4)) 
      pr_in4 = mxGetPr(prhs(4)) 
      call mxCopyPtrToReal8(pr_in4, in4, mnsize)
C corresponds to input MINB
      mnsize = mxGetM(prhs(5))*mxGetN(prhs(5)) 
      pr_in5 = mxGetPr(prhs(5)) 
      call mxCopyPtrToReal8(pr_in5, in5, mnsize)
C corresponds to input MAXB
      mnsize =mxGetM(prhs(6))*mxGetN(prhs(6)) 
      pr_in6 = mxGetPr(prhs(6)) 
      call mxCopyPtrToReal8(pr_in6, in6, mnsize)
C corresponds to input NTYPE
      mnsize = mxGetM(prhs(7))*mxGetN(prhs(7)) 
      pr_in7 = mxGetPr(prhs(7)) 
C      call mxCopyPtrToInteger4(pr_in6, in6, mnsize) %somehow it didn't work!
      call mxCopyPtrToReal8(pr_in7, rin7, mnsize)
      do I=1,14
            in7(I)=nint(rin7(I))
      enddo

C set output for mexfunction (corresponds to FITVAL) 
      plhs(1) = mxCreateDoubleMatrix(1, 1, 0)
      pl_out1 = mxGetPr(plhs(1))

C     Call the computational routine.
      call GPStleed1_wg(in1(1:strlen), in2, in3, in4,in5, in6, in7,out1)
C      call evaluate(in1, in2, %val(pr_in3),%val(pr_in4),%val(pr_in5),
C     & %val(pr_in6), %val(pr_out)) %somehow it didn't work!

       call mxCopyReal8ToPtr(out1, pl_out1, 1)

      return
      end
