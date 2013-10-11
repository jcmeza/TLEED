
C This is wrapper file for GPStleed1.f (originally named evaluate.f).
C It converts double precision real numberis (matlab's default) into 
C single precision real numbers (which GPStleed1.f and tleed codes 
C called by it used), otherwise program crashes. 
      SUBROUTINE eval_update_wg(problem_dir, PARM,NTYPE,
     & parm_out,rntype_out)
      PARAMETER (NMAX=14,NDIM=3)

      REAL*8 PARM(NMAX,NDIM), parm_out(NMAX,NDIM)
      REAL PARM_f(NMAX,NDIM), parm_out_f(NMAX,NDIM) 
C      INTEGER NTYPE(NMAX),ntype_out(NMAX)
      INTEGER NTYPE(NMAX)
       real*8 rntype_out(Nmax)
       integer ntype_out(Nmax)
       character(*) problem_dir
       DO 100 I=1,NMAX
              ntype_out(I)=ntype(I)
              do 200 J=1,3
                     parm_out(I,J)=PARM(I,J)
                     PARM_f(I,J)=PARM(I,J)
                     parm_out_f(I,J)=PARM(I,J)
 200          enddo                
 100   enddo          

       call eval_update(problem_dir, PARM_f,NTYPE,parm_out_f,ntype_out)

       DO 300 I=1,NMAX
              rntype_out(I)=ntype_out(I);
              do 400 J=1,3
                     parm_out(I,J)=parm_out_f(I,J)
 400          enddo                
 300   enddo          


       return
       end

C     The gateway subroutine
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)

      integer mxGetPr, mxCreateDoubleMatrix 
      integer mxGetM, mxGetN
      integer nlhs, nrhs, plhs(*), prhs(*)
      integer pr_in1, pr_in2, pl_out1, pl_out2

      integer in2(14), out2(14)
      real*8 in1(14,3),rin2(14), rout2(14),out1(14,3)
       
      integer mnsize
      character*100 cin1
      integer strlen

      strlen=mxGetM(prhs(1))*mxGetN(prhs(1))
      status = mxGetString(prhs(1),cin1, strlen)
C corresponds to input PARM
      mnsize =mxGetM(prhs(2))*mxGetN(prhs(2)) 
      pr_in1 = mxGetPr(prhs(2)) 
      call mxCopyPtrToReal8(pr_in1, in1, mnsize)
C corresponds to input NTYPE
      mnsize = mxGetM(prhs(3))*mxGetN(prhs(3)) 
      pr_in2 = mxGetPr(prhs(3)) 
C      call mxCopyPtrToInteger4(pr_in6, in6, mnsize) %somehow it didn't work!
      call mxCopyPtrToReal8(pr_in2, rin2, mnsize)
      do I=1,14
            in2(I)=nint(rin2(I))
      enddo

C set output for mexfunction (corresponds to FITVAL) 
      plhs(1) = mxCreateDoubleMatrix(14, 3, 0)
      pl_out1 = mxGetPr(plhs(1))
      plhs(2) = mxCreateDoubleMatrix(1,14, 0)
      pl_out2 = mxGetPr(plhs(2))

C     Call the computational routine.
      call eval_update_wg(cin1(1:strlen),in1,in2,out1,rout2)

C      do I=1,14
C            rout2(I)=out2(I)+0.000
C            rout2(I)=real(out2(I))
C      enddo

       call mxCopyReal8ToPtr(out1, pl_out1, 14*3)
       call mxCopyReal8ToPtr(rout2, pl_out2, 14)

C      do I=1,14
c            out2(I)=nint(rout2(I))
C      enddo

      return
      end
