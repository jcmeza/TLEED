      subroutine bb_tleed(xparm,fx,success)
      real  xparm(42)
      real  fx(1)
      integer i
      logical success

c
c     tleed variables
c
      PARAMETER (NMAX=14,NDIM=3)
      REAL PARM(NMAX,NDIM),MINB(NMAX,NDIM),MAXB(NMAX,NDIM),FITVAL
      REAL PARM_f(NMAX,NDIM),MINB_f(NMAX,NDIM),MAXB_f(NMAX,NDIM)
      INTEGER DIR,RANK,NTYPE(NMAX)
      character(100) problem_dir

C     on my mac
      problem_dir = "/Users/meza/MyProjects/TLEED"
c     on merced cluster
c      problem_dir = "/home/jcmeza/TLEED"
C
C
C
      success = .false.
      DIR   = 0
      DELTA = 0.4
      RANK  = 0
c
c  JCM Iterate for which rfactor = 0.2697
c
c iterate(1).x=[ -1.8757 -1.7941 -1.8067 -0.3861 0.2472 -0.0461 0.0690 0.1874 1.7112 1.7350 
c                 1.7378 1.7467 1.7751 1.7897 0 3.1141 3.0047 6.2250 -4.0621 1.2552 3.6738 
c                -4.2907 5.0398 0 5.0355 2.4703 2.5445 2.4371 0 0 -3.0047 1.2913 6.2250 1.2552
c                 1.2125 3.7093 0 0 5.0355 5.0402 0 2.4371]';
citerate(1).p={1 1 1 1 1 2 2 2 2 2 2 2 2 2};%best known solution!
c
c Set atom type: First five are Li = 1, last 9 atoms are Ni = 2
c
      do i = 1, nmax
         ntype(i) = 2
      enddo
c     first 5 atoms are Li = 1
      do i = 1, 5
         ntype(i) = 1
      enddo

c
c     Apparently KLEED/TLEED need the coordinates in the format [z, x, y]
c     Our convention was to place them in a 1-d array with all z coordinates, followed by x and y
c      write(*,*) 'bb_tleed: '
      DO j=1,NDIM
         do i=1,NMAX
            kindex = i + (j-1)*NMAX
            parm_f(I,J)=xparm(kindex)
c            write(*,*) 'k, x(k):', kindex, xparm(kindex)
            minb_f(I,J)=xparm(kindex) - DELTA*xparm(kindex)
            maxb_f(I,J)=xparm(kindex) + DELTA*xparm(kindex)
         enddo                
      enddo
      
cjcm
c      write(*,*) 'bb_tleed: before evaltleed'
c      write(*,*) NTYPE
cjcm

      call evaltleed(problem_dir,dir,rank,parm_f,minb_f,maxb_f,
     &     ntype,fitval,success)
      fx(1) = fitval
c      write(*,*) 'bb_tleed: success  = ', success
c      write(*,*) 'bb_tleed: fitval from evaltleed  = ', fitval_f
      
      return
      end
