      subroutine bb_kleed(xparm,fx)
      real*8  xparm(42)
      real*8  fx(1)
      integer i

c
c     kleed variables
c
      PARAMETER (NMAX=14,NDIM=3)
      REAL*8 PARM(NMAX,NDIM),MINB(NMAX,NDIM),MAXB(NMAX,NDIM),FITVAL
      REAL PARM_f(NMAX,NDIM),MINB_f(NMAX,NDIM),MAXB_f(NMAX,NDIM),
     & FITVAL_f
      INTEGER DIR,RANK,NTYPE(NMAX)
      character*28 problem_dir
      problem_dir = "/Users/meza/MyProjects/TLEED"
C
C
C
      DIR = 0
      DELTA = 0.4
      RANK = 0
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
c      write(*,*) 'bb_kleed: '
      DO j=1,NDIM
         do i=1,NMAX
            kindex = i + (j-1)*NMAX
            parm_f(I,J)=xparm(kindex)
cc            write(*,*) 'k, x(k):', kindex, xparm(kindex)
            minb_f(I,J)=xparm(kindex) - DELTA*xparm(kindex)
            maxb_f(I,J)=xparm(kindex) + DELTA*xparm(kindex)
         enddo                
      enddo
      
c      write(*,*) 'bb_kleed: before GPSkleed'
c      write(*,*) NTYPE

      call GPSkleed(problem_dir,dir,rank,parm_f,minb_f,maxb_f,
     &     ntype,fitval_f)
      fx(1) = fitval_f
c      write(*,*) 'bb_kleed: after GPSkleed: fitval = ', fitval_f
      
      return
      end
c
c     JCM: stub routine
c      
c      SUBROUTINE GPSkleed(problem_dir,DIR,RANK,PARM,MINB,MAXB,
c     &     NTYPE,FITVAL)
c      
c      PARAMETER (NMAX=14,NSUB=6,NIDEN=5,NDIM=3,PENALTY=1.6)
c      character*(*) problem_dir
c
c      REAL PARM(NMAX,NDIM),MINB(NMAX,NDIM),MAXB(NMAX,NDIM),FITVAL
c      INTEGER DIR,RANK,NTYPE(NMAX)
c      
c      fitval = PENALTY
c
c      return
c      end
