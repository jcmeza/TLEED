      program convertparms
      implicit none
c
c     New program created to convert old parameter files to new forma
c     Parameters are read from a new input file called tleedparms.dat
c
c     tleed variables
c
      integer nmax, ndim
      parameter (nmax=14,ndim=3)

      character(100) problem_dir, data_dir, parmsfile, newparmsfile
      real  xparm(42),minb(42), maxb(42)
      real  fx
      integer ntype(nmax), dir, rank, num_fcn
      integer maxfcn
       
      data_dir = "./tleed_data"
      problem_dir = "."
      dir   = 0
      rank  = 0
      maxfcn = 10
c
c     open parms file
c
      parmsfile = trim(data_dir)//'/tleedparms.dat'
      newparmsfile = trim(data_dir)//'/tleednewparms.dat'
      open(unit=3,file=parmsfile,status='old')
      open(unit=4,file=newparmsfile,status='old')
      
      do num_fcn=1,maxfcn
         write(*,*) '******************************'
         write(*,*) 'runtleed: num_fcn = ', num_fcn
         call setuptleed(ntype,xparm,minb,maxb)
      enddo
      
      close(unit=3)
      close(unit=4)
      end

      subroutine setuptleed(ntype,xparm,minb,maxb)
      implicit none
      integer nmax, ndim
      parameter (nmax=14,ndim=3)

      real xparm(*), minb(nmax,ndim),maxb(nmax,ndim)
      real delta

      integer i,j, kindex
      integer ntype(nmax)

      logical :: debug
c     
      delta = 0.4

      read(3,*) (ntype(i), i=1, nmax)
c
c     kleed/tleed need the coordinates in the format [z, x, y]
c     our convention was to place them in a 1-d array with all z coordinates, followed by x and y

      read(3,*) (xparm(i), i=1, nmax*ndim)

      do j=1,ndim
         do i=1,nmax
            kindex = i + (j-1)*nmax
            minb(i,j)=xparm(kindex) - delta*xparm(kindex)
            maxb(i,j)=xparm(kindex) + delta*xparm(kindex)
         enddo                
      enddo
      
c
c     write out coordinates if desired
c
      debug = .TRUE.
      if (debug) then
         
         do j=1,ndim
            do i=1,nmax
               kindex = i + (j-1)*nmax
               write(*,*) 'k, x(k),lb, ub:', kindex, 
     &              xparm(kindex), minb(i,j), maxb(i,j)
            enddo
         enddo
                  
c     now in tleed format
c         write(4,1000) nmax, ndim
 1000    format(2i6)
         do i=1,nmax
            write(4,1001) ntype(i), xparm(i), xparm(i+14), xparm(i+28)
 1001       format(i6, 3F16.8)
         enddo         
      endif
      
      return
      end

cjcm     dummy routine for debugging purposes
!      
!      subroutine evaltleedstub(problem_dir,dir,rank,xparm,minb,maxb,
!     &     ntype,fx)
!      
!      PARAMETER (NMAX=14,NSUB=6,NIDEN=5,NDIM=3,PENALTY=1.6)
!      REAL XPARM(NMAX,NDIM),MINB(NMAX,NDIM),MAXB(NMAX,NDIM),FX
!      INTEGER DIR,RANK,NTYPE(NMAX)
!      character*(*) problem_dir
!     
!      return
!      end

