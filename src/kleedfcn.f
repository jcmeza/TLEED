      program kleedfcn
      implicit none
c
c     New program derived from runkleed and to be used as
c     as interface between Matlab and KLEED
c     Parameters are read from a new inputfile called kleedinputs.dat and
c     the output function value is to be written to kleedfval.dat

      integer :: nin, nout
      integer :: my_iostat
      character(256) :: my_iomsg
c     
c     tleed/kleed variables
c
      integer nmax, ndim
      parameter (nmax=14,ndim=3)

      character(100) problem_dir, data_dir
      character(100) inputfile, outfile
      real  xparm(42),minb(42), maxb(42)
      real  fval
      integer ntype(nmax), dir, rank
c
c     data_dir is the directory used by KLEED to get all its other input parameters
c     problem_dir is the directory
      
      data_dir = "./tleed_data"
      problem_dir = "."
      dir   = 0
      rank  = 0
      fval  = 1.6
c
c     Set up file I/O 
c
      inputfile = trim(problem_dir)//'/kleedinputs.dat'
      open(newunit=nin,file=inputfile,iostat=my_iostat, iomsg=my_iomsg)
      IF(my_iostat /= 0) THEN
         WRITE(*,*) 'Failed to open kleedinputs file'
         WRITE(*,*) ' iomsg='//trim(my_iomsg)
         STOP
      ENDIF
      
      outfile = trim(problem_dir)//'/kleedfval.dat'
      open(newunit=nout,file=outfile,iostat=my_iostat, iomsg=my_iomsg)
      IF(my_iostat /= 0) THEN
         WRITE(*,*) 'Failed to open kleedfval file'
         WRITE(*,*) ' iomsg='//trim(my_iomsg)
         STOP
      ENDIF
      
      call setuptleed(nin,ntype,xparm,minb,maxb)

      call evalkleed(problem_dir,dir,rank,xparm,minb,maxb,
     &        ntype,fval)
      write(nout,*) fval

      close(nin)
      close(nout)
      
      end
      subroutine setuptleed(nin,ntype,xparm,minb,maxb)
      implicit none
      integer nin
      integer nmax, ndim
      parameter (nmax=14,ndim=3)

      real xparm(*), minb(nmax,ndim),maxb(nmax,ndim)
      real delta

      integer i,j, kindex
      integer ntype(nmax)

      logical :: debug
c
      delta = 0.1
c
      read(nin,*) (ntype(i), i=1, nmax)

c
c     kleed/tleed need the coordinates in the format [z, x, y]
c     our convention was to place them in a 1-d array with all z coordinates, followed by x and y

      read(nin,*) (xparm(i), i=1, nmax*ndim)

      do j=1,ndim
         do i=1,nmax
            kindex = i + (j-1)*nmax
            minb(i,j)=xparm(kindex) - delta*xparm(kindex)
            maxb(i,j)=xparm(kindex) + delta*xparm(kindex)
         enddo                
      enddo
c
c     write out coordinates if desired
c     mostly for debuggging
      debug = .FALSE.
      if (debug) then
         do j=1,ndim
            do i=1,nmax
               kindex = i + (j-1)*nmax
               write(*,*) 'k, x(k),lb, ub:', kindex, 
     &              xparm(kindex), minb(i,j), maxb(i,j)
            enddo
         enddo
      endif
      
      return
      end
