      program runkleed
      implicit none
c
c     New program created to compute 1 function evaluation of tleed and kleed
c     Parameters are read from a new input file called tleedparms.dat
c
c     tleed variables
c
      integer nmax, ndim
      parameter (nmax=14,ndim=3)

      character(100) problem_dir, data_dir, parmsfile
      real  xparm(42),minb(42), maxb(42)
      real  fx, t1, t2
      integer ntype(nmax), dir, rank, num_fcn
      integer maxfcn, icode
      logical success

      data_dir = "./tleed_data"
      problem_dir = "."
      dir   = 0
      rank  = 0
      maxfcn =100 
c
c     open parms file
c
      parmsfile = trim(data_dir)//'/tleedx0.dat'
      open(unit=3,file=parmsfile,status='old')

      do num_fcn=1,maxfcn
         
         write(*,*) '******************************'
         write(*,*) 'runkleed: num_fcn = ', num_fcn
         call setuptleed(ntype,xparm,minb,maxb,icode)
         if (icode < 0 ) then
c     end of file reached
            write(*,*) 'End of file reached on tleedx0.dat'
            goto 9999
         endif
         
         call cpu_time(t1)
         call evalkleed(problem_dir,dir,rank,xparm,minb,maxb,
     &        ntype,fx)
         call cpu_time(t2)

         write(*,*) 'runkleed: ', fx, t2-t1
      enddo
      
 9999 continue
      close(unit=3)
      end

      subroutine setuptleed(ntype,xparm,minb,maxb,icode)
      implicit none
      integer nmax, ndim
      parameter (nmax=14,ndim=3)

      real xparm(*), minb(nmax,ndim),maxb(nmax,ndim)
      real z(100), x(100), y(100)
      real delta

      integer i,j, kindex
      integer ntype(nmax)
      integer icode
      
      logical :: debug
c
      delta = 0.4
c
c     read(3,*) (ntype(i), i=1, nmax)
c     read(3,*) (xparm(i), i=1, nmax*ndim)
      
      icode = 0
      do i=1,nmax
         read(3,*,iostat=icode) ntype(i), z(i), x(i), y(i)
         if (icode < 0) then
c           end of file reached
            goto 9999
         endif
      enddo
      
c
c     kleed/tleed need the coordinates in the format [z, x, y]
c     our convention was to place them in a 1-d array with all z coordinates, followed by x and y
      do i=1,nmax
         xparm(i)        = z(i)
         xparm(i+nmax)   = x(i)
         xparm(i+2*nmax) = y(i)
      enddo
      
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
         do i=1,nmax
            write(*,*) i, xparm(i), xparm(i+14), xparm(i+28)
         enddo
       
      endif
 
 9999 continue
      return
      end