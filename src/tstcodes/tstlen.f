      program tstlen
      character*100 :: problem_dir
      integer dir, rank
      dir = 0
      rank = 1
      problem_dir = '.'
      call evaltleed(problem_dir,dir,rank)

      end

      subroutine evaltleed(problem_dir,dir,rank)
      character*(*) problem_dir
      integer dir,rank

      character (LEN=100) :: tleed4,tleed5
      character (LEN=100) ::tleed4doti, tleed5doti
      character(3) :: procid
      character(3) :: workid
      
      CALL INT2CHAR(DIR,WORKID,3)
      CALL INT2CHAR(RANK,PROCID,3)

      write(6,*) 'workid = ', workid
      write(6,*) 'procid = ', procid
      
      TLEED4 = trim(problem_dir)//"/work"//WORKID//"/tleed4.i"
      tleed4doti=trim(problem_dir)//"/tleed4.i"

      write(6,*) 'problem_dir = ', problem_dir
      write(6,*) 'tleed4 = ', tleed4
      write(6,*) 'tleed4doti = ', tleed4doti
      
      OPEN(UNIT=8,FILE=tleed4doti,status='old')
      OPEN(UNIT=9,FILE=TLEED4,STATUS='UNKNOWN') 

      return
      end
C======================================================================
C     SUBROUTINE INT2CHAR
C     Converts integer "x" to char string "c" of length "length"
C======================================================================

      SUBROUTINE INT2CHAR(x,c,length)

      integer x, length
      character c*(*)
      
      character*10 alph
      integer xt,i,y
      
      alph = '0123456789'
      
      xt = abs(x)
      i=length

 100  if (i .gt. 0) then
         if ( xt .ne. 0) then
            y = mod(xt,10)
            c(i:i) =  alph(y+1:y+1)
            xt = int(xt/10)
         else
            c(i:i) = '0'
         endif
         i=i-1
         goto 100
      endif
      
      return
      end
