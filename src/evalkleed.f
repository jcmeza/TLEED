C********************************************************
C FILE: evaluate.f   
C
C DESCRIPTION:
C LEED Structure validation and modification subroutine
C
C USAGE COMMENTS: To be linked with genetic algorithm code.
C      note: Link with gagen.f for sequential testing, i.e.,
C             f77 -c evaluate.f
C             f77 gagen.f evalute.o tleed1.o tleed2.o tleedlib.o 
C  
C CONSTANTS DEFINED: 
c    NMAX    -- THE MAXIMUM NUMBER OF INEQUIVALENT ATOMS IN THE SURFACE, 
C    NSUB    -- NUMBER OF INEQUIVALENT ATOMS IN SUBSTRATE LAYER. 
c    NIDEN   -- THE MAXIMUM NUMBER OF ELEMENTS. 
c    NDIM    -- THE DIMENSION OF COORDINATES.
c    NNMAX   -- THE MAXMUM NUMBER OF ATOMS OF OVERLAYER.
C    NLAY    -- THE MAXIMUM NUMBER OF COMPOSITE LAYERS.
C                                                           
C SUBROUTINES DEFINED:
C     GRAVSTRUC  -  (one line description here)
C     EVALU      -  (one line description here)
C     VALUATE    -  modify struct according to symmetry and evaluate it 
C     MODSTRUC   -  modify a random struct to minimum feasible strutcure
C     WRITECOOR  -  write coordinates to kleed5.i file
C     SORTAUold  -  (one line description here)
C     SORTAU     -  (one line description here)
C     FINDNTAU   -  (one line description here)
C     PIKSRT     -  sort an array
C     SEPSURF    -  (one line description here)
C     MODSTRUCT  -  (one line description here)
C     VALUAT     -  evaluate the validity of a structure
C     COORDAT    -  (one line description here)
C     SORT       -  (one line description here)
C     INDEXX     -  (one line description here)
C     INDEXXN    -  (one line description here)
C     INT2CHAR   -  Convert intergers to characters
C
C INPUT:
C   CHARACTER
C     DIR        - work directory for output
C     PROCID     - id for the calling process (used as a file sufffix)
C   REAL
C     PARM(NMAX,3) - position parameters
C     MINB(NMAX,3) - position parameters lower boundary
C     MAXB(NMAX,3) - position parameters upper boundary
C   INTEGER  
C     NTYPE(NMAX) - chemical element identities
C
C OUTPUT:
C     REAL    FITVAL - position parameters
C
C VARIABLES: (list of all major variables)
C
C    REAL  
C          RMIN(5)      - ?
C          RMAX(5)      - ?
C          COORSUB(NSUB,3) - ?
C    INTEGER
C          NTYPESUB(NSUB,3) - ?     
C    CHARACTER 
C          SYM - symmetry code
C          RESULT - SUCCESS or FAILURE of struct validaion             
C 
C FILE DESCRIPTOR USED:
C     UNIT=4: write to file KLEED5
C     UNIT=7: read from file "kleed5.i"
C     UNIT=8: read from file "kleed4.i"
C     UNIT=9: write to file KLEED4
C     UNIT=99: trace file for debugging messages
C
C GENERAL COMMENTS:                            
C  ATOM WITHIN A GIVEN RANGE AND VALUATED TO RETURN 'PENALTY' OR 'VALIDITY'. 
C  IF IT IS VALID, WRITE THE STRUCTURE INTO A INPUT FILE  CALLED KLEED5A.I. 
C
C LOG:
C     CREATED: 7/28/98 By Fajian Shi
C     MODIFIED: 
C     07/31/98 (Fajian Shi) EVALU fixed see C073198
C     09/08/98 (Greg Stone) Added MINB and MAXB parameters
C              to be used for range checking

CGPS  10/05/04 (Zhengji Zhao) changed file to GPStleed1.f to use NOMADm
CGPS      to solve optimization problem instead of genetic optimization
C********************************************************

      SUBROUTINE evalkleed(problem_dir,DIR,RANK,PARM,MINB,MAXB,NTYPE,
     & FITVAL)
      
      PARAMETER (NMAX=14,NSUB=6,NIDEN=5,NDIM=3,PENALTY=1.6)
      
      REAL PARM(NMAX,NDIM),MINB(NMAX,NDIM),MAXB(NMAX,NDIM),FITVAL
      real parm_t(NMAX,NDIM)
      integer ntype_t(NMAX)

      real PARM_out(NMAX,NDIM)


      INTEGER DIR,RANK,NTYPE(NMAX)

      REAL COORSUB(NSUB,NDIM),COORSUB0(NSUB,NDIM)
      REAL RMIN(NIDEN),RMAX(NIDEN),RMID(NIDEN)
      INTEGER NTYPSUB(NSUB)
      
      CHARACTER(len=16) SYMC,RESULT,NCODE(NMAX)

      CHARACTER(len=100) KLEED4,KLEED5,SEARCHS,TRACE 
      character*(*) problem_dir
      character(len=100) kleed4doti, kleed5doti, workdir
      character(len=100) XPIN,XPIN3
      CHARACTER(len=3) PROCID
      CHARACTER(len=3) WORKID
     
      common /doti/kleed4doti, kleed5doti
      COMMON KLEED4,KLEED5,SEARCHS,TRACE
      
      DATA RMIN/1.0,1.0,1.0,1.0,1.0/
      DATA RMAX/1.4,1.7,1.5,1.5,1.5/
      DATA RMID/1.23,1.245,1.2,1.2,1.2/

      FITVAL=0.0

C Define the symmetry code.
      SYMC='P4M' 
      RESULT=' '
      AA=12.45
      DO I=1,NMAX
         NCODE(I)='A'
      ENDDO
      
C Input substrate coordinates, COORSUB0, AA is the length of unit cell.

      DATA COORSUB0/0.0,0.0,0.0,0.0,0.0,0.0,
     &     0.1,0.3,0.5,0.3,0.5,0.5,
     &     0.1,0.1,0.1,0.3,0.3,0.5/
      DO I= 1,6
         COORSUB(I,1)=3.5214
         COORSUB(I,2)=AA*COORSUB0(I,2)
         COORSUB(I,3)=AA*COORSUB0(I,3)
         NTYPSUB(I)=2
c     write(*,*) 'COORSUB=',(coorsub(i,j),j=1,3)
      ENDDO

C     Setup input files and write structure to the trace file
      CALL INT2CHAR(DIR,WORKID,3)
      CALL INT2CHAR(RANK,PROCID,3)
c
      workdir = trim(problem_dir)//'/kwork'//WORKID

      KLEED4  = trim(workdir)// '/kleed4i' // PROCID
      KLEED5  = trim(workdir)// '/kleed5i' // PROCID
      SEARCHS = trim(workdir)// '/searchs' // PROCID
      TRACE   = trim(workdir)//'/trace' // PROCID
      XPIN    = trim(workdir)//'/xpin' // PROCID

c        write(*,*) 'problem_dir',problem_dir
c        write(*,>*) kleed4
c        write(*,*) kleed5
c        write(*,*) searchs
c        write(*,*) trace

      kleed4doti=trim(workdir)//'/kleed4.i'
      kleed5doti=trim(workdir)//'/kleed5.i'

      OPEN(UNIT=99,FILE=TRACE,STATUS='UNKNOWN')
      WRITE (99,*) "ID: ",PROCID," original Parms passed to evaluat
     &e.f "
c           DO I=1,NMAX
c              WRITE (99,"(I4,3F8.4)") NTYPE(I),PARM(I,1),
c     &                  PARM(I,2),PARM(I,3)
c           ENDDO
           CLOSE (UNIT=99)

           XPIN3 = trim(workdir)// '/XPIN2'
           OPEN(UNIT=98,FILE=XPIN3,STATUS='UNKNOWN',POSITION='APPEND')
           WRITE (98,*) "ID: ",PROCID," original Parms passed to evaluat
     &e.f "
c           DO I=1,NMAX
c              WRITE (98,"(I4,3F8.4)") NTYPE(I),PARM(I,1),
c     &             PARM(I,2),PARM(I,3)
c           PARM_out(I,1)=PARM(I,1)
c           PARM_out(I,2)=PARM(I,2)
c           PARM_out(I,3)=PARM(I,3)
c           ENDDO
              WRITE (98,"(14I4)") (NTYPE(I),I=1,NMAX)
              WRITE (98,"(42F16.10)") (PARM(I,1),I=1,NMAX),
     &     (PARM(I,2),I=1,NMAX),(PARM(I,3),I=1,NMAX)
           CLOSE (UNIT=98)

           DO I=1,NMAX
              ntype_t(I)= NTYPE(I)
              do 9876 J=1,3
                  parm_t(I,J)=PARM(I,J)
9876          enddo            
           ENDDO

C VALUATE turns the random surface structure(usually invalid) into a valid 
C structure with defined symmetry 'SYMC' by moving atoms laterally. 
C Then it valuates the structure and returns RESULT as 'SUCCESS'for a valid 
C structure or 'PENALTY' if any two atoms are too close or any atom is 
C alone (not touched by any atom).

      itimes = 0
 20   continue
      if (itimes.gt.0) then
cjcm         write(*,*) 'evalkleed: invalid structure'
         fitval = penalty
         return
      endif

      CALL VALUATE(PARM,NTYPE,NCODE,NMAX,SYMC,RESULT,RMIN,RMAX,NIDEN,
     &     COORSUB,NTYPSUB,NSUB,AA,nerror)

      IF(RESULT.NE.'SUCCESS'.or.nerror.eq.1) THEN
c     CALL GRAVSTUC(PARM,NTYPE,NMAX,RMIN,RMAX,RMID,NIDEN,
c     &          COORSUB,NTYPSUB,NSUB)
c     write(*,*) PROCID,' :GRAVSTUC finished',itimes 
         itimes = itimes + 1
         GOTO 20
      ELSE
CGPS
cjcm         write(*,*) 'evalkleed: valid structure'
C         return
C         DO I=1,NMAX
C            IF ((PARM(I,1).LT.MINB(I,1)).OR.
C     &           (PARM(I,1).GT.MAXB(I,1)).OR.
C     &           (PARM(I,2).LT.MINB(I,2)).OR.
C     &           (PARM(I,2).GT.MAXB(I,2)).OR.
C     &           (PARM(I,3).LT.MINB(I,3)).OR.
C     &           (PARM(I,3).GT.MAXB(I,3))) THEN
C               FITVAL = PENALTY
C            ENDIF
C         ENDDO
C         IF (FITVAL.EQ.PENALTY) THEN
Cc            write(0,*) PROCID,' :failed bounds check,returning'
C            RETURN
C         ENDIF

cjcm        print*, 4

c     write(99,*) PROCID,' :The structure coordinates:'
c     write(0,*) PROCID,' :The structure coordinates:'
c     do i=1,NMAX
c     WRITE (99,"(I4,3F8.4)") NTYPE(I),PARM(I,1),PARM(I,2),
c     &             PARM(I,3)
C     WRITE (0,"(I4,3F8.4)") NTYPE(I),PARM(I,1),PARM(I,2),
C     &             PARM(I,3)
c     enddo

         OPEN(UNIT=99,FILE=TRACE,STATUS='OLD',POSITION='APPEND')
         WRITE (99,*) "ID: ",PROCID," Parms passed to tleed codes: 
     &(may be changed by VALUATE)"
           DO I=1,NMAX
              WRITE (99,"(I4,3F8.4)") NTYPE(I),PARM(I,1),
     &             PARM(I,2),PARM(I,3)
           ENDDO
C write into the differecne of parms
           WRITE (99,*) "parm difference:" 
           DO I=1,NMAX
              WRITE (99,"(I4,3F8.4)") NTYPE(I)-ntype_t(I),
     &             PARM(I,1)-parm_t(I,1),
     &             PARM(I,2)-parm(I,2),PARM(I,3)-parm_t(I,3)
           ENDDO

           WRITE (99,*) PROCID," Calling kleed: ",PROCID
           CLOSE (UNIT=99)

           XPIN3 = trim(workdir)// '/XPIN2'
           OPEN(UNIT=98,FILE=XPIN3,STATUS='OLD',POSITION='APPEND')
           WRITE (98,*) "ID: ",PROCID," Parms passed to tleed codes: 
     &(may be changed by VALUATE)"
c           DO I=1,NMAX
c              WRITE (98,"(I4,3F8.4)") NTYPE(I),PARM(I,1),
c     &             PARM(I,2),PARM(I,3)
c           ENDDO

           WRITE (98,"(14I4)") (NTYPE(I),I=1,NMAX)
           WRITE (98,"(42F16.10)") (PARM(I,1),I=1,NMAX),
     &          (PARM(I,2),I=1,NMAX),(PARM(I,3),I=1,NMAX)
           CLOSE (UNIT=98) 

cAGL Now open trace file before the kleed subroutine
           OPEN(UNIT=99,FILE=TRACE,STATUS='OLD',POSITION='APPEND')
c
c In KLEED only one subroutine is called: kleed (instead of
c  tleed1.f and tleed2.f in TLEED)
           OPEN (UNIT=2,FILE=SEARCHS,STATUS='UNKNOWN')

c           write(*,*) 'evalkleed: calling kleed'
           CALL kleed(workdir, WORKID,PROCID,fitval)
           OPEN(UNIT=98,FILE=XPIN3,STATUS='OLD',POSITION='APPEND')
           WRITE (98,*) 'evalkeed: R-factor (kleed)',fitval
           CLOSE(98)
c           write(*,*) 'evalkeed: fitval =',fitval

           CLOSE(4)
           CLOSE(5)
           CLOSE(2)       
           CLOSE(7)

           WRITE (99,*) WORKID,":",PROCID,
     &          "...returned kleed..."
           CLOSE (UNIT=99)
           
C          Assign penalty when r-factor is less than .1
C          IF (FITVAL.LT.0.1) FITVAL = PENALTY JCM: ?????

           IF (FITVAL.LT.0.1) then
               FITVAL = PENALTY+1
           else if (fitval.eq. 100.) then
                write(*,*) 'fitval handia',fitval
                Fitval=1.62
           else
           end if
           OPEN(UNIT=99,FILE=TRACE,STATUS='OLD',POSITION='APPEND')
           WRITE (99,*) WORKID,":",PROCID,
     &          " ...Returned tleed2 fitval=",fitval
           CLOSE (UNIT=99)
        ENDIF

        RETURN
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General        
C********************************************************           
        SUBROUTINE VALUATE(COORD,NTYPE,NCODE,NMAX,SYMC,RESULT,RMIN,RMAX,
     &NIDEN,COORSUB,NTYPSUB,NSUB,AA,nerror)
C Generate a sound structure from the input random structure. 
c The 17 symmetries, contraction to the adge, corner, center,
c and restraint on minimum distance and at least one contact.
        PARAMETER (NLAY=5)

        CHARACTER(len=16) RESULT,SYMC
        CHARACTER(len=16) NCODE(NMAX)

        REAL COORD(NMAX,3),COORSUB(NSUB,3)
        REAL RNZ(87),RMIN(NIDEN),RMAX(NIDEN) 
        REAL SPAC(NLAY)
        INTEGER NTYPE(NMAX),NTYPSUB(NSUB)
        INTEGER KLAY(NLAY)

c
c Number of types of atoms(chemical elements, chemical variants)
        DATA RNZ/0.0,.373,1.,1.52,1.113,.795,.772,.549,.604,
     *0.709,1.000,1.858,1.599,1.432,1.176,1.105,1.035,.994,1.,
     *2.272,1.974,1.606,1.448,1.311,1.249,1.367,1.241,1.253,
     *1.246,1.278,1.335,1.221,1.225,1.245,1.160,1.145,1.0,2.475,
     *2.151,1.776,1.590,1.429,1.363,1.352,1.325,1.345,1.376,
     *1.445,1.489,1.626,1.405,1.450,1.432,1.331,1.000,2.655,2.174,
     *1.870,1.825,1.820,1.814,1.630,1.620,1.995,1.787,1.763,1.752,
     *1.743,1.734,1.724,1.940,1.718,1.564,1.430,1.370,1.371,
     *1.338,1.357,1.373,1.442,1.503,1.700,1.750,1.545,1.673,1.45,
     *1.000/
       
c
c ===========================================================
c Modify the rondom structure to a symmetrical structure. 
       
c        write(99,*) 'Call MODSTRUCT'
        CALL MODSTRUCT(SYMC,COORD,NTYPE,NCODE,NMAX,RMIN,RMAX,
     %NIDEN,AA)

c VALUATE the structure, returns penalty or success.
c        write(99,*) 'Call VALUAT'
        CALL VALUAT(COORD,NTYPE,NMAX,RMIN,RMAX,NIDEN,RESULT,
     &COORSUB,NTYPSUB,NSUB)
c Sort the coordinates in a increasing order. 
        CALL SORTLOCAL(NMAX,COORD,NTYPE,NCODE)
c        write(99,*) 'Put the parameter in order, finished'
        IF(RESULT.EQ.'PENALTY') THEN 
cjcm          write(*,*) 'Invalid structure, return for GRAVSTUC'
           RETURN
        ENDIF
c       write(*,*) 'Valid structure'
c Separate the surface into several composite layers,
c DSPC IS THE MINIMUM SPACING BETWEEN TWO COMPOSITE LAYERS.
        DO I=1,NLAY
           SPAC(I)=0.0
           KLAY(I)=1
        ENDDO   
c Interlayer spacing parameter: 1.5<=DSPC<=2.0  
        DSPC=1.90
        KLAYER=1
        CALL SEPSURF(COORD,NMAX,KLAY,SPAC,NLAY,KLAYER,DSPC)
c       write(*,*) 'Separate the structure finished' 
c       write(*,*) klayer,klay,spac
c Write the coordinates into input file kleed5a.i.
        CALL WRITCOOR(COORD,NTYPE,NCODE,NMAX,KLAY,SPAC,NLAY,
     *KLAYER,nerror)
        RETURN       
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General        
C********************************************************
c COORDAT RETURNS ALL COORDINATES IN THE SURFACE NEEDED IN KLEED5.I.
        
        SUBROUTINE WRITCOOR(COORD,NTYPE,NCODE,NMAX,KLAY,SPAC,NLAY,
     *KLAYER,nerror)
        PARAMETER (NNMAX=80,NLAYE=5)
        CHARACTER(len=16) NCODE(NMAX)
        CHARACTER(len=60) LINE
        REAL COORD(NMAX,3),COORA(NNMAX,3),SPAC(NLAY)
        REAL COORAL(NLAYE,NNMAX,3)
        REAL SPAC1(NLAYE,3),FR(NLAYE),VO(NLAYE),VI(NLAYE)
        REAL COORSUB(2,3),ASB(3),EIEFDE(3)
        INTEGER NTYPE(NMAX),NNTYPA(NLAYE,NNMAX)
        INTEGER NTYPA(NNMAX),KLAY(NLAY),NKK(NLAYE)
        INTEGER NTAU(NLAYE)

CGPS        CHARACTER*32 KLEED4,KLEED5,SEARCHS,TRACE 
        CHARACTER(len=100) KLEED4,KLEED5,SEARCHS,TRACE 
CGPS+
        CHARACTER(len=100) kleed4doti, kleed5doti
        COMMON /doti/kleed4doti, kleed5doti
        COMMON KLEED4,KLEED5,SEARCHS,TRACE

        DATA ASB,EIEFDE/1.7317,1.2450,1.2450,15.00,205.00,4.00/
        OPEN(UNIT=4,FILE=KLEED5,STATUS='UNKNOWN') 
CGPS        OPEN(UNIT=7,FILE='kleed5.i',status='old')
        OPEN(UNIT=7,FILE=kleed5doti,status='old')
CORIG        DO I=1,190
CORIG           READ(7,2001) LINE
cORIG           WRITE(4,2001) LINE
cORIG        ENDDO
cORIG        CLOSE(7)
        DO I=1,NLAY
           SPAC1(I,1)=SPAC(I)
           SPAC1(I,2)=0.0
           SPAC1(I,3)=0.0
           FR(I)=0.5
           VO(I)=0.0
           VI(I)=-5.0
        ENDDO
        
        do i=1,nnmax
           NTYPA(I)=1
           do j=1,3
              COORA(I,J)=0.0
           ENDDO
        ENDDO
        KK=1
      
        DO I=1,KLAYER
           K1=KLAY(I)
           K2=KLAY(I+1)-1
           CALL COORDAT(COORD,NTYPE,NCODE,K1,K2,KK,
     %COORA,NTYPA)
cAGL Hurrengo 4 lerroak nik jarri ditut
           nerror=0
           if (klayer.eq.1.and.kk.ne.58) then
                nerror=1
c                write(*,*) 'nerror',nerror
c                goto 344 
           else
                nerror=0
cjcm                write(*,*) 'nerror',nerror
           endif 
           DO j=1,kk
              DO L=1,3
                 COORAL(I,J,L)=COORA(J,L)
              ENDDO
              NNTYPA(I,J)=NTYPA(J)
           ENDDO
           NKK(I)=KK
        ENDDO
       DO I=1,190
           READ(7,2001) LINE
           WRITE(4,2001) LINE
        ENDDO
        CLOSE(7)
        WRITE(4,1007) KLAYER,KLAYER
        WRITE(4,1007) (NKK(I),I=1,KLAYER)
        IFLAG=0
        WRITE(4,1007) IFLAG

        DO I=1,NLAY
           NTAU(I)=1
        ENDDO

        CALL SORTNTAU(NNTYPA,NTYPA,NLAY,NNMAX,KLAYER,NTAU,NKK)
        WRITE(4,1007) (NTAU(I),I=1,KLAYER)
        DO I=1,KLAYER
           kk=nkk(i)
           WRITE(4,1007) (NNTYPA(I,J),J=1,KK)
        ENDDO
        DO I=1,KLAYER
           do j=1,Nkk(I)
           write(4,1005) (COORAL(I,J,L),L=1,3)
c          write(4,1005) (COORAL(I,J,L),L=1,3),NNTYPA(I,J)
           enddo
cAGL I think the next IF loop is not necessary for KLEED, since we
cAGL only have one asb vector in that case
c           IF(I.NE.KLAYER) THEN
c              write(4,1005) (spac1(i,l),l=1,3)
c              write(4,1006) FR(I),VO(I),VI(I)
c           ELSE
c              SPAC1(I,1)=COORSUB(1,1)-COORAL(I,NKK(I),1)
c              WRITE(4,1005) (SPAC1(I,L),L=1,3)
c              write(4,1006) FR(I),VO(I),VI(I)
c           ENDIF
        ENDDO
        write(4,1005) ASB 
        write(4,1006) FR(1),VO(1),VI(1)
        write(4,1006) EIEFDE
        close(4)
CGPS        OPEN(UNIT=8,FILE='kleed4.i',status='old')
C        write(0,*) kleed4, kleed5, kleed4doti, kleed5doti

        OPEN(UNIT=8,FILE=kleed4doti,status='old')
        OPEN(UNIT=9,FILE=KLEED4,STATUS='UNKNOWN') 
        DO I=1,25
           READ(8,2001) LINE
           WRITE(9,2001) LINE
        ENDDO
        DO I=1,58
           READ(8,2001) LINE
        ENDDO
        NK=0
        DO I=1,KLAYER
           NK=NK+NKK(I)
        ENDDO
        DO I=1,NK
           WRITE(9,1008) 0.0,0.0,0.0,1
        ENDDO
        DO I=1,10
           READ(8,2001) LINE
           WRITE(9,2001) LINE
        ENDDO
        close(8)
        close(9)
344        continue
1005    FORMAT(3F8.4,I3)
1006    FORMAT(3F7.2,I3)
1007    FORMAT(60I3)
1008    FORMAT(3F7.4,I3)
2001    FORMAT(A60)
        RETURN       
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C******************************************************** 
        SUBROUTINE SORTNTAU(NNTYPA,NTYPA,NLAY,NNMAX,NLAYER,NTAU,NKK)
        INTEGER NNTYPA(NLAY,NNMAX),NTAU(NLAY),NKK(NLAY)
        INTEGER NTYPA(NNMAX)
        NNTAU=0

C        write(1,*) 'NKK(i)=',NKK
        DO I=1,NLAYER
           DO J=1,NKK(I)
              NTYPA(J)=NNTYPA(I,J)
           ENDDO
c          write(1,*) 'NTYPA(i)=', (NTYPA(j),j=1,NKK(i))
           KK=NKK(I)
           NNTAU=0
           CALL FINDNTAU(KK,NTYPA,NNMAX,NNTAU)
           NTAU(I)=NNTAU 
        ENDDO
        RETURN
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C******************************************************** 
        SUBROUTINE FINDNTAU(KK,NTYPA,NNMAX,NNTAU)
        INTEGER NTYPA(NNMAX),NTA(5)
        DO I=1,5
           NTA(I)=0
           DO J=1,KK
              IF(NTYPA(J).EQ.I) NTA(I)=1
           ENDDO
           NNTAU=NNTAU+NTA(I)
        ENDDO
C        write(1,*) 'NTYPA(I)=',(ntypa(j),j=1,kk)
C        write(1,*) 'NTA(I)=',(nta(j),j=1,5)
C        write(1,*) 'NNTAU=',NNTAU
        RETURN
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************        

        SUBROUTINE SORTNTAUold(NNTYPA,NLAYER,NTAU,NKK)
        PARAMETER (NLAY=5,NNMAX=80)
        INTEGER NNTYPA(NLAY,NNMAX),NTAU(NLAY),NKK(NLAY)
        INTEGER NTYPA(NNMAX),NWKSP(NNMAX),IWKSP(NNMAX)
        DO I=1,NNMAX
           IWKSP(I)=0
        ENDDO
        NNTAU=0

        DO I=1,NLAYER
           NTAU(I)=1
           DO J=1,NKK(I)
              NTYPA(J)=NNTYPA(I,J)
           ENDDO
           KK=NKK(I)
           CALL INDEXXN(KK,NTYPA,IWKSP)
           DO 11 J=1,NKK(I)
              NWKSP(J)=NTYPA(J)
11         CONTINUE
           DO 12 J=1,NKK(I)
              NTYPA(J)=NWKSP(IWKSP(J))
12         CONTINUE
C SUBROUTINE FINDNTAU RETURNS NTAU=NNTAU. +++!!!!!
           CALL FINDNTAU(KK,NTYPA,NNMAX,NNTAU)
           NTAU(I)=NNTAU
        ENDDO
        RETURN
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************
        SUBROUTINE FINDNTAUold(N1,NTPA,NN,NU)
        PARAMETER (NNMAX=80)
        INTEGER NTPA(NN),NTA(NNMAX)
        DO I=1,N1
           NTA(I)=NTPA(I)
        ENDDO
        CALL PIKSRT(N1,NTA)
        NU=1
        DO J=2,N1
           IF(NTA(J).NE.NTA(J-1)) THEN
             NU=NU+1
           ENDIF
        ENDDO
        write(1,*) 'NU,N1,NTYP=',NU,N1,(NTA(I),I=1,N1)
C        PAUSE
        RETURN
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************
        SUBROUTINE PIKSRT(N,NRR)
C Sort an arry NRR of length N into ascending numerical order.
        INTEGER NRR(N)
        DO J=2,N
           M=NRR(J)
           DO I=J-1,1,-1
              IF(NRR(I).LE.M) GOTO 10
              NRR(I+1)=NRR(I)
           ENDDO
           I=0
10         NRR(I+1)=M
        ENDDO
        RETURN
        END

C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************
        SUBROUTINE SEPSURF(COORD,NMAX,KLAY,SPAC,NLAY,KLAYER,DSPC)
        REAL COORD(NMAX,3),SPAC(NLAY)
        INTEGER KLAY(NLAY)
        KLAYER=1
        KLAY(1)=1
        NPARM=NMAX
        DO I=2,NPARM
           SPACING=COORD(I,1)-COORD(I-1,1)
           IF(SPACING.GT.DSPC) THEN
              KLAYER=KLAYER+1
              KLAY(KLAYER)=I
              SPAC(KLAYER-1)=SPACING
           ENDIF
        ENDDO
        KLAY(KLAYER+1)=NPARM+1
        RETURN
        END

C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************
        SUBROUTINE MODSTRUCT(SYMC,COORD,NTYPE,NCODE,NMAX,RMIN,RMAX,
     %NIDEN,UNITX)
        CHARACTER(len=16) SYMC,NCODE(NMAX)
        INTEGER NTYPE(NMAX)
        REAL COORD(NMAX,3),RMIN(NIDEN),RMAX(NIDEN)

c Judge the symmetry code.
        IF(SYMC.EQ.'P4M') THEN
           AA=UNITX
           NPARM=NMAX
           do i=1,NIDEN
              RMIN(I)=RMIN(I)/AA
              RMAX(I)=RMAX(I)/AA
           enddo
c
c Loop over NPARM atoms. 
           DO I=1,NPARM
c Exchange to spacing unit.
              ZCOR=COORD(I,1)/AA
              XCOR=ABS(COORD(I,2))/AA
              YCOR=ABS(COORD(I,3))/AA
              RDMIN=RMIN(NTYPE(I))
c             RDMAX=RMAX(NTYPE(I))
C Translate all coordinates into irreduceable region.
              IF(XCOR.GT.0.5.or.ycor.gt.0.5) THEN
c                 write(99,*) 'Coordinates out side the unit cell'
                 XCOR=AMOD(XCOR,1.)
                 YCOR=AMOD(YCOR,1.)
c XCOR and YCOR may still larger .5, so judge one moew time.
                 IF(XCOR.GT.0.5) XCOR=ABS(XCOR-1.0)
                 IF(YCOR.GT.0.5) YCOR=ABS(YCOR-1.0)
c                 write(99,*) 'Coordinates translated into unit cell'
              ENDIF
c Bring coordnates into an irreducable triangle reigion.
              IF(YCOR.GT.XCOR) THEN
                 XX=XCOR
                 XCOR=YCOR
                 YCOR=XX
              ENDIF
c Converge atoms belong to edge, corner, center.
              RR=SQRT(XCOR*XCOR+YCOR*YCOR)
              PI=3.1415926
              IF(XCOR.NE.0.0) THEN
                 THETA=ATAN(YCOR/XCOR)
              ELSE
                 THETA=PI/4.
              ENDIF
              RA=YCOR
              RB=0.5-XCOR
              RC=RR*SIN(PI/4.-THETA)
              IF(RA.GT.RDMIN.AND.RB.GT.RDMIN.AND.RC.GT.RDMIN) THEN
                 XCOR=XCOR
                 YCOR=YCOR
                 NCODE(I)='G'
                 GOTO 11
              ENDIF
              IF(RA.GE.RDMIN.AND.RB.GE.RDMIN.AND.RC.LE.RDMIN) THEN
                 CPROJX=RR*COS(PI/4.-THETA)
                 XCOR=CPROJX*SQRT(2.)/2.
                 YCOR=XCOR
                 NCODE(I)='F'
                 GOTO 11
              ENDIF
              IF(RA.GT.RDMIN.AND.RB.LE.RDMIN.AND.RC.GE.RDMIN) THEN
                 XCOR=.5
                 YCOR=YCOR
                 NCODE(I)='E'
                 GOTO 11
              ENDIF
              IF(RA.LE.RDMIN.AND.RB.GT.RDMIN.AND.RC.GT.RDMIN) THEN
                 XCOR=XCOR
                 YCOR=0.0
                 NCODE(I)='D'
                 GOTO 11
              ENDIF
              IF(RA.LE.RDMIN.AND.RB.LT.RDMIN.AND.RC.GT.RDMIN) THEN
                 XCOR=0.5
                 YCOR=0.0
                 NCODE(I)='C'
                 GOTO 11
              ENDIF
              IF(RA.LE.RDMIN.AND.RB.GT.RDMIN.AND.RC.LT.RDMIN) THEN
                 XCOR=0.0
                 YCOR=0.0
                 NCODE(I)='A'
                 GOTO 11
              ENDIF
              IF(RA.GE.RDMIN.AND.RB.LT.RDMIN.AND.RC.LT.RDMIN) THEN
                 XCOR=.5
                 YCOR=.5
                 NCODE(I)='B'
                 GOTO 11
              ENDIF
11            COORD(I,1)=ZCOR
              COORD(I,2)=XCOR
              COORD(I,3)=YCOR
           ENDDO
        ENDIF
           DO I=1,NPARM
              DO J=1,3
                 COORD(I,J)=COORD(I,J)*AA
              ENDDO
cjcm        write(*,*) 'COORD=',(coord(i,j),j=1,3),ncode(i)
           ENDDO
c        stop
           DO I=1,5
              RMIN(I)=RMIN(I)*AA
              RMAX(I)=RMAX(I)*AA
c          write(*,*) 'RMIN,RMAX=',rmin(i),rmax(i)
           ENDDO
        RETURN
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************c
C Valuate the validity of structure by examning the minimum 
c distance and touching priciple.
        SUBROUTINE VALUAT(COOR,NTYP,NMA,RMIN,RMAX,NIDEN,RESULT,
     &COORSUB,NTYPSUB,NSUB)
        CHARACTER(len=16) RESULT
        PARAMETER (NMAX=40)
        REAL COOR(NMA,3)
        REAL COORSUB(NSUB,3)
        REAL RMIN(NIDEN),RMAX(NIDEN) 
        REAL COORD(NMAX,3)
        INTEGER NTYP(NMA),NTYPSUB(NSUB),NTYPE(NMAX) 
        NATOM=NMA+NSUB
        DO I=1,NATOM
           IF(I.LE.NSUB) THEN
              DO J=1,3
                 COORD(I,J)=COORSUB(I,J)
              ENDDO
              NTYPE(I)=NTYPSUB(I)
           ELSE
              K=I-NSUB
              DO J=1,3
                 COORD(I,J)=COOR(K,J)
              ENDDO
              NTYPE(I)=NTYP(K)
           ENDIF
        ENDDO 
        DO I=NSUB+1,NATOM
           DO J=1,NATOM
              IF(I.EQ.J) GOTO 10 
              RIJS=RMIN(NTYPE(I))+RMIN(NTYPE(J))
              A1=(COORD(I,1)-COORD(J,1))**2
              A2=(COORD(I,2)-COORD(J,2))**2
              A3=(COORD(I,3)-COORD(J,3))**2
              DIST=SQRT(A1+A2+A3)
              IF(DIST.LT.RIJS) THEN
c                 WRITE(99,*) 'DIST,RIJS',DIST,RIJS
c                 WRITE(99,*) 'PENALTY FOR ATOM I AND J',I,J
                 RESULT='PENALTY'
                 RETURN
              ENDIF
10         ENDDO
        ENDDO
c
        DO I=NSUB+1,NATOM
           KBONUS=0
           DO J=1,NATOM
              IF(I.EQ.J) GOTO 20 
              RIJL=RMAX(NTYPE(I))+RMAX(NTYPE(J))
              A1=(COORD(I,1)-COORD(J,1))**2
              A2=(COORD(I,2)-COORD(J,2))**2
              A3=(COORD(I,3)-COORD(J,3))**2
              DIST=SQRT(A1+A2+A3)
              IF(DIST.LE.RIJL) THEN
                 KBONUS=KBONUS+1
              ENDIF
20         ENDDO
           IF(KBONUS.LT.1) THEN
c              WRITE(99,*) 'ATOM DOESNT TOUCH ANY OTHER ATOMS',I
              RESULT='PENALTY'
              RETURN
           ENDIF
        ENDDO
        RESULT='SUCCESS'
c        WRITE(99,*) '          Success!!'
        RETURN
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************
        SUBROUTINE COORDAT(COORD,NTYPE,NCODE,K1,K2,KK,
     %COORA,NTYPA)
        PARAMETER (NMAX=14,NNMAX=80)
        CHARACTER(len=16) NCODE(NMAX) 
        REAL COORD(NMAX,3),COORA(NNMAX,3)
        INTEGER NTYPE(NMAX),NTYPA(NNMAX) 
        KK=0  
        DO K=K1,K2
           IF(NCODE(K).EQ.'G') THEN
              DO I=1,2
                 DO J=1,2
                    KK=KK+1
                    COORA(KK,1)=COORD(K,1)
                    COORA(KK,2)=COORd(K,2)*(-1)**I
                    COORA(KK,3)=COORD(K,3)*(-1)**J 
                    NTYPA(KK)=NTYPE(K)  
                 ENDDO
              ENDDO
              DO I=1,2
                 DO J=1,2
                    KK=KK+1
                    COORA(KK,1)=COORD(K,1)
                    COORA(KK,2)=COORd(K,3)*(-1)**I
                    COORA(KK,3)=COORD(K,2)*(-1)**J   
                    NTYPA(KK)=NTYPE(K)  
                 ENDDO
              ENDDO
           ENDIF
           IF(NCODE(K).EQ.'F') THEN
              DO I=1,2
                 DO J=1,2
                    KK=KK+1
                    COORA(KK,1)=COORD(K,1)
                    COORA(KK,2)=COORd(K,2)*(-1)**I
                    COORA(KK,3)=COORD(K,3)*(-1)**J   
                    NTYPA(KK)=NTYPE(K)  
                 ENDDO
              ENDDO
           ENDIF
           IF(NCODE(K).EQ.'E') THEN
              DO J=1,2
                 KK=KK+1
                 COORA(KK,1)=COORD(K,1)
                 COORA(KK,2)=COORd(K,2)
                 COORA(KK,3)=COORD(K,3)*(-1)**J   
                    NTYPA(KK)=NTYPE(K)  
              ENDDO
              DO I=1,2
                    KK=KK+1
                    COORA(KK,1)=COORD(K,1)
                    COORA(KK,2)=COORd(K,3)*(-1)**I
                    COORA(KK,3)=COORD(K,2)
                    NTYPA(KK)=NTYPE(K)  
              ENDDO
           ENDIF
           IF(NCODE(K).EQ.'D') THEN
              DO J=1,2
                 KK=KK+1
                 COORA(KK,1)=COORD(K,1)
                 COORA(KK,2)=COORd(K,2)*(-1.)**J
                 COORA(KK,3)=COORD(K,3)   
                    NTYPA(KK)=NTYPE(K)  
              ENDDO
              DO I=1,2
                    KK=KK+1
                    COORA(KK,1)=COORD(K,1)
                    COORA(KK,2)=COORd(K,3)
                    COORA(KK,3)=COORD(K,2)*(-1.)**I   
                    NTYPA(KK)=NTYPE(K)  
              ENDDO
           ENDIF
           IF(NCODE(K).EQ.'C') THEN
                    KK=KK+1
                    COORA(KK,1)=COORD(K,1)
                    COORA(KK,2)=COORd(K,2)
                    COORA(KK,3)=COORD(K,3)   
                    NTYPA(KK)=NTYPE(K)  
                    KK=KK+1
                    COORA(KK,1)=COORD(K,1)
                    COORA(KK,2)=COORd(K,3)
                    COORA(KK,3)=COORD(K,2)   
                    NTYPA(KK)=NTYPE(K)  
           ENDIF
           IF(NCODE(K).EQ.'B'.OR.NCODE(K).EQ.'A') THEN
              KK=KK+1
              COORA(KK,1)=COORD(K,1)
              COORA(KK,2)=COORd(K,2)
              COORA(KK,3)=COORD(K,3)   
                    NTYPA(KK)=NTYPE(K)  
           ENDIF
        ENDDO
        RETURN       
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************
        SUBROUTINE SORTLOCAL(N,COOR,NTYP,NCOD)
C Sort an array RA of length N into ascending numberical order
c while making rearrangements of the arrays RB, RC and RZ. An 
c index table is constructed via the routine INDEXX.
        PARAMETER (NN=20)
        CHARACTER(len=16) NCOD(N),MWKSP(NN)
        REAL COOR(N,3),RA(NN),WKSP(NN)
        INTEGER IWKSP(NN),NWKSP(NN),NTYP(N)
        DO I=1,N
           RA(I)=COOR(I,1)
           IWKSP(I)=1
        ENDDO
        CALL INDEXX(N,RA,IWKSP)
        DO 11 J=1,N
            WKSP(J)=COOR(J,1)
11      CONTINUE
        DO 12 J=1,N
            COOR(J,1)=WKSP(IWKSP(J))
12      CONTINUE
        DO 13 J=1,N
            WKSP(J)=COOR(J,2)
13      CONTINUE
        DO 14 J=1,N
            COOR(J,2)=WKSP(IWKSP(J))
14      CONTINUE
        DO 15 J=1,N
            WKSP(J)=COOR(J,3)
15      CONTINUE
        DO 16 J=1,N
            COOR(J,3)=WKSP(IWKSP(J))
16      CONTINUE
        DO 17 J=1,N
            NWKSP(J)=NTYP(J)
17      CONTINUE
        DO 18 J=1,N
            NTYP(J)=NWKSP(IWKSP(J))
18      CONTINUE
        DO 19 J=1,N
            MWKSP(J)=NCOD(J)
19      CONTINUE
        DO 20 J=1,N
            NCOD(J)=MWKSP(IWKSP(J))
20      CONTINUE
        RETURN
        END

        SUBROUTINE INDEXX(N,ARRIN,INDX)
        REAL ARRIN(N)
        INTEGER INDX(N)
        DO 11 J=1,N
            INDX(J)=J
11      CONTINUE
        IF(N.EQ.1) RETURN
        L=N/2+1
        IR=N
10      CONTINUE
           IF(L.GT.1) THEN
              L=L-1
              INDXT=INDX(L)
              Q=ARRIN(INDXT)
           ELSE
              INDXT=INDX(IR)
              Q=ARRIN(INDXT)
              INDX(IR)=INDX(1)
              IR=IR-1
              IF(IR.EQ.1) THEN
                 INDX(1)=INDXT
                 RETURN
              ENDIF
           ENDIF
           I=L
           J=L+L
20         IF(J.LE.IR) THEN
              IF(J.LT.IR) THEN
                 IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1))) J=J+1
              ENDIF
              IF(Q.LT.ARRIN(INDX(J))) THEN
                 INDX(I)=INDX(J)
                 I=J
                 J=J+1
              ELSE
                 J=IR+1
              ENDIF
           GOTO 20
           ENDIF
           INDX(I)=INDXT
        GOTO 10
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C********************************************************       
        SUBROUTINE INDEXXN(N,NRRIN,INDX)
        INTEGER NRRIN(N),INDX(N)
        DO 11 J=1,N
            INDX(J)=J
11      CONTINUE
        IF(N.EQ.1) RETURN
        L=N/2+1
        IR=N
10      CONTINUE
           IF(L.GT.1) THEN
              L=L-1
              INDXT=INDX(L)
              NQ=NRRIN(INDXT)
           ELSE
              INDXT=INDX(IR)
              NQ=NRRIN(INDXT)
              INDX(IR)=INDX(1)
              IR=IR-1
              IF(IR.EQ.1) THEN
                 INDX(1)=INDXT
                 RETURN
              ENDIF
           ENDIF
           I=L
           J=L+L
20         IF(J.LE.IR) THEN
              IF(J.LT.IR) THEN
                 IF(NRRIN(INDX(J)).LT.NRRIN(INDX(J+1))) J=J+1
              ENDIF
              IF(NQ.LT.NRRIN(INDX(J))) THEN
                 INDX(I)=INDX(J)
                 I=J
                 J=J+1
              ELSE
                 J=IR+1
              ENDIF
           GOTO 20
           ENDIF
           INDX(I)=INDXT
        GOTO 10
        END
C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS: 
c GRAVSTRUC moves atoms in the direction perpendicular to the surface
c until a valid structure is obtained.
C********************************************************
        SUBROUTINE GRAVSTUC(COORD,NTYPE,NMAX,RMIN,RMAX,RMID,NIDEN,
     &COORSUB,NTYPSUB,NSUB)
        PARAMETER (NMA=14,NSUR=40)
        REAL COORD(NMAX,3),ATOM(3),COORD1(NMA,3)
        REAL COORSUB(NSUB,3)
        REAL RMIN(NIDEN),RMAX(NIDEN),RMID(NIDEN)
        REAL COORSUR(NSUR,3)
        INTEGER NTYPE(NMAX),NTYPSUB(NSUB),NTYPSUR(NSUR)
        INTEGER NATOM,MSUB
c Initialize a array COORSUR which contains both substrate and surface 
c atoms.
c       write(*,*) 'NSUB=', NSUB
        DO I=1,NSUR 
           DO J=1,3
              IF(I.GT.NSUB) THEN
                 COORSUR(I,J)=0.0
              ELSE
                 COORSUR(I,J)=COORSUB(I,J)
                 NTYPSUR(I)=NTYPSUB(I)
              ENDIF
           ENDDO
c          write(*,*) 'COORSUR=',(COORSUR(I,J),J=1,3) 
        ENDDO
        MSUB=NSUB

c Shift the surface atom vertically to each reasonable position.        
        DO I=NMAX,1,-1
           DO J=1,3
              ATOM(J)=COORD(I,J)
           ENDDO
c          write(*,*) 'ATOM=',(ATOM(J),J=1,3) 
           NATOM=NTYPE(I)
c Examine the distance of ATOM with the COORSUR atoms, and modify 
c ATOM(1) to a proper position.
c
       CALL EVALU(ATOM,NATOM,COORSUR,NTYPSUR,NSUR,MSUB,
     &RMIN,RMAX,RMID,NIDEN)
           MSUB=MSUB+1
           DO J=1,3
              COORSUR(MSUB,J)=ATOM(J)
              NTYPSUR(MSUB)=NTYPE(I)
           ENDDO
c Write the new set of coordinates to COORD1 array.
           DO J=1,3
              COORD1(I,J)=ATOM(J)
           ENDDO
        ENDDO
        DO I=1,NMAX
           DO J=1,3
              COORD(I,J)=COORD1(I,J)
           ENDDO
        ENDDO
        do i=1,nmax
c       write(*,*) (coord1(i,j),j=1,3)
        enddo
c       pause
        RETURN
        END

C********************************************************
C SUBROUTINE:  Name of function
C
C INPUT:     Description of input parameters
C 
C OUTPUT:    Descriptions of output parameters
C 
C CALLS      Functions called by this function
C
C RETURN:    Value returned by the function/method
C
C VARIABLES: List of variable names and short description
C
C COMMENTS:  General
C******************************************************** 
        SUBROUTINE EVALU(ATOM1,NAT,COOR,NTYPE,NSUR,MSUB,
     &RMIN,RMAX,RMID,NIDEN)
        PARAMETER (NSU=40)
        REAL ATOM1(3),COOR(NSUR,3)
        REAL DIS(NSU),RRS(NSU),RRL(NSU),RRD(NSU)  
        REAL RMIN(NIDEN),RMAX(NIDEN),RMID(NIDEN)
        INTEGER NTYPE(NSUR),NAT,N,NSUR
        N=MSUB   
        NL=0
        DO I=1,N
           DIS(I)=0.0
           DO J=1,3
              DIS(I)=DIS(I)+(ATOM1(J)-COOR(I,J))**2
           ENDDO
           DIS(I)=SQRT(DIS(I))
           RRS(I)=RMIN(NTYPE(I))+RMIN(NAT)
           RRL(I)=RMAX(NTYPE(I))+RMAX(NAT)
           RRD(I)=RMID(NTYPE(I))+RMID(NAT)
c      
           IF(DIS(I).LT.RRS(I)) THEN
c              write(99,*) 'Atom is too close to substrate'
c              write(99,*) 'coord(3)=',(atom1(j),j=1,3)
              ZMID=RRD(I)**2-(COOR(I,2)-ATOM1(2))**2-
     &(COOR(I,3)-ATOM1(3))**2
              ZMID=SQRT(ZMID)
              ATOM1(1)=COOR(I,1)-ZMID
c             KK=KK+1
           ENDIF
           IF(DIS(I).GT.RRL(I)) THEN
              NL=NL+1
           ENDIF
        ENDDO
c          write(*,*) 'evalu---1',i
        IF(NL.EQ.N) THEN
c           write(99,*) 'Atom is too far away from substrate'
           A0=100.0
           DO I=1,N
              IF(A0.GT.DIS(I)) THEN
                 A0=DIS(I)
                 IN=I
              ENDIF
           ENDDO
c           write(99,*) 'A0,IN=',A0,IN
      ZMID=RRD(IN)**2-(COOR(IN,2)-ATOM1(2))**2-(COOR(IN,3)-ATOM1(3))**2
C073198>>>>
c           ZMID=SQRT(ZMID)
c           ATOM1(1)=COOR(I,1)-ZMID
           IF(ZMID.LE.0.0) THEN 
               ATOM1(1)=COOR(IN,1)
           ELSE 
               ZMID=SQRT(ZMID)
               ATOM1(1)=COOR(I,1)-ZMID
           ENDIF
C073198<<<<
c          KK=KK+1
        ENDIF
        RETURN
        END
C======================================================================
C     SUBROUTINE INT2CHAR
C     Converts integer "x" to char string "c" of length "length"
C======================================================================

      SUBROUTINE INT2CHAR(x,c,length)

      integer x, length
      character c*(*)
      
      character(len=10) alph
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



