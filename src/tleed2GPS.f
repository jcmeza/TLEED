C  file LEEDSATL.M2  Feb. 29, 1996
C
C***************************************************************************
C  Symmetrized Automated Tensor LEED (SATLEED):  program TLEED2.FOR
C  Version 4.1 of Automated Tensor LEED
C
Cga      PROGRAM TLEED2V4
CGPS      subroutine tleed2(dir,rank,rfactor)
      subroutine tleed2(problem_dir, dir,rank,rfactor)
      character*(*) problem_dir
      character*(*) dir
      character*(*) rank
      real rfactor

CGPS      character*32, tleed4i,tleed5i,shortt,gleedo,searchs,prdird,
      character*100 tleed4i,tleed5i,shortt,gleedo,searchs,prdird,
     &restartd
      character*100 expdotd, rfacdotd

c      subroutine tleed2(tleed5i,shortt,rfactor)
c      character*(*), tleed5i, shortt

C 
C PROGRAM DESCRIPTION: 
C -------------------
C 
C        Second pass of Symmetric Tensor LEED program. 
C This code performs the automated
C search for an R-factor minimum. The IV spectra at a given point in 
C structural parameter space
C are evaluated as a pertubation expansion using the tensor elements generated
C by program TLEED1V4 and input via the transfer file SHORT.T.
C
C Changes due to error messages of Fortran Lint:
C *Rename array TST by TST1; call simpl, powell, powell2, grid with tst1 
C (instead of tst).
C *complex bj ==> complex*16 bj
C *ninset(10) ==> ninset(20)
C (Martin Gierer, 7/97)
C
C AUTHORS: 
C -------
C 
C        A. Barbieri. Modified version of TLEED2V3 by Barbieri.
C                     Original version TLEED2V1.1 by Rous and Wander.
C                     Further modifications by R. Doell.
C 
C CREATION DATE:        February 1996
C -------------
C 
C REFERENCE
C ---------
C
C J. A. Nelder and R. Mead   Computer Journal 7,(1965),308.
C
C ============================================================================
C
C Parameter statements for array dimensions carried over from the 
C Van Hove/Tong package.
C
C ============================================================================
C
C       IPLMAX          The largest l value to be used (LMAX).
C
C ============================================================================
C
      PARAMETER (IPLMAX=8)
C
C ============================================================================
C
C Additional parameter statements for composite layers.
C
C ============================================================================
C
C     INLAY         Maximum number of subplanes in any single composite layer.
C     INLTOT        Total number of subplanes in all composite layer
C                   for which the tensor is collected + # nostructural
C                   parameters.
C     INLIN         Total number of inequivalent subplanes in all composite 
C                   layer (with tensor)
C                   for which the tensor is collected.
C     INTAU         No of different elements for which phase-shifts are
C                   provided
C     INST1         Number of composite Layers for which data is input.
C
C ============================================================================
C
      PARAMETER (INLAY=66,INLTOT=68,INTAU=4,INST1=5)
      PARAMETER (INLIN=25)
C
C ============================================================================
C
C Parameter statemenst specific to TLEED
C
C ============================================================================
C
C INT0              Number of beams for which the tensor is to be calculated
C                   (Number of exit beams)
C INBED             Maximum number of beams included in the theoretical
C                   or experimental data set
C IPCUT             Number of elements in Q (from output from TLEED1)
C IINERG            Number of energy points in theoretical calculation
C JSMAX             Maximum value of L used in single centre expansion
C IEERG             Maximum number of experimental data points after
C                   interpolation.
C
C ============================================================================
C
      PARAMETER (INT0=20,IPCUT=45,IINERG=100)
      PARAMETER (INBED=32,IEERG=1600,JSMAX=4)
C
C ============================================================================
C
C The following parameters are constructed from the Van Hove/Tong parameters.
C
C ============================================================================
C
C       IPNL1          }Superlattice characterization codes
C       IPNL2          }Dimensions are chosen sufficiently large to treat
C                      }a 10x10 superstructure
      PARAMETER (IPNL1=10,IPNL2=10)
      PARAMETER (JPNL=IPNL1*IPNL2)
      PARAMETER (IINCUT=IPCUT*INLIN*IINERG*INT0)
      PARAMETER (JNDIM=INLTOT*3+1,INDIM=JNDIM+1)
      PARAMETER (JSMMAX=(JSMAX+1)*(JSMAX+1),JSMAX1=JSMAX+1)
C
C ============================================================================
C
C Dimension Arrays.
C
C ============================================================================
C
      DIMENSION P(JNDIM),PT(JNDIM),PTT(JNDIM)
      DIMENSION XIT(JNDIM),XI(JNDIM,JNDIM)
      DIMENSION COORD(INDIM,JNDIM),VAL(INDIM),PR(JNDIM),PRR(JNDIM)
      DIMENSION PBAR(JNDIM),INDOM(18)
      DIMENSION ILOOK(12,INLTOT),ILOOKS(12,INLTOT)
      DIMENSION ILKBD(12,INLAY,INST1)
      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),ASA(10,3),ARB1(2)
      DIMENSION ARB2(2),RBR1(2),RBR2(2),ACOORD(12,INLTOT,3)
      DIMENSION V(JPNL,2),JJS(JPNL,6)
      DIMENSION NTAUAW(INST1),MICUT(IPCUT),MJCUT(IPCUT)
      DIMENSION ES(90),PHSS(90,80),PHSS2(90,80),ETH(INT0,IEERG)
      DIMENSION FPOS(INLAY,3),VPOS(INST1,INLAY,3)
      DIMENSION LAFLAG(INST1),LLFLAG(INLTOT)
      DIMENSION LAFLAG2(INLTOT,INST1),LPOINT(INLIN),LPBD(INLAY,INST1)
      DIMENSION NINEQ(INST1),LPSAW(INST1,INLAY)
      DIMENSION POSS(INLAY,3),POSS2(INLAY,3),ADISP(INLTOT,3)
      DIMENSION PQFEX(2,INT0),NINSET(20),IBP(INT0)
      DIMENSION WR(10),WB(INT0),AT(INT0,IEERG)
      DIMENSION DISP(INLTOT,3),PSQ(2,INT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),NBEA(INBED)
      DIMENSION BENAME(5,INBED),XPL(IEERG),YPL(IEERG),NNN(IEERG)
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),TSEPP2(INBED),TSEY2(INBED)
      DIMENSION ATP(INT0,IEERG),ATPP(INT0,IEERG),TST1(INT0),TSTY2(INT0)
      DIMENSION NST3(INT0),NST2(INT0),RAV(INT0),IBK(INT0),EET(INT0)
      DIMENSION ROS(INT0),R1(INT0),R2(INT0),RP1(INT0),RP2(INT0)
      DIMENSION RPP2(INT0),RRZJ(INT0),RMZJ(INT0),RPE(INT0),NET(INT0)
      DIMENSION AR(11),RPP1(INT0),YT(INBED,IEERG)
      DIMENSION CPVPOS(INST1,INLAY,3),CPARB1(2),CPARB2(2)
      DIMENSION LSFLAG(INLTOT),NDIML(INLTOT),DIREC(INLTOT,2,2)
      DIMENSION LLF2(INLTOT),LLF3(INLTOT)
      DIMENSION NL9(INLTOT,INLIN),NEQ(INLIN)
      DIMENSION VICL(INST1),VCL(INST1),FRCL(INST1),ASB(INST1,3)
C
C Now complex arrays
C
      COMPLEX VL(JPNL,2),QS(IINCUT),XISTS(INT0,IINERG)
C     COMPLEX JYLM(JSMMAX) ,BJ(JSMAX1)
      COMPLEX JYLM(JSMMAX)
      COMPLEX*16 BJ(JSMAX1)
      COMPLEX YLM(JSMMAX),XIST(INT0,IINERG)
C
C Final character arrays
C
      CHARACTER*4 TITLE(20)
C
C ============================================================================
C
C Common blocks
C
C ============================================================================
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
      COMMON /MS/LMAX
      COMMON /LO/CPARB1,CPARB2
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /ZMAT/IANZ(40),IZ(40,4),BL(40),ALPHA(40),BETA(40),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(40),ALPHAS(40),BETAS
     & (40),PHIR,PHIM1,PHIM2
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA,FI
      COMMON /ENY/EI,EF,DE,NERG,NSYM,NDOM,VINER,VIMAG
      COMMON /DIM/NDOMS,NSYMS
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,NCUT
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /POW/IFUNC,MFLAG,SCAL
      common /pdir/tleed4i
C
C==========================================================================
C
C Data block. INDOM is the number of symmetry equivalent domains for
C symmetry code ISYAW.
C
C==========================================================================
C
      DATA INDOM(1),INDOM(2),INDOM(3),INDOM(4),INDOM(5),INDOM(6),INDOM
     & (7),INDOM(8),INDOM(9),INDOM(10),INDOM(11),INDOM(12),INDOM(13),
     & INDOM(14),INDOM(15),INDOM(16),INDOM(17),INDOM(18)/1,2,2,2,2,2,4,
     & 4,4,8,3,6,6,2,4,4,6,12/
C
C ============================================================================
C
C  Format statements.
C
C ============================================================================
C
101   FORMAT (20I3)
102   FORMAT (/8H IDEG = ,1I3)
103   FORMAT (/,10X,25HCOORDINATES AFTER SORTING,/)
104   FORMAT (I3,F7.4)
105   FORMAT (/' 2ND PASS TLEED...PERTUBATIVE LEED CALCULATION ')
106   FORMAT (20H CORRECT TERMINATION)
107   FORMAT (20A4)
108   FORMAT (2I3,2F7.4)
109   FORMAT (3F7.4)
110   FORMAT (' CPU TIME = ',F9.2,' SECONDS ')
114   FORMAT (' ERROR: The parameter INLIN must be at least',I4 )
115   FORMAT (' ERROR: The parameter INLTOT must be at least',I4 )
116   FORMAT (' ERROR: The parameter INT0 must be at least',I4 )
119   FORMAT (' ERROR: The parameter INTAU must be at least',I4 )
120   FORMAT (' ERROR: The parameter INST1 must be at least',I4 )
121   FORMAT (' ERROR: The parameter INLAY must be at least',I4 )
123   FORMAT (' ERROR: The parameter IPCUT must be at least',I4 )
124   FORMAT (' ERROR: The parameter JSMAX must be at least',I4 )
C
      write(*,*) 'tleed2GPS: Entering tleed2'
C ============================================================================
C
C Open I/O channels.
C
C ============================================================================
C
CC
cjcm
      write(*,*) 'tleed2GPS: problem_dir = ', problem_dir
      write(*,*) 'tleed2GPS: dir         = ', dir 
      write(*,*) 'tleed2GPS: rank        = ', rank
cjcm 
c First input channels
C
cjcm
c      tleed4i  = problem_dir//'/work'//dir//'/tleed4i'//rank
c      tleed5i  = problem_dir//'/work'//dir//'/tleed5i'//rank
c      expdotd  = problem_dir//'/exp.d'
c      rfacdotd = problem_dir//'/rfac.d'
c      shortt   = problem_dir//'/work'//dir//'/shortt'//rank

      tleed4i  = './workbb1/tleed4i000'
      tleed5i  = './workbb1/tleed5i000'
      expdotd  = './workbb1/exp.d'
      rfacdotd = './workbb1/rfac.d'
      shortt   = './workbb1/shortt000'

      OPEN (UNIT=4,FILE=tleed4i,STATUS='OLD')
      OPEN (UNIT=5,FILE=tleed5i,STATUS='OLD')
      OPEN (UNIT=11,FILE=expdotd,STATUS='OLD')
      OPEN (UNIT=12,FILE=rfacdotd,STATUS='OLD')
      OPEN (UNIT=22,FILE=shortt,STATUS='OLD',form='unformatted')
      rewind(4)
      rewind(5)
      rewind(11)
      rewind(12)
      rewind(22)

C
C Now output channels
C
Cga+      OPEN (UNIT=1,FILE='gleed.o',STATUS='NEW')
C      OPEN (UNIT=2,FILE='search.s',STATUS='NEW')
C      OPEN (UNIT=15,FILE='prdir.d',STATUS='UNKNOWN')
Cga-      OPEN (UNIT=10,FILE='restart.d',STATUS='UNKNOWN')
c      gleedo = dir // '/gleedo' // rank
c      searchs = dir // '/searchs' // rank
c      prdird = dir // '/prdird' // rank
c      restartd = dir // '/restartd' // rank
CGPS      gleedo = 'work'//dir//'/gleedo'//rank
CGPS      searchs = 'work'//dir//'/searchs'//rank
CGPS      prdird = 'work'//dir//'/prdird'//rank
CGPS      restartd = 'work'//dir//'/restartd'//rank
c      gleedo = problem_dir//'/work'//dir//'/gleedo'//rank
c      searchs = problem_dir//'/work'//dir//'/searchs'//rank
c      prdird = problem_dir//'/work'//dir//'/prdird'//rank
c      restartd = problem_dir//'/work'//dir//'/restartd'//rank

      gleedo   = './workbb1/gleedo000'
      searchs  = './workbb1/searchs000'
      prdird   = './workbb1/prdird000'
      restartd = './workbb1/restartd000'

      OPEN (UNIT=1,FILE=gleedo,STATUS='unknown')
C      OPEN (UNIT=2,FILE=searchs,STATUS='UNKNOWN')
      OPEN (UNIT=15,FILE=prdird,STATUS='UNKNOWN')
      OPEN (UNIT=10,FILE=restartd,STATUS='UNKNOWN')
      rewind(1)
      rewind(2)
      rewind(15)
      rewind(10)

c      OPEN (UNIT=1,FILE=,STATUS='unknown')
c      OPEN (UNIT=2,FILE='/dev/null',STATUS='UNKNOWN')
c      OPEN (UNIT=15,FILE='/dev/null',STATUS='UNKNOWN')
c      OPEN (UNIT=10,FILE='/dev/null',STATUS='UNKNOWN')


      write(*,*) 'tleed2GPS: start of executable code'
C
C ============================================================================
C
C Start of Executable Code
C
C ============================================================================
CGSS: VOPT must be initiallized. ZZ 9/29/04 
       VOPT=0.0 
C
C      I1=MCLOCK()
Cgps:	
C      write(99, *) 'this is just a test',10000
      READ (5,107) (TITLE(I),I=1,20)
      WRITE (1,107) (TITLE(I),I=1,20)
      WRITE (1,105)
C
C Read in pass number, print control parameter and dummy parameter for
C consistency with TLEED2
C
      READ (4,101) IPR,ISTART,LRFLAG
C
C NSYM is symmetry code of surface, ASTEP is initial size of simplex
C       
      READ (4,108) NSYM,NSYMS,ASTEP,VSTEP
      NDOM=INDOM(NSYM)
      NDOMS=INDOM(NSYMS)
C
C Read in number of beams NT0, number of beam sets NSET, and the cut off
C radii for the pertubation expansion LSMAX and LLCUT
C
      READ (4,101) NT0,NSET,LSMAX,LLCUT
      IF(LSMAX.GT.JSMAX) THEN
	 WRITE(1,124) LSMAX
	 GOTO 1000
      ENDIF
      IF(NT0.GT.INT0) THEN
	 WRITE(1,116) NT0
	 GOTO 1000
      ENDIF
C
C  IDEG rotational symmetry of each layer.
C
C  NL1, NL2  Superlattice Characterization (calculated in READT1)
C
      READ (5,101) IDEG
      IF (IPR.GT.0) THEN
	 WRITE (1,102) IDEG
      ENDIF
C
C  NPSI   = No. of energies at which phase shifts are read in.
C
      READ (5,101) NPSI
C
C Read in geometry, physical parameters and convergence criteria.
C
      CALL READT1(TVA,RAR1,RAR2,INVECT,TVB,IDEG,NL,TST,
     & THETA,FI,LMMAX,NPSI,ES,PHSS,L1,IPR,NEL,
     & IPNL1,IPNL2)
C
      CALL READT2(ASA,INVECT,IDEG,NL,V,VL,JJS,IPR)
C
      IF(NEL.GT.INTAU) THEN
	 write(1,119) NEL
	 GOTO 1000
      ENDIF
       VO=.0
C
C NST1   = Number of composite layers in input.
C LAFLAG = Number of layers in each composite layer.
C NLAY   = Max. no. of subplanes in any composite layer.(ie MAX(LAFLAG)
C NNST   = number of non structural parameter 
C 
      NNST=1
      NLAY=0
      NLAYTOT=0
Cga this variable must be initialized for repeated calls
      NLTOT2=0
      READ (5,101) NST1,NSTEF,NDTENS
      IF(NST1.GT.INST1) THEN
	 WRITE(1,120) NST1
	 GOTO 1000
      ENDIF
      NBULK=NST1-NSTEF
      READ (5,101) (LAFLAG(I),I=1,NST1)
      DO 348 I=1,NST1
	 IF (LAFLAG(I).GT.NLAY) NLAY=LAFLAG(I)
	 NLTOT2=NLTOT2+LAFLAG(I)
	 IF (I.LE.NSTEF)NLAYTOT=NLAYTOT+LAFLAG(I)
348   CONTINUE
      IF(NLAYTOT+NNST.GT.INLTOT) THEN
	 WRITE(1,115) NLAYTOT+NNST
	 GOTO 1000
      ENDIF
      NLMAX=NLAY
      IF(NLMAX.GT.INLAY) THEN
	 WRITE(1,121) NLMAX
	 GOTO 1000
      ENDIF
C
C Read in additional data for composite layer.
C
      CALL READCT(NLAY,VPOS,FPOS,POSS,CPVPOS,NTAUAW,LPSAW,IPR,
     & LAFLAG,NST1,ASB,VICL,VCL,FRCL,TST,TSTS,ASA,INVECT)
C
C Set up lookup table for domain averaging
C
      NSHIFT=0
      WRITE(2,103)
CGPS
C      write(0,*) 'NSTEF=', NSTEF
      DO 351 NCL=1,NSTEF
C      write(0,*) 'NCL=', NCL
	NLAY2=LAFLAG(NCL)
	CALL LOOKUP2(CPVPOS,ILOOK,ILKBD,NLAY2,NST1,NSTEF,NSYM,
     &  NLTOT2,NSHIFT,NCL,NLMAX,IROT,IMIR,1,FPOS,POSS,POSS2,
     &  nerror_report)
C        write(0,*)'FPOS=',FPOS
C
C Set up lookup table for restricting the search
C
	CALL LOOKUP2(CPVPOS,ILOOKS,ILKBD,NLAY2,NST1,NSTEF,NSYMS,
     &  NLTOT2,NSHIFT,NCL,NLMAX,IROT,IMIR,0,FPOS,POSS,POSS2,
     &  nerror_report)
	NSHIFT=NSHIFT+NLAY2
C        write(0,*)'FPOS=',FPOS
351   CONTINUE
C
C Read in energy range and step
C
      READ (5,109) EI,EF,DE
      NERG=INT((EF-EI)/DE+1.01)
      VINER=VV
      VIMAG=VPIS
C
C Read in information relevant to the pertubative LEED calculation
C
      CALL READPL(NT0,NSET,PQFEX,NINSET,NDIM,DISP,
     & NLAYTOT,IPR,AALPHA,BBETA,GGAMMA,ITMAX,
     & FTOL1,FTOL2,MFLAG,LLFLAG,NGRID)
C
C Read in experimental IV curves and information relevant to the R-factor
C calculation.
C
      CALL RFIN(IBP,NT0,WB,WR,IPR)
      CALL EXPAN(INBED,IEERG,AE,EE,NEE,NBEA,BENAME,IPR,XPL,YPL,NNN,
     & AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,VPIS)
C
C Set up lpoint and NLTIN for the reduced tensor calculation
C
      CALL DIMSCH2(ILOOKS,NLAYTOT,NDOMS,NNDIM,LSFLAG,LLFLAG,
     & NDIML,DIREC,ADISP,ACOORD,NDIM,1,LLF2,LLF3,NL9,NEQ,
     & INLIN)
C
      CALL INEQCD(LSFLAG,LAFLAG,NST1,NLMAX,NLTOT2,
     &LAFLAG2,LPOINT,NLTIN,NSTEF,NINEQ,LPBD)
C
      IF(NLTIN.GT.INLIN) THEN
	 WRITE(1,114) NLTIN
	 GOTO 1000
      ENDIF
C
C Set up dimensionality of search acoording  to the symmetry parameter
C NSYMS. DIMSCH deals with structural parameters only 
C
      CALL DIMSCH2(ILOOKS,NLAYTOT,NDOMS,NNDIM,LSFLAG,LLFLAG,
     & NDIML,DIREC,ADISP,ACOORD,NDIM,2,LLF2,LLF3,NL9,NEQ,
     & INLIN)

C
C the last parameters are non structural (inner potential, etc.) 
C
      NNSTEF=0
      IF (NNST.GE.1) THEN
	 DO 1115 I=NLAYTOT+1,NLAYTOT+NNST
	    IF (LLFLAG(I).NE.0) THEN 
	      NNSTEF=NNSTEF+1
	    ENDIF
1115      CONTINUE
      ENDIF
C
C NNSTEF is the number of nonstructural parameters included in the search
C
      NNDIM=NNDIM+NNSTEF
      NNDIM1=NNDIM+1
C
C Call MSET to obtain the reduction information for Q
C
      LSMMAX=(LSMAX+1)*(LSMAX+1)
      LSMAX1=LSMAX+1
      CALL MSET2(LSMAX,LLCUT,ICUT,MICUT,MJCUT,IPCUT,NDTENS)
C      IF(ICUT.GT.IPCUT) THEN
C         WRITE(1,123) ICUT
C         GOTO 1000
C      ENDIF
      NCUT=ICUT*NLTIN*NERG*NT0
C
C Retrieve the tensor and plane-wave amplitudes from the transfer file
C
      CALL RETR2(NCUT,ICUT,NLTIN,NT0,XISTS,QS,THETA,FI,PSQ,PQFEX,
     & RAR1,RAR2)
C
C Perform a full optimization. What method do we use?
C
C MFLAG=1   Down-hill Simplex method of Nelder and Mead
C MFLAG=2   Powell direction-set method 
C MFLAG=3   Powell modified to compute principal directions 
C MFLAG=4   Evaluation of the rfactor along the principal direction
C           determined by 3 to check error
      write(*,*) 'tleed2GPS: MFLAG = ', MFLAG
      IF (MFLAG.EQ.1) THEN
C
C Start off by performing a simplex search until convergence by
C by FTOL1 is achieved
C
Cga make a function call for the GA
Cga      CALL SIMPL(NLAYTOT,NDIM,COORD,DISP,NNDIM,NNDIM1,ASTEP,VSTEP,PR,
      rfactor = SIMPL(NLAYTOT,NDIM,COORD,DISP,NNDIM,NNDIM1,
     & ASTEP,VSTEP,PR,
     & ADISP,FTOL1,ITMAX,PBAR,VAL,AALPHA,BBETA,GGAMMA,PRR,ISTART,IPR,
     & ILO,ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST1,TSTY2,NST3,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,
     & NDIML,DIREC,NLTIN,LPOINT)
      ELSEIF (MFLAG.EQ.2) THEN
C
C Use direction set
C

Cga make a function call for the GA
C	open (unit=29, file='tleed2.out',status='unknown')
C	write(29,*)'indside tleed2:before powell'
C	write(*,*)'indside tleed2:before powell'
C        write(29,*) 'now output the input position parameter'
C        write(*,*) 'now output the input position parameter'
C        do i=1,NNDIM
C                write(29,*) P(i)
C                write(*,*) P(i)
C        enddo
Cga       CALL POWELL(P,XI,NLAYTOT,NDIM,DISP,NNDIM,ADISP,FTOL2,ASTEP,
         write(*,*) 'tleed2GPS: Calling Powell'
       rfactor = POWELL(P,XI,NLAYTOT,NDIM,DISP,NNDIM,ADISP,
     & FTOL2,ASTEP,  
     & VSTEP,ITMAX,ISTART,IPR,PTT,PT,XIT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST1,TSTY2,NST3,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,
     & NDIML,DIREC,NLTIN,LPOINT,LRFLAG)

C	write(29,*)'indside tleed2:after powell'
C	write(*,*)'indside tleed2:after powell'
C        write(29,*) 'now output the final position parameter'
C        write(*,*) 'now output the final position parameter'
C        do i=1,NNDIM
C                write(29,*) P(i)
C                write(*,*) P(i)
C        enddo
C        write(29,*) rfactor
C        write(*,*) rfactor
C	close(29)

      ELSEIF (MFLAG.EQ.3) THEN
C
C Use principal direction method
C
Cga make a function call for the GA
Cga       CALL POWELL2(P,XI,NLAYTOT,NDIM,DISP,NNDIM,ADISP,FTOL2,ASTEP,
         write(*,*) 'tleed2GPS: Calling Powell2'
         rfactor = POWELL2(P,XI,NLAYTOT,NDIM,DISP,NNDIM,ADISP,
     & FTOL2,ASTEP,
     & VSTEP,ITMAX,ISTART,IPR,PTT,XIT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST1,TSTY2,NST3,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,
     & NDIML,DIREC,NLTIN,LPOINT)

      ELSE
	READ(4,104) MSTEP,ZSTEP
	SCAL=VSTEP/ASTEP
Cga make a function call for the GA
Cga      CALL GRID(ZSTEP,MSTEP,NLAYTOT,NDIM,DISP,XI,
      rfactor = GRID(ZSTEP,MSTEP,NLAYTOT,NDIM,DISP,XI,
     & NNDIM,ADISP,FTOL2,ITMAX,ISTART,IPR,PTT,PT,XIT,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST1,TSTY2,NST3,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,NGRID,
     & LSFLAG,NDIML,DIREC,NLTIN,LPOINT)
      ENDIF
      WRITE(1,106)
C      I2=MCLOCK()
C      FTIME=FLOAT(I2-I1)
C      CPU=FTIME/100.0
C      WRITE(2,110) CPU
Cga1000      STOP
 1000 continue
c      call ieee_flags("clear", "exception", "all", out)
      close(1)
      close(2)
      close(15)
      close(10)
      close(4)
      close(5)
      close(11)
      close(12)
      close(22)

      write(*,*) 'tleed2GPS: Leaving tleed2'

      return
      END
