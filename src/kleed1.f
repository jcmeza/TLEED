
C  file TOYF_230605.FOR
C-----------------------------------------------------------------------
c      PROGRAM KINLEED
c       subroutine kleed(problem_dir,dir,rank,nerror_report,rfactor)
       subroutine kleed(problem_dir,dir,rank,rfactor)


      character*(*) problem_dir
      character*(*) dir
      character*(*) rank
      character*100 kleed4i, kleed5i,  tleedo,gleedo,searchs,prdird,
     &restartd
      character*100 expdotd,rfacdotd
      character*100 xpin
      integer nerror_report
      real rfactor


C
C PROGRAM DESCRIPTION:
C -------------------
C
C       Kinematic LEED program. This code runs the calculation for
C the distorted structure and returns the IV spectra.
C
C AUTHORS:
C -------
C        A. Garcia-Lekue. Based on TLEED1V2 version of TLEED2 program
C                         by Wander.
C                         Original version of TLEED2 program (TLEED1V1.1)
C                         by Rous and Wander.
C
C CREATION DATE:         Spring 2005
C -------------
C ============================================================================
C
C Parameter statements for array dimensions carried over from the 
C Van Hove/Tong package.
C
C ============================================================================
C
C      IPNL1          }Superlattice characterization codes
C      IPNL2          }(NL1,NL2).
C      IPIDEG          Rotational symmetry (IDEG).
C      IPLMAX          The largest l value to be used (LMAX).
C
C ============================================================================
C
      PARAMETER (IPNL1=5,IPNL2=5,IPIDEG=4,IPLMAX=5,IPCLM=1925)
C
C ============================================================================
C
C Additional parameter statements for composite layers.
C
C ============================================================================
C
C     INLAY         Maximum number of subplanes in any composite layer.
C     INTAU         Max No of chemical elements in any composite layer.
C     INST1         Number of composite Layers for which data is input.
C     INLTOT	    Number of subplanes in composite layers that are
C 		    displaced from reference position, i.e., number
c  	            of subplanes in NEFFST1 composite layers
C
C ============================================================================
C
      PARAMETER (INLAY=58,INTAU=2,IPCAA=1820,INST1=1)
      PARAMETER (INLTOT=58,INLDISP=58)
C
C ============================================================================
C
C Parameter statemenst specific to TLEED
C
C ============================================================================
C
C     INT0          Number of beams for which the tensor is to be calculated
C                   (Number of exit beams)
C     JSMAX         Maximum L value used in expansion [=< sqrt(2*e)*mod(r)]
C     NROM          Maximum dimension of overlayer matrices(ROM,TOM,ROP,TOP)
C                   This is the maximum number of beams entering the 
C                   calculation at each energy.
C     IINERG        Number of energy points in theoretical calculation
C     INBED         Maximum number of beams included in the theoretical
C                   or experimental data set
C     IEERG         Maximum number of experimental data points after
C                   interpolation.
C     ITEMP	    Switches off (ITEMP=0) or on (ITEMP different from zero)
C		    temperature effects 
C ============================================================================
C
      PARAMETER (INT0=20,JSMAX=2)
      PARAMETER (NROM=3000,NROM2=NROM)
      PARAMETER (IINERG=100)
      PARAMETER (INBED=20,IEERG=1600)      
      PARAMETER (ITEMP=1)
C
C ============================================================================
C
C The following parameters are constructed from the Van Hove/Tong parameters.
C
C ============================================================================
C
      PARAMETER (JPN=2*IPLMAX+1,JPNN=JPN*JPN)
      PARAMETER (JPNN2=IPLMAX+1,JPNN3=IPLMAX+1,JPNN1=JPNN2+JPNN3-1)
      PARAMETER (JLMMAX=JPNN2*JPNN2,JPL1=JPNN2,JPLEV=JPNN2*(JPNN2+1)/2)
      PARAMETER (JPLOD=JLMMAX-JPLEV,JPNL=IPNL1*IPNL2)
C
C ============================================================================
C
C Dimension Arrays.
C
C ============================================================================
C
      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),ASA(10,3),ARB1(2)
      DIMENSION ARB2(2),RBR1(2),RBR2(2),ASB(INST1,3)
      DIMENSION V(JPNL,2),JJS(JPNL,IPIDEG)
      DIMENSION KNB(20),NB(20),SPQF(2,NROM),SPQ(2,NROM),PQF(2,NROM)
      DIMENSION PQ(2,NROM),NTAUAW(INST1)
      DIMENSION ES(90),PHSS(90,80),PHSS2(90,80)
      DIMENSION NLMS(9),CLM(IPCLM),YLM(JPNN),FAC2(JPNN),FAC1(JPN)
      DIMENSION PPP(JPNN1,JPNN2,JPNN3)
      DIMENSION CAA(IPCAA),NCA(9),LAFLAG(INST1)
      DIMENSION AK2M(INT0),AK3M(INT0)
      DIMENSION LX(JLMMAX),LXI(JLMMAX),LXM(JLMMAX),LT(JLMMAX)
      DIMENSION FPOS(60,3),VPOS(INST1,INLAY,3),LPS(60),LPSS(60)
      DIMENSION WPOS(INST1,INLAY,3),WPOSTF(100,3)
      DIMENSION SPOSTF(20,1500,3),SPOSTF1(30,3)
      DIMENSION CPVPOS(INST1,INLAY,3),LLFLAG(60)
      DIMENSION LPSAW(INST1,INLAY),PSQ(2,INT0)
      DIMENSION PQFEX(2,INT0),NINSET(20)
      DIMENSION IT1(5),DRPAR1(5),DR01(5),DRPER1(5)
      DIMENSION ADISP(INLTOT,3),DISP(INLDISP,3)
      DIMENSION G(2,NROM)
      DIMENSION VICL(INST1),VCL(INST1),FRCL(INST1)
      DIMENSION IELEMOL(INST1,INLTOT),IELEMOL2(INLTOT)
      DIMENSION WR(10),WB(INT0),IBP(INT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),NBEA(INBED)
      DIMENSION BENAME(5,INBED),XPL(IEERG),YPL(IEERG),NNN(IEERG)
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),TSEPP2(INBED),TSEY2(INBED)
C Now complex arrays
      COMPLEX AF(JPL1),CAF(JPL1),TSF0(6,16),TSF(6,16)
      COMPLEX VL(JPNL,2)
      COMPLEX PHSSEL0(IINERG,INTAU),PHS(16)
      COMPLEX DEL(16),PHSSEL(IINERG,INTAU,16)

C Final character arrays
C
      CHARACTER*4 TITLE(20)
C
C ============================================================================
C
C Common blocks carried over from Van Hove/Tong package. 
C
C ============================================================================
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
      COMMON /MS/LMAX
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /ZMAT/IANZ(40),IZ(40,4),BL(40),ALPHA(40),BETA(40),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(40),ALPHAS(40),BETAS
     & (40),PHIR,PHIM1,PHIM2
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /REXP/EEINCR
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
      COMMON /TEMP/IT1,TI,T01,DRPER1,DRPAR1,DR01
C ============================================================================
C
C Data statements carried over from Van Hove/Tong package.
C
C ============================================================================
C
C  NLMS is dimension of CLM as a function of LMAX.
C
      DATA NLMS(1),NLMS(2),NLMS(3),NLMS(4),NLMS(5),NLMS(6),NLMS(7),
     &NLMS(8),NLMS(9)/70,76,284,809,1925,4032,7680,13593,22693/
C
C  NCA is dimension of CAA as a function of LMAX.
C
      DATA NCA(1),NCA(2),NCA(3),NCA(4),NCA(5),NCA(6),NCA(7),NCA(8),
     &NCA(9)/70,70,264,759,1820,3836,7344,13053,21868/
C
C ============================================================================
C  Format statements.
C
C ============================================================================
C
340   FORMAT (16F7.4)
100   FORMAT (3F7.2)
101   FORMAT (20I3)
102   FORMAT (/8H IDEG = ,1I3,7H NL1 = ,1I3,7H NL2 = ,1I3)
103   FORMAT (8H VPIS = ,F9.4,8H VPIO = ,F9.4,9H DCUTS = ,F9.4,
     & 9H DCUTO = ,F9.4)
1103   FORMAT (/,10X,25HCOORDINATES AFTER SORTING,/)
104   FORMAT (8H TEMP = ,F9.4)
105   FORMAT (/' 1ST PASS TLEED...REFERENCE STRUCTURE CALCULATION ')
106   FORMAT (20H CORRECT TERMINATION)
107   FORMAT (20A4)
109   FORMAT (//,' STARTING LOOP OVER ENERGIES ',/)
110   FORMAT (/,'====================================================',
     & /)
111   FORMAT (' ERROR LSMAX BIGGER THAN DIMENSIONED')
112   FORMAT (2I3,2F7.4)
1010  FORMAT (7X,F12.8,2(5X,F12.8))


C
C ============================================================================
C
C Open I/O channels.
C
C ============================================================================
C
C First input channels
C
c      write(*,*) problem_dir,dir
      kleed4i  = '/Users/meza/MyProjects/TLEED/workbbk001//kleed4i000'
      kleed5i  = '/Users/meza/MyProjects/TLEED/workbbk001/kleed5i000'
      expdotd  = '/Users/meza/MyProjects/TLEED/workbbk001/exp.d'
      rfacdotd = '/Users/meza/MyProjects/TLEED/workbbk001/rfac.d'

      OPEN (UNIT=4,FILE=kleed4i,STATUS='OLD')
      OPEN (UNIT=5,FILE=kleed5i,STATUS='OLD')
      OPEN (UNIT=11,FILE=expdotd,STATUS='OLD')
      OPEN (UNIT=12,FILE=rfacdotd,STATUS='OLD')

      rewind(4)
      rewind(5)
      rewind(11)
      rewind(12)
C
C Now output channels
C
      tleedo = '/Users/meza/MyProjects/TLEED/workbbk001/tleedo000'
      searchs ='/Users/meza/MyProjects/TLEED/workbbk001/searchs000'

      OPEN (UNIT=1,FILE=tleedo,STATUS='unknown')
      OPEN (UNIT=2,FILE=searchs,STATUS='unknown')

      rewind(1)
      rewind(2)

c      rewind(7)

C
C ============================================================================
C
C Start of Executable Code
C
C ============================================================================
C
cjcm       write(*,*) 'starting kleed'
      READ (5,107) (TITLE(I),I=1,20)
      WRITE (1,107) (TITLE(I),I=1,20)
      WRITE (1,105)
C
C  EMACH is machine accuracy.
C
      EMACH=1.0E-6
C
C Read in pass number, print control parameter and dummy parameter for
C consistency with TLEED2
C
      READ (4,101) IPR,ISTART
C
C NSYM is symmetry code of surface
C       
      READ (4,112) NSYM,NSYMS,ASTEP,VSTEP
C
C Read in number of beams NT0, number of beam sets NSET, and the cut off
C radii for the pertubation expansion LSMAX and LLCUT
C
      READ (4,101) NT0,NSET,LSMAX,LLCUT
C
C  IDEG rotational symmetry of each layer.
C
C  NL1, NL2  Superlattice Characterization  NL1  NL2
C                      P(1*1)                1    1
C                      C(2*2)                2    1
C                      P(2*1)                2    1
C                      P(1*2)                1    2
C                      P(2*2)                2    2
C
      READ (5,101) IDEG,NL1,NL2
      IF (IPR.GT.0) WRITE (1,102) IDEG,NL1,NL2
      NL=NL1*NL2
C
C  NPSI   = No. of energies at which phase shifts are read in.
C
      READ (5,101) NPSI
C
C Read in geometry, physical parameters and convergence criteria.
C
      CALL READT(TVA,RAR1,RAR2,ASA,INVECT,TVB,IDEG,NL,V,VL,JJS,
     & TST,TSTS,THETA,FI,LMMAX,NPSI,ES,PHSS,PHSS2,L1,IPR,NEL)
C
C prepare for phase shift interpolation to be performed in TSCATF
C through a cubic spline interpolation routine from Numerical Recipes
C
      DO 660 I=1,NPSI
	DO 661 II=1,NEL
           IO=(II-1)*L1
661      CONTINUE
660   CONTINUE

      CALL FORSPLINE(NEL,L1,NPSI,PHSS,ES,PHSS2)


C
C PHSS contains the tabulated phase shifts, PHSS2 contains the second derivative
C of PHSS (as a function of the energy) which will be needed in TSCATF
C
C NST1   = Number of composite layers in input.
C LAFLAG = Number of layers in each composite layer.
C NLAY   = Max. no. of subplanes in any composite layer.(ie MAX(LAFLAG)
C
      NLAY=0
      NLAYTOT=0
c      READ (5,101) NST1,NST1EFF
      READ (5,*) NST1,NST1EFF
c      write (*,*) NST1,NST1EFF
      READ (5,101) (LAFLAG(I), I=1,NST1)
c      write (*,*) (LAFLAG(I), I=1,NST1)
      DO 348 I=1,NST1
	 IF (LAFLAG(I).GT.NLAY) NLAY=LAFLAG(I)
	 IF (I.LE.NST1EFF) NLAYTOT=NLAYTOT+LAFLAG(I)
348   CONTINUE	
C
C Read in additional data for composite layer.
C
cjcm       write(*,*) 'before READCT'
       CALL READCT(NLAY,VPOS,CPVPOS,NTAUAW,LPSAW,LMMAX,IPR,
     & LAFLAG,NST1,ASB,VICL,VCL,FRCL,TST,TSTS,ASA,INVECT)

c	do 1101 i=1,nst1
c	do 1102 j=1,laflag(i)
c	WRITE(7,*) vpos(i,j,1)
c1102	continue
c	do 2103 j=1,laflag(i)
c	WRITE(7,*) vpos(i,j,2)
c2103	continue
c	do 1104 j=1,laflag(i)
c	WRITE(7,*) vpos(i,j,3)
c1104	continue
c	do 1105 j=1,laflag(i)
c	WRITE(7,*) lpsaw(i,j) 
c1105	continue
c1101	continue

C
C Read in information relevant to the pertubative LEED calculation
C
      CALL READPL(NT0,NSET,PQFEX,NINSET,NDIM,DISP,ICOORD,
     & NSTEP,ANSTEP,NLAYTOT,IPR,ALPHA,BETA,GAMMA,ITMAX,FTOL1,
     & FTOL2,MFLAG,LLFLAG,NGRID)
cjcm        write(*,*) 'after READPL'
C
C Calculate Clebsch-Gordan coefficients
C
      T0=T01
      KLM=(2*LMAX+1)*(2*LMAX+2)/2
      LEV=(LMAX+1)*(LMAX+2)/2
      LOD=LMMAX-LEV
      LEV2=2*LEV
      NCAA=NCA(LMAX)
cjcm	write(*,*) 'before CAAA'
c	do 1589 i=1,nlaytot
c	write(*,*) vpos(1,i,1),vpos(1,i,2),vpos(1,i,3)
c1589	continue
cjcm
      CALL CAAA(CAA,NCAA,LMMAX)
cjcm	write(*,*) 'after CAAA'
c	do 1489 i=1,nlaytot
c	write(*,*) vpos(1,i,1),vpos(1,i,2),vpos(1,i,3)
cjcm1489	continue
      N=2*LMAX+1
      NN=N*N
      NLM=NLMS(LMAX)
      CALL CELMG(CLM,NLM,YLM,FAC2,NN,FAC1,N,LMAX)
cjcm	write(*,*) 'after CELMG'
c	do 1389 i=1,nlaytot
c	write(*,*) vpos(1,i,1),vpos(1,i,2),vpos(1,i,3)
cjcm1389	continue
C
C  Calculate permutations of (L-M) sequence.
C
      CALL LXGENT(LX,LXI,LT,LXM,LMAX,LMMAX)
cjcm	write(*,*) 'after LXGENT'
c	do 1289 i=1,nlaytot
c	write(*,*) vpos(1,i,1),vpos(1,i,2),vpos(1,i,3)
cjcm1289	continue
C
C Do temperature dependant phase shifts need calculating?
C
      IMARK=0
      DO 48 I=1,NEL
         IMARK=IMARK+IT1(I)
48    CONTINUE
      NN3=LMAX+1
      NN2=LMAX+1
      NN1=NN2+NN3-1
      IF (IMARK.GT.0) THEN
C
C  PPP= Clebsch Gordan coefficients for computation of temperature
C       dependant phase shifts. (Skipped if not needed).
C
cjcm	write(*,*) 'before CPPP'
         CALL CPPP(PPP,NN1,NN2,NN3)
      ENDIF

C
C Check size of JSMAX.
C
      IF (LSMAX.GT.JSMAX) THEN
         WRITE (1,111)
      ELSE

	T=TI
        IF (IPR.GT.0) WRITE (1,104) T
C =============================================================================
C
C Set up overlayer (WPOSTF) and first substrate layer (SPOSTF) atomic positions 
C
C =============================================================================
C
C Set up overlayer atom positions and phase shift assignments LPS
cjcm	write(*,*) 'before 547'
C
 	  ITOT=1	
	DO 547 NCL=1,NST1
            NLAY=LAFLAG(NCL)
            NTAU=NTAUAW(NCL)
	    LPSMAX=1
	      DO 548 I=1,NLAY
	          LPS(I)=LPSAW(NCL,I)
		  IF(LPS(I).GT.LPSMAX)LPSMAX=LPS(I)
		  DO 549 J=1,3
		     FPOS(I,J)=VPOS(NCL,I,J)
549               CONTINUE
548            CONTINUE

c Reorder the subplanes of each composite layer according to
c increasing position along the +X axis

cjcm	write(*,*) 'before SORT'


	 CALL SORT(FPOS,NLAY)

      WRITE(2,1103)
	DO 948 I=1,NLAY
	DO 947 J=1,3
	    VPOS(NCL,I,J)=FPOS(I,J)
947	CONTINUE
	  WRITE(2,1010) (FPOS(I,K),K=1,3)
948	continue


c Redefine all the atomic positions with respect to the surface plane
c located at a distance ASE from the outermost overlayer subplane.



	DO 648 I=1,NLAY
            WPOS(NCL,I,1)=VPOS(NCL,I,1)+ASE+ABS(VPOS(1,1,1))
            WPOS(NCL,I,2)=VPOS(NCL,I,2)
	    WPOS(NCL,I,3)=VPOS(NCL,I,3)
C  Assign element type 
	    IELEMOL(NCL,I)=LPSAW(NCL,I)
648	CONTINUE

C  Store the OL atomic position in WPOSTF (includes atoms
C  of all CL-s)
C ---------------------------------------------------------
        DO 220 I=1,NLAY
	   WPOSTF(ITOT,1)=WPOS(NCL,I,1)
	   WPOSTF(ITOT,2)=WPOS(NCL,I,2)
	   WPOSTF(ITOT,3)=WPOS(NCL,I,3)
           IELEMOL2(ITOT)=LPSAW(NCL,I)
	   ITOT=ITOT+1
220	CONTINUE

547	CONTINUE


C Set up atomic positions for first substrate layer
C--------------------------------------------------



cjcm	write(*,*) 'before SLPOS'
c	do 189 i=1,nlaytot
c	write(*,*) vpos(1,i,1),vpos(1,i,2),vpos(1,i,3)
c189	continue
cjcm
	CALL SLPOS(ASB,VPOS,LAFLAG,ASE,NST1,NLAYTOT,SPOSTF)



	DO 288 I=1,NL
	DO 289 J=1,3
	SPOSTF1(i,j) = SPOSTF(1,I,J)
289	CONTINUE
288	CONTINUE	
C
C
C =============================================================================
C
C Begin loop over energy range.
C
C =============================================================================
C
C  Read energy range and step.
C
	 READ (5,100) EI,EF,DE
         NERG=INT((EF-EI)/DE+1.01)
	 IF (EI.LT.0) THEN
	    WRITE (1,*) ' EI MUST BE > 0 '
	 ELSE
C
C Generate required beamsets
C
	    DFLAG=0
	    CALL BEMGEN(TST,EF,SPQF,SPQ,KNBS,KNB,RAR1,RAR2,KNT,
     &       IPR,TVA,DFLAG,NROM,G)
C
C  Start loop over given energy range.
C
	    NGAW=INT((EF-EI)/DE)+1
	    WRITE (1,109)
	    WRITE (1,110)
	    DO 1300 IEEV=1,NGAW
	       EEV=EI+(IEEV-1)*DE
	       E=EEV/27.21+VV
	       E3=E
	       WRITE (1,*) ' CALCULATING FOR E= ',EEV
	       WRITE (1,110)
C
C  Set imaginary part of muffin tin potential. (Usually -4 or -5 eV)
C
	       VPIO=VPIS
C
C =============================================================================
C
C Calculate atomic T matrix elements for all atom types.
C
C =============================================================================
C

	       DO 247 INNEL=1,NEL
		  CALL TSCATF_TOY(INNEL,L1,ES,PHSS,PHSS2,NPSI,IT1,E,0.,
     &             PPP,NN1,NN2,NN3,DR01,DRPER1,DRPAR1,T0,T,TSF0,TSF,
     &             AF,CAF,NFLAGINT,PHS,DEL,NERG,IEEV,NEL,PHSSEL)
247            CONTINUE


	       IF(NFLAGINT.eq.1) THEN
		 WRITE(1,*) 'BE CAREFUL! AT HIGH ENERGY YOU ARE DOING
     &           EXTRAPOLATION NOT INTERPOLATION OF THE PHASE SHIFTS'
	       ENDIF
	       DRPER=DRPER1(NEL)
	       DRPAR=DRPAR1(NEL)
	       DR0=DR01(NEL)
	       IT=0
C
C =============================================================================
C
C Compute components of incident wavevector parallel to surface and required
C beams at this energy.
C
C =============================================================================
C
	       NEXIT=0
C
	       CALL WAVE2(AK2,AK3,THETA,FI,E,VV,AK21,AK31,AK2M,AK3M,
     &          NT0,RAR1,RAR2,PQFEX,PSQ,NEXIT,SPQF,1,NBIN)

C
C =============================================================================
C
C Increment energy range.
C
C =============================================================================
C
1300        CONTINUE
	    WRITE (1,106)
	 ENDIF
      ENDIF


C Read in experimental IV curves and information relevant to the R-factor
C calculation.


       CALL RFIN(IBP,NT0,WB,WR,1,IPR)
       CALL EXPAN(INBED,IEERG,AE,EE,NEE,NBEA,BENAME,IPR,XPL,YPL,NNN,
     & AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,VPIS)
c  Calculate the amplitude of the diffracted beams in 
c  the kinematic limit.
c  Theoretical results stored in ivth# files
c  Experimental results stored in ivexp# files

	DO 99 I=1,INLTOT
	    DO 98 J=1,3
	  IF (I.LE.INLDISP) THEN
		ADISP(I,J)=DISP(I,J)
	  ELSE
		ADISP(I,J)=0.0
	  ENDIF
98	CONTINUE
99	CONTINUE
	
       rfactor = VINTENTF(INLTOT,DISP,PSQ,INTAU,NT0,PHSSEL,EI,EF,DE,
     & NL1,NL2,IELEMOL2,WPOSTF,TVA,SPOSTF1,PQFEX,ASA,INVECT,
     & INBED,IEERG,AE,EE,NEE,NBEA,BENAME,IPR,AP,APP,YE,
     & SE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,NERG,L1,ITEMP)

	write(*,*) 'kleed1: rfactor',rfactor
c	CLOSE(4)
c	CLOSE(5)

      RETURN 
      END
