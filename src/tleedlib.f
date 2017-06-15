C  file LEEDSATL.SB1  Feb. 29, 1996
C
C**************************************************************************
C  Symmetrized Automated Tensor LEED (SATLEED):  subroutines, part 1
C  Version 4.1 of Automated Tensor LEED
C
C =========================================================================
C Subroutine APROP takes the plane wave amplitude incident on either side
C of a composite layer and calculates the spherical wave amplitude incident
C on the origin of each subplane given the interlayer (1-X) matrix.
C
C Changed in order to correct inconsistencies detected by fortran forcheck
C program "Fortran lint". All changes are marked with
C "C++++" (M. Gierer, 7/97).
C Changes in subroutine RFIN when reading input for r factor calculation
C (marked by "C++++")
C
C Input Parameters;
C =================
C
C TSTORE        =   INTERLAYER (1-X) MATRIX (FROM MTINVT)
C AMPPLW        =   PLANE WAVE AMPLITUDES INCIDENT ON EITHER SIDE OF THE
C                   COMPOSITE LAYER (FROM RFS)
C NT            =   NUMBER OF BEAMS IN CALCULATION AT CURRENT ENERGY
C ALM           =   OUTPUT PLANE WAVES INCIDENT ON EITHER SIDE OF EACH
C                   SUBPLANE OF THE COMPOSITE LAYER
C E,VPI         =   CURRENT (COMPLEX) ENERGY
C NLAY          =   NUMBER OF SUBPLANES IN COMPOSITE LAYER
C TSF           =   PHASE SHIFTS
C NST1          =   Number of composite layers
C NCL           =   Current CL index
C
C Modified version of ROUS's routine FINALC. 
C Modifications by WANDER and BARBIERI
C
C =========================================================================
C
      SUBROUTINE APROP3(IDXN,WBDS,TSTORE,AMPPLW,NCL,NT,ALM,E,VPI,
     & NLAY,LMMAX,LMAX,LXM,NST1,LMNMAX,JLMNI,NINQ,LL1,NLMX)
C
      DIMENSION LXM(LMMAX),IDXN(LMNMAX),NINQ(NST1)
      DIMENSION LL1(NLMX,NST1)
      COMPLEX WBDS(LMNMAX)
      COMPLEX TSTORE(2,JLMNI,NT,NST1),AMPPLW(NT,2,NST1)
      COMPLEX ALM(NLAY,LMMAX),CSUM,AK,CI
C
      CI=CMPLX(0.0,1.0)
      AK=-0.5/CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      DO 100 I2=1,NINQ(NCL)
         I=LL1(I2,NCL)
         IN=(I-1)*LMMAX
         K=0
         DO 110 L=0,LMAX
            DO 120 M=-L,L
               K=K+1
               KLM=LXM(K)
               K2=IDXN(IN+KLM)
               CSUM=CMPLX(0.0,0.0)
               IF(K2.EQ.0) GOTO 119
               DO 130 JGP=1,NT
                CSUM=CSUM+TSTORE(1,K2,JGP,NCL)*AMPPLW(JGP,1,NCL)
     &          +TSTORE(2,K2,JGP,NCL)*AMPPLW(JGP,2,NCL)
130            CONTINUE
119            ALM(I2,K)=CSUM*WBDS(IN+KLM)*(CI**L)
120         CONTINUE
110      CONTINUE
100   CONTINUE
      RETURN
      END
C==================================================================== 
C                                                                     
C Subroutine CVEC evaluates the matrix G needed to construct the 
C tensor Q in subroutine QGEB.
C Two types of G matrix can be constructed by this routine for the
C incident beam (NEXIT=0), orthe time reversed beams (NEXIT>0).
C The G matrix is the product of the vector of spherical wave amplitudes
C 'TLM' and a Gaunt coefficient 'C'.
C
C  G(L1;M1,L3;M3)=
C
C  NEXIT=0 => C(L1,-M1; L3,L3; LP,MP) * TLM(LP,MP) * CI**(+LP)
C  NEXIT=1 => C(L1,+M1; L3,M3; LP,MP) * TLM(LP,MP) * CI**(-LP)
C 
C  N.B.   0 <=  L1 <= LMAX  &  0 <= L3 <= LSMAX
C
C       The Gaunt coefficients obey the usual selection rules.
C 
C Input Parameter;
C ================
C
C LLAY                 =     INDEX OF THE CURRENT LAYER.
C LMAX                 =     L VALUE OF LARGEST ANGULAR MOMENTUM COMPONENT
C                            USED TO DESCRIBE ATOMIC SCATTERING.
C TLM                  =     SPHERICAL WAVE AMPLITUDES INCIDENT UPON EACH
C                            REQUESTED LAYER FOR THE CURRENT EXIT BEAM.
C ALM                  =     WORKING SPACE
C LMMAX                =     (LMAX+1)**2.
C NEXIT                =     FLAG INDICATING WHICH G IS TO BE CALCULATED
C                            (SEE ABOVE)
C LPSS                 =     PHASE SHIFT ASSIGNMENT IN COMPOSITE LAYER
C TSF                  =     PHASE SHIFTS
C LSMMAX               =     (LSMAX+1)**2
C BELM                 =     VECTOR OF GAUNT COEFFICIENTS FROM SUBROUTINE GAUNT.
C G(NG1;LMMAX,LSMMAX)  =     G MATRIX DEFINED AS ABOVE.
C NGMAX                =     FIRST DIMENSION OF G (EITHER 1 OR NLAYER)
C NLMB                 =     DIMENSION OF BELM
C NG1                     =     STORAGE INDEX OF G (EITHER 1 OR LLAY)
C
C============================================================================
C
      SUBROUTINE CVEC3(LLAY,LMAX,TLM,NLAY,ALM,LMMAX,LSM,NEXIT,LPSS,
     & TSF,LSMMAX,BELM,G,NGMAX,NLMB,NG1,LLL1,NLMX,NST1,NCL)
C
      COMPLEX TLM(NLAY,LSM),ALM(LSM),TSF(6,16)
      COMPLEX G(NGMAX,LMMAX,LSMMAX),CI,PRE,CSUM
      DIMENSION BELM(NLMB),LPSS(NLAY),LLL1(NLMX,NST1)
C
C SET CONSTANTS:
C
      PI=4.0*ATAN(1.0)
      CI=CMPLX(0.0,1.0)
      LLAY2=LLL1(LLAY,NCL)
C
C TRANSFER SPHERICAL WAVE AMPLITUDES FOR CURRENT LAYER
C FROM TLM TO THE WORKING ARRAY ALM MULTIPLYING
C BY THE APPROPRIATE POWERS OF CI
C
      I=0
      PRE=4.0*PI/CI
      DO 110 L=0,LMAX
         PRE=PRE*CI
         DO 120 M=-L,L
            I=I+1
            ALM(I)=PRE*TLM(LLAY,I)
120      CONTINUE
110   CONTINUE
C
C START CONSTRUCTION OF G MATRIX RETRIEVING THE GAUNT COEFFICIENTS
C IN THE ORDER IN WHICH THEY WERE STORED BY GAUNT.
C
      K=0
      DO 300 I1=1,LMMAX
         L1=INT(SQRT(FLOAT(I1-1)))
         M1=I1-L1*L1-L1-1
C
C CHOOSE PREFACTOR DEPENDING ON THE VALUE OF NEXIT
C
         IF (NEXIT.GT.0) THEN
            I1A=L1*L1+L1+1-M1
            PRE=CMPLX(1.0,0.0)
         ELSE
            I1A=I1
            PRE=(-1)**(M1+L1)*CI*TSF(LPSS(LLAY2),L1+1)
         ENDIF
C
C START INNER LOOP OVER ANGULAR MOMENTUM COMPONENTS SELECTED FOR THE
C SINGLE CENTRE EXPANSION NOTE THAT ONLY VALUES OF (LP,MP) ARE
C SELECTED WHICH ARE CONSISTENT WITH THE EXISTENCE ON NON-ZERO GAUNT
C COEFFS.
C
         DO 200 I3=1,LSMMAX
            L3=INT(SQRT(FLOAT(I3-1)))
            M3=I3-L3*L3-L3-1
            MP=M1+M3
            LL1=MAX0(IABS(MP),IABS(L1-L3))
            LL2=L1+L3
C            LL2=MIN0(LMAX,L1+L3)
C            IF (MOD(LL2+L1+L3,2).NE.0) LL2=LL2-1
            CSUM=CMPLX(0.0,0.0)
            IP=LL2*LL2+LL2+1+MP
            DO 100 LP=LL2,LL1,-2
               K=K+1
               CSUM=CSUM+BELM(K)*ALM(IP)
               IP=IP-LP-LP-LP-LP+2
100         CONTINUE
            G(NG1,I1A,I3)=CSUM*PRE
200      CONTINUE
300   CONTINUE
      RETURN
      END
C============================================================================
C
C  Subroutine BEAMT2 selects those beams from the input list that are
C  needed at the current energy, based on the parameter TST which limits
C  the decay of plane waves from one layer to the next (the interlayer
C  spacing must already be incorporated in TST, which is done by subroutine
C  READT). Only symmetry inequivalent beams are considered in the
C  analysis
C
C Parameter List;
C ===============
C
C  KNBS          =   NO. OF INPUT BEAM SETS
C  KNB           =   NO. OF BEAMS IN EACH INPUT BEAM SET
C  SPQ,SPQF      =   LIST OF INPUT BEAMS (RECIPROCAL LATTICE VECTORS G)
C  KNT           =   TOTAL NO. OF INPUT BEAMS
C  AK2,AK3       =   COMPONENTS OF INCIDENT WAVEVECTOR PARALLEL TO SURFACE
C  E             =   CURRENT ENERGY
C  TST           =   CRITERION FOR BEAM SELECTION
C  NB,PQ,PQF,NT  =   KNB,SPQ,SPQF,KNT AFTER BEAM SELECTION
C  IIDXK2,IIDXK  =   IDXK2, IDXK after beam selection (see BEMGEN)
C  IIKRED        =   IKRED after beam selection (see BEMGEN)
C  NP            =   LARGEST NO. OF BEAMS SELECTED IN ANY ONE BEAM SET
C  IPR           =   PRINT CONTROL PARAMETER
C
C Modified version of routine BEAMS from the VAN HOVE/TONG LEED package.
C Modifications by BARBIERI
C
C============================================================================
C
      SUBROUTINE BEAMT3(KNBS,KNB,SPQ,SPQF,KNT,AK2,AK3,E,TST,NB,
     & PQ,PQF,NT,NP,IPR,NPRT,IDXK,
     & IIDXK,IIDXK2,IKRED,IIKRED,NROM,NSHBS,IMIN,IMAX,NKR,NK3,INDK3,
     & PH3,PQ3,INI)
C
      DIMENSION KNB(KNBS),SPQ(2,KNT),SPQF(2,KNT),NB(KNBS)
      DIMENSION PQ(2,KNT),PQF(2,KNT)
      DIMENSION NSHBS(KNBS),IMIN(KNBS),IMAX(KNBS)
      DIMENSION IDXK(NROM),IIDXK(NROM),IIDXK2(NKR)
      DIMENSION IKRED(NKR),IIKRED(NKR)
      DIMENSION INDK3(KNT),PH3(36,KNT),PQ3(2,KNT)
C
705   FORMAT (' ',1I3,' SYMMETRIC BEAMS USED  ',8I4)
706   FORMAT (' ',10(9(2X,2F6.3),/))
C
C generate PH3, PQ3, INDK3 ect in the off normal incidence case.   
C
      IF(INI.EQ.0.AND.NPRT.EQ.0) THEN
        PI2=8.*ATAN(1.)
        DO 3100 I=1,KNT
           PQ3(1,I)=SPQ(1,I)+AK2
           PQ3(2,I)=SPQ(2,I)+AK3
           INDK3(I)=IKRED(I)
           PH3(1,I)=0.
3100    CONTINUE
        NK3=KNT
C
C the mirror is along the Y=0 line (see bemgen3), if there is a mirror 
C The only beams with possibly the same magnitudes are the symmetry
C related ones
C
        DO 3200 I=1,KNT
              ZZ1=PQ3(1,I)
              ZZ2=PQ3(2,I)
              ZMOD=SQRT(ZZ1**2+ZZ2**2)
              IF(INDK3(I).EQ.2) THEN
C
C beam is not along Y=0 line
C
C first measure angle with respect to Y=0 line. The angle sought
C is double the one wrt Y=0 
C
                DOT=PQ3(2,I)
                CROS=PQ3(1,I)
C CSV can be < -1 (numerical errors)
                CSV=DOT/(ZMOD)
                IF(CSV.LT.-1.) CSV=-1.
                IF(CROS.GE.0.) THEN
                   PH3(2,I)=PI2-2.*ACOS(CSV)
                ELSE
                   PH3(2,I)= 2.*ACOS(CSV)
                ENDIF
              ENDIF
3200       CONTINUE
      ENDIF
      KNBJ=0
      NT=0
      NP=0
      NKK=0
C initialize IIDXK2 in the first call to BEAMT2
      IF (NPRT.EQ.1) THEN
        DO 800 J=1,NKR
           IIDXK2(J)=0
800     CONTINUE
      ENDIF
      DO 704 J=1,KNBS
         IMAX2=0
         IMIN(J)=10000
         IMAX(J)=0
         NSHBS(J)=0
         JPREV=J
         N=KNB(J)
         NB(J)=0
C
C identify the beamset < J (if it exists) equivalent to J 
C
         IREPJ=IDXK(KNBJ+1)
         IMAX3=IMAX(1)
         IMAX4=0
C
C identify previous inequivalent beamsets to monitor shift
C NSHBS for J2=1,J-1 is the number of inequivalent beams that are
C skipped at this energy for each beamset. NSHBS(J) is the total number 
C of inequivalent beams that are skipped before reaching beamset J
C
         DO 600 J2=1,J-1
            IMAX3=MAX0(IMAX(J2),IMAX3)
            IF(IREPJ.GT.IMAX3.AND.IMAX(J2).GT.IMAX4) THEN
               IMAX4=MAX0(IMAX(J2),IMAX4)
               NSHBS(J)=NSHBS(J)+NSHBS(J2)
            ENDIF
600      CONTINUE
         DO 703 K=1,N
            KK=K+KNBJ
            FACT1=(AK2+SPQ(1,KK))*(AK2+SPQ(1,KK))
            FACT2=(AK3+SPQ(2,KK))*(AK3+SPQ(2,KK))
C minimum inequivalent label in beamset J
            IMIN(J)=MIN0(IMIN(J),IDXK(KK))
            IF ((2.0*E-FACT1-FACT2)+TST.GE.0) THEN
               NB(J)=NB(J)+1
               NT=NT+1
C
C record maximum label of acceptable beam
C
               IMAX(J)=MAX0(IMAX(J),IDXK(KK))
               DO 702 I=1,2
                  PQ(I,NT)=SPQ(I,KK)
                  PQF(I,NT)=SPQF(I,KK)
702            CONTINUE
C
C fix labeling IDXK --> IIDXK
C
               IF (NPRT.EQ.1) THEN
                  IIDXK(NT)=IDXK(KK)-NSHBS(J)
                  IF(IIDXK2(IIDXK(NT)).EQ.0) IIDXK2(IIDXK(NT))=NT
               ENDIF
               IF (NPRT.EQ.0) IIKRED(NT)=IKRED(KK)
            ENDIF
            IMAX2=MAX0(IMAX2,IDXK(KK))
703      CONTINUE
         NSHBS(J)=IMAX2-IMAX(J)
         KNBJ=KNBJ+KNB(J)
         NP=MAX0(NP,NB(J))
704   CONTINUE
      IF (IPR.GT.0.AND.NPRT.EQ.0) THEN
         WRITE (1,705) NT,(NB(J),J=1,KNBS)
         WRITE (1,706) ((PQF(I,K),I=1,2),K=1,NT)
      ENDIF
      RETURN
      END
C======================================================================
C
C  Subroutine BEMGEN2 generates the list of beams for the TENSOR LEED
C  program. The beams are ordered by increasing ABS(G). 
C
C  Notice that the list of beams produced by BEMGEN does not respect, in
C  general, the symmetry of the surface because of the way the list is
C  generated. This is true for rather extreme unit cell cases (c2x10)
C  e.g., but can be fixed by changing the line marked by CCC . (Barbieri)
C 
C  The list of beams is also analyzed according to the symmetry of the
C  surface and NREDK ,the number of symmetry inequivalent beams, is
C  generated together with the map IDXK(I)=J giving for each beam  I the
C  correspondent reduced index J. Hence IDXK(I)=IDXK(J) means that
C  I and J are equivalent beams. Also IDXK2(I) I=1,NKRED is generated
C  giving the label of one of the inequivalent beams (in the large list) 
C  corresponding to the reduced label J
C  The # of elements in each equivalence class K (=1,NKRED) is given by
C  IKRED(K).
C
C  NK3, the number of beams with different length, is also generated
C  together with a list of such beams PQ(2,I),  INDK3(I) I=1,NK3 equal 
C  to the number of beams with same lenght as PQ(2,I), and PH3(J,I)
C  the angle between the Jth beam in class I and PQ(2,I). 
C
C  The list of external beams is modified to an equivalent list such
C  that the variable NEXIT in WAVE2 can be defined by using inequivalent
C  beams only. A new list PQFEX2 equivalent to PQFEX is also produced.
C  PQFEX(I) and PQFEX2(I) are equivalent beams, equivalent to the
C  original PQFEX(I), but now they are such that -PQFEX(I) is in the
C  symmetry reduced beam list, and PQFEX2 is in the same list
C XXXX possible problems for off-normal incidence
C Parameter List;
C ===============
C
C  TST         =   CUTOFF CRITERION FOR SELECTION OF PLANE WAVES
C  EF          =   HIGHEST ENERGY USED IN THIS CALCULATION
C  SPQ         =   LIST OF BEAMS
C  SPQF        =   SAME AS SPQ BUT IN UNITS OF RECIPROCAL LATTICE VECTORS
C  KNBS        =   NUMBER OF BEAM SETS INCLUDED IN CALCULATION
C  KNB         =   NUMBER OF BEAMS IN EACH BEAMSET
C  RAR1,RAR2   =   SUBSTRATE RECIPROCAL LATTICE VECTORS
C  KNT         =   TOTAL NUMBER OF BEAMS
C  IPR         =   PRINT CONTROL PARAMETER
C  TVA         =   AREA OF SUBSTRATE UNIT CELL
C  DFLAG       =   IF(DFLAF.EQ.1) ONLY CALCULATE IO BEAMLIST
C  NL          =   AREA OF OVERLAYER UNIT CELL IN TERMS OF SUBSTRATE UNIT
C                  CELL
C
C In Common Blocks;
C =================
C
C  ARA1,ARA2   =   SUBSTRATE LATTICE VECTORS
C  ARB1,ARB2   =   OVERLAYER LATTICE VECTORS
C  RBR1,RBR2   =   OVERLAYER RECIPROCAL LATTICE VECTORS
C  NL1,NL2     =   SUPERLATTICE CHARACTERIZATION CODES
C
C AUTHOR: VAN-HOVE, Major mods: BARBIERI
C
C======================================================================
C
      SUBROUTINE BEMGEN3(TST,EF,SPQF,SPQ,KNBS,KNB,RAR1,RAR2,KNT,IPR,TVA,
     & DFLAG,KNBMAX,G,IDXK,NSYM,NKRED,SPQFS,SPQS,IKRED,PQFEX,PQFEX2,NT0,
     & IDXK2,IWK,SPWK,KNB2M,KSNBS,KSNB,NK3,INDK3,PH3,PQ3,INI,
     & NROM2R)
C
      INTEGER LATMAT(2,2),KNB(60),DFLAG
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2),SPQF(2,KNBMAX)
      DIMENSION SPQ(2,KNBMAX),RAR1(2),RAR2(2),RBR1(2),RBR2(2)
      DIMENSION ALMR(2,2),G(2,KNBMAX),ST(4),DG(2),IDXK(KNBMAX)
      DIMENSION SPQS(2,KNB2M),SPQFS(2,KNB2M)
      DIMENSION IKRED(KNB2M),PQFEX(2,NT0),KSNB(20),PQFEX2(2,NT0)
      DIMENSION IDXK2(KNBMAX),IWK(KNB2M),SPWK(2,KNB2M)
      DIMENSION INDK3(KNB2M),PH3(36,KNB2M),PQ3(2,KNBMAX)
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
1001  FORMAT (' BEAMGENERATION CODE')
1002  FORMAT (' ===================')
165   FORMAT (' LATMAT HAS ZERO DETERMINANT')
266   FORMAT (' ERROR in parameters. NROM2R should equal ',
     & I7)
268   FORMAT (' The mirror plane in the off normal incidence case
     &  must be put along Y=0 (NSYMS=6). ')
269   FORMAT (' maximum number of beams with same length is > 36; the
     &  first dimension of PH3 wherever it appears must be equal
     &  to at least  ',I4)
270   FORMAT ('NROM too small. It should equal',I4)
271   FORMAT ('NROMR too small. It should equal',I4)
C
      IF (IPR.GT.0) THEN
         WRITE (1,1001)
         WRITE (1,1002)
      ENDIF
C
C GENERATE LATICE MATRIX
C
      IF (DFLAG.EQ.0) THEN
         CALL MATGN2(ARA1,ARA2,ARB1,ARB2,LATMAT)
      ELSE
         CALL MATGN2(ARA1,ARA2,ARA1,ARA2,LATMAT)
      ENDIF
C
C SET UP CONSTANTS 
C
      PI=3.1415926535
      PI2=2.*PI
C
C  ALMR IS MATRIX RELATING SUBSTRATE TO OVERLAYER LATTICES
C
      DET=LATMAT(1,1)*LATMAT(2,2)-LATMAT(1,2)*LATMAT(2,1)
      IF (ABS(DET).GE.1.E-5) THEN
         ALMR(1,1)=FLOAT(LATMAT(2,2))/DET
         ALMR(2,2)=FLOAT(LATMAT(1,1))/DET
         ALMR(1,2)=-FLOAT(LATMAT(2,1))/DET
         ALMR(2,1)=-FLOAT(LATMAT(1,2))/DET
         GMAX2=2.*EF/27.21+TST
         AKP2=GMAX2*TVA/(4.*3.1415926)
         IF(AKP2.LT.80.) GMAX2=80./TVA*(4.*3.1415926)
C
C  GENERATE ALL BEAMS WITHIN A BEAM CIRCLE OF RADIUS SQRT(GMAX2),
C  LIMITED TO KNBMAX BEAMS (but at least 80 beams)
C
         KNT=0
         NI1=-1
180      NI1=NI1+1
         NOP=0
         NI2=-1
190      NI2=NI2+1
         DO 280 K=1,4
            IF (K.EQ.2) THEN
               IF (NI1.EQ.0.AND.NI2.EQ.0) GOTO 290
               IF (NI1.EQ.0) GOTO 280
               II1=-NI1
               II2=NI2
            ELSEIF (K.EQ.3) THEN
               IF (NI1.EQ.0.OR.NI2.EQ.0) GOTO 280
               II1=-NI1
               II2=-NI2
            ELSEIF (K.NE.4) THEN
               II1=NI1
               II2=NI2
            ELSEIF (NI2.EQ.0) THEN
               GOTO 280
            ELSE
               II1=NI1
               II2=-NI2
            ENDIF
            IF (DFLAG.EQ.0) THEN
               GT1=FLOAT(II1)*RBR1(1)+FLOAT(II2)*RBR2(1)
               GT2=FLOAT(II1)*RBR1(2)+FLOAT(II2)*RBR2(2)
            ELSE
               GT1=FLOAT(II1)*RAR1(1)+FLOAT(II2)*RAR2(1)
               GT2=FLOAT(II1)*RAR1(2)+FLOAT(II2)*RAR2(2)
            ENDIF
            IF ((GT1*GT1+GT2*GT2).LE.GMAX2) THEN
               KNT=KNT+1
               NOP=1
               IF (NI1.EQ.0) NIT=NI2+4
CCC               IF (NI1.EQ.0) NIT=NI2+NNEW NNEW maybe 8, 10 
C               IF (KNT.GT.KNBMAX) GOTO 280
               G(1,KNT)=GT1
               G(2,KNT)=GT2
               SPQF(1,KNT)=ALMR(1,1)*II1+ALMR(2,1)*II2
               SPQF(2,KNT)=ALMR(1,2)*II1+ALMR(2,2)*II2
            ENDIF
280      CONTINUE
290      IF (NI2.LE.NIT) GOTO 190
         IF (NOP.EQ.1) GOTO 180
C
C  ORDER BEAMS BY INCREASING ABS(G)
C
         KNT1=KNT-1
         IF (KNT.GT.KNBMAX) GOTO 260
         DO 320 I=1,KNT1
            AM=G(1,I)*G(1,I)+G(2,I)*G(2,I)
            I1=I+1
            DO 310 J=I1,KNT
               AM1=G(1,J)*G(1,J)+G(2,J)*G(2,J)
               IF (AM1.LT.AM) THEN
                  ST(1)=G(1,J)
                  ST(2)=G(2,J)
                  ST(3)=SPQF(1,J)
                  ST(4)=SPQF(2,J)
                  DO 300 KK=I1,J
                     K=J+I1-KK
                     SPQF(1,K)=SPQF(1,K-1)
                     SPQF(2,K)=SPQF(2,K-1)
                     G(1,K)=G(1,K-1)
                     G(2,K)=G(2,K-1)
300               CONTINUE
                  G(1,I)=ST(1)
                  G(2,I)=ST(2)
                  SPQF(1,I)=ST(3)
                  SPQF(2,I)=ST(4)
                  AM=AM1
               ENDIF
310         CONTINUE
320      CONTINUE
C
C extract information about beams with same length if we are at normal
C incidence. Otherwise store SPQ in PQ3 to postpone the analysis in
C BEAMT2 
C
         DO 10 J=1,KNT
            DO 20 I=1,2
               SPQ(I,J)=SPQF(1,J)*RAR1(I)+SPQF(2,J)*RAR2(I)
20          CONTINUE
10       CONTINUE
         IF(INI.EQ.1) THEN
C
C normal incidence. Generate information about beams with same magnitude
C In this case an extended symmetry can be used
           ZMOD=-1.
           NK3=0
           KMAX=1
           DO 3200 I=1,KNT
              ZMOD1=SQRT(SPQ(1,I)*SPQ(1,I)+SPQ(2,I)*SPQ(2,I))
              IF(ABS(ZMOD1-ZMOD).GT.0.0001) THEN
C
C first beam of length zmod1. 
C
                NK3=NK3+1
                INDK3(NK3)=1
C
C remember shift for off normal incidence
C
                PQ3(1,NK3)=SPQ(1,I)
                PQ3(2,NK3)=SPQ(2,I)
                ZMOD=ZMOD1
                PH3(1,NK3)=0.
C
C KMAX is the maximum number of beams with same magnitude
C
                IF(NK3.GT.1) KMAX=MAX0(KMAX,INDK3(NK3-1))
              ELSE
C
C beam I has same length as PQ3(NK3). Extract info.
C
                DOT=PQ3(1,NK3)*SPQ(1,I)+PQ3(2,NK3)*SPQ(2,I)
                CROS=PQ3(1,NK3)*SPQ(2,I)-PQ3(2,NK3)*SPQ(1,I)
                INDK3(NK3)=INDK3(NK3)+1
                IF(INDK3(NK3).GT.36) GOTO 3200
C CSV can be < -1 (numerical errors)
                CSV=DOT/(ZMOD*ZMOD)
                IF(CSV.LT.-1.) CSV=-1.
                IF(CROS.GE.0.) THEN
                   PH3(INDK3(NK3),NK3)=ACOS(CSV)
                ELSE
                   PH3(INDK3(NK3),NK3)= PI2-ACOS(CSV)
                ENDIF
              ENDIF
3200       CONTINUE
           KMAX=MAX0(KMAX,INDK3(NK3))
           IF(KMAX.GT.36) GOTO 250
         ELSE
C
C the number of beams with same length is assumed to be no more than 36.
C In the generic case only symmetry equivalent beams can have the same
C magnitude. The angles between symmetry equivalent beams however
C depend on the energy. Redo calculation in BEAMT2
C
C check that the mirror is along Y=0
C
           IF(NSYM.NE.1.AND.NSYM.NE.6) GOTO 240
           KMAX=36
         ENDIF
         NKRED=1
         DO 30 J=1,KNT
            IDXK(J)=1
30       CONTINUE
         IDXK2(1)=1
         IDXK(1)=0
C
C generate Symmetric list. IDXK is used here as working space
C Notice that even off normal incidence the labeling of equivalent
C beams is independent on the energy 
C
         DO 40 J=2,KNT
            IF(IDXK(J).EQ.1) THEN
              NKRED=NKRED+1
              IDXK(J)=0
              IDXK2(J)=NKRED
C              IF(NKRED.GT.KNB2M) NKRED=KNB2M
              DO 80 I=1,2
                 SPQS(I,NKRED)=SPQ(I,J)
                 SPQFS(I,NKRED)=SPQF(I,J)
80            CONTINUE
              IF(J+KMAX.GT.KNT) THEN
                 KMAX2=KNT-J
              ELSE
                 KMAX2=KMAX
              ENDIF
C
C generate symmetry equivalent beams, set IDXK(K)=0 for all k's equivalent
C to J 
C
              CALL LOOKK(SPQ,SPQS,KNT,NSYM,
     &             IDXK2,NKRED,KMAX2,J,IDXK)
            ENDIF
40       CONTINUE
         IF (NKRED.GT.KNB2M) GOTO 261
C
C count # beams in each equivalence class to define normalization factor
C IKRED
C
         J=0
         DO 50 K=1,KNT
            IK=0
            IF(IDXK2(K).GT.J) THEN
               IF(K+KMAX.GT.KNT) THEN
                  KMAX2=KNT-K
               ELSE
                  KMAX2=KMAX
               ENDIF
               DO 60 J=K,K+KMAX2
                 IF(IDXK2(K).EQ.IDXK2(J)) IK=IK+1
60             CONTINUE
               J=IDXK2(K)
               IKRED(IDXK2(K))=IK
            ENDIF
50       CONTINUE
C
C  ORDER BEAMS BY BEAM SET
C
         TWPI=2.*3.1415926535
         I=1
         KNBS=1
330      CONTINUE
         KNB(KNBS)=1
         IDXK(I)=IDXK2(I)
         JL=I+1
         IF (JL.LE.KNT) THEN
            DO 360 J=JL,KNT
               DG(1)=G(1,J)-G(1,I)
               DG(2)=G(2,J)-G(2,I)
C
C  TEST WHETHER BEAMS I AND J BELONG TO SAME SUBSET. IF THEY DO,
C  BRING THEM TOGETHER IN THE LIST
C
               B=ABS(DG(1)*ARA1(1)+DG(2)*ARA1(2))+0.001
               B=AMOD(B,TWPI)/TWPI
               IF (ABS(B).LT.0.002) THEN
                  B=ABS(DG(1)*ARA2(1)+DG(2)*ARA2(2))+0.001
                  B=AMOD(B,TWPI)/TWPI
                  IF (ABS(B).LT.0.002) THEN
C
C  IF J=I+1 NO REORDERING NEEDED
C
                     IF (J.NE.I+1) THEN
                        ST(1)=G(1,J)
                        ST(2)=G(2,J)
                        ST(3)=SPQF(1,J)
                        ST(4)=SPQF(2,J)
                        III=IDXK2(J)
                        I2=I+2
                        DO 340 KK=I2,J
                           K=J+I2-KK
                           SPQF(1,K)=SPQF(1,K-1)
                           SPQF(2,K)=SPQF(2,K-1)
                           IDXK2(K)=IDXK2(K-1)
                           G(1,K)=G(1,K-1)
                           G(2,K)=G(2,K-1)
340                     CONTINUE
                        G(1,I+1)=ST(1)
                        G(2,I+1)=ST(2)
                        SPQF(1,I+1)=ST(3)
                        SPQF(2,I+1)=ST(4)
                        IDXK(I+1)=III
                     ELSE
                        IDXK(I+1)=IDXK2(J)
                     ENDIF
                     I=I+1
                     KNB(KNBS)=KNB(KNBS)+1
                     IF (I.EQ.KNT) GOTO 1003
                  ENDIF
               ENDIF
360         CONTINUE
            I=I+1
            IF (I.NE.KNT) THEN
               KNBS=KNBS+1
               GOTO 330
            ENDIF
         ENDIF
1003     CONTINUE
C
C IDXK(I)=K means that the Ith beam in the list SPQF is equivalent to
C the Kth beam in the reduced list.  
C
         DO 101 J=1,KNT
            DO 201 I=1,2
               SPQ(I,J)=SPQF(1,J)*RAR1(I)+SPQF(2,J)*RAR2(I)
201          CONTINUE
101      CONTINUE
C
C order list of symmetric beams according to beamsets
C First consider only inequivalent beams in IDXK2
C
         IDXK2(1)=1
         JK=1
         DO 123 K=2,KNT
            KK=IDXK(K)
            DO 124 J=1,JK
               IF (KK.EQ.IDXK2(J)) GOTO 123
124         CONTINUE     
            JK=JK+1
            IDXK2(JK)=KK
123      CONTINUE
C
C Now order IDXK2 
C
         DO 122 K=1,NKRED
            DO 126 J=1,NKRED
               IF(IDXK2(J).EQ.K) THEN
                 IWK(K)=J
                 GOTO 122
               ENDIF
126         CONTINUE
122      CONTINUE
         DO 128 K=1,NKRED
            SPWK(1,K)=SPQFS(1,IDXK2(K))
            SPWK(2,K)=SPQFS(2,IDXK2(K))
128      CONTINUE
         DO 177 K=1,NKRED
            IDXK2(K)=IKRED(IDXK2(K))
177      CONTINUE
C
C reorder IDXK, IKRED and SPQFS
C
         DO 125 K=1,KNT
            IDXK(K)=IWK(IDXK(K))
125      CONTINUE
         DO 129 K=1,NKRED
            SPQFS(1,K)=SPWK(1,K)
            SPQFS(2,K)=SPWK(2,K)
            IKRED(K)=IDXK2(K)
129      CONTINUE
         DO 131 J=1,NKRED
            DO 132 I=1,2
               SPQS(I,J)=SPQFS(1,J)*RAR1(I)+SPQFS(2,J)*RAR2(I)
132          CONTINUE
131      CONTINUE
           
C
C regroup beamsets taking symmetry into account. Generate KSNBS, # of
C symmetric beamsets, KSNB(I) I=1,KSNBS number of inequivalent beams
C in each symmetric beamset.
C
         KSNBS=0
         J=0
         JL=0
         DO 140 K=1,KNBS
            I=JL
            DO 145 K2=1,KNB(K)
               J=J+1
               JL=MAX0(IDXK(J),JL)
145         CONTINUE
            IF(JL.GT.I) THEN
               KSNBS=KSNBS+1
               KSNB(KSNBS)=JL-I
            ENDIF
140      CONTINUE
C
C generate IDXK2
C
         DO 90 K=1,NKRED
            DO 95 KK=1,KNT
C
C select the first beam with inequivalent label K 
C
               IF(IDXK(KK).EQ.K) THEN
C           BCH=ABS(SPQFS(1,K)-SPQF(1,KK))+ABS(SPQFS(2,K)-SPQF(2,KK))
C                 IF(BCH.LE.001) THEN
                  IDXK2(K)=KK
                  GOTO 90
C                 ENDIF
               ENDIF             
95       CONTINUE
90       CONTINUE
C
C reorder SPQFS and SPQS according to IDXK2
C
         DO 150 K=1,NKRED
            SPQFS(1,K)=SPQF(1,IDXK2(K))
            SPQFS(2,K)=SPQF(2,IDXK2(K))
150      CONTINUE
         DO 151 J=1,NKRED
            DO 152 I=1,2
               SPQS(I,J)=SPQFS(1,J)*RAR1(I)+SPQFS(2,J)*RAR2(I)
152          CONTINUE
151      CONTINUE

C
C modify external beam list
C
         DO 70 K=1,NT0
C
C is the time reversed beam one of the inequivalent beams selected?
C
            ZV1=PQFEX(1,K)
            ZV2=PQFEX(2,K)
            NKK=0
            DO 75 KK=1,KNT
               ZK=ABS(ZV1+SPQF(1,KK))+ABS(ZV2+SPQF(2,KK))
               ZK2=ABS(ZV1-SPQF(1,KK))+ABS(ZV2-SPQF(2,KK))
               IF(ZK.LT.0.001) THEN
C
C yes, then replace it by equivalent beam, while keep the old list
C in PQFEX2
C
                  PQFEX(1,K)=-SPQFS(1,IDXK(KK))
                  PQFEX(2,K)=-SPQFS(2,IDXK(KK))
                  IF (KK.EQ.1) THEN
C
C The vector is trivial
C
                    PQFEX2(1,K)=PQFEX(1,K)
                    PQFEX2(2,K)=PQFEX(2,K)
                  ENDIF
                  NKK=NKK+1
                  IF(NKK.EQ.2) GOTO 70
                ELSEIF(ZK2.LT.0.001) THEN
                  PQFEX2(1,K)=SPQFS(1,IDXK(KK))
                  PQFEX2(2,K)=SPQFS(2,IDXK(KK))
                  IF(NKK.EQ.2) GOTO 70
               ENDIF
75       CONTINUE
70       CONTINUE
C
C check some parameter statements in main
C
         NCH1=0
         DO 102 K=1,KSNBS
           NCH1=MAX0(NCH1,KSNB(K))
102      CONTINUE
         IF(NCH1.GT.NROM2R) GOTO 238
         RETURN
238      WRITE (1,266) NCH1 
         GOTO 1000
240      WRITE (1,268) 
         GOTO 1000
250      WRITE (1,269) KMAX
         GOTO 1000
260      WRITE (1,270) KNT
         GOTO 1000
261      WRITE (1,271) NKRED
      ELSE
         WRITE (1,165)
      ENDIF
1000  STOP
      END
C=========================================================================
C
C  Subroutine CAAA2 computes Clebsch-Gordan coefficients for use by
C  routine GHD in the same order as they are used.  
C
C Parameter List;
C ===============
C
C  CAA      =  CLEBSCH-GORDAN COEFFICIENTS
C  NCAA     =  MAXIMUM NUMBER OF CLEBSCH-GORDAN COEFFICIENTS NEEDED IN
C              THIS CALCULATION
C  LMMAX    =  (LMAX+1)**2
C
C  This is a modified version of the routine CAAA from the Van Hove/Tong
C  LEED package. Modifications by WANDER.
C
C============================================================================
C
      SUBROUTINE CAAA2(CAA,NCAA,LMMAX,KOUNT)
C
      DIMENSION CAA(NCAA),KOUNT(LMMAX,LMMAX)
C
      II=1
      DO 60 I=1,LMMAX
         L1=INT(SQRT(FLOAT(I-1)+0.00001))
         M1=I-L1-L1*L1-1
         DO 50 J=1,LMMAX
            L2=INT(SQRT(FLOAT(J-1)+0.00001))
            M2=J-L2-L2*L2-1
            M3=M2-M1
            IL=IABS(L1-L2)
            IM=IABS(M3)
            LMIN=MAX0(IL,IM+MOD(IL+IM,2))
            LMAX=L1+L2
            LMIN=LMIN+1
            LMAX=LMAX+1
            KOUNT(I,J)=(LMAX-LMIN)/2 +1
            DO 40 ILA=LMIN,LMAX,2
               LA=ILA-1
               CAA(II)=CA(L1,M1,L2,M2,LA,M3)
               II=II+1
40          CONTINUE
50       CONTINUE
60    CONTINUE
      II=II-1
      RETURN
      END
C ========================================================================
C
C Subroutine DIMSCH generates the dimensionality of the search
C using the information contained in the input NSYMS
C
C Input Parameters;
C =================
C
C NLAY       =   NUMBER OF LAYERS INCALCULATION
C CPVPOS     =   sorted cartesian coordinates of the composite layers
C ARB1,ARB2  =   OVERLAYER LATTICE VECTORS
C FPOS       =   CARTESIAN COORDINATES OF ATOMS IN TOP COMP. LAYER
C ILOOK      =   LOOKUP TABLE OF RELATED ATOMS CORRESPONDING TO NSYM
C NSYM       =   SYMMETRY CODE SPECIFYING THE TYPE OF DISPLACEMENTS
C                ALLOWED BY THE SEARCH
C NDOM       =   NUMBER OF SYMMETRY OPERATIONS COORESPONDING TO NSYM
C NDIM       =   DIMENSIONALITY OF THE SEARCH
C LLFLAG     =   ARRAYS SPECIFYING WHETHER THE ATOM HAS TO BE INCLUDED
C                (LLFLAG(i)=1) IN THE SEARCH OR NOT (LLFLAG(i)=0)
C
C OUTPUT
C =============
C
C NNDIM      =  DIMENSIONALITY OF THE STRUCTURAL PARAMETER SEARCH
C LSFLAG     =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C               NSYM) IN THE COMPOSITE LAYER. LSFLAG(i)=LSFLAG(j)
C               INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C               AS EQUIVALENT IN THE SEARCH
C NDIML      =  ARRAY GIVING THE EFFECTIVE DIMENSIONALITY OF ATOM
C               j (ACTUALLY TO BE USED IN CONJUNCTION WITH LLFLAG) 
C DIREC        =  SET OF DIRECTIONS
C
C
C AUTHOR: BARBIERI
C
C =========================================================================
C
C
      SUBROUTINE DIMSCH2(ILOOK,NLAY,NDOM,NNDIM,LSFLAG,
     % LLFLAG,NDIML,DIREC,ADISP,ACOORD,NDIM,NPASS,LLFLAG2,
     & LLFLAGC,NL,NEQ,INLIN)
C
      DIMENSION ILOOK(12,NLAY),LLFLAG(NLAY),NORD(12)
      DIMENSION LSFLAG(NLAY),LLFLAG2(NLAY),NDIML(NLAY)
      DIMENSION DIREC(NLAY,2,2),NL(NLAY,INLIN),LLFLAGC(NLAY)
      DIMENSION ADISP(NLAY,3),ACOORD(12,NLAY,3),NEQ(INLIN)
C
!177   FORMAT(2F7.4)
C
C COPY LLFLAG
C
      DO 50 I=1,NLAY
             LLFLAG2(I)=LLFLAG(I)
             LLFLAGC(I)=LLFLAG(I)
50    CONTINUE
      IF(NPASS.EQ.1) THEN
         DO 55 I=1,NLAY
             LLFLAG2(I)=1.
             LLFLAG(I)=1.
55       CONTINUE
      ENDIF
C
C First subdivide atoms in subgroups and determine LSFLAG
C NEF is the number of effective layers
C
      NEF=0
      DO 70 I=1,NLAY
C
C check whether the atom is to be included in the search
C
            IF (LLFLAG(I).EQ.0) THEN
                LSFLAG(I)=0
C
C yes
C
            ELSEIF(LLFLAG2(I).NE.0) THEN 
                NEF=NEF+1
C
C copy the string of atoms equivalent to atom I
C
                DO 75 K=1,NDOM
                   NORD(K)=ILOOK(K,I)
                   DO 76 KK=1,NLAY
                      IF(KK.EQ.NORD(K)) THEN
                         LSFLAG(KK)=NEF
                         LLFLAG2(KK)=0
                      ENDIF
76                 CONTINUE
75              CONTINUE
            ENDIF
70    CONTINUE
C
C++++
      IF (NEF.GT.INLIN) GOTO 1999
C
C  If no strucural parameter is included in the search return with NNDIM=0
C
      IF(NEF.EQ.0) THEN
        NNDIM=0
        GOTO 1000
      ENDIF
C
C Now that we have the information about which atoms have to be symmetry
C related (in LSFLAG) we determine the effective dimensionality associated
C to each atom. Start by defining a set of generic test displacements
C ADISP(defined on inequivalent atoms only), and use DISDOM to 
C generate all the other set of displacements
C related by symmetry
C
      DISY=.6
      DISZ=.8
      DMOD=DISY**2+DISZ**2
      KK=0
      DO 90 I=1,NLAY
         ADISP(I,1)=.0
         ADISP(I,2)=.0
         ADISP(I,3)=.0
         IF(LSFLAG(I).GT.KK) THEN
            KK=KK+1
            ADISP(I,1)=.0
            ADISP(I,2)=DISY
            ADISP(I,3)=DISZ
         ENDIF
90    CONTINUE
      CALL DISDOM(ADISP,ACOORD,ILOOK,NLAY,1)
C
C
C loop over atoms with different NEF. 
C
      DO 110 K=1,NEF
C
C identify equivalent atoms.  
C
         KK=0
         DO 100 I=1,NLAY
            IF(LSFLAG(I).EQ.K) THEN
              KK=KK+1
              NL(KK,K)=I
            ENDIF
100      CONTINUE
         NEQ(K)=KK
C      
C KK is the number of equivalent atoms(stored in NEQ), NL points 
C to their layer index. 
C Starting with a generic in plane displacement to be assigned
C to atom NL(1,K) we compute all the symmetry related displacements to be
C assigned to the same atom  due to the fact that the atom is
C left invariant by some symmetry operation. If no symmetry operation
C leaves the atom invariant NDIML=3. Otherwhise an exam of these
C displacements shows whether NDIML=1 or =2.
C
C
C Check the different symmetry operations
C
         L=NL(1,K)
         DISFY=.0
         DISFZ=.0
         NCOUNT=0
         DO 105 I=1,NDOM
            IF(ILOOK(I,L).EQ.L) THEN
               DISFY=DISFY + ACOORD(I,L,2)
               DISFZ=DISFZ + ACOORD(I,L,3)
               NCOUNT=NCOUNT+1
            ENDIF
105         CONTINUE
            DISFY=DISFY/FLOAT(NCOUNT)
            DISFZ=DISFZ/FLOAT(NCOUNT)
            DFMOD=DISFY**2+DISFZ**2
            IF (DFMOD.LT.1.E-4) THEN
                NDIML(K)=1
            ELSEIF (ABS(DFMOD-DMOD).LE.1.E-4) THEN
                NDIML(K)=3
                IF (NDIM.EQ.1) NDIML(K)=1
                DIREC(L,1,1)=1.
                DIREC(L,1,2)=0.
                DIREC(L,2,1)=0.
                DIREC(L,2,2)=1.
            ELSE
                NDIML(K)=2
                IF (NDIM.EQ.1) NDIML(K)=1
                DIREC(L,1,1)=DISFY/SQRT(DFMOD)
                DIREC(L,1,2)=DISFZ/SQRT(DFMOD)
            ENDIF
110    CONTINUE      
C
C Now complete the definition of DIREC(NLAY,2,2) containing,
C for each layer, a set of
C (at most 2) 2-dimensional vectors which respect the symmetry
C specified for the search. DIREC will be used in SETCOR and
C it is defined only for layers with non trivial dimensionality(>1)
C
      DO 115 II=1,2
C
C II=1 corresponds to the first of the possible two vectors(it will be the only
C one if NDIML(K)=2).In each case define a new set of displaced vectors
C
        KK=0
        DO 120 I=1,NLAY
         ADISP(I,1)=.0
         ADISP(I,2)=.0
         ADISP(I,3)=.0
C
C define the displacement for nonequivalent atoms only
C
         IF(LSFLAG(I).GT.KK) THEN
            KK=KK+1
            IF(NDIML(KK).GT.1) THEN
              ADISP(I,1)=.0
              ADISP(I,2)=DIREC(I,II,1)
              ADISP(I,3)=DIREC(I,II,2)
            ENDIF
          ENDIF
120      CONTINUE
        CALL DISDOM(ADISP,ACOORD,ILOOK,NLAY,1)
C
C Complete the definition of DIREC
C
        DO 125 K=1,NEF
           IF(NDIML(K).GT.1.AND.NEQ(K).GT.1) THEN
            DO 127 KK=1,NEQ(K)
              L1=NL(1,K)
              L=NL(KK,K)
C
C Identify the symmetry operation moving L1 into L
C
              DO 129 IS=1,NDOM
                  IF(ILOOK(IS,L).EQ.L1) THEN
C                     L2=ILOOK(IS,L1)
C
C If there is more than one symmetry operation the DIREC should be
C neverthless identical
C
                     DIREC(L,II,1)=ACOORD(IS,L,2)
                     DIREC(L,II,2)=ACOORD(IS,L,3)
                  ENDIF
129           CONTINUE
127         CONTINUE
           ENDIF
125     CONTINUE
115   CONTINUE
C    
C Now compute NNDIM (the dimensionality of the search)
C
       NNDIM=0
       DO 145 I=1,NEF
          NNDIM=NNDIM+NDIML(I)
145    CONTINUE
      IF(NPASS.EQ.1) THEN
         DO 150 I=1,NLAY
             LLFLAG(I)=LLFLAGC(I)
150       CONTINUE
      ENDIF
C    
C Now compute CINEQ (the dimensionality of the search)
C
C      DO 150 I=1,NLAY
C        WRITE(*,*) 'NLAY=',I 
C      DO 151 J=1,2
C        WRITE(*,177) (DIREC(I,J,K),K=1,2)
C151   CONTINUE
C150   CONTINUE
C      DO 152 I=1,NEF
C        WRITE(*,*) NDIML(I),NEQ(I)
C152   CONTINUE
1000  RETURN
C++++
1999  write (6,1998) nef
1998  format (' ERROR DIMSCH2: PARAMETER INLIN MUST BE AT LEAST',I5)
      stop
      END
C ========================================================================
C
C Subroutine DISDOM takes the input coordinates ADISP and uses them to
C generate the coordinates for all the symmetry related domains.
C
C Input Parameters;
C =================
C
C ADISP      =   CURRENT DISPLACEMENT
C NLAY       =   NUMBER OF LAYERS INCALCULATION
C ACOORD     =   OUTPUT STORAGE OF COORDINATES FOR EACH DOMAIN
C ILOOK      =   LOOKUP TABLE OF RELATED ATOMS
C NSYM       =   SYMMTERY CODE OF SURFACE  (in common block)
C NDOM       =   NUMBER OF SYMMETRY EQUIVALENT DOMAINS  (in common block)
C EI,EF,DE   =   ENERGY RANGE (in common block-not used)
C NERG       =   NUMBER OF ENERGY POINTS (=INT((EF-EI)/DE)+1)(in common
C               block-not used)
C NSYMS       =   SYMMETRY CODE SPECIFYING THE TYPE OF DISPLACEMENTS
C                ALLOWED BY THE SEARCH
C NDOMS       =   NUMBER OF SYMMETRY OPERATIONS COORESPONDING TO NSYM
C NSCH       =   flag specifying whether DISDOM is being used for domain
C                averaging(0), or for search purposes (1)
C VV         =   REAL PART OF INNER POTENTIAL (in common block-not used)
C
C AUTHOR: WANDER
C MODIFICATIONS: BARBIERI
C
C =========================================================================
C
      SUBROUTINE DISDOM(ADISP,ACOORD,ILOOK,NLAY,NSCH)
C
      DIMENSION ADISP(NLAY,3),ILOOK(12,NLAY),ACOORD(12,NLAY,3)
C
      COMMON /ENY/EI,EF,DE,NERG,NSYM,NDOM,VV,VPIS
      COMMON /DIM/NDOMS,NSYMS
C
1500  FORMAT (' SYMMETRY CODES >15 NOT CURRENTLY IMPLEMENTED')
C
C GENERATE PREFACTORS
C
      NDOM1=NDOMS
      NSYM1=NSYMS
      IF(NSCH.EQ.0) THEN
         NDOM1=NDOM
         NSYM1=NSYM
      ENDIF
      PI=0.0
      PI=2.0*ACOS(PI)
C
C SET FIRST DOMAIN FOR EACH LAYER
C
      DO 100 IJK=1,NLAY
         DO 110 J=1,3
            ACOORD(1,IJK,J)=ADISP(ILOOK(1,IJK),J)
110      CONTINUE
         IF (NDOM1.GT.1) THEN
C
C X-COORDINATE DOES NOT ROTATE SO SET THIS FIRST
C
            DO 120 LLL=2,NDOM1
               ACOORD(LLL,IJK,1)=ADISP(ILOOK(LLL,IJK),1)
120         CONTINUE
         ENDIF
100   CONTINUE
      DO 130 IJK=1,NLAY
C
C BRANCH ACCORDING TO VALUE OF ISYAW
C FIRST NO SYMMETRY 
C
         IF (NSYM1.NE.1) THEN
C
C TWO FOLD ROTATION AXIS
C
            IF (NSYM1.EQ.2) THEN
               DO 200 I=2,3
C
C MIIROR THROUGH 45 DEG (X=-Y)
C
                  ACOORD(2,IJK,I)=-ADISP(ILOOK(2,IJK),I)
200            CONTINUE
            ELSEIF (NSYM1.EQ.3) THEN
               ACOORD(2,IJK,2)=-ADISP(ILOOK(2,IJK),3)
               ACOORD(2,IJK,3)=-ADISP(ILOOK(2,IJK),2)
C
C MIIROR THROUGH 45 DEG (X=Y)
C
            ELSEIF (NSYM1.EQ.4) THEN
               ACOORD(2,IJK,2)=ADISP(ILOOK(2,IJK),3)
               ACOORD(2,IJK,3)=ADISP(ILOOK(2,IJK),2)
C
C MIRROR THROUGH Z=0
C
            ELSEIF (NSYM1.EQ.5) THEN
               ACOORD(2,IJK,2)=ADISP(ILOOK(2,IJK),2)
               ACOORD(2,IJK,3)=-ADISP(ILOOK(2,IJK),3)
C
C MIRROR OR GLIDE THROUGH Y=0
C
            ELSEIF (NSYM1.EQ.6.OR.NSYM1.EQ.14) THEN
               ACOORD(2,IJK,2)=-ADISP(ILOOK(2,IJK),2)
               ACOORD(2,IJK,3)=ADISP(ILOOK(2,IJK),3)
C
C 4-FOLD AXIS
C
            ELSEIF (NSYM1.EQ.7) THEN
C
C SET UP DOMAIN 2 AS INVERSION
C
               ACOORD(2,IJK,2)=-ADISP(ILOOK(2,IJK),2)
               ACOORD(2,IJK,3)=-ADISP(ILOOK(2,IJK),3)
C
C NOW ROTATE
C
               DO 700 I=3,4
C
C MIRROR PLANES Y=Z,Y=-Z
C
                  ACOORD(I,IJK,2)=ADISP(ILOOK(I,IJK),3)*(-1)**I
                  ACOORD(I,IJK,3)=ADISP(ILOOK(I,IJK),2)*(-1)**(I+1)
700            CONTINUE
            ELSEIF (NSYM1.EQ.8) THEN
C
C FIRST INVERT
C
               ACOORD(2,IJK,2)=-ADISP(ILOOK(2,IJK),2)
               ACOORD(2,IJK,3)=-ADISP(ILOOK(2,IJK),3)
C
C NOW MIRROR
C
               DO 800 I=3,4
                  ACOORD(I,IJK,2)=ADISP(ILOOK(I,IJK),3)*(-1)**(I+1)
                  ACOORD(I,IJK,3)=ADISP(ILOOK(I,IJK),2)*(-1)**(I+1)
800            CONTINUE
C
C MIRROR PLANES Y=0,Z=0
C
            ELSEIF (NSYM1.EQ.9.OR.NSYM1.EQ.15.OR.NSYM1.EQ.16) THEN
               DO 900 I=2,4
                  ACOORD(I,IJK,2)=ADISP(ILOOK(I,IJK),2)*(-1)**I
                  ACOORD(I,IJK,3)=ADISP(ILOOK(I,IJK),3)*(-1)**(I+1)
900            CONTINUE
               ACOORD(4,IJK,2)=-ACOORD(4,IJK,2)
C
C 4 MIRROR PLANES Y=0,Z=0,Y=Z,Y=-Z
C
            ELSEIF (NSYM1.EQ.10) THEN
               DO 1000 I=2,4
                  ACOORD(I,IJK,2)=ADISP(ILOOK(I,IJK),2)*(-1)**I
                  ACOORD(I,IJK,3)=ADISP(ILOOK(I,IJK),3)*(-1)**(I+1)
1000           CONTINUE
               ACOORD(4,IJK,2)=-ACOORD(4,IJK,2)
               ACOORD(5,IJK,2)=ADISP(ILOOK(5,IJK),3)
               ACOORD(5,IJK,3)=ADISP(ILOOK(5,IJK),2)
               DO 1001 I=6,8
                  ACOORD(I,IJK,2)=ADISP(ILOOK(I,IJK),3)*(-1)**(I+1)
                  ACOORD(I,IJK,3)=ADISP(ILOOK(I,IJK),2)*(-1)**I
1001           CONTINUE
               ACOORD(8,IJK,3)=-ACOORD(8,IJK,3)
            ELSE
C
C SYMMETRY CODES GT 10 ALL REQUIRE AT LEAST A THREE FOLD ROTATION, SO FIRST
C ROTATE.
C
               FACT1=COS(2*PI/3)
               FACT2=SIN(2*PI/3)
               FACT3=COS(4*PI/3)
               FACT4=SIN(4*PI/3)
               ACOORD(2,IJK,2)=ADISP(ILOOK(2,IJK),2)*FACT1
               ACOORD(2,IJK,2)=ACOORD(2,IJK,2)-ADISP(ILOOK(2,IJK),3)
     &          *FACT2
               ACOORD(2,IJK,3)=ADISP(ILOOK(2,IJK),2)*FACT2
               ACOORD(2,IJK,3)=ACOORD(2,IJK,3)+ADISP(ILOOK(2,IJK),3)
     &          *FACT1
               ACOORD(3,IJK,2)=ADISP(ILOOK(3,IJK),2)*FACT3
               ACOORD(3,IJK,2)=ACOORD(3,IJK,2)-ADISP(ILOOK(3,IJK),3)
     &          *FACT4
               ACOORD(3,IJK,3)=ADISP(ILOOK(3,IJK),2)*FACT4
               ACOORD(3,IJK,3)=ACOORD(3,IJK,3)+ADISP(ILOOK(3,IJK),3)
     &          *FACT3
               IF (NSYM1.NE.11) THEN
C
C PRODUCE REFLECTIONS
C
                  IF (NSYM1.EQ.12.OR.NSYM1.EQ.13) THEN
                     IAWFAC=1
                     IF (NSYM1.EQ.13) IAWFAC=-1
                     ACOORD(4,IJK,2)=ADISP(ILOOK(4,IJK),2)*IAWFAC
                     ACOORD(4,IJK,3)=-ADISP(ILOOK(4,IJK),3)*IAWFAC
                     ACOORD(5,IJK,2)=ADISP(ILOOK(5,IJK),2)*FACT1
                     ACOORD(5,IJK,2)=ACOORD(5,IJK,2)-ADISP(ILOOK(5,IJK)
     &                ,3)*FACT2
                     ACOORD(5,IJK,2)=ACOORD(5,IJK,2)*IAWFAC
                     ACOORD(5,IJK,3)=ADISP(ILOOK(5,IJK),2)*FACT2
                     ACOORD(5,IJK,3)=ACOORD(5,IJK,3)+ADISP(ILOOK(5,IJK)
     &                ,3)*FACT1
                     ACOORD(5,IJK,3)=-ACOORD(5,IJK,3)*IAWFAC
                     ACOORD(6,IJK,2)=ADISP(ILOOK(6,IJK),2)*FACT3
                     ACOORD(6,IJK,2)=ACOORD(6,IJK,2)-ADISP(ILOOK(6,IJK)
     &                ,3)*FACT4
                     ACOORD(6,IJK,2)=ACOORD(6,IJK,2)*IAWFAC
                     ACOORD(6,IJK,3)=ADISP(ILOOK(6,IJK),2)*FACT4
                     ACOORD(6,IJK,3)=ACOORD(6,IJK,3)+ADISP(ILOOK(6,IJK)
     &                ,3)*FACT3
                     ACOORD(6,IJK,3)=-ACOORD(6,IJK,3)*IAWFAC
                  ELSEIF (NSYM1.GT.15) THEN
                     GOTO 1002
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
130   CONTINUE
      RETURN
1002  WRITE (1,1500)
      STOP
      END
C ========================================================================
C
C Subroutine FOLDCOMP (called from LOOKUP) compare symmetry transformed 
C coordinates of the surface composite layer to produce entries
C for the lookup table. For details see LOOKUP
C
C Input Parameters;
C =================
C
C NLAY    =  number of layers in surface composite layer
C IJK     =  one of the NLAY atoms
C FPOSTR  =  symmetry transformed coordinates of atoms in the top comp. layer
C FPOSFD  =  folded coordinates of atoms in the top composite layer
C Q1X etc =  coordinates (cartesian) of ARB1 and ARB2(in common block)
C A1MOD   =  modulus of ARB1(in common)
C A2MOD   =  modulus of ARB2(in common)
C DET     =  area spanned bye ARB1 and ARB2(in common) 
C 
C The variable EPS sets the tolerance with which we identify two points
C on the lattice. This depends a lot on the accuracy of th input geometry.
C If in single precision eps=1.0e-3 should always work. 
C AUTHOR: BARBIERI
C
C ===========================================================================
C
      SUBROUTINE FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
C
      DIMENSION FPOSFD(NLAY,3),FPOSTR(NLAY,3),TRFD(3)
C
      COMMON /FOL/Q1X,Q1Y,Q2X,Q2Y,A1MOD,A2MOD,DET
Cgss
      integer nerror_stat
C
1000  FORMAT('either EPS is too small, or (more likely) NSYM is not
     %consistent with the geometry for the sorted coordinate #  ',I4)
      EPS=5.0E-2
      Q1Q2=2.0*(Q1X*Q2X + Q1Y*Q2Y)
C
C no change in X coordinate
C
      TRFD(1)=FPOSTR(IJK,1)
C new coordinates in y-z plane
C
      TRFD(2)=(Q2Y*FPOSTR(IJK,2) - Q2X*FPOSTR(IJK,3))/DET
      TRFD(3)=(-Q1Y*FPOSTR(IJK,2) + Q1X*FPOSTR(IJK,3))/DET
C
C now fold coordinates. Be careful of the definition of AMOD for
C negative variables to produce only inequivalent points
C
      TRFD(2)=AMOD(TRFD(2),1.0)
      TRFD(3)=AMOD(TRFD(3),1.0)
      IF (TRFD(2).LT.0.D0) TRFD(2)=TRFD(2)+1.0
      IF (TRFD(3).LT.0.D0) TRFD(3)=TRFD(3)+1.0
C
C compare TRFD and FPOSFD to determine the entry in the lookup table
C This is done by computing the distance between the TRFD vector
C and all the FPOSFD vectors. Notice that the distance is computed by
C taking into account the non-orthogonality of coordinates.
C
      DO 10 I=1,NLAY
            X=TRFD(1) - FPOSFD(I,1)
            Y=TRFD(2) - FPOSFD(I,2)
            YY=ABS(Y)-1.
            IF(ABS(YY).LE.1.E-3) Y=.0
            Z=TRFD(3) - FPOSFD(I,3)
            ZZ=ABS(Z)-1.
            IF(ABS(ZZ).LE.1.E-3) Z=.0
            DIST= X**2 + Y**2*A1MOD**2 + Z**2*A2MOD**2 + Z*Y*Q1Q2
            IF (DIST.LE.EPS) THEN
               NNN=I
               nerror_stat=0
               GOTO 20
            ENDIF
10    CONTINUE
      WRITE(1,*) 'PROBLEM IN FOLDCOMP'
      WRITE(1,1000)IJK 
               nerror_stat=1
Cgss      STOP
20    RETURN
      END
C======================================================================
C                                                                       
C Subroutine GAUNT generated a vector of Gaunt coefficients BLM in the
C same order as they are used in subroutine CVEC. Only non-zero elements
C are stored, and are selected by observing the usual selection rules;
C
C C(l1,m1; l3,m3; lp,mp)     = 0    Unless;
C                                 m1+m3-mp=0  (ie m1=mp-m3)
C                                l1+l3+lp is even
C
C Parameter List;
C ===============
C
C NLMB          =   DIMENSION OF BELM (SEE BELOW)
C LMAX          =   LARGEST L-VALUE USED TO DESCRIBE ATOMIC
C                   SCATTERING WITHIN THE REFERENCE SURFACE.
C LSMAX         =   LARGEST L-VALUE FOR SINGLE CENTRE EXPANSION
C                   (N.B. LSMAX<=2*LMAX.)
C BELM(NLMB)    =   VECTOR OF GAUNT COEFFS.
C
C
C =========================================================================
C
      SUBROUTINE GAUNT2(BELM,NLMB,LMAX,LSMAX)
C
      DIMENSION BELM(NLMB)
C      
      K=0
      LSMMAX=(LSMAX+1)*(LSMAX+1)
      LMMAX=(LMAX+1)*(LMAX+1)
      DO 12 I1=1,LMMAX
         L1=INT(SQRT(FLOAT(I1-1)))
         M1=I1-L1*L1-L1-1
         DO 1 I3=1,LSMMAX
            L3=INT(SQRT(FLOAT(I3-1)))
            M3=I3-L3*L3-L3-1
            MP=M1+M3
            LL1=MAX0(IABS(MP),IABS(L1-L3))
C            LL2=MIN0(LMAX,L1+L3)
            LL2=L1+L3
C            IF (MOD(LL2+L1+L3,2).NE.0) LL2=LL2-1
            DO 2 LP=LL2,LL1,-2
               K=K+1
               BELM(K)=BLMT2(L3,M3,L1,M1,LP,-MP,LMAX,LL2)
C               WRITE(*,*) l3,m3,l1,m1,lp,mp
C               WRITE(*,*) BELM(K)
2           CONTINUE
1        CONTINUE
12    CONTINUE
      RETURN
      END
      FUNCTION BLMT2(L1,M1,L2,M2,L3,M3,LMAX,LMAX3)
C
cjcm      DOUBLEPRECISION a,b,bn,c,cn
      REAL a,b,bn,c,cn
C                                ,blmt
C
40    FORMAT (' INVALID ARGUMENTS FOR BLMT ',6(I3,','))
C
      PI=3.14159265
      IF (M1+M2+M3.EQ.0) THEN
         IF (L1-LMAX-LMAX.LE.0) THEN
            IF (L2.LE.LMAX) THEN
               IF (L3.LE.LMAX3) THEN
                  IF (L1.GE.IABS(M1)) THEN
                     IF (L2.GE.IABS(M2)) THEN
                        IF (L3.GE.IABS(M3)) THEN
                           IF (MOD(L1+L2+L3,2).NE.0) GOTO 420
                           NL1=L1
                           NL2=L2
                           NL3=L3
                           NM1=IABS(M1)
                           NM2=IABS(M2)
                           NM3=IABS(M3)
                           IC=(NM1+NM2+NM3)/2
                           IF (MAX0(NM1,NM2,NM3).GT.NM1) THEN
                              IF (MAX0(NM2,NM3).GT.NM2) THEN
                                 NL1=L3
                                 NL3=L1
                                 NM1=NM3
                                 NM3=IABS(M1)
                              ELSE
                                 NL1=L2
                                 NL2=L1
                                 NM1=NM2
                                 NM2=IABS(M1)
                              ENDIF
                           ENDIF
                           IF (NL2.LT.NL3) THEN
                              NTEMP=NL2
                              NL2=NL3
                              NL3=NTEMP
                              NTEMP=NM2
                              NM2=NM3
                              NM3=NTEMP
                           ENDIF
                           IF (NL3.GE.IABS(NL2-NL1)) THEN
C
C      CALCULATION OF FACTOR \A\
C
                              IS=(NL1+NL2+NL3)/2
                              IA1=IS-NL2-NM3
                              IA2=NL2+NM2
                              IA3=NL2-NM2
                              IA4=NL3+NM3
                              IA5=NL1+NL2-NL3
                              IA6=IS-NL1
                              IA7=IS-NL2
                              IA8=IS-NL3
                              IA9=NL1+NL2+NL3+1
                              A=((-1.0)**IA1)/FACT(IA3)*FACT(IA2)
     &                         /FACT(IA6)*FACT(IA4)
                              A=A/FACT(IA7)*FACT(IA5)/FACT(IA8)
     &                         *FACT(IS)/FACT(IA9)
                              A=A*(10.0**(IA2-IA3+IA4+IA5-IA6-IA7
     &                         -IA8-IA9+IS))
C
C      CALCULATION OF SUM \B\
C
                              IB1=NL1+NM1
                              IB2=NL2+NL3-NM1
                              IB3=NL1-NM1
                              IB4=NL2-NL3+NM1
                              IB5=NL3-NM3
                              IT1=MAX0(0,-IB4)+1
                              IT2=MIN0(IB2,IB3,IB5)+1
                              B=0.
                              SIGN=(-1.0)**(IT1)
                              IB1=IB1+IT1-2
                              IB2=IB2-IT1+2
                              IB3=IB3-IT1+2
                              IB4=IB4+IT1-2
                              IB5=IB5-IT1+2
                              DO 520 IT=IT1,IT2
                                 SIGN=-SIGN
                                 IB1=IB1+1
                                 IB2=IB2-1
                                 IB3=IB3-1
                                 IB4=IB4+1
                                 IB5=IB5-1
                                 BN=SIGN/FACT(IT-1)*FACT(IB1)
     &                            /FACT(IB3)*FACT(IB2)
                                 BN=BN/FACT(IB4)/FACT(IB5)
                                 BN=BN*(10.0**(IB1+IB2-IB3-IB4-IB5
     &                            -IT+1))
                                 B=B+BN
520                           CONTINUE
C
C      CALCULATION OF FACTOR \C\
C
                              IC1=NL1-NM1
                              IC2=NL1+NM1
                              IC3=NL2-NM2
                              IC4=NL2+NM2
                              IC5=NL3-NM3
                              IC6=NL3+NM3
                              CN=FLOAT((2*NL1+1)*(2*NL2+1)*(2*NL3+1))
     &                         /PI
                              C=CN/FACT(IC2)*FACT(IC1)/FACT(IC4)
     &                         *FACT(IC3)/FACT(IC6)*FACT(IC5)
                              C=C*(10.0**(IC1-IC2+IC3-IC4+IC5-IC6))
                              C=(SQRT(C))/2.
                              BLMT2=((-1.0)**IC)*A*B*C
                           ELSE
                              BLMT2=0.0
                           ENDIF
                           GOTO 531
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         WRITE (1,40) L1,M1,L2,M2,L3,M3
         GOTO 531
      ENDIF
420   BLMT2=0.0
531   RETURN
      END
C ========================================================================
C
C Subroutine INEQCD generates the coordinates  of symmetry inequivalent
C atoms. The symmetry is specified by the input parameter NSYMS.
C
C Input Parameters;
C =================
C
C NST1       =   NUMBER OF COMPOSITE LAYERS IN CALCULATION
C NCL        =   COMPOSITE LAYER UNDER CONSIDERATION
C NLAY2      =   MAXIMUM NUMBER OF ATOMS IN COMPOSITE LAYERS
C NLAY3      =   TOTAL # of ATOMS IN ALL CL'S
C LAFLAG     =  ARRAY SPECIFYING # OF ATOMS per CL
C LSFLAG     =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C               NSYMS) IN THE COMPOSITE LAYER. LSFLAG(i)=LSFLAG(j)
C               INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C               AS EQUIVALENT IN THE SEARCH
C
C OUTPUT
C =============
C
C LAFLAG2    = # of inequivalent atoms for each ineq. sublayer
C NLTIN      = total # of inequivalent layers in the NSTEF CL's
C LPOINT(i)  = for i=1,total # ineq atoms is the label of the ith
C              inequivalent atoms (1<LPOINT(i)<NLAY3
C LPBD(i,NCL)  i=1,total # inequiv atoms in plane, point to the corresponding
C              atom, in the NCL composite layer.
C AUTHOR: BARBIERI
C
C =========================================================================
C
C
      SUBROUTINE INEQCD(LSFLAG,LAFLAG,NST1,NLAY2,NLAY3
     %,LAFLAG2,LPOINT,NLTIN,NSTEF,NINEQ,LPBD)
C
      DIMENSION LAFLAG(NST1),LSFLAG(NLAY3),LAFLAG2(NLAY2,NST1)
      DIMENSION LPOINT(NLAY3),NINEQ(NST1)
      DIMENSION LPBD(NLAY2,NST1)
C
      NAT2=0
      NLTIN=0
      DO 5 NCL=1,NSTEF
        NSH=0
        NMIN=1000
        NMAX=0
        DO 50 I=1,NCL-1
          NSH=NSH+LAFLAG(I)
50      CONTINUE
        NLAY=LAFLAG(NCL)
        DO 55 I=1,NLAY
          IF (LSFLAG(I+NSH).LT.NMIN)NMIN=LSFLAG(I+NSH)
          IF (LSFLAG(I+NSH).GT.NMAX)NMAX=LSFLAG(I+NSH)
55      CONTINUE
        NT=NMAX-NMIN+1
C
C select RX,RY,RZ
C
        NAT=0
        DO 70 I=1,NLAY
C
C take only the first inequivalent atom and set pointer to that atom
C
          IF (LSFLAG(I+NSH).GE.NMIN) THEN
             NAT=NAT+1
             NAT2=NAT2+1
             LAFLAG2(NAT,NCL)=1
             DO 77 II=1,NLAY
               IF(LSFLAG(II+NSH).EQ.LSFLAG(I+NSH))
     &          LAFLAG2(NAT,NCL)=LAFLAG2(NAT,NCL)+1
77           CONTINUE
             LAFLAG2(NAT,NCL)=LAFLAG2(NAT,NCL)-1
             NLTIN=NLTIN+1
             LPOINT(NAT2)=I+NSH
             LPBD(NAT,NCL)=I
             NMIN=NMIN+1
          ENDIF
70      CONTINUE
        NINEQ(NCL)=NAT
5     CONTINUE

      RETURN
      END
C ========================================================================
C
C Subroutine LOOKUP generates the lookup table (previously in look.i)
C used in DISDOM
C
C Input Parameters;
C =================
C
C NST1       =   NUMBER OF COMPOSITE LAYERS. INPUT DATA SHOULD
C                BE SUCH THAT THE FIRST COMPOSITE LAYER CORRESPONDS
C                TO THE SURFACE LAYER. IF INSTEAD THE LAST LAYER IS AT
C                THE SURFACE SEARCH FOR FLAGNST1 AND MODIFY THE LINE
C NLAY       =   NUMBER OF LAYERS INCALCULATION
C CPVPOS     =   presorted cartesian coordinates of the composite layers
C ARB1,ARB2  =   OVERLAYER LATTICE VECTORS
C FPOS       =   CARTESIAN COORDINATES OF ATOMS IN TOP COMP. LAYER
C FPOSFD     =   CARTESIAN COORDINATES OF ATOMS IN TOP COMP. LAYER
C                AFTER FOLDING IN THE REGION (x*ARB1,y*ARB2) WITH
C                0<=x<1, 0<=y<1
C FPOSTR     =   COORDINATES OF ATOMS IN TOP COMP. LAYER AFTER SYMMETRY
C                TRANFORMATION
C FPOSTRFD   =   FOLDED FPOSTR
C ILOOK      =   OUTPUT LOOKUP TABLE OF RELATED ATOMS (ILOOK(IS,L)=M
C                MEANS THAT IS[M]=L
C ILKBD      =   SIMILAR to ILOOK with domains ordered in a different
C                way. First rotations by 2PI/IROT (i-1) i=1,IROT,
C                then rotations in the same order followed by the mirror
C IROT       =   output IROT=N if the surface has an N fold rotation
C                axis
C IMIR       =   output mirror plane symmetry code of the surface.
C                IMIR =1 if (surface has the Y=0 mirror plane)
C                IMIR=2 (Z=0 mirror), IMIR=3 (Y=Z),IMIR=4 (Y=-Z),
C                IMIR=5 (Z=0 and Y=0), IMIR=6 (Z=Y and Z=-Y) and
C                nore can be added. IMIR=0 if no mirror is present.
C                IROT and IMIR are used to reduce the tensor. For the
C                reduction to work properly it is sufficient that the
C                N-fold axis and the mirror symmetry be independent
C                and generate the whole symmetry. (e.g. P3M1 has (3,1),
C                an FCC(100) surface has (4,1) or (4,2) but (4,5) is
C                redundant because if you have a 4fold axis and Y=0,
C                we have also Z=0 ect.
C NSYM       =   SYMMETRY CODE OF SURFACE  (in common block)
C
C In common blocks
C CPVPOS is obtained from the subroutine READCT and corresponds
C to the coordinates (cartesian) of the atoms in the composite layers.
C This should correspond to the input data with the first composite
C layer being at the surface. CPARB1 and CPARB2 are copies of ARB1 and ARB2
C in units of angstroms. They come from subroutine READT.
C
C AUTHOR: BARBIERI
C
C =========================================================================
C
C FPOS contains the coordinates of the atoms in the surface composite layer.
C as well as the bulk layers needed for the substrate.
C We label these atoms and, for each symmetry operation, we transform
C their coordinates. Because the transformation corresponds to a symmetry
C the new points obtained this way are equivalent(modulus lattice vectors
C translation) to the original set of points modulus a permutation.
C This permutation defines the new entry for the lookup table corresponding
C to the symmetry under consideration.
C
      SUBROUTINE LOOKUP2(CPVPOS,ILOOK,ILKBD,NLAY,NST1,NSTEF,NSYM,NTOT,
     & NSH,NCL,NLAY2,IROT,IMIR,NWC,FPOS,FPOSFD,FPOSTR,nerror_report)
C
      DIMENSION FPOS(NLAY,3),FPOSFD(NLAY,3),FPOSTR(NLAY,3)
      DIMENSION CPARB1(2),CPARB2(2),CPVPOS(NST1,NLAY2,3)
      DIMENSION ILKBD(12,NLAY2,NST1),ILOOK(12,NTOT)
      DIMENSION RBR1(2),RBR2(2),ARB1(2),ARB2(2),ARA1(2),ARA2(2)
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C      COMMON /LO/CPARB1,CPARB2
      COMMON /FOL/Q1X,Q1Y,Q2X,Q2Y,A1MOD,A2MOD,DET
C
CGss
      Integer nerror_report, nerror_stat
Cgss default nerror_report and initialize nerror_stat
      nerror_report=0 
      nerror_stat=0
1500  FORMAT (' SYMMETRY CODES >16 NOT CURRENTLY IMPLEMENTED')
1600  FORMAT (' SYMMETRY CODE = 14 or 16  HAS a Y=0 GLIDE. ARB1 or ARB2 
     & SHOULD BE ALONG Y=0 ')
1610  FORMAT (' SYMMETRY CODE = 15 REQUIRES ARB1 AND ARB2 ALONG Y=0 and
     & Z=0 ')
1010  FORMAT (7X,F12.8,2(5X,F12.8))
C1010  FORMAT (7X,F10.4,2(5X,F10.4))
C
C GENERATE PREFACTORS
C
      PI=0.0
      PI=2.0*ACOS(PI)
C
C generate CPVPOS. For tensor calculation use overlayer vectors, for
C substrate use substrate lattice vectors
C
      IF(NCL.LE.NSTEF) THEN
         DO 40 I=1,2
           CPARB1(I)=ARB1(I)*.529
           CPARB2(I)=ARB2(I)*.529
40       CONTINUE
      ELSE
         DO 41 I=1,2
           CPARB1(I)=ARA1(I)*.529
           CPARB2(I)=ARA2(I)*.529
41       CONTINUE
      ENDIF
C GENERATE FPOS
C FLAGNST1: IF THE LAST COMPOSITE LAYER IS AT THE SURFACE CPVPOS(1,I,J)
C           SHOULD BE REPLACED BY CPVPOS(NST1,I,J)
      DO 50 I=1,NLAY
             DO 70 J=1,3
                  FPOS(I,J)=CPVPOS(NCL,I,J)
70           CONTINUE
      IF (NWC.EQ.1) WRITE(2,1010) (FPOS(I,K),K=1,3)
CGPS      IF (NWC.EQ.1) WRITE(0,1010) (FPOS(I,K),K=1,3)
50    CONTINUE
C
C NOW FOLD THE FPOS COORDINATES IN THE DOMAIN (x*ARB1,y*ARB2) 0<=x<1,0<=y<1
C This is done more conveniently by expressing the cartesian coordinates in
C the the y-z plane using ARB1 and ARB2 as basis vectors. Hence while FPOS
C is expressed in orthogonal coordinates, FPOSFD will be expressed in
C "X,ARB1,ARB2" coordinates.
C
      A1MOD=SQRT(CPARB1(1)**2 + CPARB1(2)**2)
      A2MOD=SQRT(CPARB2(1)**2 + CPARB2(2)**2)
      Q1X=CPARB1(1)
      Q2X=CPARB2(1)
      Q1Y=CPARB1(2)
      Q2Y=CPARB2(2)
      DET=Q1X*Q2Y - Q2X*Q1Y
      IF (ABS(DET).LE.1.0d-7) THEN
         WRITE (1,*) 'HEY MAN, THE OVERLAYER LATTICE VECTORS MUST BE
     &    INDEPENDENT'
      ENDIF
C
C new coordinates in y-z plane
C
      DO 51 I=1,NLAY
             FPOSFD(I,2)=(Q2Y*FPOS(I,2) - Q2X*FPOS(I,3))/DET
             FPOSFD(I,3)=(-Q1Y*FPOS(I,2) + Q1X*FPOS(I,3))/DET
51    CONTINUE
C
C now fold coordinates. Be careful of the definition of AMOD for
C negative variables. We only want inequivalent points
C First the X coordinate is unchanged
C
      DO 45 J=1,NLAY
            FPOSFD(J,1)=FPOS(J,1)
            FPOSTR(J,1)=FPOS(J,1)
45    CONTINUE
      DO 52 I=1,NLAY
             FPOSFD(I,2)=AMOD(FPOSFD(I,2),1.0) 
             FPOSFD(I,3)=AMOD(FPOSFD(I,3),1.0)
             IF (FPOSFD(I,2).LT.0.D0) FPOSFD(I,2)=FPOSFD(I,2)+1.0
             IF (FPOSFD(I,3).LT.0.D0) FPOSFD(I,3)=FPOSFD(I,3)+1.0
52    CONTINUE
C
C SET FIRST ENTRY IN LOOKUP TABLE
C
      DO 100 IJK=1,NLAY
            ILOOK(1,IJK+NSH)=IJK+NSH
            ILKBD(1,IJK,NCL)=IJK
100   CONTINUE
      DO 130 IJK=1,NLAY
C
C BRANCH ACCORDING TO VALUE OF ISYAW
C FIRST NO SYMMETRY 
C NOTICE THAT THE 1(X) COORDINATE IS NOT CHANGED BY THE TRANSFORMATION IN
C THE YZ PLANE
C
         IF (NSYM.EQ.1) THEN
            IROT=1
            IMIR=0
         ENDIF
         IF (NSYM.NE.1) THEN
C
C TWO FOLD ROTATION AXIS
C
            IF (NSYM.EQ.2) THEN
               IROT=2
               IMIR=0
               DO 200 I=2,3
                  FPOSTR(IJK,I)=-FPOS(IJK,I)
200            CONTINUE
C
C now fold transformed coordinates and compare them witn FPOSFD to produce
C the entry ILOOK(.,IJK) given by the output NNN
C
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
               if (nerror_stat.eq.1) then
                   nerror_report=1
                   goto 99999
               end if

               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C MIRROR THROUGH 45 DEG (X=-Y)
C
            ELSEIF (NSYM.EQ.3) THEN
               IROT=1
               IMIR=4
               FPOSTR(IJK,2)=-FPOS(IJK,3)
               FPOSTR(IJK,3)=-FPOS(IJK,2)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C MIIROR THROUGH 45 DEG (X=Y)
C
            ELSEIF (NSYM.EQ.4) THEN
               IROT=1
               IMIR=3
               FPOSTR(IJK,2)=FPOS(IJK,3)
               FPOSTR(IJK,3)=FPOS(IJK,2)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C MIRROR THROUGH Z=0
C
            ELSEIF (NSYM.EQ.5) THEN
               IROT=1
               IMIR=2
               FPOSTR(IJK,2)=FPOS(IJK,2)
               FPOSTR(IJK,3)=-FPOS(IJK,3)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C MIRROR THROUGH Y=0
C
            ELSEIF (NSYM.EQ.6) THEN
               IROT=1
               IMIR=1
               FPOSTR(IJK,2)=-FPOS(IJK,2)
               FPOSTR(IJK,3)=FPOS(IJK,3)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C GLIDE THROUGH Y=0
C
            ELSEIF (NSYM.EQ.14) THEN
               IROT=1
               IMIR=1
C
C check that one lattice vector is along Y=0 and generate shift
C GSHY (GSHZ) is the shift along the Y (Z) axis
C
               IF(ABS(CPARB1(1)).LE.0.00001)THEN
                  GSHZ=CPARB1(2)/2.
                  GSHY=0.
               ELSEIF (ABS(CPARB2(1)).LE.0.00001)THEN
                  GSHZ=CPARB2(2)/2.
                  GSHY=0.
               ELSE
                  WRITE (1,1600)
               ENDIF
               FPOSTR(IJK,2)=-FPOS(IJK,2)
               FPOSTR(IJK,3)=FPOS(IJK,3)+GSHZ
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C 4-FOLD AXIS
C
            ELSEIF (NSYM.EQ.7) THEN
               IROT=4
               IMIR=0
C
C SET UP DOMAIN 2 AS INVERSION 
C
               FPOSTR(IJK,2)=-FPOS(IJK,2)
               FPOSTR(IJK,3)=-FPOS(IJK,3)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(3,NNN,NCL)=IJK
C
C NOW ROTATE
C
               DO 700 I=3,4
                  FPOSTR(IJK,2)=FPOS(IJK,3)*(-1)**I
                  FPOSTR(IJK,3)=FPOS(IJK,2)*(-1)**(I+1)
                  CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                  ILOOK(I,NNN+NSH)=IJK+NSH
                  IF(I.EQ.3) THEN
                    ILKBD(2,NNN,NCL)=IJK
                  ELSE
                    ILKBD(4,NNN,NCL)=IJK
                  ENDIF
700            CONTINUE
            ELSEIF (NSYM.EQ.8) THEN
               IROT=2
               IMIR=3
C
C FIRST INVERT
C
               FPOSTR(IJK,2)=-FPOS(IJK,2)
               FPOSTR(IJK,3)=-FPOS(IJK,3)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C MIRROR PLANES Y=Z,Y=-Z
C
               DO 800 I=3,4
                  FPOSTR(IJK,2)=FPOS(IJK,3)*(-1)**(I+1)
                  FPOSTR(IJK,3)=FPOS(IJK,2)*(-1)**(I+1)
                  CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                  ILOOK(I,NNN+NSH)=IJK+NSH
                  ILKBD(I,NNN,NCL)=IJK
800            CONTINUE
C
C MIRROR PLANES Y=0,Z=0 AND INVERSION
C
            ELSEIF (NSYM.EQ.9) THEN
               IROT=2
               IMIR=1
               DO 900 I=2,4
                  FPOSTR(IJK,2)=FPOS(IJK,2)*(-1)**I
                  FPOSTR(IJK,3)=FPOS(IJK,3)*(-1)**(I+1)
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                    CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                    ILOOK(I,NNN+NSH)=IJK+NSH
                    IF(I.EQ.2) ILKBD(4,NNN,NCL)=IJK
                    IF(I.EQ.3) ILKBD(3,NNN,NCL)=IJK
                  ENDIF
900            CONTINUE
               FPOSTR(IJK,2)=-FPOSTR(IJK,2)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(4,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
C
C GLIDE PLANES Y=0,Z=0 (ILKBD should be checked for glides) 
C
            ELSEIF (NSYM.EQ.15) THEN
               IROT=1
               IMIR=5
C
C check that lattice vectors are along Y=0,Z=0 and generate shifts
C
               IF(ABS(CPARB1(1)).LE.0.001.AND.ABS(CPARB2(2)).LE.
     &           0.001) THEN
                  GSHY=CPARB2(1)/2.
                  GSHZ=CPARB1(2)/2.
               ELSEIF(ABS(CPARB2(1)).LE.0.001.AND.ABS(CPARB1(2)).LE.
     &           0.001) THEN
                  GSHY=CPARB1(1)/2.
                  GSHZ=CPARB2(2)/2.
               ELSE
                  WRITE (1,1610)
               ENDIF
               DO 901 I=2,4
                  IF(I.EQ.2) THEN
                     GSH1=GSHY
                     GSH2=0.
                  ELSEIF(I.EQ.3) THEN
                     GSH1=0.
                     GSH2=GSHZ
                  ELSE
                     GSH1=-GSHY
                     GSH2=-GSHZ
                  ENDIF
                  FPOSTR(IJK,2)=FPOS(IJK,2)*(-1)**I+GSH1
                  FPOSTR(IJK,3)=FPOS(IJK,3)*(-1)**(I+1)+GSH2
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                    CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                    ILOOK(I,NNN+NSH)=IJK+NSH
                    ILKBD(I,NNN,NCL)=IJK
                  ENDIF
901            CONTINUE
               FPOSTR(IJK,2)=-FPOSTR(IJK,2)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(4,NNN+NSH)=IJK+NSH
               ILKBD(4,NNN,NCL)=IJK
C
C GLIDE PLANE Y=0, MIRROR Z=0 
C
            ELSEIF (NSYM.EQ.16) THEN
               IROT=1
               IMIR=5
C
C check which lattice vector is along Y=0 and generate shifts
C
               IF(ABS(CPARB1(1)).LE.0.001.AND.ABS(CPARB2(2)).LE.
     &           0.001) THEN
                  GSHY=.0
                  GSHZ=CPARB1(2)/2.
               ELSEIF(ABS(CPARB2(1)).LE.0.001.AND.ABS(CPARB1(2)).LE.
     &           0.001) THEN
                  GSHY=.0
                  GSHZ=CPARB2(2)/2.
               ELSE
                  WRITE (1,1600)
               ENDIF
               DO 903 I=2,4
                  IF(I.EQ.2) THEN
                     GSH1=GSHY
                     GSH2=0.
                  ELSEIF(I.EQ.3) THEN
                     GSH1=0.
                     GSH2=GSHZ
                  ELSE
                     GSH1=-GSHY
                     GSH2=-GSHZ
                  ENDIF
                  FPOSTR(IJK,2)=FPOS(IJK,2)*(-1)**I+GSH1
                  FPOSTR(IJK,3)=FPOS(IJK,3)*(-1)**(I+1)+GSH2
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                    CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                    ILOOK(I,NNN+NSH)=IJK+NSH
                    ILKBD(I,NNN,NCL)=IJK
                  ENDIF
903            CONTINUE
               FPOSTR(IJK,2)=-FPOSTR(IJK,2)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(4,NNN+NSH)=IJK+NSH
               ILKBD(4,NNN,NCL)=IJK
C
C 4 MIRROR PLANES Y=0,Z=0,Y=Z,Y=-Z
C
            ELSEIF (NSYM.EQ.10) THEN
               IROT=4
               IMIR=1
               DO 1000 I=2,4
                  FPOSTR(IJK,2)=FPOS(IJK,2)*(-1)**I
                  FPOSTR(IJK,3)=FPOS(IJK,3)*(-1)**(I+1)
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                    CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                    ILOOK(I,NNN+NSH)=IJK+NSH
                    IF(I.EQ.2) ILKBD(7,NNN,NCL)=IJK
                    IF(I.EQ.3) ILKBD(5,NNN,NCL)=IJK
                  ENDIF
1000           CONTINUE
               FPOSTR(IJK,2)=-FPOSTR(IJK,2)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(4,NNN+NSH)=IJK+NSH
               ILKBD(3,NNN,NCL)=IJK
               FPOSTR(IJK,2)=FPOS(IJK,3)
               FPOSTR(IJK,3)=FPOS(IJK,2)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(5,NNN+NSH)=IJK+NSH
               ILKBD(6,NNN,NCL)=IJK
               DO 1001 I=6,8
                  FPOSTR(IJK,2)=FPOS(IJK,3)*(-1)**(I+1)
                  FPOSTR(IJK,3)=FPOS(IJK,2)*(-1)**I
                  IF (I.EQ.6.OR.I.EQ.7) THEN
                    CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                    ILOOK(I,NNN+NSH)=IJK+NSH
                    IF(I.EQ.6) ILKBD(2,NNN,NCL)=IJK
                    IF(I.EQ.7) ILKBD(4,NNN,NCL)=IJK
                  ENDIF
1001           CONTINUE
               FPOSTR(IJK,3)=-FPOSTR(IJK,3)
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(8,NNN+NSH)=IJK+NSH
               ILKBD(8,NNN,NCL)=IJK
            ELSE
               IROT=3
               IMIR=0
C
C SYMMETRY CODES GT 10 (except 14 and 15,16) ALL REQUIRE AT LEAST A THREE FOLD 
C ROTATION, SO FIRST ROTATE.
C
               FACT1=COS(2*PI/3)
               FACT2=SIN(2*PI/3)
               FACT3=COS(4*PI/3)
               FACT4=SIN(4*PI/3)
               FPOSTR(IJK,2)=FPOS(IJK,2)*FACT1
               FPOSTR(IJK,2)=FPOSTR(IJK,2)-FPOS(IJK,3)
     &          *FACT2
               FPOSTR(IJK,3)=FPOS(IJK,2)*FACT2
               FPOSTR(IJK,3)=FPOSTR(IJK,3)+FPOS(IJK,3)
     &          *FACT1
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(2,NNN+NSH)=IJK+NSH
               ILKBD(2,NNN,NCL)=IJK
               FPOSTR(IJK,2)=FPOS(IJK,2)*FACT3
               FPOSTR(IJK,2)=FPOSTR(IJK,2)-FPOS(IJK,3)
     &          *FACT4
               FPOSTR(IJK,3)=FPOS(IJK,2)*FACT4
               FPOSTR(IJK,3)=FPOSTR(IJK,3)+FPOS(IJK,3)
     &          *FACT3
               CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
               ILOOK(3,NNN+NSH)=IJK+NSH
               ILKBD(3,NNN,NCL)=IJK
               IF (NSYM.NE.11) THEN
C
C PRODUCE REFLECTIONS
C
                  IF (NSYM.EQ.12.OR.NSYM.EQ.13) THEN
                     IAWFAC=1
                     IF (NSYM.EQ.13) THEN
                        IAWFAC=-1
                        IMIR=1
                     ELSE
                        IMIR=2
                     ENDIF
                     FPOSTR(IJK,2)=FPOS(IJK,2)*IAWFAC
                     FPOSTR(IJK,3)=-FPOS(IJK,3)*IAWFAC
                     CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                     ILOOK(4,NNN+NSH)=IJK+NSH
                     ILKBD(4,NNN,NCL)=IJK
                     FPOSTR(IJK,2)=FPOS(IJK,2)*FACT1
                     FPOSTR(IJK,2)=FPOSTR(IJK,2)-FPOS(IJK
     &                ,3)*FACT2
                     FPOSTR(IJK,2)=FPOSTR(IJK,2)*IAWFAC
                     FPOSTR(IJK,3)=FPOS(IJK,2)*FACT2
                     FPOSTR(IJK,3)=FPOSTR(IJK,3)+FPOS(IJK
     &                ,3)*FACT1
                     FPOSTR(IJK,3)=-FPOSTR(IJK,3)*IAWFAC
                     CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                     ILOOK(5,NNN+NSH)=IJK+NSH
                     ILKBD(5,NNN,NCL)=IJK
                     FPOSTR(IJK,2)=FPOS(IJK,2)*FACT3
                     FPOSTR(IJK,2)=FPOSTR(IJK,2)-FPOS(IJK
     &                ,3)*FACT4
                     FPOSTR(IJK,2)=FPOSTR(IJK,2)*IAWFAC
                     FPOSTR(IJK,3)=FPOS(IJK,2)*FACT4
                     FPOSTR(IJK,3)=FPOSTR(IJK,3)+FPOS(IJK
     &                ,3)*FACT3
                     FPOSTR(IJK,3)=-FPOSTR(IJK,3)*IAWFAC
                     CALL FOLDCOMP(IJK,NLAY,FPOSTR,FPOSFD,NNN,
     & nerror_stat)
                     if (nerror_stat.eq.1) then
                        nerror_report=1
                        goto 99999
                     endif
                     ILOOK(6,NNN+NSH)=IJK+NSH
                     ILKBD(6,NNN,NCL)=IJK
                  ELSEIF (NSYM.GT.16) THEN
                     WRITE (1,1500)
                     nerror_report=1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
130   CONTINUE
99999 continue
      RETURN
      END
C ========================================================================
C
C Subroutine LOOKK generates list of symmetry equivalent beams
C
C Input Parameters;
C =================
C
C
C AUTHOR: BARBIERI
C
C =========================================================================
C
C
      SUBROUTINE LOOKK(SPQ,SPQS,KNT,NSYM,IDXK,NK,KMAX,J,NC)
C
      DIMENSION SPQ(2,KNT),IDXK(KNT)
      DIMENSION SPQS(2,KNT),ZVEC(2),ZTR(2),NC(KNT)
C
C
1500  FORMAT (' SYMMETRY CODES >16 NOT CURRENTLY IMPLEMENTED')
C
C vector to be symmetry transformed
C
      ZVEC(1)=SPQS(1,NK)
      ZVEC(2)=SPQS(2,NK)
      PI=0.0
      PI=2.0*ACOS(PI)
      EPS=0.001
C
C BRANCH ACCORDING TO VALUE OF ISYAW
C FIRST NO SYMMETRY 
C
         IF (NSYM.EQ.1) THEN
            IDXK(J)=NK
         ENDIF
         IF (NSYM.NE.1) THEN
C
C TWO FOLD ROTATION AXIS
C
            IF (NSYM.EQ.2) THEN
               DO 200 I=1,2
                  ZTR(I)=-ZVEC(I)
200            CONTINUE
C
C check which beams are equivalent to the transformed one         
C
               DO 210 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
210            CONTINUE
C
C MIRROR THROUGH 45 DEG (X=-Y)
C
            ELSEIF (NSYM.EQ.3) THEN
               ZTR(1)=-ZVEC(2)
               ZTR(2)=-ZVEC(1)
               DO 211 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
211            CONTINUE
C
C MIRROR THROUGH 45 DEG (X=Y)
C
            ELSEIF (NSYM.EQ.4) THEN
               ZTR(1)=ZVEC(2)
               ZTR(2)=ZVEC(1)
               DO 212 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
212            CONTINUE
C
C MIRROR THROUGH Z=0
C
            ELSEIF (NSYM.EQ.5) THEN
               ZTR(1)=ZVEC(1)
               ZTR(2)=-ZVEC(2)
               DO 213 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
213            CONTINUE
C
C MIRROR THROUGH Y=0
C
            ELSEIF (NSYM.EQ.6) THEN
               ZTR(1)=-ZVEC(1)
               ZTR(2)=ZVEC(2)
               DO 214 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
214            CONTINUE
C
C GLIDE THROUGH Y=0 (glides are equivalent to mirrors in K space)
C
            ELSEIF (NSYM.EQ.14) THEN
C
C
               ZTR(1)=-ZVEC(1)
               ZTR(2)=ZVEC(2)
               DO 215 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
215            CONTINUE
C
C 4-FOLD AXIS
C
            ELSEIF (NSYM.EQ.7) THEN
C
C SET UP DOMAIN 2 AS INVERSION 
C
               ZTR(1)=-ZVEC(1)
               ZTR(2)=-ZVEC(2)
               DO 216 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
216            CONTINUE
C
C NOW ROTATE
C
               DO 700 I=3,4
                  ZTR(1)=ZVEC(2)*(-1)**I
                  ZTR(2)=ZVEC(1)*(-1)**(I+1)
                  DO 217 K=J,J+KMAX
                    ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                    IF(ZMOD.LT.EPS) NC(K)=0
                    IF(ZMOD.LT.EPS) IDXK(K)=NK
217               CONTINUE
700            CONTINUE
            ELSEIF (NSYM.EQ.8) THEN
C
C FIRST INVERT
C
               ZTR(1)=-ZVEC(1)
               ZTR(2)=-ZVEC(2)
               DO 218 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
218            CONTINUE
C
C MIRROR PLANES Y=Z,Y=-Z
C
               DO 800 I=3,4
                  ZTR(1)=ZVEC(2)*(-1)**(I+1)
                  ZTR(2)=ZVEC(1)*(-1)**(I+1)
                  DO 219 K=J,J+KMAX
                    ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                    IF(ZMOD.LT.EPS) NC(K)=0
                    IF(ZMOD.LT.EPS) IDXK(K)=NK
219               CONTINUE
800            CONTINUE
C
C MIRROR PLANES Y=0,Z=0 AND INVERSION
C
            ELSEIF (NSYM.EQ.9) THEN
               DO 900 I=2,4
                  ZTR(1)=ZVEC(1)*(-1)**I
                  ZTR(2)=ZVEC(2)*(-1)**(I+1)
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                     DO 220 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
220                  CONTINUE
                  ENDIF
900            CONTINUE
               ZTR(1)=-ZTR(1)
               DO 221 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
221            CONTINUE
C
C GLIDE PLANES Y=0,Z=0 (ILKBD should be checked for glides) 
C
            ELSEIF (NSYM.EQ.15) THEN
               DO 901 I=2,4
                  ZTR(1)=ZVEC(1)*(-1)**I
                  ZTR(2)=ZVEC(2)*(-1)**(I+1)
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                     DO 222 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
222                  CONTINUE
                  ENDIF
901            CONTINUE
               ZTR(1)=-ZTR(1)
               DO 223 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
223            CONTINUE
C
C GLIDE PLANE Y=0, MIRROR Z=0 
C
            ELSEIF (NSYM.EQ.16) THEN
               DO 903 I=2,4
                  ZTR(1)=ZVEC(1)*(-1)**I
                  ZTR(2)=ZVEC(2)*(-1)**(I+1)
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                     DO 224 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
224                  CONTINUE
                  ENDIF
903            CONTINUE
               ZTR(1)=-ZTR(1)
               DO 225 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
225            CONTINUE
C
C 4 MIRROR PLANES Y=0,Z=0,Y=Z,Y=-Z
C
            ELSEIF (NSYM.EQ.10) THEN
               DO 1000 I=2,4
                  ZTR(1)=ZVEC(1)*(-1)**I
                  ZTR(2)=ZVEC(2)*(-1)**(I+1)
                  IF (I.EQ.2.OR.I.EQ.3) THEN
                     DO 226 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
226                  CONTINUE
                  ENDIF
1000           CONTINUE
               ZTR(1)=-ZTR(1)
               DO 227 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
227            CONTINUE
               ZTR(1)=ZVEC(2)
               ZTR(2)=ZVEC(1)
               DO 228 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
228            CONTINUE
               DO 1001 I=6,8
                  ZTR(1)=ZVEC(2)*(-1)**(I+1)
                  ZTR(2)=ZVEC(1)*(-1)**I
                  IF (I.EQ.6.OR.I.EQ.7) THEN
                     DO 229 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
229                  CONTINUE
                  ENDIF
1001           CONTINUE
               ZTR(2)=-ZTR(2)
               DO 230 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
230            CONTINUE
            ELSE
C
C SYMMETRY CODES GT 10 (except 14 and 15,16) ALL REQUIRE AT LEAST A THREE FOLD 
C ROTATION, SO FIRST ROTATE.
C
               FACT1=COS(2*PI/3)
               FACT2=SIN(2*PI/3)
               FACT3=COS(4*PI/3)
               FACT4=SIN(4*PI/3)
               ZTR(1)=ZVEC(1)*FACT1
               ZTR(1)=ZTR(1)-ZVEC(2)
     &          *FACT2
               ZTR(2)=ZVEC(1)*FACT2
               ZTR(2)=ZTR(2)+ZVEC(2)
     &          *FACT1
               DO 231 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
231            CONTINUE

               ZTR(1)=ZVEC(1)*FACT3
               ZTR(1)=ZTR(1)-ZVEC(2)
     &          *FACT4
               ZTR(2)=ZVEC(1)*FACT4
               ZTR(2)=ZTR(2)+ZVEC(2)
     &          *FACT3
               DO 232 K=J,J+KMAX
                  ZMOD=(ZTR(1)-SPQ(1,K))**2 +(ZTR(2)-SPQ(2,K))**2  
                  IF(ZMOD.LT.EPS) NC(K)=0
                  IF(ZMOD.LT.EPS) IDXK(K)=NK
232            CONTINUE
               IF (NSYM.NE.11) THEN
C
C PRODUCE REFLECTIONS
C
                  IF (NSYM.EQ.12.OR.NSYM.EQ.13) THEN
                     IAWFAC=1
                     IF (NSYM.EQ.13) IAWFAC=-1
                     ZTR(1)=ZVEC(1)*IAWFAC
                     ZTR(2)=-ZVEC(2)*IAWFAC
                     DO 233 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
233                  CONTINUE
                     ZTR(1)=ZVEC(1)*FACT1
                     ZTR(1)=ZTR(1)-ZVEC(2)*FACT2
                     ZTR(1)=ZTR(1)*IAWFAC
                     ZTR(2)=ZVEC(1)*FACT2
                     ZTR(2)=ZTR(2)+ZVEC(2)*FACT1
                     ZTR(2)=-ZTR(2)*IAWFAC
                     DO 234 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
234                  CONTINUE
                     ZTR(1)=ZVEC(1)*FACT3
                     ZTR(1)=ZTR(1)-ZVEC(2)*FACT4
                     ZTR(1)=ZTR(1)*IAWFAC
                     ZTR(2)=ZVEC(1)*FACT4
                     ZTR(2)=ZTR(2)+ZVEC(2)*FACT3
                     ZTR(2)=-ZTR(2)*IAWFAC
                     DO 235 K=J,J+KMAX
                        ZMOD=(ZTR(1)-SPQ(1,K))**2 +
     &                      (ZTR(2)-SPQ(2,K))**2  
                        IF(ZMOD.LT.EPS) NC(K)=0
                        IF(ZMOD.LT.EPS) IDXK(K)=NK
235                  CONTINUE
                  ELSEIF (NSYM.GT.16) THEN
                     WRITE (1,1500)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      RETURN
      END
C===================================================================
C
C Subroutine MSET calculates the information for the truncation of the
C Q matrix.
C
C Input Parameter;
C ================
C
C LSMAX                   =    MAXIMUM L VALUE FOR THE TENSOR Q
C LLCUT                   =    LARGEST VALUE OF L+L' TO BE INCLUDED
C IPCUT                   =    DIMENSION OF MICUT, MJCUT
C Q(MICUT(K),MJCUT(K))    =    THE KTH INDEPENDENT ELEMENT OF Q
C ICUT                    =    NUMBER OF INDEPENDENT ELEMENTS OF Q
C
C =========================================================================
C
      SUBROUTINE MSET2(LSMAX,LLCUT,ICUT,MICUT,MJCUT,IPCUT,NDT)
C
      DIMENSION MICUT(IPCUT),MJCUT(IPCUT)
C    
100   FORMAT (' NUMBER OF ELEMENTS IN Q IS ',I3)
C
      ICUT=0
      I=0
      DO 101 L=0,LSMAX
         DO 102 M=-L,L
            I=I+1
            IP=0
            DO 103 LP=0,LSMAX
               DO 10 MP=-LP,LP
                  IP=IP+1
                  IF ((L+LP).LE.LLCUT) THEN
                    IF (IP.GE.I) THEN
                      IF((M.EQ.0.AND.MP.EQ.0.).OR.NDT.EQ.3) THEN
                        ICUT=ICUT+1
                        MICUT(ICUT)=I
                        MJCUT(ICUT)=IP
                      ENDIF
                    ENDIF
                  ENDIF
10             CONTINUE
103         CONTINUE
102      CONTINUE
101   CONTINUE
      WRITE (1,100) ICUT
      RETURN
      END
C  file LEEDSATL.SB2  Feb. 29, 1996
C
C**************************************************************************
C  Symmetrized Automated Tensor LEED (SATLEED):  subroutines, part 2
C  Version 4.1 of Automated Tensor LEED
C
C =========================================================================
C
C  Subroutine MTNVTSYM computes reflection and transmission matrices for
C  an atomic layer consisting of NLAY subplanes. A Beeby-type matrix
C  inversion is used.
C
C Parameter List;
C ===============
C
C RA1,TA1,....      =     OUTPUT REFLECTION AND TRANSMISSION MATRICES
C                         (LETTERS R AND T STAND FOR REFLECTION AND 
C                       TRANSMISSION, RESP., NUMBERS 1 AND 2 REFER TO 
C                       INCIDENCE TOWARDS +X AND -X, RESP.).
C N                   =   NO. OF BEAMS IN CURRENT BEAM SET.
C NM,NP               =   DIMENSIONS OF MATRICES (.GE.N)
C AMULT,CYLM          =   WORKING SPACE.
C PQ                  =   LIST OF RECIPROCAL LATTICE VECTORS G (BEAMS).
C NT                  =   TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
C FLMS                =   LATTICE SUMS FROM SUBROUTINE FMAT.
C NL                  =   NO. OF SUBLATTICES CONSIDERED IN SLIND AND FMAT.
C KLM                 =   (2*LMAX+1)*(2*LMAX+2)/2
C LX,LXI,LT,LXM       =   PERMUTATIONS OF (L,M) SEQUENCE, FROM SUBROUTINE
C                             LXGENT.
C LMMAX               =   (LMAX+1)**2.
C KLM                 =   (2*LMAX+1)*(2*LMAX+2)/2.
C XEV,TAU,TAUG,TAUGM  =   WORKING SPACE.
C LEV                 =   (LMAX+1)*(LMAX+2)/2.
C LEV2                =   2*LEV.
C LMT                 =   NTAU*LMMAX.
C LTAUG               =   (NTAU+NINV)*LMMAX.
C CLM                 =   CLEBSCH-GORDAN COEFFICIENTS, FROM SUBROUTINE CELMG.
C NLM                 =   DIMENSION OF CLM (SEE TLEED1).
C POS                 =   INPUT ATOMIC POSITIONS IN UNIT CELL (ONE ATOM PER 
C                        SUBPLANE).
C POSS,MGH,DRL        =   WORKING SPACE
C NUGH,TEST           =   WORKIG SPACE
C GH,RG,TS,TG,VT      =   WORKING SPACES
C NLAY                =   NO. OF SUBPLANES IN LAYER.
C NLAY2               =   NLAY*(NLAY-1)/2.
C LM2N                =   2*LMNI.
C CAA                 =   CLEBSCH-GORDAN COEFFICIENTS, FROM SUBROUTINE CAAA.
C NCAA                =   DIMENSION OF CAA (SEE TLEED1).
C TH                  =   WORKING SPACE.
C LMNI                =   NLAY*LMMAX.
C LMAX                =   LARGEST VALUE OF L.
C LPS                 =   CHEMICAL ELEMENT ASSIGNMENT FOR EACH SUBPLANE 
C                       LPS(I)=J MEANS LAYER NO. I USES ATOMIC T-MATRIX 
C                         ELEMENTS STORED IN TSF(J,K),K=1,LMAX+1.
C LPSS                =   WORKING SPACE.
C TSF                 =   ATOMIC T-MATRIX ELEMENTS
C NEL                 =    NUMBER OF CHEMICAL ELEMENTS IN COMPOSITE LAYER
C NST1                =  number of composite layers
C NCL                 =  current CL index
C NLYMX             =  Maximum LMN over the different CL (dimensioning TSTORE)
C CTR                  phase factors (+1 or -1) to be used in TFLSYM
C CTT
C
C
C In Common Blocks;
C =================
C
C E,VPI               =   CURRENT COMPLEX ENERGY.
C AK2,AK3             =   PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C NA                  =   OFFSET OF CURRENT BEAM SET IN LIST PQ.
C NS                  =   OFFSET FOR FIRST INDEX OF MATRIX ELEMENTS IN RA1,...
C                           (NORMALLY 0).
C LAY                 =   1 IF CURRENT CALCULATION REFERS TO OVERLAYER
C                        2 IF CURRENT CALCULATION REFERS TO SUBSTRATE
C LM                  =   LMAX+1.
C NTAU                =   NO. OF CHEMICAL ELEMENTS TO BE USED.
C TST                 =   CUTOFF PARAMETER FOR RECIPROCAL LATTICE SUM IN 
C                       SUBROUTINE GHMAT
C TV                  =   AREA OF UNIT CELL.
C DCUT                =   CUTOFF RADIUS FOR DIRECT LATTICE SUM IN GHD.
C NOPT                =   1  COMPUTE ALL OUTPUT MATRICES.
C                     =   2  LAYER IS SYMMETRICAL IN +-X  COMPUTE ONLY OUTPUT 
C                       MATRICES WITH NUMBER 1 IN THEIR NAME. BY SYMMETRY 
C                       RA2=RA1,TA2=TA1, ETC. (RA2, TA2, ETC. ARE NOT 
C                       WRITTEN ONTO).
C NOPT                =   3  COMPUTE ONLY A REFLECTION VECTOR (NOT MATRIX) FOR 
C                        (00)-BEAM INCIDENCE (FROM BOTH SIDES +-X).
C NEW                 =   -1  INPUT FOR SKIPPING COMPUTATION OF THOSE 
C                       QUANTITIES THAT HAVE NOT CHANGED SINCE THE PREVIOUS 
C                         CALL TO MTINVT (I.E. THE PREVIOUS CALL TO MTINVT WAS 
C                         FOR THE SAME ENERGY, BUT A DIFFERENT GEOMETRY, BUT 
C                         ALWAYS USE NEW=+1 IN THE CASE OF MULTIPLE CALLS TO 
C                         MTINVT FOR DIFFERENT BEAM SUBSETS, AS FOR SUBSTRATE
C                         LAYERS WHEN A SUPER-LATTICE SURFACE IS PRESENT).
C                         (NOTE- IF USING NEW=-1, DO NOT OVERWRITE CYLM,TAU,
C                        AMULT,TH OR GH IN MAIN PROGRAM, SINCE THEIR OLD 
C                         VALUES WILL BE REUSED.
C NEW                 =   +1  NO REUSE OF OLD VALUES OF CYLM,TAU,AMULT,TH OR 
C                         GH WILL OCCUR
C
C OTHER COMMON BLOCKS;
C ====================
C
C /SL/  DATA PASSED TO SUBROUTINE GHMAT (SEE SUBROUTINE SLIND).
C /MFB/ PASSES DATA TO SUBROUTINE MFOLT.
C
C Modified version of routine MTINV from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE MTSYM6(RA1,TA1,RA2,TA2,N,NM,NP,AMULT,CYLM,PQ,NT,FLMS,
     & NL,LXI,LT,LXM,LX2,LT2,LXM2,LMMAX,KLM,XEV,LEV,LEV2,TAU,LMT,
     & LMT2,TAUG,TAUG2,TAUGM,TAUGM2,CLM,NLM,POS,POSS,MGH,NLAY,
     & DRL,NUGH,LEVV,NLAY2,TEST,RG,TS,
     & TS2,TG,LM2N,VT,CAA,NCAA,TH,LMNI,LMN,IPL,TSTORE,NST1,NSTEF,NST2,
     & NCL,NLYMX,NTAUSH,CTR,CTT,CYTR,CYTT,VVST,VT2,LSM,
     & IND1,IND2,NRCP,NINEQ,NLMAX,LL1,
     & WBDS,IDXS,IDXN,ZRED,LMNBD,WBDS2,IDXS2,IDXN2,LMNBD2,
     & TAUINV,NLTU,LAN,KOUNT,IKRED,IKBS,LMNMX,LMSMX,
     & NBLK,INDK3,PQ3,PH3,NK3,KSNB,NB,NSU,PH4,VEC,LPS,LPSS,TSF)
C
      COMPLEX XEV,RA1,TA1,RA2,TA2,CYLM,AMULT,XA,YA,ST,CF,CT,CI,RU
      COMPLEX XB,YB,AM,AM2,FLMS,CTJ,CZ,AK
      COMPLEX CTR(LMMAX),CTT(LMMAX),CYTR(LMN,NT),CYTT(LMN,NT)
      COMPLEX CSQRT,CEXP,TSTORE(2,NLYMX,NT,NSTEF),VVST(LMN)
      COMPLEX TH(LMN,LMN),VT2(LSM)
      COMPLEX TAU(LMT,LEV),TSF(6,16),RG(4,NLAY,N)
      COMPLEX TS(LMN),TG(2,LM2N),VT(LM2N),TAUG(LMT),TAUGM(LMT)
      COMPLEX TAUG2(LMT2),TAUGM2(LMT2),TS2(NLYMX)
      DIMENSION RA1(NM,NP,NST1),TA1(NM,NP,NST1),RA2(NM,NP,NST1)
      DIMENSION TA2(NM,NP,NST1),AMULT(N)
      DIMENSION XEV(LEV,LEV2),CYLM(NT,LSM),PQ(2,NT)
      DIMENSION CLM(NLM),GP(2),CAA(NCAA)
      DIMENSION FLMS(NL,KLM),LT2(LSM),LXM2(LSM),LX2(LSM)
      DIMENSION LXI(LMMAX),LT(LMMAX),LXM(LMMAX)
      DIMENSION POS(NLMAX,3),POSS(NLMAX,3),LPS(NLMAX),LPSS(NLMAX)
      DIMENSION MGH(NLAY,NLAY),DRL(NLAY2,3)
      DIMENSION NUGH(NLAY2)
      DIMENSION RBR1(2),RBR2(2),TEST(NLAY2),IPL(LMN)
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2)
      DIMENSION IND1(NT),IND2(NT),NB(KSNB)
      DIMENSION INDK3(NK3),PH3(36,NK3),PQ3(2,NK3)
      DIMENSION PH4(12,NT),VEC(2,NST1)
C
C For Symmetry calculations
C a) purely structural

      COMPLEX WBDS(LMNMX),WBDS2(LMSMX)
      DIMENSION IDXN(LMNMX),IDXS(LMNMX),LL1(NLMAX,NST2)
      DIMENSION IDXN2(LMSMX),IDXS2(LMSMX)
      DIMENSION ZRED(LMN),NINEQ(NST2),IKRED(NT),IKBS(NT)
      DIMENSION LMNBD(2),LMNBD2(2)
      LOGICAL  LOGTR,LOGP,SUBNS

      COMPLEX TAUINV(NLTU,LEV)
      DIMENSION LAN(LMMAX),KOUNT(LMMAX,LMMAX)
C
!1610  FORMAT (' NO SYMMETRY in MTNVSYM for GLIDE SYMMETRY. Lower
!     & the symmetry of NSYMS (for the TLEED1.f run) to exclude the
!     & glide ')
      COMMON E,AK2,AK3,VPI
      COMMON /MFB/GP,L1,LAY1
      COMMON /MS/LMAX
      COMMON /MPT/NA,NS,LAY,LLM,NTAU,TST,TV,DCUT,NPERT,NOPT,NEW
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
      AK=-0.5/CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      NA=0  
      L1=LLM
C LM2=LMAX+JSMAX+1
      LM2=INT(SQRT(FLOAT(LSM-1)+0.00001))
      LAY1=LAY
      L2M=2*LMAX+1
      LMSS=L2M*L2M
      LOD=LMMAX-LEV
      LODD=LSM-LEVV
      LEEE=(LM2/2+1)*(LM2/2+1)
      LEE=(LMAX/2+1)*(LMAX/2+1)
      LOE=((LMAX-1)/2+1)*((LMAX-1)/2+1)
      LOEE=((LM2-1)/2+1)*((LM2-1)/2+1)
      CI=CMPLX(0.0,1.0)
      CZ=CMPLX(0.0,0.0)
      RU=CMPLX(1.0,0.0)
      PI=4.0*ATAN(1.0)
      PI4=16.0*ATAN(1.0)
      PI3=PI/3.
      NCLC=NCL
C loop parameters used only in substrate calculation
      KSNB2=1
      KBL1=1
      KBL2=1
      KSBL=1
      SUBNS=.FALSE.
C SUBNS is true if we are doing a substrate layer calculation and
C overlayer and substrate vectors are identical
      IF(NSU.EQ.0.AND.NCL.GT.NSTEF) THEN
         SUBNS=.TRUE.
         KSNB2=KSNB
         KBL1=NCL-NSTEF
         KBL2=NST1
         KSBL=NBLK
      ENDIF
C
C Time reversal invariance? If yes set logical variable to reduce
C # of equations to be solved. (Not used for substrate calculation)
C
      LOGTR=.TRUE.
C      NDTR=2
      IF(NRCP.EQ.NT.OR.SUBNS) THEN 
         LOGTR=.FALSE.
         NRCP=NT
C         NDTR=1
      ENDIF
C
C inversion symmetry of the composite layer is not yet implemented
C
      LOGP=.FALSE.
      NDP=1
      NNN1=0
      NNN2=0
C
C define CTR and CTU
C
      JS1=1
      DO 87 I=1,LLM
         CTJ=-RU
         LL=I+I-1
         DO 88 JM=1,LL
C
C given L and M CTR=(-1)^(L+M)
C
            CTR(JS1)=-CTJ
            CTT(JS1)=RU
            CTJ=-CTJ
            JS1=JS1+1
88       CONTINUE
87       CONTINUE
      NLL=NL
      IF (LAY.EQ.1) NLL=1
C
C  SORT SUBPLANES ACCORDING TO INCREASING POSITION ALONG +X AXIS.
C  FIGURE OUT WHICH OLD GH^S CAN BE KEPT, WHICH MUST BE RECOMPUTED
C  AND WHICH ARE MUTUALLY IDENTICAL
C
      IF(NLAY.GT.1) THEN
          CALL SRLAY2(POS,POSS,LPS,LPSS,MGH,NLAY,DRL,NLAY2,NUGH,
     &    DCUT,NLAY,NINEQ(NCL),LL1(1,NCL),NLMAX)
      ELSE
          DO 89 I=1,3
             POSS(1,I)=POS(1,I)
89        CONTINUE
          LPSS(1)=LPS(1)
      ENDIF
      YA=CMPLX(2.0*E,-2.0*VPI+0.000001)
      YA=CSQRT(YA)
C
      DO 90 IG=1,N
         JG=IG+NA
         BK2=PQ(1,JG)+AK2
         BK3=PQ(2,JG)+AK3
         C=BK2*BK2+BK3*BK3
         XA=CMPLX(2.0*E-C,-2.0*VPI+0.000001)
         XA=CSQRT(XA)
         RG(3,1,IG)=CEXP(CI*XA*(POSS(NLAY,1)-POSS(1,1)))
         DO 60 I=1,NLAY
C
C  GENERATE PLANE-WAVE PROPAGATORS. 
C 
           IF(SUBNS) THEN 
              RG(1,I,IG)=CEXP(CI*XA*(POSS(I,1)-POSS(1,1)))
              RG(2,I,IG)=CEXP(CI*XA*(POSS(I,1)-POSS(1,1)))
              RG(4,I,IG)=CEXP(CI*XA*(POSS(NLAY,1)-POSS(I,1)))
           ELSE
            X=BK2*POSS(I,2)+BK3*POSS(I,3)
            RG(1,I,IG)=CEXP(CI*(XA*(POSS(I,1)-POSS(1,1))+X))
            RG(2,I,IG)=CEXP(CI*(XA*(POSS(I,1)-POSS(1,1))-X))
            RG(4,I,IG)=CEXP(CI*(XA*(POSS(NLAY,1)-POSS(I,1))-X))
           ENDIF
60       CONTINUE
C
C  IF NEW=1 (I.E. NEW ENERGY) COMPUTE NEW SPHERICAL HARMONICS FOR THE
C  BEAM DIRECTIONS IN SPHRM
C
         IF (NEW.NE.-1) THEN
            B=0.0
            CF=RU
            IF (C.GT.1.0E-7) THEN
               B=SQRT(C)
               CF=CMPLX(BK2/B,BK3/B)
            ENDIF
            CT=XA/YA
            ST=B/YA
C
C  GENERATE PREFACTOR OF REFLECTION AND TRANSMISSION MATRIX ELEMENTS
C 
            AMULT(IG)=-16.0*PI*PI*CI/(TV*XA)
C            CALL SPHRM(LMAX,VT,LMMAX,CT,ST,CF)
            CALL SPHRM(LM2,VT2,LSM,CT,ST,CF)
C
C  STORE THE SPHERICAL HARMONICS
C
C            DO 85 K=1,LMMAX
            DO 85 K=1,LSM
               CYLM(IG,K)=VT2(K)
85          CONTINUE
         ENDIF
90    CONTINUE
C
C  start loop over the possibly 2 blocks (because of inversion) that
C  the symmetric rep. decomposes into. 
C  In case of substrate calculation split calculation into beamsets 
C
      DO 901 NBL=1,NDP
         LBD=LMNBD(NBL)
         LBD2=LMNBD2(NBL)
         IF(NCL.GT.NSTEF)LBD2=LBD
      DO 902 NBSET=1,KSNB2
C
C  COMPUTE TAU and TAUINV FOR EACH CHEMICAL ELEMENT
C  NEW could be used if last call generated the needed stuff
C
         CALL TAUMAT3(TAU,LMT,NTAU,XEV,LEV,LEV2,LOD,TSF,
     &   LMMAX,LMAX,FLMS,NL,KLM,LLM,CLM,NLM,LXI,NT,PQ,NA,
     &   NLL,TG,NTAUSH,TAUINV,NLTU,LAN)
C
C  COMPUTE (OR COPY) THE GH MATRICES (INTERLAYER PROPAGATORS IN (LM)-
C  SPACE) and create TH
C
C      I2=MCLOCK()
         CALL GHSYM6(LMMAX,MGH,NLAY,NUGH,NLAY2,
     &   TEST,VT,L2M,VT,LLM,TG,LMSS,DRL,TV,LXM,LEV,DCUT,CAA,NCAA,
     &   LAY,PQ(1,1+NA),LBD,TH,KOUNT,
     &   LPSS,IDXS,IDXN,WBDS,TAUINV,NLTU,LMNI,ZRED,TSF,LAN,
     &   NINEQ(NCL),LL1(1,NCL),NLMAX,NK3)
C      I3=MCLOCK()
C      FTIME=FLOAT(I3-I2)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU
C
C  PREPARE TH FOR INVERSION (GAUSSIAN ELIMINATION)
C
          CALL LUDCMP(TH,LBD,LBD,IPL,D,VVST)
C      I4=MCLOCK()
C      FTIME=FLOAT(I4-I3)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU
C
C loop over different substrate layers to generate distinct sustrate
C matrices
C
      DO 839 ISU=KBL1,KBL2 ,KSBL
         NCP=N
         IF(SUBNS) NCP=NB(NBSET)

C
C  Store combinations used in TFLSYM using symmetry projected states
C
         DO 840 IG2=1,NCP
             IG=IG2+NA
            IF(SUBNS) THEN 
              BK2=PQ(1,IG)+AK2
              BK3=PQ(2,IG)+AK3
              C=BK2*BK2+BK3*BK3
              GKS2=SQRT(C)
            ENDIF
            DO 841 K=1,LBD
               CYTR(K,IG)=CZ
               CYTT(K,IG)=CZ
841         CONTINUE
            DO 84 I=1,NLAY
               II=(I-1)*LMMAX
               IF(SUBNS) THEN 
                  PS2=POSS(I,2)+VEC(1,ISU)
                  PS3=POSS(I,3)+VEC(2,ISU)
                  DRLS=SQRT(PS2**2+PS3**2)
                  DOT=BK2*PS2+BK3*PS3
                  CROS=-BK3*PS2+BK2*PS3
                  IF(GKS2.LE.0.001.OR.DRLS.LE.0.0001) THEN
                     PH0=0.
                  ELSE
                     CCC=DOT/(GKS2*DRLS)
                     IF(CCC.GT.1.) CCC=1.
                     IF(CCC.LT.-1.) CCC=-1.
                     PH0=ACOS(CCC)
                     IF(CROS.LT.0) PH0=2*PI-PH0
                  ENDIF
               ENDIF
C
C  PH0 is the angle between Ripar and Kpar. Use angle between equivalent
C  beams in same beamset to symmetrize over equivalent beams 
C
               DO 86 K=1,LMMAX
                  IP=LXM(K)
C M value 
                  KK=IDXS(II+IP)
                  IF (KK.EQ.0) GOTO 86
                  IF(SUBNS) THEN 
                     L2=INT(SQRT(FLOAT(K-1)+0.00001))
                     M=K-L2-L2*L2-1
                     DO 61 KEQ=1,IKBS(IG) 
                        PHKR=-DRLS*GKS2*COS(PH0-PH4(KEQ,IG))
                        PHLM=FLOAT(M)*PH4(KEQ,IG)
                        CYTR(KK,IG)=CYTR(KK,IG)+CEXP(CI*(PHKR+PHLM))
                        CYTT(KK,IG)=CYTT(KK,IG)+CEXP(CI*(PHKR+PHLM))
61                   CONTINUE
                         CYTR(KK,IG)=CYTR(KK,IG)*CYLM(IG,K)*CTR(K)
     &                  *RG(2,I,IG)/FLOAT(IKRED(IG))
                        CYTT(KK,IG)=CYTT(KK,IG)*CYLM(IG,K)
     &                  *RG(4,I,IG)/FLOAT(IKRED(IG))
                  ELSE
                     CYTR(KK,IG)=CYTR(KK,IG)+CYLM(IG,K)*CTR(K)
     &               *WBDS(II+IP)*RG(2,I,IG)
                     CYTT(KK,IG)=CYTT(KK,IG)+CYLM(IG,K)
     &               *WBDS(II+IP)*RG(4,I,IG)
                  ENDIF
86          CONTINUE
84          CONTINUE
840      CONTINUE
C
C  START LOOP OVER INCIDENT BEAMS
C
C Possible symmetries to take care:
C
C Time reversal: applicable to normal incidence data; relates TA1 to TA2
C                half components of RA1 and RA2 to the other halves
C
C Inversion symmetry: check the coordinates of the composite layer
C                     for this. Relates TA1 to TA2, RA1 to RA2
C                     In particular this could be useful for substrate
C                     matrices.
C
C 2-D point group: symmetry code is input. It relates various
C                  components of each diffraction matrix
C
C First symmetry to be implemented is general for normal incidence: time 
C reversal invariance which leads to reciprocity. 
C     
      NCPTR=1
      NRCP2=NRCP
      IF(SUBNS) NRCP2=NB(NBSET)
4300  CONTINUE
      DO 430 JGPTR2=1,NRCP2
         JGPTR=JGPTR2+NA
         IF(NCPTR.EQ.1) THEN
           JGP=IND1(JGPTR)
           JGP2=IND2(JGPTR)
         ELSE
            JGP=IND2(JGPTR)
            JGP2=IND1(JGPTR)
C skip if calculation already performed
            IF(JGP.EQ.JGP2) GOTO 430
         ENDIF
C
C if Off-normal incidence forget about reciprocity
C
         IF(.NOT.LOGTR) JGP=JGPTR
C
C  IF ONLY DIFFRACTION FROM (00) BEAM REQUIRED, SKIP FURTHER LOOPING
C
         IF (NOPT.EQ.3.AND.JGP.EQ.2) GOTO 531
         GP(1)=PQ(1,JGP)
         GP(2)=PQ(2,JGP)
         C=GP(1)*GP(1)+GP(2)*GP(2)
         GKS2=SQRT(C)
C
C  GENERATE QUANTITIES TAUG,TAUGM INTO WHICH INVERSE OF TH WILL BE
C  MULTIPLIED
C
         CALL TAUT3(TAUG,TAUGM,LMT,LEV,CYLM,NT,LMMAX,LT,
     &    NTAU,LOD,LEE,LOE,JGP,NTAUSH)
C
C  GENERATE QUANTITIES TAUG2,TAUGM2 for tensor extension 
C
         IF(NCL.LE.NSTEF) THEN
           CALL TAUT3(TAUG2,TAUGM2,LMT2,LEVV,CYLM,NT,LSM,LT2,
     &     NTAU,LODD,LEEE,LOEE,JGP,NTAUSH)
         ENDIF
C
C  FIRST CONSIDER INCIDENCE TOWARDS +X
C
         INC=1
         IXC=1
285      CONTINUE
         IF(NCPTR.EQ.1)THEN
           JGA=1
           NCTR=1
         ELSE
           JGA=JGPTR
           NCTR=2
         ENDIF
         IF(.NOT.LOGTR) JGA=1
C
C  INCLUDE APPROPRIATE PLANE-WAVE PROPAGATING FACTORS
C  and consider only the symmetric component.
C  Notice that the momenta in the star (i.e. equivalent to) of a generic JGP 
C  (JGP.NE.1) will contribute in different ways to the symmetric
C  component. However the sum over symmetry related momenta makes the
C  non-symmetric components equal to zero. As a result
C  the solution from LUBKSB corresponds to the average over
C  symmetry related input wavevectors.
C
C  First initialize TS. Then symmetrize by multiplying by phase factors. 
C
         DO 314 I=1,LBD2
           IF(I.LE.LBD) TS(I)=(0.,0.)
           TS2(I)=(0.,0.)
314      CONTINUE
         DO 315 I=1,NLAY
            IF(SUBNS) THEN 
             PS2=POSS(I,2)+VEC(1,ISU)
             PS3=POSS(I,3)+VEC(2,ISU)
             DRLS=SQRT(PS2**2+PS3**2)
             DOT=GP(1)*PS2+GP(2)*PS3
             CROS=-GP(2)*PS2+GP(1)*PS3
             IF(GKS2.LE.0.001.OR.DRLS.LE.0.0001) THEN
               PH0=0.
             ELSE
                 CCC=DOT/(GKS2*DRLS)
                 IF(CCC.GT.1.) CCC=1.
                 IF(CCC.LT.-1.) CCC=-1.
                 PH0=ACOS(CCC)
                 IF(CROS.LT.0) PH0=2*PI-PH0
             ENDIF
            ENDIF
            IN=(I-1)*LMMAX
            IN2=(I-1)*LSM
            LP=(LPSS(I)-1-NTAUSH)*LMMAX
            LP2=(LPSS(I)-1-NTAUSH)*LSM
            IF (INC.EQ.-1) THEN
               ST=RU/RG(2,I,JGP)
            ELSE
               ST=RG(1,I,JGP)
            ENDIF
C
C symmetrization is different for substrate and overlayer because for the
C the substrate we might not use symmetry for the diagonalization of TH
C (For substrate it is convenient to symmetrize on the left (CYTR))
C               DO 310 K=1,LMMAX
               LSM2=LSM
               IF(NCL.GT.NSTEF)LSM2=LMMAX
               DO 310 K=1,LSM2
                  IF(K.LE.LMMAX) K2=IDXS(IN+K)
                  K3=IDXS2(IN2+K)
C                  IF(K2.NE.0) THEN
C
C the substrate calculation does not need a tensor extension
C
                     IF(SUBNS.AND.K2.NE.0) THEN 
                         PHKR=DRLS*GKS2*COS(PH0)
                         IF (INC.EQ.-1) THEN
                           TS(K2)=TS(K2)+TAUGM(LP+K)*ST*CEXP(CI*PHKR)
                         ELSE
                           TS(K2)=TS(K2)+TAUG(LP+K)*ST*CEXP(CI*PHKR)
                         ENDIF
                     ELSEIF(.NOT.SUBNS) THEN
                       IF (INC.EQ.-1) THEN
                          IF(K.LE.LMMAX.AND.K2.NE.0) THEN
                             TS(K2)=TS(K2)+
     &                       TAUGM(LP+K)*ST*CONJG(WBDS(IN+K))
                          ENDIF
                          IF(K3.NE.0) TS2(K3)=TS2(K3)+
     &                    TAUGM2(LP2+K)*ST*CONJG(WBDS2(IN2+K))
                       ELSE
                          IF(K.LE.LMMAX.AND.K2.NE.0) THEN
                            TS(K2)=TS(K2)+
     &                      TAUG(LP+K)*ST*CONJG(WBDS(IN+K))
                          ENDIF
                          IF(K3.NE.0) TS2(K3)=TS2(K3)+
     &                    TAUG2(LP2+K)*ST*CONJG(WBDS2(IN2+K))
                       ENDIF
                     ENDIF
C                  ENDIF
310            CONTINUE
315      CONTINUE
C
C  DO INVERSION WITH MULTIPLICATION
C
         CALL LUBKSB(TH,LBD,LBD,IPL,TS)
C
C extra normalization factor. TS corresponds to the average
C of the unreduced TS's over symmetry related momenta.
C TS*N(JGP) correspond to the sum over symmetry related momenta.
C
         DO 3144 I=1,LBD
            TS(I)=TS(I)*FLOAT(IKRED(JGP))
3144      CONTINUE
         NNN1=NNN1+1
C
C  INCLUDE FURTHER PLANE-WAVE PROPAGATING FACTORS
CCCC 
C We begin with a tensor without symmetry reductions.
C Hence we need to transform back the reduced TS to the unreduced one.
C Note that the input TS vector (the vector Y in the equation
C Y=TH*X that LUBKSB is solving) can always be taken to be fully
C symmetric and this has been used already to arrive at TS. 
C If this is done (by summing over symmetrically related
C momenta) the solution X (output TS in LUB) is also symmetric,
C and the unreduced solution can be obtained very simply 
C
C Take care of symmetry reductions
C a) collects only data from inequivalent layers
C b) collects only nonzero elements
C
         IF(LAY.EQ.1) THEN
C
C collect TSTORE for overlayers calculations
C apart from a phase factor, TSTORE correspond to the amplitude BEFORE
C the last t scattering. (if we multiply by 1-Gt=TAUINV we get the
C amplitude before the last tau scattering)
C
C we also collect an extension of the tensor corresponding to the
C spherical components greater than LMAX which we assume have not
C been scattered. These components are relevant in the center
C expansion
C
          DO 335 I2=1,NINEQ(NCL)
            I=LL1(I2,NCL)
            IN=(I-1)*LMMAX
            IN2=(I-1)*LSM
C            DO 330 K=1,LMMAX
            DO 330 K=1,LSM
                  KNAT=LX2(K)
                  K3=IDXN2(IN2+LXM2(KNAT))
                  IF(KNAT.LE.LMMAX) THEN
                    K2=IDXN(IN+LXM(KNAT))
                  ENDIF
C
C LXM2(KNAT) is the symmetric label corresponding to K in the LSM list. 
C K2=0 iff K3=0
                  IF(K3.EQ.0) GOTO 330
C                  IF(K2.EQ.0) GOTO 330
                  IF (INC.EQ.-1) THEN
                   IF(KNAT.LE.LMMAX) THEN
                      TSTORE(IXC,K3,JGP,NCL)=PI4*TS(K2)*
C                      TSTORE(IXC,K2,JGP,NCL)=PI4*TS(K2)*
     &                RG(3,1,JGP)
                   ELSE
                      TSTORE(IXC,K3,JGP,NCL)=PI4*TS2(K3)*
     &                FLOAT(IKRED(JGP))*RG(3,1,JGP)
                   ENDIF
                  ELSE
                   IF(KNAT.LE.LMMAX) THEN
                     TSTORE(IXC,K3,JGP,NCL)=PI4*TS(K2)
C                     TSTORE(IXC,K2,JGP,NCL)=PI4*TS(K2)
                   ELSE
                      TSTORE(IXC,K3,JGP,NCL)=PI4*TS2(K3)*
     &                FLOAT(IKRED(JGP))
                   ENDIF
                  ENDIF
330         CONTINUE
335       CONTINUE
         ENDIF
C
C  START LOOP OVER SCATTERED BEAMS but first redefine TS as the
C  amplitude AFTER scattering (multiply by t or tau depending on
C  how TS was computed).
      DO 1500 NLR=1,NLAY 
         NII2=(NLR-1)*LMMAX
         DO 1510 I=1,LMMAX 
            K=IDXN(NII2+I)
            IF (K.EQ.0) GOTO 1510
            LL=LAN(I)+1
            IT2=LPSS(NLR)
C
C normalization factor?
C
            TS(K)=TS(K)*AK*TSF(IT2,LL)
C            DO 1550 KP=1,LMNI 
C               TS2(K)=TAU*TS
C1550        CONTINUE
1510     CONTINUE
1500  CONTINUE

C  The vector TS at this point correspond to the symmetric components
C  of a symmetric vector. The projection over symmetry related
C  JG momenta will give the same result. However we only need to collect
C  such result for a single momentum in each equivalence class. 
C
4100    CONTINUE
        DO 410 JGTR2=JGA,NRCP2
            JGTR=JGTR2+NA
            IF(NCTR.EQ.1.OR.NCPTR.EQ.2) THEN
               JG=IND1(JGTR)
               JG2=IND2(JGTR)
            ELSE
               JG=IND2(JGTR)
               JG2=IND1(JGTR)
            ENDIF
            IF(.NOT.LOGTR) JG=JGTR
            IF(LAY.NE.1.AND.SUBNS) NCLC=ISU
            IF(LAY.NE.1.AND..NOT.SUBNS) NCLC=NCL-NSTEF
C
C  COMPLETE COMPUTATION IN TFLSYM INCLUDING REGISTRY SHIFTS 
C
            IF(INC.EQ.1) THEN
             CALL TFLSYM(TS,LBD,
     &       YB,XB,CYTR(1,JG),CYTT(1,JG))
             NNN2=NNN2+1
            ELSE
             CALL TFLSYM(TS,LBD,
     &       YB,XB,CYTT(1,JG),CYTR(1,JG))
             YB=YB*RG(3,1,JGP)
             XB=XB*RG(3,1,JGP)
             NNN2=NNN2+1
            ENDIF
            AM=AMULT(JG)
            AM2=AMULT(JGP)*FLOAT(IKRED(JG))/FLOAT(IKRED(JGP))
            IF (INC.EQ.-1) THEN
               RA2(JG,JGP-NA,NCLC)=YB*AM
               TA2(JG,JGP-NA,NCLC)=XB*AM
C with reciprocity
               IF(LOGTR) THEN
                 RA2(JGP2,JG2-NA,NCLC)=YB*AM2
                 TA1(JGP2,JG2-NA,NCLC)=XB*AM2
               ENDIF
            ELSE
               RA1(JG,JGP-NA,NCLC)=YB*AM
               TA1(JG,JGP-NA,NCLC)=XB*AM
C with reciprocity
               IF(LOGTR) THEN
                 RA1(JGP2,JG2-NA,NCLC)=YB*AM2
                 TA2(JGP2,JG2-NA,NCLC)=XB*AM2
               ENDIF
             ENDIF
410      CONTINUE
         IF (INC.NE.-1) THEN
             IF (NOPT.EQ.2) GOTO 430
C
C The beams have been divided into two groups + and - whose label
C in the large list of beams is given by IND1 and IND2 
C +,+ has allowed the computation of -,- using reciprocity
C Now we compute +,-   (INC=1)
C
            IF(NCTR.EQ.1.AND.LOGTR) THEN
               NCTR=2
               JGA=JGPTR
C               IF(JGPTR.EQ.1) JGA=2
               GOTO 4100
            ENDIF
            INC=-1
            IXC=2
            GOTO 285
         ENDIF
C
C Now we compute +,-   (INC=-1)
C
            IF(NCTR.EQ.1.AND.LOGTR) THEN
               NCTR=2
               JGA=JGPTR
C               IF(JGPTR.EQ.1) JGA=2
               GOTO 4100
            ENDIF
430   CONTINUE
C
C Finally we need -,+
C
         IF(NCPTR.EQ.1.AND.LOGTR) THEN
            NCPTR=2
            GOTO 4300
         ENDIF
C
C Add unit matrix
C
839   CONTINUE
      NA=NA+NCP  
902   CONTINUE
901   CONTINUE
      JBM=0
      DO 439 K=1,KSNB
        DO 440 JG1=1,NB(K) 
          JG=JG1+JBM
          DO 438 ISU=KBL1,KBL2 ,KSBL
            IF(LAY.EQ.1) THEN
               JG2=JG
               NCLC=NCL
            ELSEIF(SUBNS) THEN
               JG2=JG1
               NCLC=ISU
            ELSE
               JG2=JG1
               NCLC=NCL-NSTEF
            ENDIF
            TA1(JG,JG2,NCLC)=TA1(JG,JG2,NCLC)+RG(3,1,JG)
            TA2(JG,JG2,NCLC)=TA2(JG,JG2,NCLC)+RG(3,1,JG)
438       CONTINUE
440     CONTINUE
        JBM=JBM+NB(K)
439   CONTINUE
C      WRITE(*,*) NRCP
C      WRITE(*,*) NNN1,NNN2
531   RETURN
      END
C=======================================================================
C
C  Subroutine PRPGAT performs one propagation through one layer for RFS
C  subroutines like RFSO2 and RFSO3. A layer asymmetrical in +-X can be
C  handled, as well as independent beam sets.
C
C Parameter List;
C ===============
C
C  RA,TA,RB,TB    =  Input diffraction matrices for current layer (R for
C                    reflection, T for transmission, A,B for incidence towards
C                    +-X, resp.).
C  N              =  Total No. of beams at current energy.
C  NM             =  Largest No. of beams in any current beam set (But =N for 
C                    overlayer)
C  AW             =  Working space.
C  I              =  Index of interlayer spacing to which current call to 
C                   PRPGAT will lead.
C  L1,L2          =  Current choice of plane-wave propagators, referring to 
C                   second index of matrix PK. L1 is to describe propagation 
C                   in direction towards interlayer spacing I, L2 in opposite 
C                   direction.
C  CRIT           =  Criterion for penetration convergence (Set in calling 
C                   routine).
C  IR             =  Output flag for penetration convergence.
C  IA             =  +-1 indicates propagation towards +-X.
C  BNORM          =  Output measure of current wavefield amplitude at current 
C                   layer.
C
C Modified version of routine RFSO2 from the VAN HOVE/TONG LEED package.
C Modifications by ROUS and WANDER.
C Modifications by BARBIERI
C
C =========================================================================
C
      SUBROUTINE PRPGAT2(RA,TA,RB,TB,N,NM,ANEW,ND,NSL,NL,AW,I,PK,L1,L2,
     & CRIT,IR,IA,BNORM,NROM,IKRED)
C
      COMPLEX RA(NROM,NM),TA(NROM,NM),ANEW(N,ND),AW(N,2),PK(N,8),CZ
      COMPLEX RB(NROM,NM),TB(NROM,NM)
      DIMENSION NSL(NL),IKRED(N)
C
6     FORMAT ('***THIS ORDER TOO DEEP')
C
      CZ=CMPLX(0.0,0.0)
      BNORM=0.0
      NA=0
      DO 4 NN=1,NL
         NB=NSL(NN)
         DO 1 K=1,NB
            KNA=K+NA
C
C  PROPAGATE WAVEFIELD TO CURRENT LAYER FROM THE TWO NEAREST LAYERS
C
            AW(K,1)=PK(KNA,L1)*ANEW(KNA,I-IA)
            AW(K,2)=PK(KNA,L2)*ANEW(KNA,I)
1        CONTINUE
         DO 10 J=1,NB
            JNA=J+NA
            ANEW(JNA,I)=CZ
C
C  SELECT FORMULA ACCORDING TO PROPAGATION DIRECTION
C
            IF (IA.EQ.-1) THEN
               DO 9 K=1,NB
                  ANEW(JNA,I)=ANEW(JNA,I)+TB(JNA,K)*AW(K,1)+RA(JNA,K)
     &             *AW(K,2)
9              CONTINUE
            ELSE
               DO 2 K=1,NB
C
C  TRANSMIT AND REFLECT AND ADD WAVES
C
                  ANEW(JNA,I)=ANEW(JNA,I)+TA(JNA,K)*AW(K,1)+RB(JNA,K)
     &             *AW(K,2)
2              CONTINUE
            ENDIF
C
C  GET MEASURE OF NEW WAVEFIELD MAGNITUDE
C
            FACT1=REAL(ANEW(JNA,I))*REAL(ANEW(JNA,I))
            FACT2=AIMAG(ANEW(JNA,I))*AIMAG(ANEW(JNA,I))
            BNORM=BNORM+(FACT1+FACT2)*FLOAT(IKRED(JNA))
10       CONTINUE
         NA=NA+NB
4     CONTINUE
C
C  CHECK AGAINST PENETRATION LIMIT, EXCEPT WHEN EMERGING (IN WHICH CASE
C  IR=1)
C
      IF (.NOT.((I.LT.ND-1).OR.(IR.EQ.1))) THEN
         WRITE (1,6)
         IR=1
      ENDIF
C
C  CHANGE I TO INDICATE TO WHICH NEW INTERLAYER SPACING THE NEXT CALL
C  TO PRPGAT WILL HAVE TO LEAD
C
      I=I+IA
C
C  CHECK ON PENETRATION CONVERGENCE
C
      IF (BNORM.LE.CRIT) IR=1
      RETURN
      END
C==================================================================== 
C                                                                     
C Routine QGEN calculates the Q matrix for the current layer and Energy.
C 
C Input Parameters;
C =================
C
C GA,GB              = G VECTORS FROM CVEC
C LSMMAX             = (LSMAX+1)**2
C LMMAX              = (LMAX+1)**2
C NLLAY              = INDEX OF CURRENT LAYER
C NLAY               = NUMBER OF SUBPLANES IN COMPOSITE LAYER
C E,VPI              = CURRENT (COMPLEX) ENERGY
C NEXIT              = INDEX OF CURRENT BEAM
C NT0                = TOTAL NUMBER OF EXIT BEAMS
C AK2M,AK3M          = PARALLEL COMPONENTS OF MOMENTUM OF EACH EXIT BEAM
C TV                 = AREA OF UNRECONSTRUCTED UNIT CELL
C Q                  = THE Q MATRIX FOR THE CURRENT LAYER AND EXIT BEAM
C                     DIRECTION AT THE CURRENT ENERGY
C
C =========================================================================
C
      SUBROUTINE QGEN2(GA,GB,Q,LMMAX,LSMMAX,NLLAY,NLAY,E,VPI,NEXIT,NT0,
     & AK2M,AK3M,TV,NLMX,LAFLG,NSTEF,NIND,NCL)
C
      COMPLEX GA(NLAY,LMMAX,LSMMAX),GB(LMMAX,LSMMAX)
      COMPLEX Q(LSMMAX,LSMMAX)
      COMPLEX CSUM,CI,PREG,CAK,XA
      DIMENSION AK2M(NT0),AK3M(NT0),LAFLG(NLMX,NSTEF)
C
      CI=CMPLX(0.0,1.0)
      D=0.0
      ZN=FLOAT(LAFLG(NIND,NCL))
      IF (NEXIT.GT.0) D=AK2M(NEXIT)**2+AK3M(NEXIT)**2
      CAK=CMPLX(2.0*E,-2.0*VPI+0.0000001)
      CAK=CSQRT(CAK)
      XA=CMPLX(2.0*E-D,-2.0*VPI+0.0000001)
      XA=CSQRT(XA)
      PREG=ZN/(2.0*TV*CAK*XA)
      DO 100 I2=1,LSMMAX
         DO 110 I3=1,LSMMAX
            CSUM=CMPLX(0.0,0.0)
            DO 120 I1=1,LMMAX
               CSUM=CSUM+GB(I1,I2)*GA(NLLAY,I1,I3)
120         CONTINUE
            Q(I2,I3)=PREG*CSUM
110      CONTINUE
100   CONTINUE
      RETURN
      END
C============================================================================
      SUBROUTINE REDMTR(RA1,TA1,RA1R,TA1R,NT,NP,NSS,NBM, 
     +IDXB,IDXB2,NBAL,KNBS,NB,KSNBS,VEC,NTAL,NPAL,PQAL)
C
C   Project the substrate diffraction matrices by taking into account the
C   symmetry of the beams  
C
      COMPLEX RA1(NTAL,NPAL),TA1(NTAL,NPAL), RA1R(NT,NP,NBM)
      COMPLEX TA1R(NT,NP,NBM),CZ,CI,CPH
      DIMENSION VEC(2,NBM),NB(KSNBS)
      DIMENSION IDXB(NTAL),IDXB2(NSS),NBAL(KNBS),PQAL(2,NTAL)
      LOGICAL TEST
C
C IDXB(i)=IDXB(j), i,j=1,NTAL (NTAL beams altogether) indicates beams i and
C j are equivalent. The reduced matrix is obtained by summing over
C equivalent indices: MR(I,J)=SUM M(I,EQUIV(J))
C
      EPS=.001
      CZ=(0.,0.)
      CI=(0.,1.)
C         DO 11 K=1,NSS
C         DO 12 KK=1,NSS
C
C go through the reduced K's
C
      KSNBJ=0
      DO 11 K1=1,KSNBS
        NN1=NB(K1)
        DO 111 K2=1,NN1
          K=K2+KSNBJ
          JJJ=IDXB2(K)
C
C loop over reduced KK's in the same reduced beamset as K
C
          DO 12 KK=1,NN1
            DO 13 NBK=1,NBM
               RA1R(K,KK,NBK)=CZ
               TA1R(K,KK,NBK)=CZ
13          CONTINUE
C
C Identify beams equivalent to KK and sum over
C
            KNBJ=0
            DO 1 I=1,KNBS
               NN=NBAL(I)
C
C if TEST,  JJ-K is in an integer beam and JJ can be tested with KK for
C equivalence
C
               TEST=JJJ.GT.KNBJ.AND.JJJ.LE.KNBJ+NN
               IF(TEST) THEN
                DO 2 II=1,NN
                  JJ=II+KNBJ
C
C is JJ equivalent to KK?
C
                  IF(IDXB(JJ).EQ.KK+KSNBJ) THEN
C
C Yes! Then add appropriate phase factor for each bulk matrix.
C
                    ZK1=(PQAL(1,JJ)-PQAL(1,JJJ))
                    ZK2=(PQAL(2,JJ)-PQAL(2,JJJ))
C                    ZK1=(0.,0.)
C                    ZK2=(0.,0.)
                    DO 20 NBK=1,NBM
                     CPH=CEXP(CI*(ZK1*VEC(1,NBK)+ZK2*VEC(2,NBK)))
                  RA1R(K,KK,NBK)=RA1R(K,KK,NBK)+RA1(IDXB2(K),II)*CPH
                  TA1R(K,KK,NBK)=TA1R(K,KK,NBK)+TA1(IDXB2(K),II)*CPH
20                  CONTINUE
                  ENDIF
2              CONTINUE
             ENDIF  
             KNBJ=KNBJ+NN
1           CONTINUE
12       CONTINUE
111     CONTINUE
        KSNBJ=KSNBJ+NN1
11    CONTINUE
      RETURN
      END
C====================================================================== 
C                                                                       
C  Subroutine RFSSYM applies to substrates with a 1- or 2- layer periodicity 
C  of the layer diffraction matrices, and a 1- or 2- layer periodicity
C  of the interlayer vectors (two alternating interlayer vectors are used.)
C  It deals with matrices which have been reduced according to symmetry.
C
C  LAYER NO.    I   INTERL.VECTOR   INTERL.PROPAG\S   LAYER DIFFR.MATR\S
C
C  0(SURF.)   ---------------------------------------
C               1        ASE           1,2
C  a(OVERLAY.)--------------------------------------- ROP,TOP,ROM,TOM(1)
C               2a       ASB(1)        3,4    1,2 
C  b(OVERLAY.)--------------------------------------- ROP,TOP,ROM,TOM(1)
C               2b       ASB(2)               3,4
C
C                  ...............................
C
C  2(SUBSTR.) --------------------------------------- RA1,TA1
C               3        AS1              5,6
C  3          --------------------------------------- RA2,TA2
C               4        AS2              7,8
C  4          --------------------------------------- RA1,TA1
C               5        AS1              5,6
C  5          --------------------------------------- RA2,TA2
C               6        AS2              7,8
C  ETC.
C
C Parameter List;
C ===============
C
C ROP,TOP,ROM,TOM        =  OVERLAYER DIFFRACTION MATRICES (R FOR REFLECTION,
C                           T FOR TRANSMISSION, P FOR INCIDENCE TOWARDS +X, 
C                         M FOR INCIDENCE TOWARDS -X).
C                         These matrices have an index for the composite 
C                         layer. If NBULK.NE.0 They are used also to make 
C                         the bulk. ROP etc corresponding to the CL NST1,NST1-1,
C                         .... NST1-NBULK+1 (NST1 is the last computed and 
C                         deepest CL) are used to make the bulk.
C NBULK                  = 0 if the bulk is built out of diffraction matrices
C                          computed in MSMFT (RA and TA), .NE.0 if composite
C                          layer diffraction matrices are used instead.
C NROM                   =  Dimension of the above matrices
C RA1R,TA1R,RA2R,TA2R    =  SUBSTRATE LAYER DIFFRACTION MATRICES (R FOR 
C                         REFLECTION, T FOR TRANSMISSION, 1 FOR LAYERS 2,4,
C                               6,.., 2 FOR LAYERS 3,5,7,..)
C                           reduced with symmetry
C N                      =  NO. OF BEAMS USED AT CURRENT ENERGY (INCL. ALL BEAM 
C                         SETS).
C NSL                    =  NO. OF BEAMS IN EACH SET OF BEAMS.
C NL                     =  NO. OF BEAM SETS.
C NM                     =  2ND DIMENSION OF SUBSTRATE DIFFRACTION MATRICES 
C                         (= LARGEST VALUE OF NSL).
C WV                     =  OUTPUT REFLECTED BEAM AMPLITUDES.
C PQ                     =  LIST OF BEAMS.
C PKCL                   =  Interlayer propagator (between Composite layers)
C PK                    =  Interlayer propagator (within substrate)
C VICL                   =  Imaginary part of inner potential in the CL
C VCL                    =  Real part of inner potential in the C
C                           (Currently implemented for VCL(i)=0.)
C FRCL                   =  FRACTION OF SPACING BETWEEN composite layers NCL
C                           and NCL+1 which is allotted to CL NCL.
C                           (NST1+1 is the substrate)
C,AW,ANEW                =  WORKING SPACES.
C ND                     =  NO.OF THE LAYER TO WHICH PENETRATION INTO SURFACE 
C                         IS ALLOWED (DIMENSIONS ANEW).
C ASB                    =  INTERLAYER VECTOR BETWEEN OVERLAYER AND TOP 
C                         SUBSTRATE LAYER.
C ASA                    =  SUBSTRATE INTERLAYER VECTOR.
C AMPPLW                 =  PLANE WAVE AMPLITUDES INCIDENT ON EITHER SIDE
C                           OF THE COMPOSITE LAYERS
C NBIN                       =  LABEL OF CURRENT INCIDENT BEAM
C 
C In Common Blocks;
C =================
C
C ASL                    =  NOT USED.
C ASE                    =  SPACING BETWEEN SURFACE AND OVERLAYER NUCLEI.
C VPIS,VPIO              =  OPTICAL POTENTIAL IN SUBSTRATE AND OVERLAYER, RESP.
C VO,VV                  =  LEVEL OF CONSTANT POTENTIAL IN OVERLAYER AND 
C                         VACUUM, RESP., REFERRED TO LEVEL OF CONSTANT 
C                         POTENTIAL IN SUBSTRATE.
C E                      =  CURRENT ENERGY.
C VPI                    =  NOT USED.
C BK2,BK3                =  PARALLEL COMPONENTS OF INCIDENT K-VECTOR.
C AS                     =  NOT USED.
C
C Routine RFSO2T is a modified version of the routine RFSO2 from the
C VAN HOVE/TONG LEED package. Modifications by ROUS and WANDER.
C Modifications to include multiple composite layers by BARBIERI
C
C =========================================================================
C
      SUBROUTINE RFSSYM(ROP,TOP,ROM,TOM,RA1,TA1,RA2,TA2,N,NSL,NL,NM,
     &WV,PQ,PK,AW,ANEW,ND,ASB,AS1,AS2,IPR,AMPPLW,NBIN,NST1,NST2,NBULK,
     &VICL,VCL,PKCL,FRCL,NROM,NSTB,INVECT,IKRED,NEXIT)
C
      COMPLEX RA1,TA1,RA2,TA2,WV,ANEW,AW,PK,EK1,EK2,PKCL
      COMPLEX ROP,TOP,ROM,TOM,EL,EK,EKP,CZ,CI
      COMPLEX CSQRT,CEXP,AMPPLW(N,2,NST1)
      DIMENSION RA1(N,NM,NSTB),TA1(N,NM,NSTB)
      DIMENSION RA2(N,NM,NSTB),TA2(N,NM,NSTB)
      DIMENSION ROP(NROM,NROM,NST1),TOP(NROM,NROM,NST1)
      DIMENSION ROM(NROM,NROM,NST1),TOM(NROM,NROM,NST1)
      DIMENSION NSL(NL),WV(N),MSL(1),ASB(NST2,3),AS1(3),AS2(3)
      DIMENSION PQ(2,N),ANEW(N,ND),AW(N,2),PK(N,8)
      DIMENSION VICL(NST2),VCL(NST2),PKCL(N,8,NST2)
      DIMENSION FRCL(NST2),IKRED(N)
C      COMPLEX RA1,TA1,WV,ANEW,AW,PK,EK1,EK2,PKCL
C
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /X4/E,VPI,BK2,BK3
C
89    FORMAT (/,' RFS CONVERGENCE TEST ')
 32   FORMAT (' IMAX = ',I2,5X,'BNORM = ',E12.4,' ANORM = ',E12.4)
38    FORMAT (/,' ***NO CONVERGENCE AFTER',1I4,' ITERATIONS')
C
      IF (IPR.GT.0) WRITE (1,89)
C
C Note that the convergence criterium is slightly different depending on
C the symmetry, because in PRPGAT  BNORM is obtained by summing over
C inequivalent beams only.
C
      CRIT=0.002
      CZ=CMPLX(0.0,0.0)
      CI=CMPLX(0.0,1.0)
      DO 1 J=1,N
         AK2=BK2+PQ(1,J)
         AK3=BK3+PQ(2,J)
         A=2.0*E-AK2*AK2-AK3*AK3
C
C compute propagators from the surface to the first composite layer
C
         EL=CSQRT(CMPLX(A-2.0*VCL(1),-2.0*VICL(1)+0.000001))
         PK(J,1)=CEXP(CI*EL*ASE)
         PK(J,2)=PK(J,1)
C
C Now from one composite layer to the next
C
         DO 5 NCL=1,NST2
           EL=CSQRT(CMPLX(A-2.0*VCL(NCL),-2.0*VICL(NCL)+0.000001))
           IF (NCL.EQ.1) THEN
             PKCL(J,1,NCL)=CEXP(CI*EL*ASE)
             PKCL(J,2,NCL)=PKCL(J,1,NCL)
           ELSE
             PKCL(J,1,NCL)=PKCL(J,3,NCL-1)
             PKCL(J,2,NCL)=PKCL(J,4,NCL-1)
           ENDIF
           M1=NCL
           M2=NCL+1
           IF (NCL.LT.NST2) THEN
             EK1=CSQRT(CMPLX(A-2.0*VCL(M1),-2.0*VICL(M1)+0.000001))
             EK2=CSQRT(CMPLX(A-2.0*VCL(M2),-2.0*VICL(M2)+0.000001))
           ELSE
             M1=NCL
             EK1=CSQRT(CMPLX(A-2.0*VCL(M1),-2.0*VICL(M1)+0.000001))
             EK2=CSQRT(CMPLX(A,-2.0*VPIS+0.000001))
           ENDIF
           EKP=ASB(NCL,1)*(FRCL(NCL)*EK1+(1.0-FRCL(NCL))*EK2)
           PKCL(J,3,NCL)=CEXP(CI*(EKP+AK2*ASB(NCL,2)+AK3*ASB(NCL,3)))
           PKCL(J,4,NCL)=CEXP(CI*(EKP-AK2*ASB(NCL,2)-AK3*ASB(NCL,3)))
CCC
           IF (NCL.EQ.NST2) THEN
             PKCL(J,3,NCL)=CEXP(CI*EKP)
             PKCL(J,4,NCL)=CEXP(CI*EKP)
           ENDIF
CCC
5        CONTINUE
C
C  Finally for the substrate
C  TWO SETS OF DIFFERENT INTERLAYER PROPAGATORS ARE PRODUCED FOR THE
C  SUBSTRATE, USING THE BULK INTERLAYER VECTORS AS1 AND AS2
C
         EK=CSQRT(CMPLX(A,-2.0*VPIS+0.000001))
         EKP=AK2*AS1(2)+AK3*AS1(3)
         PK(J,5)=CEXP(CI*(EK*AS1(1)))
         PK(J,6)=CEXP(CI*(EK*AS1(1)))
         PK(J,3)=PKCL(J,3,NST2)
         PK(J,4)=PKCL(J,4,NST2)
         EKP=AK2*AS2(2)+AK3*AS2(3)
         PK(J,7)=CEXP(CI*(EK*AS2(1)))
         PK(J,8)=CEXP(CI*(EK*AS2(1)))
1     CONTINUE
      DO 202 I=1,N
         WV(I)=CZ
         DO 2 J=1,ND
             ANEW(I,J)=CZ
2        CONTINUE
202   CONTINUE
      ANEW(NBIN,1)=CMPLX(1.0,0.0)
C
C initialize the wavefunction amplitude
C
      DO 203 I=1,N
         DO 22 J=1,2
           DO 20 NCL=1,NST1
            AMPPLW(I,J,NCL)=CZ
20         CONTINUE
22       CONTINUE
203   CONTINUE
      MSL(1)=N
      IO=0
      ANORM1=1.0E-6
3     ANORM2=0.0
      IO=IO+1
      IR=0
      I=2
C
C Propagate downward
C
      DO 7 NCL=1,NST2
         IF(NCL.LE.NST1) THEN
           CALL PRPGAT2(ROP(1,1,NCL),TOP(1,1,NCL),ROM(1,1,NCL),
     &     TOM(1,1,NCL),N,NROM,ANEW,ND,MSL,1,AW,I,PKCL(1,1,NCL),
     &     1,4,CRIT,IR,1,BNORM,NROM,IKRED)
         ELSE
C           IF (IR.NE.1) THEN 
             NCLC=NCL-NST1
             CALL PRPGAT2(RA1(1,1,NCLC),TA1(1,1,NCLC),RA2(1,1,NCLC),
     &       TA2(1,1,NCLC),N,NM,ANEW,ND,NSL,NL,AW,I,PKCL(1,1,NCL),
     &       1,4,CRIT,IR,1,BNORM,N,IKRED)
C           ENDIF
         ENDIF
7     CONTINUE 
      L1=3
C
C  PROPAGATE THE WAVEFIELD THROUGH THE NEXT TWO LAYERS AND REPEAT UNTIL
C  CONVERGENCE IN PENETRATION
C
      IF (IR.NE.1) THEN
          NBLK=NBULK
50        NBLK=MOD(NBLK+1,NSTB)
          IF(NBLK.EQ.0)NBLK=NSTB
          CALL PRPGAT2(RA1(1,1,NBLK),TA1(1,1,NBLK),RA2(1,1,NBLK),
     &    TA2(1,1,NBLK),N,NM,ANEW,ND,NSL,NL,AW,I,
     &    PK,L1,6,CRIT,IR,1,BNORM,N,IKRED)
          IF (IR.NE.1) THEN
C
C Need extra propagation
C
            L1=7
            IF(INVECT.EQ.2) THEN
C
C with a different propagating factor
C
              NBLK=MOD(NBLK+1,NSTB)
              IF(NBLK.EQ.0)NBLK=NSTB
              CALL PRPGAT2(RA1(1,1,NBLK),TA1(1,1,NBLK),RA2(1,1,NBLK),
     &        TA2(1,1,NBLK),N,NM,ANEW,ND,NSL,NL,
     &        AW,I,PK,5,8,CRIT,IR,1,BNORM,N,IKRED)
              IF (IR.LE.0) GOTO 50
            ELSE
C
C or with the same PK (but possibly different matrix)
C
              GOTO 50
            ENDIF
          ENDIF
      ENDIF
      I=I-1
      IMAX=I-1
      BNORMS=BNORM
      IMAXS=IMAX
C
C SUM UP THE POSITIVE TRAVELLING PLANE WAVES BETWEEN EACH LAYER INTO
C AMPPLW MULTIPLYING BY THE APPROPRIATE PHASE FACTOR (PK) SO THAT THE
C AMPLITUDES ARE EXPRESSED W.R.T AN ORIGIN CENTRED UPON THE ORIGIN ATOM
C OF EACH LAYER.(on the top of the corresponding composite layer)
C
      DO 100 IG=1,N
         DO 103 NCL=1,NST1
           AMPPLW(IG,1,NCL)=AMPPLW(IG,1,NCL)+
     &      ANEW(IG,NCL)*PKCL(IG,1,NCL)
103      CONTINUE
100   CONTINUE
C
C Now propagate upward
C
      DO 21 IG=1,N
         ANEW(IG,I+1)=CZ
21    CONTINUE
      L3=7
      IF (I.EQ.NST2) GOTO 61
      II=MOD(I-NST2+1,2)+1
      NBLK=MOD(NBLK+1,NSTB)
      IF(NBLK.EQ.0)NBLK=NSTB
         IF (II.EQ.1.OR.INVECT.EQ.1) GOTO 58
57       CALL PRPGAT2(RA1(1,1,NBLK),TA1(1,1,NBLK),RA2(1,1,NBLK),
     &    TA2(1,1,NBLK),N,NM,ANEW,ND,NSL,NL,AW,I,PK,
     &    8,5,CRIT,IR,-1,BNORM,N,IKRED)
         NBLK=MOD(NBLK-1,NSTB)
         IF(NBLK.EQ.0)NBLK=NSTB
58       IF (I.EQ.NST2+1) L3=3
         CALL PRPGAT2(RA1(1,1,NBLK),TA1(1,1,NBLK),RA2(1,1,NBLK),
     &    TA2(1,1,NBLK),N,NM,ANEW,ND,NSL,NL,AW,I,PK,
     &    6,L3,CRIT,IR,-1,BNORM,N,IKRED)
         NBLK=MOD(NBLK-1,NSTB)
         IF(NBLK.EQ.0)NBLK=NSTB
         IF (I.GT.NST2) THEN
           IF(INVECT.EQ.2) GOTO 57
           GOTO 58
         ENDIF
61    DO 27 NCL=NST2,1,-1
         IF(NCL.GT.NST1) THEN
           NCLC=NCL-NST1
           CALL PRPGAT2(RA1(1,1,NCLC),TA1(1,1,NCLC),RA2(1,1,NCLC),
     &     TA2(1,1,NCLC),N,NM,ANEW,ND,NSL,NL,AW,I,PKCL(1,1,NCL),
     &     4,1,CRIT,IR,-1,BNORM,N,IKRED)
         ELSE
           CALL PRPGAT2(ROP(1,1,NCL),TOP(1,1,NCL),ROM(1,1,NCL),
     &     TOM(1,1,NCL),N,NROM,ANEW,ND,MSL,1,AW,I,PKCL(1,1,NCL),
     &     4,1,CRIT,IR,-1,BNORM,NROM,IKRED)
         ENDIF
27    CONTINUE
C
C SUM UP THE NEGATIVE TRAVELLING PLANE WAVES BETWEEN EACH LAYER INTO AMPPLW
C MULTIPLIED BY THE APPROPRIATE PHASE FACTOR (PK) SO THAT THE AMPLITUDES ARE
C EXPRESSED WITH RESPECT TO AN ORGIN CENTRED UPON THE ORIGIN OF EACH LAYER.
C (at the bottom of the corresponding composite layer)
      DO 101 IG=1,N
         DO 104 NCL=1,NST1
            AMPPLW(IG,2,NCL)=AMPPLW(IG,2,NCL)+
     &       ANEW(IG,NCL+1)*PKCL(IG,4,NCL)
104      CONTINUE
101   CONTINUE
C
C Check convergence
C
      DO 29 IG=1,N
         EK=ANEW(IG,1)*PKCL(IG,2,1)
         AB=CABS(EK)
         ANORM2=ANORM2+AB*AB
         WV(IG)=WV(IG)+EK
         ANEW(IG,1)=CZ
29    CONTINUE
      ANORM1=ANORM1+ANORM2
      IF (IPR.GT.0) WRITE (1,32) IMAXS,BNORMS,ANORM1
      IF (ANORM2/ANORM1-0.001.LE.0) GOTO 204
      IF (IO.LT.6) GOTO 3
      WRITE (1,38) IO
C
C normalize amplitude in case of time reverse calculation. Right
C now the input field correspond to a symmetric star.
C
204   IF(NEXIT.GT.0) THEN
        DO 302 IG=1,N
           DO 303 NCL=1,NST1
              AMPPLW(IG,1,NCL)=AMPPLW(IG,1,NCL)/FLOAT(IKRED(NBIN))
              AMPPLW(IG,2,NCL)=AMPPLW(IG,2,NCL)/FLOAT(IKRED(NBIN))
303      CONTINUE
302     CONTINUE
      ENDIF
      RETURN
      END
C============================================================================
C
C SUB SYMBLK computes how many symmetric (but different) bulk matrices
C are necessary to form the bulk. The origin of each matrix will
C be the same and equal to the symmetric origin in tleed5.i. 
C The output is:
C NST1B # of bulk matrices
C VEC(I) I=1,NST1B giving the 2D vector displacement from the origin
C of the computed bulk matrix (either the position of an atom if the
C matrix is generated in MSMFT, or not necessarely so if NBULK.NE.0) to
C the symmetric origin, for the NST1B bulk matrices.
C
C Notice that if the symmetry is a single mirror plane, we only need to worry
C about displacements perpendicular to such plane.(not implemented)
C 

      SUBROUTINE SYMBLK (ASA,INVECT,NST1B,
     & INST1B,ARA1,ARA2,NST1,NBULK,VEC)
      DIMENSION ARA1(2),ARA2(2),ASA(10,3),VEC(2,INST1B)
      DIMENSION VECC(2)
C
C information abot the substrate vectors
C
      A1MOD=SQRT(ARA1(1)**2 + ARA1(2)**2)
      A2MOD=SQRT(ARA2(1)**2 + ARA2(2)**2)
      Q1X=ARA1(1)
      Q2X=ARA2(1)
      Q1Y=ARA1(2)
      Q2Y=ARA2(2)
      DET=Q1X*Q2Y - Q2X*Q1Y
      NST1B=NBULK
      NBL=1
C
C NBULK truly independent layer substrate matrices are needed
C The code currently coded to allow only 2 different ASA vectors
C (bulk interlayer vectors). This subroutine is more general in that
C works for INVECT vectors
C 
      DO 20 I=1,NBULK
           VEC(1,I)=0.
           VEC(2,I)=0.
20    CONTINUE
C
C Check for successive layer. VECC is the candidate for the registry
C corresponding to the next substrate layer
C
30         VECC(1)=VEC(1,NST1B)+ASA(NBL,2)
           VECC(2)=VEC(2,NST1B)+ASA(NBL,3)
           IF(NBL.LT.INVECT) THEN
              NBL=NBL+1
           ELSE
              NBL=1
           ENDIF
           VECX=VECC(1)-VEC(1,1)
           VECY=VECC(2)-VEC(2,1)
C
C Check a) whether the registry is the same as the first
C          (i.e. whether VECX,VECY is equivalent to a substrate
C           lattice vector)
C AND   b) whether the bulk matrix involved is the same
C Until these two conditions are satisfied keep increasing
C NST1B and recording VECC into VEC
C For details see LOOKUP 
C
C if there is only one mirror things can be fixed to minimize NST1B
C ( displacements perpendicular to the mirror
C are the only ones which matter; but then the PK's must be fixed in RFS)
C
           VECX2=VECX
           VXFD=(Q2Y*VECX2 - Q2X*VECY)/DET +100.
           VYFD=(-Q1Y*VECX2 + Q1X*VECY)/DET +100.
           VXFD=AMOD(VXFD,1.0)
           VYFD=AMOD(VYFD,1.0)
           IF(ABS(VXFD-1.).LT..01)VXFD=VXFD-1.
           IF(ABS(VYFD-1.).LT..01)VYFD=VYFD-1.
      IF (ABS(VXFD)+ABS(VYFD).GT.0.001.OR.NBL.NE.1) THEN
                 NST1B=NST1B+1
                 IF(NST1B.GT.INST1B) GOTO 50
C
C Store vector
C
                 VEC(1,NST1B)=VECC(1)
                 VEC(2,NST1B)=VECC(2)
                 GOTO 30
      ENDIF
C
C One should have NST1B=NBULK*n where n is an integer
C
      RETURN
50    WRITE (1,*) ' INST1B is too small, or you did not need INVECT=2'
      STOP
      WRITE (1,*) ' Put the mirror plane along the z axis, symmetry
     &  treatment in REDMTR will not be correct otherwise'
      STOP
      END
C======================================================================
C
C  Subroutine TREV select the beams that are equivalent because
C  of time-reversal invariance.  Author: Barbieri
C  Also compute angles (stored in PH4) between the irreducible K used 
C  in MTSYM and other symmetry equivalent K's
C
      SUBROUTINE TREV(IND1,IND2,PQAL,PQ,NTAL,NT,INI,NRCP,
     &     IDXK,IDXK2,IKBS,PH4,KNBS,KNB)
C
      DIMENSION IND1(NT),IND2(NT),PQAL(2,NTAL),PQ(2,NT),IDXK(NTAL)
      DIMENSION IDXK2(NT),IKBS(NT),PH4(12,NT),KNB(KNBS)
C
C determine map connecting +K and -K (will be used for normal incidence
C only) used in MTNVSYM
C The set of all inequivalent K's is divided into two subsets {K} and {-K}
C For this purpose K1 and K2 will be considered equivalent if the symmetry
C class of -K1 is equivalent to the symmetry class of K2  
C The number of truly inequivalent beams is then reduced from NT to NRCP
C J=1,NRCP IND1(J)=JJ where the JJ is the index of the beam in the big list
C IND2(J)=JJ is the index of - the beam corresponding to IND1(J) 
C
C        N1=(NT-1)/2+1
      IF(INI.EQ.1) THEN
            KII=0
            DO 122 II=1,NT
               DO 121 I=1,NTAL
                  DIFF1=ABS(PQ(1,II)+PQAL(1,I))
                  DIFF2=ABS(PQ(2,II)+PQAL(2,I))
                  IF ((DIFF1.LT.1.0E-02).AND.(DIFF2.LT.1.0E-02)) THEN
                     K2=IDXK(I)
                     IF(K2.GE.II) THEN
                        KII=KII+1
                        IND1(KII)=II
                        IND2(KII)=K2
                     ENDIF
                  ENDIF
121            CONTINUE
122         CONTINUE
            NRCP=KII
      ELSE
C
C No extra symmetry to be used
C
            NRCP=NT
      ENDIF
C
C set up angle PH4 between each irreducible beam and equivalent beams
C
      PI2=8.*ATAN(1.)
      DO 100 II=1,NT
C
C identify equivalent beams but only if they belong to the same beamset
C
         IIR=IDXK2(II)
C identify beamset
         IK=IDXK(IIR)
         JBM=0
         DO 103 IBS=1,KNBS
          IF(IIR.GT.JBM.AND.IIR.LE.JBM+KNB(IBS)) IBS2=IBS
          JBM=JBM+KNB(IBS)
103      CONTINUE
         NEQ=1
         PH4(1,II)=0.
         ZMOD=SQRT(PQAL(1,IIR)**2+PQAL(2,IIR)**2)
         JBM=0
         DO 101 IBS=1,KNBS
         DO 102 I2=1,KNB(IBS)
            I=I2+JBM
            IF(IDXK(I).EQ.IK.AND.I.NE.IIR.AND.IBS.EQ.IBS2) THEN
C
C the beam I is equivalent to IIR (different) and they belong to the
C same unsymmetrized beamset
C              
                DOT=PQAL(1,IIR)*PQAL(1,I)+PQAL(2,IIR)*PQAL(2,I)
                CROS=PQAL(1,IIR)*PQAL(2,I)-PQAL(2,IIR)*PQAL(1,I)
C CSV can be < -1 (numerical errors)
                CSV=DOT/(ZMOD*ZMOD)
                IF(CSV.LT.-1.) CSV=-1.
                NEQ=NEQ+1
                IF(CROS.GE.0.) THEN
C                   PH4(NEQ,II)=ACOS(DOT/(ZMOD*ZMOD))*360./PI2
                   PH4(NEQ,II)=ACOS(CSV)
                ELSE
                   PH4(NEQ,II)= (PI2-ACOS(CSV))
                ENDIF
           ENDIF
102      CONTINUE
         JBM=JBM+KNB(IBS)
101      CONTINUE
         IKBS(II)=NEQ
100   CONTINUE
      RETURN
      END
C =========================================================================
C
C  Subroutine GHSYM computes the reduced matrix TH corresponding to the
C  symmetric block of the unreduced TH.
C
C  The unreduced matrix is different from the usual matrix from equation
C  48 in Van-Hove Tong. Block I,J (I.NE.J) in the old matrix was given
C  by  -TAU(I)GH(I,J). Diagonal blocks were Unit matrices.
C  Now the unreduced matrix has  - t(I)*GH(I,J) as IJ block (I.NE.J) and
C  TAUINV(I)*t^-1 for the diagonal block I.
C  This implies a modification of the subroutine TAUT which has now
C  become TAUT2
C
C  From this new matrix  we extract the completely symmetric block.
C  This sub first computes
C  (LM)-space interplanar propagators GH for
C  the composite layer treated by subroutine MTNVSYM. Depending on the
C  interplanar apacing, either a reciprocal-space summation or a direct
C  space summation is performed (The latter in routine GHD). 
C  The reciprocal lattice sum has been recoded completely for more
C  efficiency
C
C  
C Parameter List;
C ===============
C
C  LMMAX        =   (LMAX+1)**2.
C  MGH          =   MATRIX CONTAINING KEY TO POSITION OF INDIVIDUAL GH^S IN THE
C                   MATRIX GH  MGH(I,J) IS SEQUENCE NUMBER OF GH(I,J) IN 
C                   COLUMNAR MATRIX GH.
C  NLAY          =  NO. OF SUBPLANES CONSIDERED.
C  NUGH          =  LIST OF THOSE GH^S THAT MUST BE COMPUTED.
C  NLAY2         =  NLAY*(NLAY-1)/2.
C  TST           =  INPUT QUANTITY FOR DETERMINING NO. OF POINTS REQUIRED IN
C                   RECIPROCAL LATTICE SUMMATION.
C  TEST,Y1,Y,S   =  WORKING SPACE.
C  L2M           =  2*LMAX+1.
C  LM            =  LMAX+1.
C  LMS           =  (2*LMAX+1)**2.
C  DRL           =  SET OF INTERPLANAR VECTORS.
C  TV            =  AREA OF UNIT CELL OF EACH SUBPLANE.
C  LXM           =  PERMUTATION OF (LM) SEQUENCE.
C  LEV           =  (LMAX+1)*(LMAX+2)/2.
C  DCUT          =  CUTOFF RADIUS FOR LATTICE SUMMATION.
C  CAA           =  CLEBSCH-GORDAN COEFFICIENTS.
C  NCAA          =  NO. OF CLEBSCH-GORDAN COEFFICIENTS.
C  LAY         =   1 IF CURRENT CALCULATION REFERS TO OVERLAYER.
C                  0 IF CURRENT CALCULATION REFERS TO SUBSTRATE.
C  PQ          =   ADDITIONAL OFFSET INTRODUCED FOR UNKNOWN PURPOSES
C 
C New for symmetric calculation
C  LMNI         reduced dimension of the block
C  IDX(I)       I=1,NLAY*LMMAX. Assign the reduced coordinate corresponding
C               to the old label (L in the SYMMETRIC order) I.
C               1<=IDX<=LMNI        
C  WBD(I)       weight to assign to the old label I to obtain the 
C               symmetric representation
C LL1(I,NCL)     Position of the Ith inequivalent layer among the NLAY
C i=1,NINEQ  layers
C
C In Common Blocks;
C =================
C
C  E             =  CURRENT ENERGY.
C  AK2,AK3       =  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C  VPI           =  IMAGINARY PART OF ENERGY.
C  ARA1,ARA2     =  BASIS VECTORS OF SUBSTRATE LAYER LATTICE.
C  ARB1,ARB2     =  BASIS VECTORS OF SUPERLATTICE.
C  RBR1,RBR2     =  RECIPROCAL LATTICE OF SUPERLATTICE.
C  NL1,NL2       =  SUPERLATTICE CHARACTERIZATION CODES.
C
C =========================================================================
C
      SUBROUTINE GHSYM6(LMMAX,MGH,NLAY,NUGH,NLAY2,
     & TEST,Y1,L2M,Y,LM,S,LMS,DRL,TV,LXM,LEV,DCUT,CAA,NCAA,LAY,PQ,LMNI,
     & TH,KOUNT,LPS,IDX,IDXN,WBD,TAUINV,NLTU,LMNO,
     & ZRED,TSF,LAN,NINEQ,LL1,NLMAX,NK3)
C
      DIMENSION KOUNT(LMMAX,LMMAX),LAN(LMMAX)
      DIMENSION MGH(NLAY,NLAY),NUGH(NLAY2)
      DIMENSION RBR1(2),RBR2(2),DRL(NLAY2,3),TEST(NLAY2),LXM(LMMAX)
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2),CAA(NCAA)
      DIMENSION RXR1(2),RXR2(2),PQ(2)
      COMPLEX Y(LM,LM),Y1(L2M,L2M),S(LMS),TH(LMNI,LMNI)
      COMPLEX CI,CZ,KPRG,K0,Z,T1,T2,T3,BS,CS,CFAC
      COMPLEX AK
      DIMENSION LPS(NLAY),IDX(LMNO),IDXN(LMNO)
      DIMENSION ZRED(LMNI),LL1(NLMAX)
      COMPLEX TAUINV(NLTU,LEV),WBD(LMNO),TSF(6,16)
C
      COMMON E,AK2,AK3,VPI
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
      PI2=8.*ATAN(1.)
      LOD=LMMAX-LEV
      CZ=(0.0,0.0)
      CI=(0.0,1.0)
      AK=-0.5/CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      IF (LAY.EQ.1) THEN
         RXR1(1)=RBR1(1)
         RXR1(2)=RBR1(2)
         RXR2(1)=RBR2(1)
         RXR2(2)=RBR2(2)
      ELSE
         PI=3.14159265
         ATV=2.0*PI/TV
         RXR1(1)=ARA2(2)*ATV
         RXR1(2)=-ARA2(1)*ATV
         RXR2(1)=-ARA1(2)*ATV
         RXR2(2)=ARA1(1)*ATV
      ENDIF
      K0=CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      CFAC=-16.0*(3.14159265)*(3.14159265)*CI/TV
C
C From each sum from 1 to NLAY2 where GH are defined we extract the 
C NLR NLC layer indeces
C NLR being the row, NLC the column index. See numbering on page 55
C Van-Hove Tong
C
C  INITIALIZE TH IN PREPARATION FOR NEW VALUES
C
            DO 1462 J=1,LMNI
               DO 7 K=1,LMNI
                  TH(J,K)=CZ
                  TH(J,K)=CZ
7              CONTINUE
1462        CONTINUE
C
C  TSTS= ESTIMATED NO. OF POINTS IN DIRECT LATTICE SUM
C
      TSTS=DCUT*DCUT*3.14159265/TV
      DO 10 NLR2=1,NINEQ
      DO 11 NLC=1,NLAY
         NLR=LL1(NLR2)
         IF(NLR.EQ.NLC) GOTO 11
         IZ=MGH(NLR,NLC)
         K=1
C         IF (NLR.GT.NLC) IZ=MGH(NLC,NLR)
         IF (NLR.GT.NLC) K=2
            IF (ABS(DRL(IZ,1)).GT.0.001) THEN
C
C  AKP2= ESTIMATED NO. OF POINTS IN RECIPROCAL LATTICE SUM
C
C               FACT1=ALOG(TST)/DRL(IZ,1)
               FACT1=ALOG(0.0002)/DRL(IZ,1)
               AKP2=(2.0*E+FACT1*FACT1)*TV/(4.*3.1415926)
C
C  SKIP DIRECT LATTICE SUM, IF RECIPROCAL LATTICE SUM FASTER (BUT NUMBER
C  OF REC. LATT. POINTS IS TO BE RESTRICTED DUE TO A CONVERGENCE
C  PROBLEM)
C
               IF ((TSTS.GE.2.0*AKP2).AND.(AKP2.LT.80.0)) GOTO 11
            ENDIF
C
C  PRODUCE GH(I,J) AND GH(J,I) WITH DIRECT LATTICE SUM FOR TWO
C  PROPAGATION DIRECTIONS
C
              CALL GHDSYM(IZ,K,LMMAX,S,LMS,Y1,L2M,DRL,NLAY2,K0,
     &        DCUT,CAA,NCAA,LXM,LAY,PQ,TH,LMNI,NLR,NLC,WBD,IDX,IDXN,
     &        LMNO,KOUNT)
C
C  GH(I,J) AND GH(J,I) NOW NO LONGER NEED TO BE COMPUTED
C
            NUGH(IZ)=0
C         ENDIF
11    CONTINUE
10    CONTINUE
      TSTS=0.0
      DO 12 NLR2=1,NINEQ
      DO 13 NLC=1,NLAY
         NLR=LL1(NLR2)
         IF(NLR.EQ.NLC) GOTO 13
         I=MGH(NLR,NLC)
C         IF (NLR.GT.NLC) I=MGH(NLC,NLR)
         IF (.NOT.((NUGH(I).EQ.0).OR.(ABS(DRL(I,1)).LE.0.001))) THEN
            DRL(I,1)=DRL(I,1)
C
C  TEST(I) WILL SERVE AS CUTOFF IN RECIPROCAL LATTICE SUM
C  Notice that the convergence of the planar sum is not
C  related to the convergence of the RFS scheme controlled by TST
C
CC            TEST(I)=ABS(ALOG(TST)/DRL(I,1))
            TEST(I)=ABS(ALOG(.0002)/DRL(I,1))
            TSTS=AMAX1(TEST(I),TSTS)
         ENDIF
13    CONTINUE
12    CONTINUE
C
C We have TH(K,KP) for all K and KP whose NLR and
C NLC  correspond to direct lattice sums 
C But for each K we do not have yet TH(K,KP) for those KP corresponding
C to reciprocal lattice sums
C
      IF (TSTS.GT.0.00001) THEN
C
C  START OF TWO 1-DIMENSIONAL SUMMATION LOOPS IN ONE QUADRANT OF
C  RECIPROCAL SPACE
C
         NCOUNT=0
         NUMG=0
         JJ1=0
1171     JJ1=JJ1+1
         JJ2=0
1172     JJ2=JJ2+1
         NOG=0
         J1=JJ1-1
         J2=JJ2-1
C
C  START OF LOOP OVER QUADRANTS
C
         DO 1370 KK=1,4
            IF (KK.EQ.2) THEN
               IF (J1.EQ.0.AND.J2.EQ.0) GOTO 1380
               IF (J2.EQ.0) GOTO 1370
               NG1=J1
               NG2=-J2
            ELSEIF (KK.EQ.3) THEN
               IF (J1.EQ.0) GOTO 1370
               NG1=-J1
               NG2=J2
            ELSEIF (KK.NE.4) THEN
               NG1=J1
               NG2=J2
            ELSEIF (J1.EQ.0.OR.J2.EQ.0) THEN
               GOTO 1370
            ELSE
               NG1=-J1
               NG2=-J2
            ENDIF
C
C  CURRENT RECIPROCAL LATTICE POINT
C
            GX=NG1*RXR1(1)+NG2*RXR2(1)
            GY=NG1*RXR1(2)+NG2*RXR2(2)
            GKX=GX+AK2+PQ(1)
            GKY=GY+AK3+PQ(2)
            GK2=GKX*GKX+GKY*GKY
C
C  TEST FOR CUTOFF
C
            KPRG=CSQRT(CMPLX(2.0*E-GK2,-2.0*VPI+0.000001))
            AKP2=AIMAG(KPRG)
            IF (AKP2.LE.(TSTS)) THEN
               NUMG=NUMG+1
               NOG=1
               Z=KPRG/K0
               FY=0.0
               IF (GK2.GT.1.0E-8) THEN
                  CFY=GKX/SQRT(GK2)
                  IF (ABS(ABS(CFY)-1.).LE.1.E-6) THEN
                     IF (CFY.LT.0.0) FY=3.14159265
                  ELSE
                     FY=ACOS(CFY)
                  ENDIF
                  IF (GKY.LT.0.0) FY=-FY
               ENDIF
C
C  FIND APPROPRIATE SPHERICAL HARMONICS
C
               CALL SH(LM,Z,FY,Y)
C
C  START OF LOOP OVER INTERPLANAR VECTORS
C
               DO 1361 NLR2=1,NINEQ
                 NLR=LL1(NLR2)
                 NII2=(NLR-1)*LMMAX
C                 IT2=(LPS(NLR)-1)+NTAUSH
                 DO 1350 I=1,LMMAX
C
C I will be used as a natural index
C
                   IP=LXM(I)
                   K=IDXN(NII2+IP)
                   IF (K.EQ.0) GOTO 1350
C
C symmetry limit calculations
C
                      DO 1362 NLC=1,NLAY
                        NNI2=(NLC-1)*LMMAX
                        IF(NLC.NE.NLR) THEN
                          IZ=MGH(NLR,NLC)
C                          IF (NLR.GT.NLC) IZ=MGH(NLC,NLR)
C
C  SKIP IF NEW GH(I,J) NOT NEEDED (it might have been calculated with
C  real space sums
C 
                          IF (NUGH(IZ).NE.0) THEN
C
C  SKIP IF THIS CONTRIBUTION OUTSIDE CUTOFF
C
                            IF (AKP2.LE.(TEST(IZ))) THEN
                              T2=GKX*DRL(IZ,2)+GKY*DRL(IZ,3)
C                              T2P=-T2
                              T3=KPRG*ABS(DRL(IZ,1))
                              T1=(CEXP(CI*(T2+T3))/KPRG)*CFAC
C                              T1P=(CEXP(CI*(T2P+T3))/KPRG)*CFAC
                              L1=INT(SQRT(FLOAT(I-1)+0.00001))
                              M1=I-L1-L1*L1-1
                              IF (M1.GT.0) THEN
                                 BS=Y(M1,L1+1)
                              ELSE
                                 BS=(-1)**(MOD(M1,2))*Y(L1+1,-M1+1)
                              ENDIF
                              DO 1340 J=1,LMMAX
                                 JP=LXM(J)
                                 KP=IDX(NNI2+JP)
                                 IF (KP.EQ.0) GOTO 1340
                                 L2=INT(SQRT(FLOAT(J-1)+0.00001))
                                 M2=J-L2-L2*L2-1
                                 IF (M2.GE.0) THEN
                                    CS=Y(L2+1,M2+1)
                                 ELSE
                                    CS=(-1)**(MOD(M2,2))*Y(-M2,L2+1)
                                 ENDIF
C
C multiply standard term by symmetry weight. Notice the extra -
C (-GH(I,J))
C
            IF(NLR.GT.NLC) THEN
C                 TH(K,KP)=TH(K,KP)-T1P*BS*CS*WBD(NNI2+JP)
                 TH(K,KP)=TH(K,KP)-T1*BS*CS*WBD(NNI2+JP)
            ELSE
                 TH(K,KP)=TH(K,KP)-(-1)**(MOD(L1+M1+L2+M2,2))
     &               *T1*BS*CS*WBD(NNI2+JP)
            ENDIF
1340                          CONTINUE
                            ENDIF
                          ENDIF
                        ENDIF
1362                  CONTINUE
1350             CONTINUE
1361           CONTINUE
            ENDIF
1370     CONTINUE
1380     IF (NOG.EQ.1) GOTO 1172
         IF (JJ2.NE.1) GOTO 1171
         IF (JJ2.EQ.1) NCOUNT=NCOUNT+1
         IF (NCOUNT.LE.3) GOTO 1172

      ENDIF
C      I3=MCLOCK()
C      FTIME=FLOAT(I3-I2)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU
C
CCCC multiply each row by t (analogous to multiply by tau in the old version)
C multiply each column by t (analogous to multiply by tau in the old version)
C In this way the output TS will be the amplitude of the wave before being
C last scattered by the atom of the corresponding sublayer
C
      DO 1500 NLR=1,NLAY 
         NII2=(NLR-1)*LMMAX
         DO 1510 I=1,LMMAX 
            K=IDXN(NII2+I)
            IF (K.EQ.0) GOTO 1510
C
C identify L of the row and element type
C
            LL=LAN(I)+1
            IT2=LPS(NLR)
            CS=AK*TSF(IT2,LL)
            DO 1550 KP=1,LMNI 
C               TH(K,KP)=TH(K,KP)*CS
               TH(KP,K)=TH(KP,K)*CS
1550        CONTINUE
1510     CONTINUE
1500  CONTINUE
C
C take care of NLR=NLC by inserting TAUINV. Notice that L's are in the
C symmetric order in TAUINV, in the same order in TH
C
      DO 1440 NLR=1,NLAY
         NII2=(NLR-1)*LMMAX
         DO 1450 I=1,LMMAX
            K=IDXN(NII2+I)
            IF (K.EQ.0) GOTO 1450
C
C symmetry limit calculations
C
              IN2=(LPS(NLR)-1)*LMMAX
              LSH=0
              IF(I.GT.LEV) THEN
                 LSH=LEV
              ENDIF
              DO 1460 J=1,LMMAX
                 KP=IDX(NII2+J)
                 IF (KP.EQ.0) GOTO 1460
                 IF((I.LE.LEV.AND.J.LE.LEV).OR.
     &           (I.GT.LEV.AND.J.GT.LEV)) THEN
                   TH(K,KP)=TH(K,KP)+TAUINV(IN2+I,J-LSH)*
     &             WBD(NII2+J)
                 ENDIF
1460          CONTINUE
1450     CONTINUE
1440   CONTINUE
C
C Fix the normalization of TH 
C
      DO 1600 K=1,LMNI 
      DO 1650 KP=1,LMNI 
         TH(K,KP)=TH(K,KP)*ZRED(K)
1650  CONTINUE
1600  CONTINUE
C      I4=MCLOCK()
C      FTIME=FLOAT(I4-I3)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU
      RETURN
      END
C =========================================================================
C
C  Subroutine GHDSYM computes direct lattice sums for (LM)-space
C  propagators GH between two subplanes (having Bravais lattices) of
C  a composite layer. The subplanes may be coplanar.
C
C Parameter List;
C ==============
C
C  IZ          =   SERIAL NO. OF CURRENT INTERPLANAR VECTOR DRL.
C  IS          =   1 FOR PROPAGATION FROM FIRST TO SECOND SUBPLANE.
C                  2 FOR PROPAGATION FROM SECOND TO FIRST SUBPLANE.
C  GH          =   OUTPUT INTERPLANAR PROPAGATOR.
C  S           =   WORKING SPACE (LATTICE SUM).
C  LMS         =   (2*LMAX+1)**2.
C  Y           =   WORKING SPACE (SPHERICAL HARMONICS).
C  L2M         =   2*LMAX+1.
C  K0          =   COMPLEX MAGNITUDE OF WAVEVECTOR.
C  DCUT        =   CUTOFF RADIUS FOR LATTICE SUM.
C  CAA         =   CLEBSCH-GORDAN COEFFICIENTS FROM SUBROUTINE CAAA.
C  NCAA        =   NO. OF CLEBSCH-GORDAN COEFFICIENTS IN CAA.
C  LXM         =   PERMUTATION OF (LM) SEQUENCE FROM SUBROUTINE LXGENT.
C  LAY         =   1 IF CURRENT CALCULATION REFERS TO OVERLAYER.
C                  0 IF CURRENT CALCULATION REFERS TO SUBSTRATE.
C  PQ          =   BEAM LIST.
C
C Note;
C =====
C
C Dimension 33 set for (LMAX.LE.15).
C
C
C Modified version of routine GHD from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE GHDSYM(IZ,IS,LMMAX,S,LMS,Y,L2M,DRL,NLAY2,K0,DCUT,
     & CAA,NCAA,LXM,LAY,PQ,TH,LMNI,NLR,NLC,WBD,IDX,IDXN,LMNO,
     & KOUNT)
C
      DIMENSION KOUNT(LMMAX,LMMAX)
      COMPLEX Y(L2M,L2M),H(33),S(LMS)
      COMPLEX RU,CI,CZ,K0,FF,Z,Z1,Z2,Z3,ST,TH(LMNI,LMNI)
      DIMENSION DRL(NLAY2,3),ARA1(2),ARA2(2),ARB1(2),ARB2(2)
      DIMENSION RBR1(2),RBR2(2),CAA(NCAA)
      DIMENSION V(3),LXM(LMMAX),PQ(2)
      DIMENSION IDX(LMNO),IDXN(LMNO)
      COMPLEX WBD(LMNO)
C
      COMMON E,AK2,AK3,VPI
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
      RU=(1.0,0.0)
      CZ=(0.0,0.0)
      CI=(0.0,1.0)
      DCUT2=DCUT*DCUT
      DO 5 I=1,LMS
         S(I)=CZ
5     CONTINUE
C
C V(I)=R(NLR)-R(NLC) if NLR<NLC
C V(I)=-R(NLR)+R(NLC) if NLR>NLC
C
      DO 10 I=1,3
         V(I)=DRL(IZ,I)
10    CONTINUE
C
C  TURN INTERPLANAR VECTOR AROUND IF IS=2
C
C      IF (IS.GE.2) THEN
C         DO 30 I=1,3
C            V(I)=-V(I)
C30       CONTINUE
C      ENDIF
C
C  START OF TWO 1-DIMENSIONAL LATTICE LOOPS FOR SUMMATION OVER 1 QUADRANT
C
C      I1=MCLOCK()
      NCOUNT=0
      NUMR=0
      JJ1=0
50    JJ1=JJ1+1
      JJ2=0
60    JJ2=JJ2+1
      NOR=0
      J1=JJ1-1
      J2=JJ2-1
C
C  START OF LOOP OVER QUADRANTS
C
      DO 140 KK=1,4
         IF (KK.EQ.2) THEN
            IF (J1.EQ.0.AND.J2.EQ.0) GOTO 150
            IF (J2.EQ.0) GOTO 140
            NR1=J1
            NR2=-J2
         ELSEIF (KK.EQ.3) THEN
            IF (J1.EQ.0) GOTO 140
            NR1=-J1
            NR2=J2
         ELSEIF (KK.NE.4) THEN
            NR1=J1
            NR2=J2
         ELSEIF (J1.EQ.0.OR.J2.EQ.0) THEN
            GOTO 140
         ELSE
            NR1=-J1
            NR2=-J2
         ENDIF
         IF (LAY.EQ.1) THEN
            PX=NR1*ARB1(1)+NR2*ARB2(1)
            PY=NR1*ARB1(2)+NR2*ARB2(2)
         ELSE
            PX=NR1*ARA1(1)+NR2*ARA2(1)
            PY=NR1*ARA1(2)+NR2*ARA2(2)
         ENDIF
         X1=(PX+V(2))*(PX+V(2))+(PY+V(3))*(PY+V(3))
C
C  CUTOFF OF LATTICE SUMMATION AT RADIUS DCUT
C
         IF (X1.LE.DCUT2) THEN
            NOR=1
            NUMR=NUMR+1
            Z1=CEXP(-CI*(PX*(AK2+PQ(1))+PY*(AK3+PQ(2))))
            X2=SQRT(X1+V(1)*V(1))
            X1=SQRT(X1)
            Z2=K0*X2
            Z=CMPLX(V(1)/X2,0.0)
            ZR=V(1)/X2
            FY=0.0
            IF (ABS(X1).GE.1.E-6) THEN
               CFY=(PX+V(2))/X1
               IF (ABS(ABS(CFY)-1.).LE.1.E-6) THEN
                  IF (CFY.LT.0.0) FY=3.14159265
               ELSE
                  FY=ACOS(CFY)
               ENDIF
               IF ((PY+V(3)).LT.0.0) FY=-FY
            ENDIF
C
C  COMPUTE REQUIRED BESSEL FUNCTIONS AND SPHERICAL HARMONICS
C
            CALL SB(Z2,H,L2M)
            CALL SHR(L2M,ZR,FY,Y)
            ST=RU
            DO 130 L=1,L2M
               ST=ST*CI
               Z3=ST*H(L)*Z1
               L1=L*L-L
               CSGN=-1.
               DO 120 M=1,L
                  CSGN=CSGN*(-1.)
                  M1=M-1
                  IPM=L1+M
                  IMM=L1-M+2
                  S(IPM)=S(IPM)+Z3*Y(L,M)
                  IF (M.NE.1) S(IMM)=S(IMM)+Z3*Y(M-1,L)*CSGN
C
C  S NOW CONTAINS THE LATTICE SUM
C
120            CONTINUE
130         CONTINUE
         ENDIF
140   CONTINUE
150   IF (NOR.EQ.1) GOTO 60
      IF (JJ2.NE.1) GOTO 50
      IF (JJ2.EQ.1) NCOUNT=NCOUNT+1
      IF (NCOUNT.LE.3) GOTO 60
C
C  PRINT NUMBER OF LATTICE POINTS USED IN SUMMATION
C
      FF=-8.0*3.14159265*K0
C
C  USE SUBROUTINE GHSC TO MULTIPLY LATTICE SUM INTO CLEBSCH-GORDAN
C  COEFFICIENTS
C
      CALL GHSCSY(LMMAX,S,LMS,CAA,NCAA,FF,LXM,
     & TH,LMNI,NLR,NLC,WBD,IDX,IDXN,LMNO,KOUNT)
      RETURN
      END
C =========================================================================
C
C  Subroutine GHSCSY computes the (LM)-space interplanar propagators GH
C  from direct lattice sums produced in subroutine GHD and Clebsch-
C  Gordan coefficients from subroutine CAAA.
C
C Parameter List;
C ===============
C
C  IZ       =    SERIAL NO. OF CURRENT INTERPLANAR VECTOR DRL.
C  IS       =    1 FOR PROPAGATION FROM FIRST TO SECOND SUBPLANE.
C                2 FOR PROPAGATION FROM SECOND TO FIRST SUBPLANE.
C  FF       =    PREFACTOR OF GH.
C  LXM      =    PERMUTATION OF (LM) SEQUENCE
C
C  For other quantities see subroutines GHD and GHMAT.
C
C Modified version of routine GHSC from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE GHSCSY(LMMAX,S,LMS,CAA,NCAA,FF,LXM,
     & TH,LMNI,NLR,NLC,WBD,IDX,IDXN,LMNO,KOUNT)
C
      DIMENSION KOUNT(LMMAX,LMMAX)
      COMPLEX S(LMS),FF,TH(LMNI,LMNI)
      DIMENSION CAA(NCAA),LXM(LMMAX)
      DIMENSION IDX(LMNO),IDXN(LMNO)
      COMPLEX WBD(LMNO),ZT
C
      II=1
      NII2=(NLR-1)*LMMAX
      NNI2=(NLC-1)*LMMAX
      DO 1350 I=1,LMMAX
         IP=LXM(I)
         K=IDXN(NII2+IP)
C
C symmetry limits calculations but need to increase II defining CC
C
           L1=INT(SQRT(FLOAT(I-1)+0.00001))
           M1=I-L1-L1*L1-1
           DO 1340 J=1,LMMAX
              JP=LXM(J)
              KP=IDX(NNI2+JP)
              KCH=K*KP
              IF (KCH.EQ.0) THEN
                II=II+KOUNT(I,J)
                GOTO 1340
              ENDIF
              L2=INT(SQRT(FLOAT(J-1)+0.00001))
              M2=J-L2-L2*L2-1
              M3=M2-M1
              IL=IABS(L1-L2)
              IM=IABS(M3)
              LMIN=MAX0(IL,IM+MOD(IL+IM,2))
              LMAX=L1+L2
              LMIN=LMIN+1
              LMAX=LMAX+1
              ZT=(0.,0.)
              DO 203 ILA=LMIN,LMAX,2
                 LA=ILA-1
                 CC=CAA(II)
                 II=II+1
C                 TH(NII2+IP,NNI2+JP)=TH(NII2+IP,NNI2+JP)+CC*
C     &            S(LA*LA+LA+M3+1)
                 ZT=ZT+CC*S(LA*LA+LA+M3+1)
203           CONTINUE
                 TH(K,KP)=TH(K,KP)-ZT*FF*WBD(NNI2+JP)
C              TH(NII2+IP,NNI2+JP)=FF*TH(NII2+IP,NNI2+JP)
1340       CONTINUE
1350  CONTINUE
      RETURN
      END
C
C =========================================================================
C
C Sub.  
      SUBROUTINE SYMCL(NINEQ,LAFLG2,NDBD,IROT,IMIR,NLMAX,LMNI,NCL,LL1,
     &LL2,WBDS,IDXS,IDXN,ZRED,NST1,LMMAX,LX,LXM,LMNBD,NSU,NSTEF,NLAY)
C
      COMPLEX OM(7),WBDS(LMNI)
      COMPLEX CI,RU
      DIMENSION IDXN(LMNI),IDXS(LMNI),LL1(NLMAX,NST1)
      DIMENSION ZRED(LMNI),NINEQ(NST1)
      DIMENSION LAFLG2(NLMAX,NST1)
      DIMENSION LMNBD(2),LL2(12,NLMAX,NST1)
      DIMENSION LX(LMMAX),LXM(LMMAX)
      LOGICAL  LOGPG,LOGP,LOGAX,LOGM,REDR,REDM
cjcm
      DIMENSION RBR1(2),RBR2(2),ARB1(2),ARB2(2),ARA1(2),ARA2(2)
cjcm

C
1610  FORMAT (' NO SYMMETRY in MTNVSYM for GLIDE SYMMETRY. Lower
     & the symmetry of NSYMS (for the TLEED1.f run) to exclude the
     & glide ')
      COMMON E,AK2,AK3,VPI
      COMMON /MFB/GP,L1,LAY1
      COMMON /MS/LMAX
      COMMON /MPT/NA,NS,LAY,LLM,NTAU,TST,TV,DCUT,NPERT,NOPT,NEW
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
      CI=CMPLX(0.0,1.0)
      RU=CMPLX(1.0,0.0)
      PI=4.0*ATAN(1.0)
      PI4=16.0*ATAN(1.0)
      PI3=PI/3.
      CNM=1./SQRT(2.)
      CNAX=1./SQRT(FLOAT(IROT))
      CNMAX=1./SQRT(2.*FLOAT(IROT))
C
C Inversion Symmetry and/or 2-d Point Group Symmetry?
C Inversion symmetry of the composite layer can be used if present.
C Remember however that, even if a composite layer has I symmetry,
C the problem does not.
C We allow for Cn + one mirror as maximal symmetry
C Determine total number of blocks (NDSYM).
C 
C
      LOGPG=.FALSE.
      LOGP=.FALSE.
      NDP=1
      IF(NDBD.EQ.1) THEN
C no symmetry
         LOGPG=.FALSE.
         LOGM=.FALSE.
         LOGAX=.FALSE.
         NDM=1
      ELSEIF(IROT.EQ.1) THEN
C mirror plane only
         LOGPG=.TRUE.
         LOGM=.TRUE.
         LOGAX=.FALSE.
         NDM=2
      ELSEIF(IROT.GT.1.AND.NDBD.EQ.IROT) THEN
C Cn axis only 
         LOGPG=.TRUE.
         LOGAX=.TRUE.
         LOGM=.FALSE.
         NDM=1
      ELSEIF(IROT.GT.1.AND.NDBD.EQ.2*IROT) THEN
C Cn axis and mirror
         LOGPG=.TRUE.
         LOGAX=.TRUE.
         LOGM=.TRUE.
         NDM=2
      ENDIF
      NDPG=NDBD
C
C Total symmetry.      
C
      NDSYM=NDPG*NDP
C
C Define weights for n-fold symmetry 
C 2-fold:
      IF(IROT.EQ.2) THEN
         OM(1)=RU
         OM(2)=-RU
         OM(3)=RU
      ELSEIF(IROT.EQ.3)THEN
C 3-fold:
         OM(1)=RU
         OM(2)=CEXP(2.*CI*PI3)
         OM(3)=OM(2)*OM(2)
         OM(4)=OM(1)
      ELSEIF(IROT.EQ.4)THEN
C 4-fold:
         OM(1)=RU
         OM(2)=CI
         OM(3)=OM(2)*OM(2)
         OM(4)=OM(3)*OM(2)
         OM(5)=OM(1)
      ELSEIF(IROT.EQ.6)THEN
C 6-fold:
         OM(1)=RU
         OM(2)=CEXP(CI*PI3)
         OM(3)=OM(2)*OM(2)
         OM(4)=OM(3)*OM(2)
         OM(5)=OM(4)*OM(2)
         OM(6)=OM(5)*OM(2)
         OM(7)=OM(1)
      ELSE
C trivial
         OM(1)=RU
         OM(2)=RU
      ENDIF
C
C We are only interested in the maximally symmetric representation
C (Maximal symmetry refer to the 2d symmetry, not inversion) 
C In the general case we could not label states by the rotational
C eigenvalue AND the mirror plane symmetry. But in the symmetrical
C subspace (which is left invariant by the mirror plane operation)
C we can.
C
C Compute dimension of totally symmetric block LMNBD() ( the new LMNI)
C split into two in case of inversion symmetry.
C The inversion symmetry is not implemented yet
C
C We need:
C LL1(I,NCL)     Position of the Ith inequivalent layer among the NLAY
C i=1,NINEQ  layers
C
C LL2(S,I,NCL) Position of the Ith layer after Symmetry operation S
C           (The order has been fixed in the following way:
C            S=1,IROT corresponds to clockwise rotation of 2PI/IROT (S-1) 
C            S=IROT+1, 2IROT (if there is an additional mirror plane)
C              to mirror followed by same rotation (in that order,
C              notice that mirror and rotation do not commute)
C
C generate:
C IDXS(I), I=1, NLAY*LMMAX give the reduced index
C corresponding to unreduced  symmetric ( that is the L labeling is in
C the symmetric order) index .
C IDXN(I), I=1, NLAY*LMMAX gives the reduced index for one selected
C element in each equivalent class. 0 for all other elements
C
C WBDS(I), I=1, NLAY*LMMAX give the weight to obtain
C the symmetric representation (used in GHSYM and to transform back
C the symmetric solution) corresponding to the unreduced index
C (symmetric). To understand WBDS consider the transformation properties
C of spherical harmonics
C
C ZRED,LMNBD are also generated
C

      DO 20 IP=1,NDP
C
         K=0
         DO 23 I=1,NINEQ(NCL)
            IDX1=(LL1(I,NCL)-1)*LMMAX
C
C First IF consider the case of trivial point group symmetry and, if
C LOGPG, the case of trivially transforming sublayers.
C If LOGPG each inequivalent sublayer is in a class with either 1, IROT, 
C 2 or 2*IROT elements.  The first IF deals with the case where the
C sublayer is invariant under rotation (REDR=.TRUE.) or under a mirror
C and as a result the number of elements in the class is < IROT 
C The second and third 
C with the case  where the sublayer transforms under rot.(REDR=.FALSE.)
C
            REDR=.FALSE.
            REDM=.FALSE.
            IF (LOGAX) THEN
              IF(LL2(2,LL1(I,NCL),NCL).EQ.LL1(I,NCL)) REDR=.TRUE.
            ENDIF
            IF (LOGM) THEN
              IF(LL2(IROT+1,LL1(I,NCL),NCL).EQ.LL1(I,NCL)) REDM=.TRUE.
            ENDIF
C
            IF(.NOT.LOGPG.OR.(LOGAX.AND.LAFLG2(I,NCL).LT.IROT).OR.
     &       (.NOT.LOGAX.AND.LAFLG2(I,NCL).EQ.1).OR.(LOGAX.AND
     &       .REDR)) THEN
C
C different m's correspond to different point group symmetry eigenvalues
C Dimension is reduced
C
                 DO 25 LM=1,LMMAX
                      L=INT(SQRT(FLOAT(LX(LM)-1)+0.0001))
                      M=LX(LM)-L-L*L-1
                      LTEST=MOD(L,2)
                      IDX1S=IDX1+LM
C
C IDX2 is the index of -M, IDX1S and IDX2S refer to symmetric indexing
C N(I)=(L,M) for I=L*L+L+M+1
C
                      LM2=LXM(L*L+L+1-M)
                      CNM2=1.
                      IDX2S=IDX1+LM2
C
C mods if sublayer transforms under mirror
C
                      IF(LOGM.AND..NOT.REDM) THEN
                        CNM2=CNM
                        IDX2S=(LL2(IROT+1,LL1(I,NCL),NCL)-1)
     &                 *LMMAX + LM2
                      ENDIF
                      IF(LOGAX.AND..NOT.REDR) THEN
                        CNAX2=CNM
                        CNMAX2=CNAX2*CNM2
                        IROT2=2
                      ENDIF
C
C If parity is a symmetry then we count only even L's for IP=1
C odd l's for IP=2 (not implemented)
C
                      IF(NDP.EQ.1.OR.IP-1.EQ.LTEST) THEN
C
C shift to make sure MTEST>0
C
                         ICN=1
                         MTEST=M-ICN+1+20*IROT
C
C Only states with positive M are of interest if there is a mirror plane
C symmetry and LAFLG2(I,NCL).EQ.1 . 
C We count the states and generate weight and index of the
C equivalent states.
C
                         IF(MOD(MTEST,IROT).EQ.0.AND.
     &                      (REDR.OR..NOT.LOGAX)) THEN
                            IF(.NOT.LOGM.OR.M.EQ.0) THEN
                               K=K+1
                               WBDS(IDX1S)=CNM2
                               IDXS(IDX1S)=K
                               IDXN(IDX1S)=K
                               ZRED(K)=1./CNM2
C
C equivalent states if LAFLG2(I,NCL).EQ.2
C
                               IF((.NOT.REDM).AND.LOGM) THEN
                                  IF(IMIR.EQ.1) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)
                                  ELSEIF(IMIR.EQ.2) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(-1)**M
                                  ELSEIF(IMIR.EQ.3) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(-CI)**M
                                  ELSEIF(IMIR.EQ.4) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(CI)**M
                                  ENDIF
                                  IDXS(IDX2S)=K
                                  IDXN(IDX2S)=0
                               ENDIF
C
C if the sublayer is invariant only M>0 are of
C interest, otherwise both M<>0 are important.
C
                            ELSEIF(M.GT.0.OR..NOT.REDM) THEN
                               K=K+1
                               WBDS(IDX1S)=CNM
                               IDXS(IDX1S)=K
                               IDXN(IDX1S)=K
                               ZRED(K)=1./CNM
C
C equivalent states
C
                               IF(IMIR.EQ.1) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)
                               ELSEIF(IMIR.EQ.2) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(-1)**M
                               ELSEIF(IMIR.EQ.3) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(-CI)**M
                               ELSEIF(IMIR.EQ.4) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(CI)**M
                               ELSE
                                  WRITE (1,1610)
                                  STOP
                               ENDIF
                               IDXS(IDX2S)=K
                               IDXN(IDX2S)=0
                            ENDIF
                         ELSEIF(.NOT.REDR.AND.LOGAX) THEN
C
C take care of states which transform under rotation but such that
C LAFLG2(I)<IROT. In this case the effective rot. group is a subgroup:
C Can only happen for IROT=4 IROTEF=2. However there is still the
C possibility that the two equivalent subl are connected also by the
C mirror
C
                            IF(MOD(MTEST,IROT2).EQ.0)THEN
C
C state acceptable under rot (no mirror to worry about)
C
                              IF(.NOT.LOGM.OR.M.EQ.0) THEN
                                 K=K+1
                                 WBDS(IDX1S)=CNAX2
                                 IDXS(IDX1S)=K
                                 IDXN(IDX1S)=K
                                 ZRED(K)=1./CNAX2
C
C equivalent states
C
                                 IDX2S=(LL2(2,LL1(I,NCL),NCL)-1)
     &                           *LMMAX + LM
                                 I3=MOD(M*(2-1)+20*IROT,IROT)+1
                                 WBDS(IDX2S)=WBDS(IDX1S)*OM(I3)
                                 IDXS(IDX2S)=IDXS(IDX1S)
                                 IDXN(IDX2S)=0
                              ELSEIF(M.GT.0) THEN
C
C possible mirror to worry about
C
                                 K=K+1
                                 WBDS(IDX1S)=CNAX2
                                 ZRED(K)=1./CNAX2
                                 IF(.NOT.REDM) THEN
C the sublayer transform
                                   WBDS(IDX1S)=CNMAX2
                                   ZRED(K)=1./CNMAX2
                                 ENDIF
                                 IDXS(IDX1S)=K
                                 IDXN(IDX1S)=K
C
C equivalent state under mirror
C
                                 IDX2S=(LL2(IROT+1,LL1(I,NCL),NCL)-1)
     &                           *LMMAX + LM2
                                 IF((.NOT.REDM).AND.LOGM) THEN
                                  IF(IMIR.EQ.1) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)
                                  ELSEIF(IMIR.EQ.2) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(-1)**M
                                  ELSEIF(IMIR.EQ.3) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(-CI)**M
                                  ELSEIF(IMIR.EQ.4) THEN
                                   WBDS(IDX2S)=WBDS(IDX1S)*(CI)**M
                                  ENDIF
                                  IDXS(IDX2S)=K
                                  IDXN(IDX2S)=0
                                 ENDIF
C
C equivalent state under rot
C

                                 IDX2S=(LL2(2,LL1(I,NCL),NCL)-1)
     &                           *LMMAX + LM
                                 I3=MOD(M*(2-1)+20*IROT,IROT)+1
                                 WBDS(IDX2S)=WBDS(IDX1S)*OM(I3)
                                 IDXS(IDX2S)=IDXS(IDX1S)
                                 IDXN(IDX2S)=0
C
C equivalent state under rot + mir
C
                                 IDX2S=(LL2(1,LL1(I,NCL),NCL)-1)
     &                           *LMMAX + LM2
                                 IF((.NOT.REDM).AND.LOGM) THEN
                                  IF(IMIR.EQ.1) THEN
                               WBDS(IDX2S)=WBDS(IDX1S)*OM(I3)
                                  ELSEIF(IMIR.EQ.2) THEN
                               WBDS(IDX2S)=WBDS(IDX1S)*(-1)**M*OM(I3)
                                  ELSEIF(IMIR.EQ.3) THEN
                               WBDS(IDX2S)=WBDS(IDX1S)*(-CI)**M*OM(I3)
                                  ELSEIF(IMIR.EQ.4) THEN
                               WBDS(IDX2S)=WBDS(IDX1S)*(CI)**M*OM(I3)
                                  ENDIF
                                  IDXS(IDX2S)=K
                                  IDXN(IDX2S)=0
                                 ENDIF
                              ENDIF
                            ENDIF
                         ELSE
C
C these states do not belong to the symmetric rep.
C
                            WBDS(IDX1S)=0.
                            IDXS(IDX1S)=0
                            IDXN(IDX1S)=0
                         ENDIF
                      ELSE
C
C when inversion is implemented these states should be related to
C equivalent states.
C
                         WBDS(IDX1S)=0.
                         IDXS(IDX1S)=0
                         IDXN(IDX1S)=0
                      ENDIF
25               CONTINUE
C
C Now go through all the other cases
C
         ELSEIF(LOGAX.AND.LAFLG2(I,NCL).EQ.IROT) THEN
C
C the counting should give        K=K+(LMAX+1)*(LMAX+2)/2 if LOGM
C but we need weights             K=K+(LMAX+1)*(LMAX+1) if .NOT.LOGM
C
                 DO 30 LMS=1,LMMAX
                      L=INT(SQRT(FLOAT(LX(LMS)-1)+0.0001))
                      M=LX(LMS)-L-L*L-1
                      LM2S=LXM(L*L+L+1-M)
                      IF(M.EQ.0) THEN
                        K=K+1
                        WBDS(IDX1+LMS)=CNAX
                        ZRED(K)=1./CNAX
                        IDXS(IDX1+LMS)=K
                        IDXN(IDX1+LMS)=K
C
C equivalent states under rotation and mirror
C
                        DO 34 I2=2,IROT
                           IDX2=(LL2(I2,LL1(I,NCL),NCL)-1)*LMMAX
                           WBDS(IDX2+LMS)=WBDS(IDX1+LMS)
                           IDXS(IDX2+LMS)=IDXS(IDX1+LMS)
                           IDXN(IDX2+LMS)=0
34                      CONTINUE 
                      ELSEIF(M.GT.0.OR..NOT.LOGM) THEN
                        K=K+1
                        WBDS(IDX1+LMS)=CNAX
                        ZRED(K)=1./CNAX
C
C distinguish the two cases LOGM and .NOT.LOGM
C
                        IF(LOGM) THEN
                           WBDS(IDX1+LMS)=CNMAX
                           ZRED(K)=1./CNMAX
                        ENDIF
                        IDXS(IDX1+LMS)=K
                        IDXN(IDX1+LMS)=K
                        DO 35 I2=2,NDM*IROT
                           IDX2=(LL2(I2,LL1(I,NCL),NCL)-1)*LMMAX
                           I3=MOD(M*(I2-1)+20*IROT,IROT)+1
                           I33=MOD(M*(2*IROT-I2+1)+20*IROT,IROT)+1
                           IF(I2.LE.IROT) THEN
C equivalent states under rotation
                             WBDS(IDX2+LMS)=WBDS(IDX1+LMS)*OM(I3)
                             IDXS(IDX2+LMS)=IDXS(IDX1+LMS)
                             IDXN(IDX2+LMS)=0
                           ELSE
C and mirror (IDX2 corresponds to Rotation*Mirror=Mirror*Rotation^(-1)
                               IF(IMIR.EQ.1) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)
                               ELSEIF(IMIR.EQ.2) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)*(-1)**M
                               ELSEIF(IMIR.EQ.3) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)*(-CI)**M
                               ELSEIF(IMIR.EQ.4) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)*(CI)**M
                               ELSE
                                  WRITE (1,1610)
                                  STOP
                               ENDIF
                             IDXS(IDX2+LM2S)=IDXS(IDX1+LMS)
                             IDXN(IDX2+LM2S)=0
                           ENDIF
35                      CONTINUE 
                      ENDIF
30               CONTINUE
         ELSEIF(LOGM.AND.LAFLG2(I,NCL).EQ.2*IROT) THEN
C
C the counting should give        K=K+(LMAX+1)*(LMAX+1)
C but we need weights
C
                 DO 40 LMS=1,LMMAX
                      L=INT(SQRT(FLOAT(LX(LMS)-1)+0.0001))
                      M=LX(LMS)-L-L*L-1
                      LM2S=LXM(L*L+L+1-M)
                      K=K+1
                      WBDS(IDX1+LMS)=CNM
                      ZRED(K)=1./CNM
                      IDXS(IDX1+LMS)=K
                      IDXN(IDX1+LMS)=K
C
C distinguish the two cases LOGAX and .NOT.LOGAX
C
                      IF(LOGAX) THEN
                           WBDS(IDX1+LMS)=CNMAX
                           ZRED(K)=1./CNMAX
                      ENDIF
                      DO 45 I2=2,2*IROT
                           IDX2=(LL2(I2,LL1(I,NCL),NCL)-1)*LMMAX
                           I3=MOD(M*(I2-1)+20*IROT,IROT)+1
                           I33=MOD(M*(2*IROT-I2+1)+20*IROT,IROT)+1
                           IF(.NOT.LOGAX) THEN
                             I3=1
                             I33=1
                           ENDIF
                           IF(I2.LE.IROT) THEN
C equivalent states under rotation (M can be negative, hence shift)
                             WBDS(IDX2+LMS)=WBDS(IDX1+LMS)*OM(I3)
                             IDXS(IDX2+LMS)=IDXS(IDX1+LMS)
                             IDXN(IDX2+LMS)=0
                           ELSE
C and mirror(as before)
                               IF(IMIR.EQ.1) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)
                               ELSEIF(IMIR.EQ.2) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)*(-1)**M
                               ELSEIF(IMIR.EQ.3) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)*(-CI)**M
                               ELSEIF(IMIR.EQ.4) THEN
                             WBDS(IDX2+LM2S)=WBDS(IDX1+LMS)
     &                          *OM(I33)*(CI)**M
                               ELSE
                                  WRITE (1,1610)
                                  STOP
                               ENDIF
                             IDXS(IDX2+LM2S)=IDXS(IDX1+LMS)
                             IDXN(IDX2+LM2S)=0
                           ENDIF
45                    CONTINUE 
40               CONTINUE
            ENDIF
23       CONTINUE
         LMNBD(IP)=K
20    CONTINUE
      IF(NSU.EQ.0.AND.NCL.GT.NSTEF) THEN
C
C no symmetry will be used for substrate calculation
C 
         K=0
         DO 10 I=1,NLAY
          DO 11 L=1,LMMAX
            K=K+1
            IDXN(K)=K
            IDXS(K)=K
            WBDS(K)=RU
            ZRED(K)=1.
11        CONTINUE
10       CONTINUE
         LMNBD(1)=LMMAX*NLAY
         RETURN
      ENDIF
      RETURN
      END
      SUBROUTINE TAUMAT3(TAU,LMT,NTAU,X,LEV,LEV2,LOD,TSF,LMMAX,LMAX,
     & FLMS,NL,KLM,LM,CLM,NLM,LXI,NT,PQ,NA,NLL,FLM,NTAUSH,TAUINV,
     & NLTU,LAN)
C
      DIMENSION CLM(NLM),LXI(LMMAX),PQ(2,NT),LAN(LMMAX)
      DIMENSION BR1(2),BR2(2),AR1(2),AR2(2),RAR1(2),RAR2(2)
      COMPLEX AK,CZ,TAU(LMT,LEV),X(LEV,LEV2),TSF(6,16),FLMS(NL,KLM)
      COMPLEX DET,FLM(KLM),CI,XA
      COMPLEX TAUINV(NLTU,LEV)
C      COMPLEX TCH(50,50)
C
      COMMON E,AK2,AK3,VPI
      COMMON /SL/BR1,BR2,AR1,AR2,RAR1,RAR2,NL1,NL2
C
      CZ=(0.0,0.0)
      CI=(0.0,1.0)
      AK=-0.5/CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      DO 100 K=1,KLM
         FLM(K)=CZ
100   CONTINUE
      BK2=PQ(1,1+NA)
      BK3=PQ(2,1+NA)
      JS=1
      S1=0.
      DO 130 J=1,NL1
         S2=0.
         DO 120 K=1,NL2
            ADR1=S1*BR1(1)+S2*BR2(1)
            ADR2=S1*BR1(2)+S2*BR2(2)
            ABR=ADR1*BK2+ADR2*BK3
            XA=CEXP(ABR*CI)
            DO 110 I=1,KLM
               FLM(I)=FLM(I)+FLMS(JS,I)*XA
110         CONTINUE
            IF (NLL.EQ.1) GOTO 144
            JS=JS+1
            S2=S2+1.
120      CONTINUE
         S1=S1+1.
130   CONTINUE
144   CONTINUE
      DO 141 IT=1,NTAU
         IT2=IT+NTAUSH
         DO 13 IL=1,2
            DO 142 I=1,LEV
               DO 1 J=1,LEV2
                  X(I,J)=CZ
1              CONTINUE
142         CONTINUE
            LL=LOD
            IF (IL.EQ.2) LL=LEV
C
C  GENERATE MATRIX 1-X FOR L+M= ODD (IL=1), LATER FOR L+M= EVEN (IL=2)
C  The output X is equal to 1-Gt (notice Gt and not tG if the
C  call is to XMT2) of Van Hove-Tong eq. 45
C
           CALL XMT2(IL,FLM,X,LEV,LL,TSF,IT2,LM,LXI,LMMAX,KLM,
     &       CLM,NLM,NST)
            LD2=(IT2-1)*LMMAX
            IF (IL.EQ.1) LD2=LD2+LEV
            DO 1422 I=1,LL
               IF(IL.EQ.1) LLL=LAN(I+LEV)+1
               IF(IL.EQ.2) LLL=LAN(I)+1
               DO 1100 J=1,LL
                     TAUINV(I+LD2,J)=X(I,J)
C                  IF (CABS(TSF(IT2,LLL)).GE.1.0E-06) THEN
C                     TAUINV(I+LD2,J)=X(I,J)/(AK*TSF(IT2,LLL))
C                  IF (CABS(TSF(IT2,LLL)).GE.1.0E-06.AND
C     &               .I.EQ.J) THEN
C                     TAUINV(I+LD2,J)=1./(AK*TSF(IT2,LLL))
C                  ELSE
C                     TAUINV(I+LD2,J)=CZ
C                  ENDIF
1100           CONTINUE
1422         CONTINUE
C
C  PREPARE QUANTITIES INTO WHICH INVERSE OF 1-X WILL BE MULTIPLIED
C
            IF (IL.LT.2) THEN
               IS=0
               LD1=0
               L=1
7              LD=LD1+1
               LD1=LD+L-1
               DO 8 I=LD,LD1
C notice that (1-tG)^-1 * t = t * (1-Gt)^-1
C for 1-tG
C                  X(I,I+LOD)=AK*TSF(IT2,L+1)
C for 1-Gt
                  X(I,I+LOD)=AK
8              CONTINUE
               L=L+2
               IF (L.LE.LMAX) GOTO 7
               IS=IS+1
               L=2
               IF (IS.LE.1) GOTO 7
            ELSE
               IS=0
               LD1=0
               L=0
3              LD=LD1+1
               LD1=LD+L
               DO 4 I=LD,LD1
C notice that (1-tG)^-1 * t = t * (1-Gt)^-1
C                  X(I,I+LEV)=AK*TSF(IT2,L+1)
                  X(I,I+LEV)=AK
4              CONTINUE
               L=L+2
               IF (L.LE.LMAX) GOTO 3
               IS=IS+1
               L=1
               IF (IS.LE.1) GOTO 3
            ENDIF
            LL2=LL+LL
C
C  PERFORM INVERSION AND MULTIPLICATION
C
            CALL CXMTXT(X,LEV,LL,LL,LL2,MARK,DET,-1)
            LD=(IT-1)*LMMAX
            IF (IL.EQ.1) LD=LD+LEV
            DO 143 I=1,LL
               DO 11 J=1,LL
C
C  PUT RESULT IN TAU 
C
C notice that (1-tG)^-1 * t = t * (1-Gt)^-1
C if X is (1-tG)^-1  TAU=X, if X=(1-Gt)^-1  TAU=tX
C                  TAU(LD+I,J)=X(I,J+LL)
                   IF(IL.EQ.1) LLL=LAN(I+LEV)+1
                   IF(IL.EQ.2) LLL=LAN(I)+1
                   TAU(LD+I,J)=X(I,J+LL)*TSF(IT2,LLL)
11             CONTINUE
143         CONTINUE
13       CONTINUE
141   CONTINUE
      RETURN
      END
      SUBROUTINE TAUT3(TAUG,TAUGM,LTAUG,LEV,CYLM,NT,LMMAX,LT,
     & NTAU,LOD,LEE,LOE,JGP,NTAUSH)
C
      COMPLEX CZ,CF,AK
      COMPLEX CYLM(NT,LMMAX)
      COMPLEX TAUG(LTAUG),TAUGM(LTAUG)
      DIMENSION LT(LMMAX)
      COMMON E,AK2,AK3,VPI
C
      CZ=(0.0,0.0)
      AK=-0.5/CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
C
C  PERFORM MATRIX PRODUCT TAU*YLM(G+-) FOR EACH CHEMICAL ELEMENT
C  for this modified version use TAU=1 instead 
C
      DO 275 I=1,NTAU
         IT2=I+NTAUSH
         IS=(I-1)*LMMAX
         DO 250 JLM=1,LEV
               KLP=LT(JLM)
               CF=CYLM(JGP,KLP)
               IF (JLM.GT.LEE) THEN
                  TAUG(IS+JLM)=-CF
                  TAUGM(IS+JLM)=-CF
               ELSE
                  TAUG(IS+JLM)=CF
                  TAUGM(IS+JLM)=CF
               ENDIF
250      CONTINUE
         IS=IS+LEV
         DO 270 JLM=1,LOD
               KLP=LT(JLM+LEV)
               CF=CYLM(JGP,KLP)
               IF (JLM.GT.LOE) THEN
                  TAUG(IS+JLM)=-CF
                  TAUGM(IS+JLM)=CF
               ELSE
                  TAUG(IS+JLM)=CF
                  TAUGM(IS+JLM)=-CF
               ENDIF
270      CONTINUE
275   CONTINUE
      RETURN
      END
C =========================================================================
C
C  Subroutine TFLSYM has the same function as TFOLD, but is designed for
C  the case of composite layers with symmetry. TFLSYM uses input from MTINVT.
C
C Parameter List;
C ===============
C
C TS       =  INPUT FROM MTINVT.
C LMN      =  NLAY*LMMAX.
C RG       =  PLANE-WAVE PROPAGATORS BETWEEN SUBPLANES.
C NLAY     =  NO. OF SUBPLANES.
C JGP      =  INDEX OF INCIDENT BEAM.
C EGS      =  WORKING SPACE.
C LXM      =  PERMUTATION OF (LM) SEQUENCE.
C INC      =  +-1  INDICATES INCIDENCE IN DIRECTION OF +-X.
C POSS     =  ATOMIC POSITIONS IN UNIT CELL.
C
C This is a modified version of routine MFOLT from the VAN HOVE/TONG
C LEED package. Modifications by BARBIERI.
C
C =========================================================================
C
      SUBROUTINE TFLSYM(TS,LMN,
     & RA,TA,CTR,CTT)
C
      COMPLEX CZ,RA,TA
      COMPLEX TS(LMN)
      COMPLEX CTR(LMN),CTT(LMN)
C
      COMMON /MFB/GP,LM,LAY
C 
      CZ=(0.0,0.0)
      RA=CZ
      TA=CZ
C
C  SUM OVER ATOMS and spherical indices IN UNIT CELL 
C
      DO 199 J=1,LMN
         RA=RA+CTR(J)*TS(J)
         TA=TA+CTT(J)*TS(J)
199   CONTINUE
      RETURN
      END
C======================================================================
C
C MATGEN generates the lattice matrix for routine BEMGEN.   
C
C Parameter List;
C ===============
C
C ARA1,ARA2       =   Substrate lattice vectors.
C ARB1,ARB2       =   Overlayer Lattice vectors.
C LATMAT          =   Matrix linking the overlayer to the substrate.
C
C AUTHOR: WANDER
C mod Barbieri
C======================================================================
C
      SUBROUTINE MATGN2(ARA1,ARA2,ARB1,ARB2,LATMAT)
C
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2)
      DIMENSION JAW(4),A(4)
      INTEGER LATMAT(2,2)
      COMPLEX WORK(4,4),CONT(4),VVST(4)
C
      EMACH=1E-06
C
C SET UP COEFFICIENT MATRIX WORK
C
      DO 149 I=1,4
         DO 101 J=1,4
            WORK(I,J)=CMPLX(0.0,0.0)
101      CONTINUE
149   CONTINUE
      DO 100 I=1,2
         WORK(I,1)=ARA1(I)
         WORK(I,2)=ARA2(I)
         WORK(I+2,3)=WORK(I,1)
         WORK(I+2,4)=WORK(I,2)
100   CONTINUE
C
C SET UP CONSTANT MATRIX
C
      DO 102 I=1,2
         CONT(I)=ARB1(I)
         CONT(I+2)=ARB2(I)
102   CONTINUE
C
C SOLVE SET OF EQUATIONS TO GENERATE ELEMENTS OF LATMAT
C
C      CALL ZGE(WORK,JAW,4,4,EMACH)
       CALL LUDCMP(WORK,4,4,JAW,D,VVST)
C      CALL ZSU(WORK,JAW,CONT,4,4,EMACH)
      CALL LUBKSB(WORK,4,4,JAW,CONT)

C
C COPY SOLUTIONS INTO CORRECT POSITION IN LATMAT. FIRST CHECK FOR
C ROUNDING ERRORS.
C
      DO 148 I=1,4
         A(I)=REAL(CONT(I))
         IF (A(I).GT.0) THEN
            A(I)=A(I)+0.1
         ELSE
            A(I)=A(I)-0.1
         ENDIF
148   CONTINUE
      LATMAT(1,1)=INT(A(1))
      LATMAT(1,2)=INT(A(2))
      LATMAT(2,1)=INT(A(3))
      LATMAT(2,2)=INT(A(4))
      RETURN
      END
C =========================================================================
C
C  Subroutine SRLAY2 is used by subroutine MTINV to reorder the subplanes
C  of a composite layer according to increasing position along the +X axis.
C  (Reordering the chemical element assignment LPS as well). If NEW=-1, it
C  finds which results fo GH will be available by simply copying old values
C  (of a previous call to MTINV), which ones need to be recomputed, and which
C  ones can be obtained by copying other recomputed ones.
C
C Parameter List;
C ===============
C
C POS         =    INPUT SUBPLANE POSITIONS (ATOMIC POSITIONS IN UNIT CELL).
C POSS        =    OUTPUT REORDERED POS.
C LPS         =    INPUT CHEMICAL ELEMENT ASSIGNMENT.
C LPSS        =    OUTPUT REORDERED LPS.
C MGH         =    OUTPUT  INDICATES ORGANIZATION OF GH IN MPERTI OR MTINV.
C MGH(I,J)    =    K MEANS GH(I,J) (I,J= SUBPLANE INDICES) IS TO BE FOUND IN
C                  K-TH POSITION IN COLUMNAR MATRIX GH.
C NLAY        =    NO. OF SUBPLANES.
C DRL         =    OUTPUT INTERPLANAR VECTORS.
C NLAY2       =    NLAY*(NLAY-1)/2.
C NUGH        =    OUTPUT  NUGH(K)=1 MEANS THE GH(I,J) (I,J= SUBPLANE INDICES)
C                  FOR WHICH MGH(I,J)=K MUST BE COMPUTED AFRESH. NUGH(K)=0 
C                  MEANS NO COMPUTATION NECESSARY FOR THAT GH(I,J).
C NEW         =    +-1 (INPUT) DETERMINES WHETHER OLD VALUES OF GH CAN BE 
C                  REUSED OR NOT.
C DCUT        =    INPUT RADIAL LIMIT OF LATTICE SUMMATION FOR GH, USED TO 
C                  DETECT VANISHING GH MATRICES WHEN DCUT IS SMALL. DCUT IS 
C                  APPLIED ONLY WHEN NEW=+1. WHEN NEW=-1 THE SAME GH MATRICES 
C                  ARE MADE TO VANISH AS IN THE PREVIOUS CALL TO MTINVT. THIS 
C                  ASSUMES THAT THE REORDERING OF SUBPLANES IS THE SAME AS IN 
C                  THE PREVIOUS CALL TO MTINVT. 
C LAY         =    1 FOR OVERLAYER SUPERLATTICE, 
C             =    2 FOR SUBSTRATE LATTICE (INPUT)
C
C In Common Blocks;
C =================
C
C ARA1,ARA2   =    SUBSTRATE LATTICE BASIS VECTORS
C ARB1,ARB2   =    SUPERLATTICE BASIS VECTORS
C
C Modified version of routine SRTLAY from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE SRLAY2(POS,POSS,LPS,LPSS,MGH,NLAY,DRL,NLAY2,NUGH,
     & DCUT,LAY,NIN,LL1,NLMAX)
C
      DIMENSION POS(NLMAX,3),POSS(NLMAX,3),POSA(3),DRL(NLAY2,3)
      DIMENSION AR1(2),AR2(2),LL1(NLMAX)
      DIMENSION LPS(NLMAX),LPSS(NLMAX),MGH(NLAY,NLAY),NUGH(NLAY2)
C
      COMMON /SL/ARA1(2),ARA2(2),ARB1(2),ARB2(2),RBR1(2),RBR2(2),NL1,
     & NL2
C
C
cjcm
c      write(*,*) 'SRLAY2 (JCM): possible trouble with common block SL'
cjcm
      D2=DCUT*DCUT
      DO 146 I=1,NLAY
         LPSS(I)=LPS(I)
         DO 1 J=1,3
            POSS(I,J)=POS(I,J)
1        CONTINUE
146   CONTINUE
C
C  ANALYSE ORDER OF SUBPLANE POSITIONS ALONG X-AXIS AND REORDER IN
C  ASCENDING POSITION ALONG +X-AXIS (EQUALLY POSITIONED SUBPLANES ARE
C  NOT PERMUTED)
C
      NLAY1=NLAY-1
      DO 7 I=1,NLAY1
         II=I+1
         KM=I
         PM=POSS(I,1)
         DO 2 K=II,NLAY
            IF (POSS(K,1).LT.PM) THEN
               PM=POSS(K,1)
               KM=K
            ENDIF
2        CONTINUE
         IF (KM.NE.I) THEN
            DO 3 J=1,3
               POSA(J)=POSS(KM,J)
3           CONTINUE
            LPSA=LPSS(KM)
            DO 5 KK=II,KM
               K=KM+II-KK
               DO 4 J=1,3
                  POSS(K,J)=POSS(K-1,J)
4              CONTINUE
C
C  REORDER CHEMICAL ASSIGNMENTS CORRESPONDINGLY
C
               LPSS(K)=LPSS(K-1)
5           CONTINUE
            DO 6 J=1,3
               POSS(I,J)=POSA(J)
6           CONTINUE
            LPSS(I)=LPSA
         ENDIF
7     CONTINUE
      IF (LAY.EQ.2) THEN
         DO 73 I=1,2
            AR1(I)=ARA1(I)
            AR2(I)=ARA2(I)
73       CONTINUE
      ELSE
         DO 71 I=1,2
            AR1(I)=ARB1(I)
            AR2(I)=ARB2(I)
71       CONTINUE
      ENDIF
C
C  GENERATE INTERPLANAR VECTORS DRL
C
      NN=1
      DO 100 NLR2=1,NIN
      DO 111 NLC=1,NLAY
         NLR=LL1(NLR2)
         IF(NLR.EQ.NLC) GOTO 111
         DO 9 K=1,3
            DRL(NN,K)=POSS(NLR,K)-POSS(NLC,K)
9        CONTINUE
         NUGH(NN)=1
         MGH(NLR,NLC)=NN
         NN=NN+1
111   CONTINUE
100   CONTINUE
      RETURN
      END
      SUBROUTINE XMT2(IL,FLM,X,LEV,LL,TSF,IT,LM,LXI,LMMAX,KLM,
     & CLM,NLM,N)
C
      COMPLEX TSF(6,16)
      complex acc
      COMPLEX X,FLM,CZERO,CI,ST,SU,RU
      DIMENSION CLM(NLM),X(LEV,LL),FLM(KLM),LXI(LMMAX)
C
      LMAX=LM-1
      RU=(1.0,0.0)
      CZERO=CMPLX(0.0,0.0)
      CI=CMPLX(0.0,-1.0)
      L2MAX=LMAX+LMAX
C
C  IF IL=1, CONSIDER L+M= ODD ONLY
C  IF IL=2, CONSIDER L+M= EVEN ONLY
C
      IF (IL.LT.2) THEN
         JSET=1
         MM=LEV
         N=1
      ELSE
         JSET=0
         MM=0
      ENDIF
      J=1
      L=JSET
470   M=-L+JSET
      JL=L+1
      MEXP=MOD(L,4)
      ST=CI**MEXP
480   K=1
      LPP=JSET
490   MPP=-LPP+JSET
      JLPP=LPP+1
      MEXP=MOD(LPP,4)
      SU=CI**MEXP/ST
500   MPA=IABS(MPP-M)
      LPA=IABS(LPP-L)
      IF (LPA.GT.MPA) MPA=LPA
      MP1=MPP-M+L2MAX+1
      LP1=L+LPP+1
      ACC=CZERO
530   JLM=(LP1*LP1+MP1-L2MAX)/2
      ACC=ACC+CLM(N)*FLM(JLM)
      N=N+1
      LP1=LP1-2
      IF (LP1-1-MPA.GE.0) GOTO 530
      JX=LXI(J+MM)
      KX=LXI(K+MM)
C
C X=1-tG
C      X(KX,JX)=-ACC*TSF(IT,JLPP)*SU
C X=1-Gt
      X(KX,JX)=-ACC*TSF(IT,JL)*SU
      IF (J.EQ.K) X(KX,JX)=X(KX,JX)+RU
      K=K+1
      MPP=MPP+2
      IF (LPP.GE.MPP) GOTO 500
      LPP=LPP+1
      IF (LMAX.GE.LPP) GOTO 490
      J=J+1
      M=M+2
      IF (L.GE.M) GOTO 480
      L=L+1
      IF (LMAX.GE.L) GOTO 470
      RETURN
      END
C  file LEEDSATL.SB3  Feb. 29, 1996
C
C**************************************************************************
C  Symmetrized Automated Tensor LEED (SATLEED):  subroutines, part 3
C  Version 4.1 of Automated Tensor LEED
C
C======================================================================  
C
C Function BLMT provides the integral of the product of three Legendre
C polynomials, each of which can be expressed as a prefactor times a 
C Legendre function. The three prefactors are grouped together as factor C
C and the integral is carried out following Gaunt's summation scheme set out
C in Slater 'Atomic Structure' Vol 1, P309-10.
C
C Parameter List;
C ===============
C
C  L1,M1,L2,M2,L3,M3   =  ANGULAR MOMENTUM VALUES FOR THE THREE LEGENDRE 
C                         FUNCTIONS
C  LMAX                =  MAXIMUM ANGULAR MOMENTUM VALUE USED IN THE 
C                         CALCULATION 
C
C This is a modified version of routine BLMT by PENDRY. Modifications
C by WANDER.
C                                                               
C======================================================================        
C
      FUNCTION BLMT(L1,M1,L2,M2,L3,M3,LMAX)
C
cjcm      DOUBLEPRECISION a,b,bn,c,cn
      REAL a,b,bn,c,cn
C                                ,blmt
C
40    FORMAT (' INVALID ARGUMENTS FOR BLMT ',6(I3,','))
C
      PI=3.14159265
      IF (M1+M2+M3.EQ.0) THEN
         IF (L1-LMAX-LMAX.LE.0) THEN
            IF (L2.LE.LMAX) THEN
               IF (L3.LE.LMAX) THEN
                  IF (L1.GE.IABS(M1)) THEN
                     IF (L2.GE.IABS(M2)) THEN
                        IF (L3.GE.IABS(M3)) THEN
                           IF (MOD(L1+L2+L3,2).NE.0) GOTO 420
                           NL1=L1
                           NL2=L2
                           NL3=L3
                           NM1=IABS(M1)
                           NM2=IABS(M2)
                           NM3=IABS(M3)
                           IC=(NM1+NM2+NM3)/2
                           IF (MAX0(NM1,NM2,NM3).GT.NM1) THEN
                              IF (MAX0(NM2,NM3).GT.NM2) THEN
                                 NL1=L3
                                 NL3=L1
                                 NM1=NM3
                                 NM3=IABS(M1)
                              ELSE
                                 NL1=L2
                                 NL2=L1
                                 NM1=NM2
                                 NM2=IABS(M1)
                              ENDIF
                           ENDIF
                           IF (NL2.LT.NL3) THEN
                              NTEMP=NL2
                              NL2=NL3
                              NL3=NTEMP
                              NTEMP=NM2
                              NM2=NM3
                              NM3=NTEMP
                           ENDIF
                           IF (NL3.GE.IABS(NL2-NL1)) THEN
C
C      CALCULATION OF FACTOR \A\
C
                              IS=(NL1+NL2+NL3)/2
                              IA1=IS-NL2-NM3
                              IA2=NL2+NM2
                              IA3=NL2-NM2
                              IA4=NL3+NM3
                              IA5=NL1+NL2-NL3
                              IA6=IS-NL1
                              IA7=IS-NL2
                              IA8=IS-NL3
                              IA9=NL1+NL2+NL3+1
                              A=((-1.0)**IA1)/FACT(IA3)*FACT(IA2)
     &                         /FACT(IA6)*FACT(IA4)
                              A=A/FACT(IA7)*FACT(IA5)/FACT(IA8)
     &                         *FACT(IS)/FACT(IA9)
                              A=A*(10.0**(IA2-IA3+IA4+IA5-IA6-IA7
     &                         -IA8-IA9+IS))
C
C      CALCULATION OF SUM \B\
C
                              IB1=NL1+NM1
                              IB2=NL2+NL3-NM1
                              IB3=NL1-NM1
                              IB4=NL2-NL3+NM1
                              IB5=NL3-NM3
                              IT1=MAX0(0,-IB4)+1
                              IT2=MIN0(IB2,IB3,IB5)+1
                              B=0.
                              SIGN=(-1.0)**(IT1)
                              IB1=IB1+IT1-2
                              IB2=IB2-IT1+2
                              IB3=IB3-IT1+2
                              IB4=IB4+IT1-2
                              IB5=IB5-IT1+2
                              DO 520 IT=IT1,IT2
                                 SIGN=-SIGN
                                 IB1=IB1+1
                                 IB2=IB2-1
                                 IB3=IB3-1
                                 IB4=IB4+1
                                 IB5=IB5-1
                                 BN=SIGN/FACT(IT-1)*FACT(IB1)
     &                            /FACT(IB3)*FACT(IB2)
                                 BN=BN/FACT(IB4)/FACT(IB5)
                                 BN=BN*(10.0**(IB1+IB2-IB3-IB4-IB5
     &                            -IT+1))
                                 B=B+BN
520                           CONTINUE
C
C      CALCULATION OF FACTOR \C\
C
                              IC1=NL1-NM1
                              IC2=NL1+NM1
                              IC3=NL2-NM2
                              IC4=NL2+NM2
                              IC5=NL3-NM3
                              IC6=NL3+NM3
                              CN=FLOAT((2*NL1+1)*(2*NL2+1)*(2*NL3+1))
     &                         /PI
                              C=CN/FACT(IC2)*FACT(IC1)/FACT(IC4)
     &                         *FACT(IC3)/FACT(IC6)*FACT(IC5)
                              C=C*(10.0**(IC1-IC2+IC3-IC4+IC5-IC6))
                              C=(SQRT(C))/2.
                              BLMT=((-1.0)**IC)*A*B*C
                           ELSE
                              BLMT=0.0
                           ENDIF
                           GOTO 531
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         WRITE (1,40) L1,M1,L2,M2,L3,M3
         GOTO 531
      ENDIF
420   BLMT=0.0
531   RETURN
      END
C============================================================================
C
C Function CA performs the same purpose for CAAA as BLMT does for CELMG.
C
C Parameter List;
C ===============
C
C  L1,MA1,L2,M2,L3,MA3   =  ANGULAR MOMENTUM VALUES
C
C Modified version of function CA from the Van Hove/Tong LEED package.
C Modifications by WANDER,BARBIERI.
C
C============================================================================
C
      FUNCTION CA(L1,MA1,L2,M2,L3,MA3)
C
cjcm      DOUBLEPRECISION yint,del,pref,prefa,prefb,preft,sum
      REAL yint,del,pref,prefa,prefb,preft,sum
C
      M1=-MA1
      M3=-MA3
      IF (.NOT.((IABS(M1).GT.L1).OR.(IABS(M2).GT.L2).OR.(IABS(M3).GT.
     & L3))) THEN
         IF ((M1+M2).EQ.(-M3)) THEN
            IF ((L3.GE.IABS(L1-L2)).AND.(L3.LE.(L1+L2))) THEN
               IF (MOD(L1+L2+L3,2).EQ.0) THEN
                  IF ((L1.EQ.0).OR.(L2.EQ.0).OR.(L3.EQ.0)) THEN
                     YINT=((2.0*FLOAT(L1)+1.0)*(2.0*FLOAT(L2)+1.0)*(2.0
     &                *FLOAT(L3)+1.0)/(12.56637))**(0.5)*(-1.0)**((IABS
     &                (M1)+IABS(M2)+IABS(M3))/2)/(FLOAT(IABS(L1)+IABS
     &                (L2)+IABS(L3))+1.0)
                  ELSE
                     DELT=FACT(L1+L2-L3)*FACT(L2+L3-L1)*FACT(L3+L1-L2)
     &                /FACT(L1+L2+L3+1)*((-1)**(L1-L2-M3+((L1+L2+L3)/2)
     &                ))
                     PREF=FACT(L1+M1)*FACT(L1-M1)*FACT(L2+M2)*FACT(L2
     &                -M2)*FACT(L3-M3)*FACT(L3+M3)
                     PREFA=SQRT(PREF)*(10.0**(L1+L2+L3))
                     PREFB=FACT((L1+L2+L3)/2)/(FACT((L1+L2-L3)/2)*FACT(
     &                (L2+L3-L1)/2)*FACT((L3+L1-L2)/2))*(0.1)
                     PREFT=DELT*PREFA*PREFB
                     IGT=MIN0(L2+M2,L1-M1,L1+L2-L3)
                     ILST=MAX0(L1+M2-L3,-M1+L2-L3,0)
                     NUM=IGT-ILST+1
                     SUM=0.0
                     DO 1560 IT=1,NUM
                        ITA=ILST+IT-1
                        SUM=SUM+((-1.0)**ITA)/(FACT(ITA)*FACT(L3-L2+ITA
     &                   +M1)*FACT(L3-L1+ITA-M2)*FACT(L1+L2-L3-ITA)
     &                   *FACT(L1-ITA-M1)*FACT(L2-ITA+M2))*(10.0**(-L1
     &                   -L2-L3))
1560                 CONTINUE
                     YINT=((2.0*FLOAT(L1)+1.0)*(2.0*FLOAT(L2)+1.0)*(2.0
     &                *FLOAT(L3)+1.0)/(12.56637))**(0.5)*PREFT*SUM
                  ENDIF
                  GOTO 1590
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      YINT=0.0
1590  CA=((-1.0)**(MA1+MA3))*YINT
      RETURN
      END
C============================================================================
C
C Routine CELMGT tabulates the Clebsch-Gordan type coefficients CLM and 
C ELM, for use by routines XM and XMT. The non-zero values are tabulated
C first for (L2+M2) and (L3+M3) odd and then for those quantities even.
C The same scheme is used by XM and XMT when the coefficients are retrieved.
C
C This is a modified version of routine CELMG by PENDRY. Modifications by
C WANDER,BARBIERI.
C
C============================================================================
C
      SUBROUTINE CELMG(CLM,NLM,YLM,FAC2,NN,FAC1,N,LMAX)
C
      DIMENSION CLM(NLM),YLM(NN),FAC2(NN),FAC1(N)
C      DOUBLEPRECISION FAC(100),blmt,a,b
cjcm       DOUBLEPRECISION FAC(100),a,b
       REAL FAC(100),a,b
C
      COMMON /F/FAC
C
      PI=3.14159265
      L2MAX=LMAX+LMAX
      NF=4*LMAX+1
      FAC(1)=1.0
      DO 340 I=1,NF
         FAC(I+1)=FLOAT(I)*FAC(I)
340   CONTINUE
C
COMMENT THE ARRAY YLM IS FIRST LOADED WITH SPHERICAL
C       HARMONICS, ARGUMENTS THETA=PI/2.0, FI=0.0
C
      LM=0
      CL=0.0
      A=1.0
      B=1.0
      ASG=1.0
      LL=L2MAX+1
C
COMMENT MULTIPLICATIVE FACTORS REQUIRED
C
      DO 240 L=1,LL
         FAC1(L)=ASG*SQRT((2.0*CL+1.0)*A/(4.0*PI*B*B))
         CM=-CL
         LN=L+L-1
         DO 230 M=1,LN
            LO=LM+M
            FAC2(LO)=SQRT((CL+1.0+CM)*(CL+1.0-CM)/((2.0*CL+3.0)*(2.0*CL
     &       +1.0)))
            CM=CM+1.0
230      CONTINUE
         CL=CL+1.0
         A=A*2.0*CL*(2.0*CL-1.0)/4.0
         B=B*CL
         ASG=-ASG
         LM=LM+LN
240   CONTINUE
C
COMMENT FIRST ALL THE YLM FOR M=+-L AND M=+-(L-1) ARE
C       ARE CALCULATED BY EXPLICIT FORMULAE
C
      LM=1
      CL=1.0
      ASG=-1.0
      YLM(1)=FAC1(1)
      DO 250 L=1,L2MAX
         LN=LM+L+L+1
         YLM(LN)=FAC1(L+1)
         YLM(LM+1)=ASG*FAC1(L+1)
         YLM(LN-1)=0.0
         YLM(LM+2)=0.0
         CL=CL+1.0
         ASG=-ASG
         LM=LN
250   CONTINUE
C
COMMENT USING YLM AND YL(M-1) IN A RECURRENCE RELATION
C       YL(M+1) IS CALCULATED
C
      LM=1
      LL=L2MAX-1
      DO 270 L=1,LL
         LN=L+L-1
         LM2=LM+LN+4
         LM3=LM-LN
         DO 260 M=1,LN
            LO=LM2+M
            LP=LM3+M
            LQ=LM+M+1
            YLM(LO)=-(FAC2(LP)*YLM(LP))/FAC2(LQ)
260      CONTINUE
         LM=LM+L+L+1
270   CONTINUE
      K=1
      II=0
280   LL=LMAX+II
      DO 341 IL2=1,LL
         L2=IL2-II
         M2=-L2+1-II
         DO 310 I2=1,IL2
            DO 342 IL3=1,LL
               L3=IL3-II
               M3=-L3+1-II
               DO 300 I3=1,IL3
                  LA1=MAX0(IABS(L2-L3),IABS(M2-M3))
                  LB1=L2+L3
                  LA11=LA1+1
                  LB11=LB1+1
                  M1=M3-M2
                  DO 290 L11=LA11,LB11,2
                     L1=LA11+LB11-L11-1
                     L=(L3-L2-L1)/2+M3
                     M=L1*(L1+1)-M1+1
                     ALM=((-1.0)**L)*4.0*PI*BLMT(L1,M1,L2,M2,L3,-M3,
     &                LMAX)
                     CLM(K)=YLM(M)*ALM
                     K=K+1
290               CONTINUE
                  M3=M3+2
300            CONTINUE
342         CONTINUE
            M2=M2+2
310      CONTINUE
341   CONTINUE
      IF (II.LE.0) THEN
         II=1
         GOTO 280
      ENDIF
      RETURN
      END
C=========================================================================
C
C Subroutine CPPP tabulates the function PPP(I1,I2,I3) for each
C element. It contains the integral of the product of three Legendre
C functions P(I1),P(I2),P(I3). The integrals are evluated using Gaunt's
C summation scheme given in Slater 'Atomic Structure' Vol 1, P309-10.
C PPP is used by routine PSTEMT in computing temperature-dependant phase
C shifts.
C
C This is a modified version of routine CPPP by PENDRY. Modifications by
C WANDER,BARBIERI
C
C
C============================================================================
C
      SUBROUTINE CPPP(PPP,N1,N2,N3)
C
      DIMENSION PPP(N1,N2,N3)
cjcm      DOUBLEPRECISION SUM,A
      REAL SUM,A
C
      DO 461 I1=1,N1
         DO 462 I2=1,N2
            DO 460 I3=1,N3
               IM1=I1
               IM2=I2
               IM3=I3
               IF (I1.LT.I2) THEN
                  IM1=I2
                  IM2=I1
               ENDIF
               IF (IM2.LT.I3) THEN
                  IM3=IM2
                  IM2=I3
                  IF (IM1.LT.IM2) THEN
                     J=IM1
                     IM1=IM2
                     IM2=J
                  ENDIF
               ENDIF
               A=0.0
               IS=I1+I2+I3-3
               IF (MOD(IS,2).NE.1) THEN
                  IF (IABS(IM2-IM1)-IM3+1.LE.0) THEN
                     SUM=0.0
                     IS=IS/2
                     SIGN=1.0
                     DO 450 IT=1,IM3
                        SIGN=-SIGN
                        IA1=IM1+IT-1
                        IA2=IM1-IT+1
                        IA3=IM3-IT+1
                        IA4=IM2+IM3-IT
                        IA5=IM2-IM3+IT
                        SUM=SUM-SIGN*(10.0**(IA1-IA2-IA3+IA4-IA5-IT
     &                   +2))*FACT(IA1-1)*FACT(IA4-1)/(FACT(IA2-1)
     &                   *FACT(IA3-1)*FACT(IA5-1)*FACT(IT-1))
450                  CONTINUE
                     IA1=2+IS-IM1
                     IA2=2+IS-IM2
                     IA3=2+IS-IM3
                     IA4=3+2*(IS-IM3)
                     A=-(-1.0)**(IS-IM2)*(10.0**(IA4+IS+IM3-IA1-IA2
     &                -IA3-2*IS))*FACT(IA4-1)*FACT(IS)*FACT(IM3-1)
     &                *SUM/(FACT(IA1-1)*FACT(IA2-1)*FACT(IA3-1)
     &                *FACT(2*IS+1))
                  ENDIF
               ENDIF
               PPP(I1,I2,I3)=A
460         CONTINUE
462      CONTINUE
461   CONTINUE
      RETURN
      END
C=======================================================================
C
C  Subroutine CXMTXT solves a set of linear equations.
C
C Parameter List;
C ===============
C
C  A     =  INPUT MATRIX. A IS DESTROYED DURING THE COMPUTATION AND
C           IS REPLACED BY THE UPPER TRIANGULAR MATRIX RESULTING FROM THE
C           GAUSSIAN ELIMINATION PROCESS (WITH PARTIAL PIVOTING).
C  NC    =  FIRST DIMENSION OF A IN CALLING ROUTINE. NC.GE.NR.
C  NR    =  ORDER OF A
C  NSYS  =  NO. OF SYSTEMS TO BE SOLVED. IF INVERSE OPTION IS CHOSEN
C           NSYS MUST BE AT LEAST AS LARGE AS NR. COEFFICIENT MATRIX MUST 
C           BE STORED IN A(I,J), I=1,NR, J=1,NR CONSTANT VECTORS MUST BE 
C           STORED IN A(I,NR+1), I=1,NR AND A(I,NR+2), I=1,NR  .... 
C           A(I,NR+NSYS), I=1,NR. RESULT OVERWRITES CONSTANT VECTORS.
C  NTOT  =  NR + NSYS
C  MARK  =  SINGULARITY INDICATOR (MARK=1 FOR SINGULAR A)
C  DET   =  DET(A)
C  INOPT =  -1 FOR SYSTEM SOLN. AND DET.
C            0 FOR DET ONLY.
C           +1 FOR INVERSE AND DET.
C
C NOTE  
C ====
C
C Dimension of X array must be at least as large as the first dimension of
C A array.
C
C Modified version of routine CXMTXT from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C=======================================================================
C
      SUBROUTINE CXMTXT(A,NC,NR,NSYS,NTOT,MARK,DET,INOPT)
C
      REAL AMAX,QZ
      COMPLEX A,X
      COMPLEX AGG,DET,CONST,SIGN,TEMP
      DIMENSION A(NC,NTOT),X(150)
C
310   FORMAT (//,1X,'SINGULAR MATRIX ',//)
C
C     PRESET PARAMETERS
C
      SIGN=(1.0E00,0.0E00)
      MARK=0
      IFLAG=INOPT
      N=NR
      NPL=N+1
      NMI=N-1
      NN=N+N
      NPLSY=N+NSYS
      IF (IFLAG.GT.0) THEN
C
C     INVERSE OPTION - PRESET AUGMENTED PART TO I
C
         DO 1161 I=1,N
            DO 930 J=NPL,NN
               A(I,J)=(0.0E00,0.0E00)
930         CONTINUE
1161     CONTINUE
         DO 940 I=1,N
            J=I+N
            A(I,J)=(1.0E00,0.0E00)
940      CONTINUE
         NPLSY=NN
      ENDIF
C
C     TRIANGULARIZE A
C
      DO 1162 I=1,NMI
         IPL=I+1
C
C     DETERMINE PIVOT ELEMENT
C
         MAX=I
         AMAX=CABS(A(I,I))
         DO 970 K=IPL,N
            QZ=CABS(A(K,I))
            IF (AMAX.LT.QZ) THEN
               MAX=K
               AMAX=CABS(A(K,I))
            ENDIF
970      CONTINUE
         IF (MAX.NE.I) THEN
C
C     PIVOTING NECESSARY - INTERCHANGE ROWS
C
            DO 990 L=I,NPLSY
               TEMP=A(I,L)
               A(I,L)=A(MAX,L)
               A(MAX,L)=TEMP
990         CONTINUE
            SIGN=-SIGN
         ENDIF
C
C     ELIMINATE A(I+1,I)---A(N,I)
C
         DO 1020 J=IPL,N
            TEMP=A(J,I)
            QZ=CABS(TEMP)
            IF (QZ.GE.1.0E-10) THEN
               CONST=-TEMP/A(I,I)
               DO 1010 L=I,NPLSY
                  A(J,L)=A(J,L)+A(I,L)*CONST
1010           CONTINUE
            ENDIF
1020     CONTINUE
1162  CONTINUE
C
C     COMPUTE VALUE OF DETERMINANT
C
      TEMP=(1.0E00,0.0E00)
      DO 1030 I=1,N
         AGG=A(I,I)
         QZ=CABS(AGG)
         IF (QZ.LE.1.0E-10) GOTO 1163
         TEMP=TEMP*AGG
1030  CONTINUE
      DET=SIGN*TEMP
      GOTO 1040
C
C     MATRIX SINGULAR
C
1163  WRITE (1,310)
      MARK=1
C
C     EXIT IF DET ONLY OPTION
C
1040  IF (IFLAG.NE.0) THEN
C
C     CHECK FOR INVERSE OPTION OR SYSTEMS OPTION
C
         IF (IFLAG.EQ.1) THEN
C
C     INVERSE OPTION - ABORT IF A IS SINGULAR
C
            IF (MARK.EQ.1) GOTO 1160
         ENDIF
C
C     BACK SUBSTITUTE TO OBTAIN INVERSE OR SYSTEM SOLUTION(S)
C
         DO 1150 I=NPL,NPLSY
            K=N
1080        X(K)=A(K,I)
            IF (K.NE.N) THEN
               DO 1100 J=KPL,N
                  X(K)=X(K)-A(K,J)*X(J)
1100           CONTINUE
            ENDIF
            X(K)=X(K)/A(K,K)
            IF (K.NE.1) THEN
               KPL=K
               K=K-1
               GOTO 1080
            ENDIF
C
C     PUT SOLN. VECT. INTO APPROPRIATE COLUMN OF A
C
            DO 1140 L=1,N
               A(L,I)=X(L)
1140        CONTINUE
1150     CONTINUE
      ENDIF
1160  RETURN
      END
C ========================================================================
C
C  Subroutine DEBWAT computes Debye-Waller factors for diffraction from
C  beam G(prime) into beam G for, in general, anisotropic atomic vibration,
C  including zero-temperature vibration.
C
C Parameter List;
C ===============
C
C  GP        =   INCIDENT BEAM G(PRIME).
C  E         =   ENERGY.
C  VPI       =   OPTICAL POTENTIAL
C  AK2,AK3   =   PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C  T         =   ACTUAL TEMPERATURE.
C  T0        =   REFERENCE TEMPERATURE FOR VIBRATION AMPLITUDES.
C  DRX       =   RMS VIBRATION AMPLITUDE PERPENDICULAR TO SURFACE.
C  DRY       =   RMS VIBRATION AMPLITUDE PARALLEL TO SURFACE.
C  D04       =   FOURTH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE
C                (ISOTROPIC).
C  EDW       =   OUTPUT DEBYE-WALLER FACTOR (EDW(1) FOR REFLECTION, EDW(2)
C                FOR TRANSMISSION).
C
C This ia a modified version of routine DEBWAL from the VAN HOVE/TONG
C LEED package. Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE DEBWAT(G,GP,E,VPI,AK2,AK3,T,T0,DRX,DRY,D04,EDW)
C
      COMPLEX CSQRT
      COMPLEX EDW
      DIMENSION GP(2),G(2),EDW(2)
C
      A1=GP(1)+AK2
      A2=GP(2)+AK3
      CC=REAL(CSQRT(CMPLX(2.0*E-A1*A1-A2*A2,-2.0*VPI)))
      A1=G(1)+AK2
      A2=G(2)+AK3
      DD=REAL(CSQRT(CMPLX(2.0*E-A1*A1-A2*A2,-2.0*VPI)))
C
C  C IS PERPENDICULAR COMPONENT OF SCATTERING VECTOR FOR REFLECTION,
C  D IS SAME FOR TRANSMISSION
C
      C=CC+DD
      D=CC-DD
      A1=GP(1)-G(1)
      A2=GP(2)-G(2)
C
C  D2 IS MEAN-SQUARE VIBRATION AMPLITUDE PARALLEL TO SURFACE AT ACTUAL
C  TEMPERATURE
C
      D2=DRY*DRY*T/T0
C
C  ZERO-TEMPERATURE VIBRATION IS NOW MIXED IN
C
      D2=SQRT(D2*D2+D04)
      A1=(A1*A1+A2*A2)*D2
C
C  D2 IS SAME AS ABOVE, BUT FOR PERPENDICULAR COMPONENTS
C
      D2=DRX*DRX*T/T0
      D2=SQRT(D2*D2+D04)
      EDW(1)=CMPLX(EXP(-0.166667*(A1+C*C*D2)),0.0)
      EDW(2)=CMPLX(EXP(-0.166667*(A1+D*D*D2)),0.0)
      RETURN
      END
C=====================================================================
C
C Subroutine DUMP dumps the complex amplitude for each exit beam
C to the transfer file.
C 
C Parameter List;
C ===============
C
C PQF     =  LIST OF ALL BEAMS IN UNITS OF RECIPROCAL LATTICE VECTORS
C PQFEX   =  LIST OF TLEED EXIT BEAMS
C NT0     =  NUMBER OF EXIT BEAMS
C NT      =  NUMBER OF BEAMS IN CALCULATION AT CURRENT ENERGY
C XI      =  AMPLITUDES OF EACH LEED BEAM IN THE SAME ORDER AS PQF
C E       =  CURRENT ENERGY
C VV      =  INNER POTENTIAL
C VPIS    =  IMAGINARY PART OF INNER POTENTIAL
C NFILE   =  OUTPUT CHANNEL NUMBER
C
C Modified version of routines XIDUMP and DUMP1 by ROUS. Modifications by
C WANDER.
C
C =========================================================================
C
      SUBROUTINE DUMP(PQF,PQFEX,NT0,NT,XI,XIST,E,NFILE)
C
      DIMENSION PQF(2,NT),PQFEX(2,NT0)
      COMPLEX XI(NT),XIST(NT0)
C
!1000  FORMAT (3F7.4)
!1010  FORMAT (30(2E13.5,/))
C
C FIRST REORDER THE AMPLITUDES STORED IN XIST FROM THE ORDER OF THE LIST
C PQF TO THE ORDER OF THE LIST PQFEX.
C
      DO 100 IG=1,NT0
         XIST(IG)=CMPLX(0.0,0.0)
         DO 200 JG=1,NT
            P1=ABS(PQF(1,JG)-PQFEX(1,IG))
            P2=ABS(PQF(2,JG)-PQFEX(2,IG))
            IF ((P1+P2).LT.0.01) XIST(IG)=XI(JG)
200      CONTINUE
100   CONTINUE
C
C NOW DUMP THE ENERGY AND AMPLITUDES TO THE TRANSFER FILE
C
C      WRITE (NFILE,1000) E
C      WRITE (NFILE,1010) (XIST(I),I=1,NT0)
      WRITE (NFILE) E
      WRITE (NFILE) (XIST(I),I=1,NT0)
      RETURN
      END
C=========================================================================
C
C  Function FACT computes the factorial using an asymptotic expansion.
C  Modified version of function FACT from the Van Hove/Tong LEED package.
C  Modifications by WANDER.
C
C=========================================================================
C
      FUNCTION FACT(L)
!     jcm why the heck are we writing and using our own factorial program!!!!
!     does this even compute the correct value? L=1,2,3 don't seem right
!     jcm      DOUBLEPRECISION DFACT,X
      
      REAL DFACT,X
C
      IF (L.GT.4) THEN
         X=L+1
         DFACT=EXP(-X)*(10.0**((X-0.5)*LOG10(X)-(X-1.0)))*
     &    (SQRT(6.283185307179586))*(1.0+(1.0/(12.0*X))+(1.0/
     &    (288.0*(X**2)))-(139.0/(51840.0*(X**3)))-(571.0/
     &    (2488320.0*(X**4))))
!         FACT=SNGL(DFACT)
         FACT=DFACT
      ELSE
         IF (L.EQ.0) FACT=1.0
         IF (L.EQ.1) FACT=0.1
         IF (L.EQ.2) FACT=0.02
         IF (L.EQ.3) FACT=6.0*0.001
         IF (L.EQ.4) FACT=24.0*0.0001
      ENDIF
!jcm      write(*,*), 'FACT: ***** L, FACT = ', L, FACT
      RETURN
      END
C=========================================================================
C
C Subroutine FMAT calculates the values of the sum FLMS(JS,LM) over lattice
C points of each sublattice JS, where LM=(0,0),(1,-1),(1,1).......
C Note that for (L+M) odd, FLMS=0. Dimensions 96 and 33 are set for LMAX.LE.8
C
C Parameter List;
C ===============
C
C  FLMS         =   OUTPUT LATTICE SUMS.
C  V,JJS        =   INPUT FROM SUBROUTINE SLIND.
C  NL           =   NO. OF SUBLATTICES
C  NLS          =   ACTUAL NO. OF SUBLATTICE SUMS DESIRED (E.G. 1 OR NL).
C  DCUT         =   CUTOFF DISTANCE FOR LATTICE SUM.
C  IDEG         =   DEGREE OF SYMMETRY OF LATTICE (IDEG-FOLD ROTATION AXIS).
C                   (NB  Do not use IDEG=1. IDEG=3 preferable over IDEG=6.)
C  LMAX         =   LARGEST VALUE OF L.
C  KLM          =   (2*LMAX+1)*(2*LMAX+2)/2.
C
C In Common Blocks;
C =================
C
C  E            =   CURRENT ENERGY.
C  AK               =   PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C  VPI          =   IMAGINARY PART OF ENERGY.
C  AR1,AR2      =   BASIS VECTORS OF SUPERLATTICE.
C  BR1,BR2      =   NOT USED.
C  RAR1,RAR2    =   NOT USED.
C  NL1,NL2      =   NOT USED.
C
C Note;
C =====
C
C Dimensions 180 (2*LMAX*IDEGREE) and 61 (4*lmax+1) require (LMAX.LE.15)
C if IEGREE=6.
C
C Modified version of PENDRY'S routine FMAT. Modifications by WANDER,BARBIERI.
C
C=========================================================================
C
      SUBROUTINE FMAT(FLMS,V,JJS,NL,NLS,DCUT,IDEG,LMAX,KLM)
C
      COMPLEX FLMS,SCC,SA,RTAB,CZERO,CI,KAPPA,SC,SD,SE,Z,ACC,RF
      complex acs
      COMPLEX CSQRT,CEXP,CCOS,CSIN
      DIMENSION FLMS(NL,KLM),V(NL,2),JJS(NL,IDEG),BR1(2)
      DIMENSION BR2(2),SCC(6,61),SA(180),ANC(6),ANS(6),RTAB(4)
      DIMENSION AK(2),AR1(2),AR2(2),RAR1(2),RAR2(2),R(2)
C
      COMMON E,AK,VPI
      COMMON /SL/BR1,BR2,AR1,AR2,RAR1,RAR2,NL1,NL2
C
      PI=3.14159265
      CZERO=CMPLX(0.0,0.0)
      CI=CMPLX(0.0,1.0)
      KAPPA=CMPLX(2.0*E,-2.0*VPI+0.000001)
      KAPPA=CSQRT(KAPPA)
      AG=SQRT(AK(1)*AK(1)+AK(2)*AK(2))
C
COMMENT ANC,ANS AND SA ARE PREPARED TO BE USED IN THE SUM
C       OVER SYMMETRICALLY RELATED SECTORS OF THE LATTICE
C
      L2MAX=LMAX+LMAX
      LIM=L2MAX*IDEG
      LIML=L2MAX+1
      ANG=2.0*PI/FLOAT(IDEG)
      D=1.0
      DO 10 J=1,IDEG
         ANC(J)=COS(D*ANG)
         ANS(J)=SIN(D*ANG)
         D=D+1.0
10    CONTINUE
      D=1.0
      DO 20 J=1,LIM
         SA(J)=CEXP(-CI*D*ANG)
         D=D+1.0
20    CONTINUE
      DO 161 J=1,NL
         DO 30 K=1,KLM
            FLMS(J,K)=CZERO
30       CONTINUE
161   CONTINUE
C
COMMENT THE LATTICE SUM STARTS.THE SUM IS DIVIDED INTO ONE
C       OVER A SINGLE SECTOR,THE OTHER (IDEG-1) SECTORS
C       ARE RELATED BY SYMMETRY EXCEPT FOR FACTORS
C       INVOLVING THE DIRECTION OF R
C  THE RANGE OF SUMMATION IS LIMITED BY DCUT
C
      D=SQRT(AR1(1)*AR1(1)+AR1(2)*AR1(2))
      LI1=INT(DCUT/D)+2
      D=SQRT(AR2(1)*AR2(1)+AR2(2)*AR2(2))
      LI2=INT(DCUT/D)+2
      DCUT2=DCUT*DCUT
C
C  ONE SUBLATTICE AT A TIME IS TREATED IN THE FIRST SECTOR
C
      DO 160 JS=1,NLS
         LI11=LI1
         LI22=LI2
         ASST=0.0
         ADD=1.0
         ANT=-1.0
         ADR1=V(JS,1)*COS(V(JS,2))
         ADR2=V(JS,1)*SIN(V(JS,2))
C
C  SHIFT POINT ADR1,2 BY MULTIPLES OF AR1 AND AR2
C  INTO SUPERLATTICE UNIT CELL NEAR ORIGIN, DEFINED BY THE LIMITS
C  (-0.001 TO 0.999)*AR1 AND (0.001 TO 1.001)*AR2
C
         DET=AR1(1)*AR2(2)-AR1(2)*AR2(1)
         B1=(ADR1*AR2(2)-ADR2*AR2(1))/DET
         B2=(AR1(1)*ADR2-AR1(2)*ADR1)/DET
         BP1=AMOD(B1,1.)
         IF (BP1.LT.-.001) BP1=BP1+1.
         IF (BP1.GT.0.999) BP1=BP1-1.
         BP2=AMOD(B2,1.)
         IF (BP2.LT.0.001) BP2=BP2+1.
         IF (BP2.GT.1.001) BP2=BP2-1.
         ADR1=BP1*AR1(1)+BP2*AR2(1)
         ADR2=BP1*AR1(2)+BP2*AR2(2)
         AST=-1.0
60       AN1=ANT
         DO 162 I1=1,LI11
            AN1=AN1+ADD
            AN2=AST
            DO 130 I2=1,LI22
               AN2=AN2+1.0
C
COMMENT R=THE CURRENT LATTICE VECTOR IN THE SUM
C       AR=MOD(R)
C       RTAB(1)=-EXP(I*FI(R))
C
               R(1)=AN1*AR1(1)+AN2*AR2(1)+ADR1
               R(2)=AN1*AR1(2)+AN2*AR2(2)+ADR2
               AR=R(1)*R(1)+R(2)*R(2)
               IF (AR.LE.DCUT2) THEN
                  AR=SQRT(AR)
                  RTAB(1)=-CMPLX(R(1)/AR,R(2)/AR)
                  ABC=1.0
                  ABB=0.0
                  IF (AG.GT.1.0E-4) THEN
                     ABC=(AK(1)*R(1)+AK(2)*R(2))/(AG*AR)
                     ABB=(-AK(2)*R(1)+AK(1)*R(2))/(AG*AR)
                  ENDIF
                  SC=CI*AG*AR
C
COMMENT SCC CONTAINS FACTORS IN THE SUMMATION DEPENDENT ON
C       THE DIRECTION OF R. CONTRIBUTIONS FROM SYMMETRICALLY
C       RELATED SECTORS CAN BE GENERATED SIMPLY AND ARE
C       ACCUMULATED FOR EACH SECTOR, INDEXED BY THE SUBSCRIPT J.
C       THE SUBSCRIPT M IS ORDERED M=(-L2MAX),(-L2MAX+1)....
C       (+L2MAX)
C
                  DO 163 J=1,IDEG
                     AD=ABC*ANC(J)-ABB*ANS(J)
                     SD=CEXP(SC*AD)
                     SCC(J,LIML)=SD
                     MJ=0
                     SE=RTAB(1)
                     DO 90 M=1,L2MAX
                        MJ=MJ+J
                        MP=LIML+M
                        MM=LIML-M
                        SCC(J,MP)=SD*SA(MJ)/SE
                        SCC(J,MM)=SD*SE/SA(MJ)
                        SE=SE*RTAB(1)
90                   CONTINUE
163               CONTINUE
                  Z=AR*KAPPA
                  ACS=CSIN(Z)
                  ACC=CCOS(Z)
C
COMMENT RTAB(3)=SPHERICAL HANKEL FUNCTION OF THE FIRST KIND,L=0
C       RTAB(4)=SPHERICAL HANKEL FUNCTION OF THE FIRST KIND,L=1
C
                  RTAB(3)=(ACS-CI*ACC)/Z
                  RTAB(4)=((ACS/Z-ACC)-CI*(ACC/Z+ACS))/Z
                  AL=0.0
C
COMMENT THE SUMMATION OVER FACTORS INDEPENDENT OF THE
C       DIRECTION OF R IS ACCUMULATED IN FLM, FOR EACH
C       SUBLATTICE INDEXED BY SUBSCRIPT JSP.  THE SECOND
C       SUBSCRIPT ORDERS L AND M AS  (0,0),(1,-1),(1,1),(2,-2),
C       (2,0),(2,2)...
C
                  JF=1
                  DO 120 JL=1,LIML
                     RF=RTAB(3)*CI
                     JM=L2MAX+2-JL
                     DO 110 KM=1,JL
C
C  CONSIDER THE CORRESPONDING LATTICE POINTS IN THE OTHER SECTORS AND
C  GIVE THEIR CONTRIBUTION TO THE APPROPRIATE SUBLATTICE
C
                        DO 100 J=1,IDEG
                           JSP=JJS(JS,J)
                           FLMS(JSP,JF)=FLMS(JSP,JF)+SCC(J,JM)*RF
100                     CONTINUE
                        JF=JF+1
                        JM=JM+2
110                  CONTINUE
C
COMMENT SPHERICAL HANKEL FUNCTIONS FOR HIGHER L ARE
C       GENERATED BY RECURRENCE RELATIONS
C
                     ACS=(2.0*AL+3.0)*RTAB(4)/Z-RTAB(3)
                     RTAB(3)=RTAB(4)
                     RTAB(4)=ACS
                     AL=AL+1.0
120               CONTINUE
               ENDIF
130         CONTINUE
162      CONTINUE
C
COMMENT SPECIAL TREATMENT IS REQUIRED IF IDEG=2
C  TWO SECTORS REMAIN TO BE SUMMED OVER
C
         IF (IDEG.LE.2) THEN
            IF (ASST.LE.0) THEN
               ASST=1.0
               ADD=-1.0
               AST=-1.0
               IF (BP2.GT.0.999) AST=-2.
               ANT=0.0
               GOTO 60
            ENDIF
         ENDIF
160   CONTINUE
      RETURN
      END
C ============================================================================
C Subroutine FORSPLINE prepare for the phase shift interpolation to be 
C performed in TSCATF through a cubic spline interpolation routine
C
C input: 
C NEL         = NUMBER OF CHEMICAL ELEMENT
C L1          = LMAX+1.
C ES          = LIST OF ENERGIES AT WHICH PHASE SHIFTS ARE TABULATED.
C PHSS        = TABULATED PHASE SHIFTS.
C NPSI        = NO. OF ENERGIES AT WHICH PHASE SHIFTS ARE GIVEN.
C output:
C PHSS2       = second derivative of the tabulated phase shifts(as a function
C               of the energy) needed for cubic spline interpolation
C
C Author: Barbieri
C==============================================================================
      SUBROUTINE FORSPLINE(NEL,L1,NPSI,PHSS,ES,PHSS2)
      DIMENSION PHSS(NPSI,80),PHSS2(NPSI,80),ES(NPSI)
      DIMENSION PHSL(90),PHSL2(90)
C
C natural spline
      YP1=1.E10
      YPN=1.E10
C
      DO 300 IEL=1,NEL
         IO=(IEL-1)*L1
         DO 330 L=1,L1
            DO 360 I=1,NPSI
               PHSL(I)=PHSS(I,L+IO)
360         CONTINUE
            CALL SPLINE(ES,PHSL,NPSI,YP1,YPN,PHSL2)
            DO 370 I=1,NPSI
               PHSS2(I,L+IO)=PHSL2(I)
370         CONTINUE
330      CONTINUE
300    CONTINUE
       RETURN
       END
C======================================================================
C                                                                       
C Subroutine GAUNT generated a vector of Gaunt coefficients BLM in the
C same order as they are used in subroutine CVEC. Only non-zero elements
C are stored, and are selected by observing the usual selection rules;
C
C C(l1,m1; l3,m3; lp,mp)     = 0    Unless;
C                                 m1+m3-mp=0  (ie m1=mp-m3)
C                                l1+l3+lp is even
C
C Parameter List;
C ===============
C
C NLMB          =   DIMENSION OF BELM (SEE BELOW)
C LMAX          =   LARGEST L-VALUE USED TO DESCRIBE ATOMIC
C                   SCATTERING WITHIN THE REFERENCE SURFACE.
C LSMAX         =   LARGEST L-VALUE FOR SINGLE CENTRE EXPANSION
C                   (N.B. LSMAX<=2*LMAX.)
C BELM(NLMB)    =   VECTOR OF GAUNT COEFFS.
C
C Table of the dimension of BELM : NLMB against the number of L
C values used for atomic scattering LMAX and the number used for
C the single centre expansion of the Q matrix.
C
C                                 LSMAX                                 
C                                                                       
C      |    0    1    2   3    4     5     6     7     8     9     10   
C LMAX +--------------------------------------------------------------- 
C      |                                                                
C  0   |    1                                                           
C  1   |    4   10   19                                                 
C  2   |    9   33   71  101   126                                      
C  3   |   16   70  167  269   379   449   498                          
C  4   |   25  121  307  537   804  1042  1256  1382  1463              
C  5   |   36  186  491  905  1419  1935  2458  2880  3230  3428  3549  
C  6   |   49  265  719 1373  2224  3150  4138  5054  5905  6559  7077  
C  7   |   64  358  991 1941  3219  4687  6322  7942  9554 10966 12217  
C  8   |   81  465 1307 2609  4404  6546  9010
C  9   |  100  586 1667 3377  5779  8599
C 10   |           2071
C 11   |           2519
C 12   |           3011
C      |                                                                
C
C =========================================================================
C
      SUBROUTINE GAUNT(BELM,NLMB,LMAX,LSMAX)
C
      DIMENSION BELM(NLMB)
C      
      K=0
      LSMMAX=(LSMAX+1)*(LSMAX+1)
      LMMAX=(LMAX+1)*(LMAX+1)
      DO 12 I1=1,LMMAX
         L1=INT(SQRT(FLOAT(I1-1)))
         M1=I1-L1*L1-L1-1
         DO 1 I3=1,LSMMAX
            L3=INT(SQRT(FLOAT(I3-1)))
            M3=I3-L3*L3-L3-1
            MP=M1+M3
            LL1=MAX0(IABS(MP),IABS(L1-L3))
            LL2=MIN0(LMAX,L1+L3)
            IF (MOD(LL2+L1+L3,2).NE.0) LL2=LL2-1
            DO 2 LP=LL2,LL1,-2
               K=K+1
               BELM(K)=BLMT(L3,M3,L1,M1,LP,-MP,LMAX)
2           CONTINUE
1        CONTINUE
12    CONTINUE
      RETURN
      END
C=========================================================================
C
C  Subroutine GEOMVT produces cartesian coordinates for a given structure.
C
C  The input is through READCT and is described below. It consists of
C  a set of basic atomic positions given by bond lengths, bond angles
C  and dihedral angles.
C
C  Description of input and method for molecular structures -
C  First, NATOMS atom positions are read in. Each atom position is
C  to be given as
C    IANZ,IZ(1),BLS,IZ(2),ALPHAS,IZ(3),BETAS,IZ(4)  which means -
C  This atom (dummy if IANZ=0) has a bond length BLS (Angstroms) to
C  atom no. IZ(1) (Atoms are numbered in the entry sequence), makes a
C  bond angle of ALPHAS(degrees) in the sequence Atom IZ(2)/Atom IZ(1)/
C  this atom, and makes a dihedral angle BETAS(degrees) in the sequence
C  Atom IZ(3)/Atom IZ(2)/Atom IZ(1)/This atom.
C  (BETAS is the angle through ehich the plane of atoms IZ(1)/IZ(2)/
C  IZ(3) is to be rotated to coincide with the plane of atoms
C  IZ(2)/IZ(3)/This atom, as viewed along the bond from Atom IZ(2) to
C  Atom IZ(3)). The input field IZ(4) is to be left blank. Dummy atoms
C  are only used to help define others.
C  Three dummy atoms are implicitly provided in the program initially;
C    1  (X,Y,Z)=(0,0,0) Local origin
C    2  (X,Y,Z)=(1,0,0)
C    3  (X,Y,Z)=(0,1,0)
C
C  Parameters not defined above-   
C   IP=1 Produces printout of molecular structure.                     
C
C Modified version of routine GEOMV from the Van Hove/Tong LEED package.
C Modifications by WANDER.
C
C=========================================================================
C
      SUBROUTINE GEOMV(IP,VPOS,NLAY,NLMX)
C
      DIMENSION VPOS(NLMX,3)
C
      COMMON /ZMAT/IANZ(40),IZ(40,4),BL(40),ALPHA(40),BETA(40),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(40),ALPHAS(40),BETAS
     & (40),PHIR,PHIM1,PHIM2
C
1000  FORMAT (//,
     & ' PARAMETER OUT OF RANGE OF DEFINED Z MATRIX  $$$ STOP $$$')
1001  FORMAT ('NATOMS(=',I4,') .NE.NLAY(=',I4,') IN GEOMV  *********')
C
      DO 111 I=1,3
         BL(I)=.0
         ALPHA(I)=.0
         IF (I.EQ.2) BL(I)=1.
         IF (I.EQ.3) THEN
           BL(I)=1.
           ALPHA(I)=90.
         ENDIF
111   CONTINUE
      DO 5 K=4,NZ
         BL(K)=BLS(K)
         ALPHA(K)=ALPHAS(K)
         BETA(K)=BETAS(K)
5     CONTINUE
      DO 1004 J=1,NUM
         NPARR=NIPAR(J)
         DO 60 I=1,NPARR
            K=IPAR(I,J)
            IF (K.LE.NZ) THEN
               BL(K)=BLS(K)
            ELSEIF (K-2*NZ.LE.0) THEN
               ALPHA(K-NZ)=ALPHAS(K-NZ)
            ELSEIF (K-3*NZ.GT.0) THEN
               GOTO 70
            ELSE
               BETA(K-2*NZ)=BETAS(K-2*NZ)
            ENDIF
60       CONTINUE
1004  CONTINUE
      CALL GEOMXY(IP,VPOS,NLMX)
      IF (NATOMS.NE.NLAY) WRITE (1,1001) NATOMS,NLAY
      GOTO 1005
70    WRITE (1,1000)
      STOP
1005  RETURN
      END
C=========================================================================
C
C  Subroutine GEOMVP produces VP=X cross Y
C
C=========================================================================
C
      SUBROUTINE GEOMVP(VP,X,Y)
C
      DIMENSION VP(3),X(3),Y(3)
C
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
C=========================================================================
C
C  Subroutine GEOMVU produces U which is a unit vector along c(j)-c(k)
C
C=========================================================================
C
      SUBROUTINE GEOMVU(U,C,J,K)
C
      DIMENSION C(40,3),R(3),U(3)
C
      DATA ZERO/0.0/
C
      R2=ZERO
      DO 10 I=1,3
         R(I)=C(J,I)-C(K,I)
         R2=R2+R(I)*R(I)
10    CONTINUE
      R2=SQRT(R2)
      DO 20 I=1,3
         U(I)=R(I)/R2
20    CONTINUE
      RETURN
      END
C=========================================================================
C
C  Subroutine GEOMXY produces cartesian coordinates for a given composite
C  layer geometry.
C
C=========================================================================
C
      SUBROUTINE GEOMXY(IP,C,NLMX)
C
      DIMENSION ALPHA(40),BETA(40)
      DIMENSION A(40),B(40),CZ(40,3),D(40),U1(3),U2(3),U3(3),U4(3)
      DIMENSION VJ(3),VP(3),V3(3)
      DIMENSION C(NLMX,3),DIS(40)
C
      COMMON /ZMAT/IANZ(40),IZ(40,4),BL(40),ALPH(40),BET(40),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(40),ALPHAS(40),BETAS
     & (40),PHIR,PHIM1,PHIM2
C
      DATA ZERO/0.0/,ONE/1.0/,TWO/2.0/
      DATA TENM5/1.0E-5/,TENM6/1.0E-6/
      DATA TORAD/1.74532925199E-02/
      DATA PI/3.141592654/
C
1000  FORMAT (//,15X,'COORDINATES',/,13X,'X',14X,'Y',14X,'Z')
1010  FORMAT (1X,1I4,2X,F10.5,2(5X,F10.5))
1020  FORMAT (//,30X,'Z MATRIX'//1X,'IN IN1 AN',3X,'Z1',4X,'BL',19X,
     & 'Z2',2X,'ALPHA',17X,'Z3',3X,'BETA',17X,'Z4',/)
1030  FORMAT (3I3)
1040  FORMAT (3I3,I5,G14.7,'(',I3,')')
1050  FORMAT (3I3,I5,G14.7,'(',I3,')   ',I5,G14.7,'(',I3,')')
1060  FORMAT (3I3,I5,G14.7,'(',I3,')   ',I5,G14.7,'(',I3,')   ',I5,
     & G14.7,'(',I3,')   ',I5)
5000  FORMAT (11X,15I7)
5003  FORMAT (5X,I3,5X,15F7.3)
 5001 FORMAT (//5X,'DISTANCE MATRIX',/)
C
      DO 5004 I=1,3
         DO 510 J=1,4
            IZ(I,J)=0
510      CONTINUE
5004  CONTINUE
      IZ(2,1)=1
      IZ(3,1)=1
      IZ(3,2)=2
C
C     PRINT Z MATRIX
C
      IF (IP.EQ.1) THEN
         WRITE (1,1020)
         IN=1
         IN1=0
         WRITE (1,1030) IN,IN1,IANZ(1)
         IF (NZ.GT.1) THEN
            I1=2
            IN=IN+1
            WRITE (1,1040) IN,IN1,IANZ(2),IZ(2,1),BL(2),I1
            IF (NZ.GT.2) THEN
               I1=I1+1
               I2=NZ+I1
               IN=IN+1
               WRITE (1,1050) IN,IN1,IANZ(3),IZ(3,1),BL(3),I1,IZ(3,2),
     &          ALPH(3),I2
               IF (NZ.GT.3) THEN
                  I3=2*NZ+I1
                  DO 4 I=4,NZ
                     I1=I1+1
                     I2=I2+1
                     I3=I3+1
                     IN=IN+1
                     IN1=IN-3
                     WRITE (1,1060) IN,IN1,IANZ(I),IZ(I,1),BL(I),I1,IZ
     &                (I,2),ALPH(I),I2,IZ(I,3),BET(I),I3,IZ(I,4)
4                 CONTINUE
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
C     ZERO TEMPORARY COORDINATE ARRAY CZ
C
      DO 5005 I=1,NZ
         DO 10 J=1,3
            CZ(I,J)=ZERO
10       CONTINUE
5005  CONTINUE
C
C     CONVERT ANGLES FROM DEGREES TO RADIANS
C
      DO 20 I=1,NZ
         ALPHA(I)=ALPH(I)*TORAD
         BETA(I)=BET(I)*TORAD
20    CONTINUE
      PHR=PHIR*TORAD
      PHM1=PHIM1*TORAD
      PHM2=PHIM2*TORAD
      CZ(2,3)=BL(2)
      IF (NZ.GE.3) THEN
         CZ(3,1)=BL(3)*SIN(ALPHA(3))
         IF (IZ(3,1).EQ.1) THEN
            CZ(3,3)=BL(3)*COS(ALPHA(3))
         ELSE
            CZ(3,3)=CZ(2,3)-BL(3)*COS(ALPHA(3))
         ENDIF
         DO 80 I=4,NZ
            IF (ABS(CZ(I-1,1)).GE.TENM5) GOTO 90
            CZ(I,1)=BL(I)*SIN(ALPHA(I))
            ITEMP=IZ(I,1)
            JTEMP=IZ(I,2)
            CZ(I,3)=CZ(ITEMP,3)-BL(I)*COS(ALPHA(I))*SIGN(ONE,CZ(ITEMP,
     &       3)-CZ(JTEMP,3))
80       CONTINUE
90       K=I
         IF (K.LE.NZ) THEN
            DO 250 J=K,NZ
               DCAJ=COS(ALPHA(J))
               DSAJ=SIN(ALPHA(J))
               DCBJ=COS(BETA(J))
               DSBJ=SIN(BETA(J))
               IF (IZ(J,4).EQ.0) THEN
                  CALL GEOMVU(U1,CZ,IZ(J,2),IZ(J,3))
                  CALL GEOMVU(U2,CZ,IZ(J,1),IZ(J,2))
                  CALL GEOMVP(VP,U1,U2)
                  FACT1=U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
                  R=SQRT(ONE-FACT1*FACT1)
                  DO 120 I=1,3
                     U3(I)=VP(I)/R
120               CONTINUE
                  CALL GEOMVP(U4,U3,U2)
                  DO 130 I=1,3
                     VJ(I)=BL(J)*(-U2(I)*DCAJ+U4(I)*DSAJ*DCBJ+U3(I)
     &                *DSAJ*DSBJ)
                     ITEMP=IZ(J,1)
                     CZ(J,I)=VJ(I)+CZ(ITEMP,I)
130               CONTINUE
               ELSEIF (IABS(IZ(J,4)).EQ.1) THEN
                  CALL GEOMVU(U1,CZ,IZ(J,1),IZ(J,3))
                  CALL GEOMVU(U2,CZ,IZ(J,2),IZ(J,1))
                  ZETA=-(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
                  A(J)=(-DCBJ+ZETA*DCAJ)/(ONE-ZETA*ZETA)
                  B(J)=(DCAJ-ZETA*DCBJ)/(ONE-ZETA*ZETA)
                  R=ZERO
                  GAMMA=PI/TWO
                  IF (ZETA.NE.0) THEN
                     IF (ZETA.LE.0) R=PI
                     GAMMA=ATAN(SQRT(ONE-ZETA*ZETA)/ZETA)+R
                  ENDIF
                  D(J)=ZERO
                  IF (ABS(GAMMA+ALPHA(J)+BETA(J)-TWO*PI).GE.TENM6) D(J)
     &             =IZ(J,4)*(SQRT(ONE+A(J)*DCBJ-B(J)*DCAJ))/SQRT(ONE
     &             -ZETA*ZETA)
                  CALL GEOMVP(V3,U1,U2)
                  DO 200 I=1,3
                     U3(I)=A(J)*U1(I)+B(J)*U2(I)+D(J)*V3(I)
                     VJ(I)=BL(J)*U3(I)
                     ITEMP=IZ(J,1)
                     CZ(J,I)=VJ(I)+CZ(ITEMP,I)
200               CONTINUE
               ELSE
                  CALL GEOMVU(U1,CZ,IZ(J,1),IZ(J,3))
                  CALL GEOMVU(U2,CZ,IZ(J,2),IZ(J,1))
                  ZETA=-(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
                  CALL GEOMVP(V3,U1,U2)
                  V3MAG=SQRT(V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
                  A(J)=V3MAG*DCBJ/(ONE-ZETA*ZETA)
                  B(J)=SQRT((ONE-DCAJ*DCAJ-A(J)*DCBJ*V3MAG)/(ONE-ZETA
     &             *ZETA))
                  IF (IZ(J,4).NE.2) B(J)=-B(J)
                  D(J)=B(J)*ZETA+DCAJ
                  DO 240 I=1,3
                     U3(I)=B(J)*U1(I)+D(J)*U2(I)+A(J)*V3(I)
                     VJ(I)=BL(J)*U3(I)
                     ITEMP=IZ(J,1)
                     CZ(J,I)=VJ(I)+CZ(ITEMP,I)
240               CONTINUE
               ENDIF
250         CONTINUE
         ENDIF
      ENDIF
C
C     ELIMINATE DUMMY ATOMS (THESE ARE CHARACTERIZED BY AN ATOMIC
C     NUMBER OF ZERO) AND CHANGE TO LEED PROGRAM COORDINATE SYSTEM,
C     WHICH HAS X POINTING INTO THE SURFACE,
C     Y AND Z ALONG THE SURFACE AND IS LEFT-HANDED.
C
      NATOMS=0
      DO 290 I=4,NZ
         IF (IANZ(I).NE.0) THEN
            NATOMS=NATOMS+1
            CZ(NATOMS,1)=CZ(I,3)
            CZ(NATOMS,2)=CZ(I,1)
            CZ(NATOMS,3)=CZ(I,2)
         ENDIF
290   CONTINUE
      NAT=NATOMS
      IF (PHIR.LE.360.) THEN
C
C PRODUCE ROTATED ATOMS
C
         IROT=INT(2.*PI/PHR+0.01)-1
         DO 410 IR=1,IROT
            CP=COS(PHR*IR)
            SP=SIN(PHR*IR)
            DO 420 I=1,NAT
               NATOMS=NATOMS+1
               CZ(NATOMS,1)=CZ(I,1)
               CZ(NATOMS,2)=CP*CZ(I,2)-SP*CZ(I,3)
               CZ(NATOMS,3)=SP*CZ(I,2)+CP*CZ(I,3)
C
C  PREVENT DUPLICATION OF ATOMS
C
               NAT1=NATOMS-1
               DO 405 I1=1,NAT1
                  DIZ=ABS(CZ(NATOMS,1)-CZ(I1,1))+ABS(CZ(NATOMS,2)-CZ
     &             (I1,2))+ABS(CZ(NATOMS,3)-CZ(I1,3))
                  IF (DIZ.LE.0.001) GOTO 5008
405            CONTINUE
               GOTO 420
5008           NATOMS=NATOMS-1
420         CONTINUE
410      CONTINUE
      ENDIF
C
C  PRODUCE MIRRORED ATOMS
C
      M2=1
      NAT=NATOMS
430   IF (PHM1.LE.2.*PI) THEN
         CP=COS(PHM1)
         SP=SIN(PHM1)
         DO 450 I=1,NAT
            NATOMS=NATOMS+1
            DOT=CZ(I,2)*CP+CZ(I,3)*SP
            CZ(NATOMS,1)=CZ(I,1)
            CZ(NATOMS,2)=2.*DOT*CP-CZ(I,2)
            CZ(NATOMS,3)=2.*DOT*SP-CZ(I,3)
C
C  PREVENT DUPLICATION
C
            NAT1=NATOMS-1
            DO 440 I1=1,NAT1
               DIZ=ABS(CZ(NATOMS,1)-CZ(I1,1))+ABS(CZ(NATOMS,2)-CZ(I1,2)
     &          )+ABS(CZ(NATOMS,3)-CZ(I1,3))
               IF (DIZ.LE.0.001) GOTO 5009
440         CONTINUE
            GOTO 450
5009        NATOMS=NATOMS-1
450      CONTINUE
      ENDIF
      IF (M2.NE.0) THEN
         PHM1=PHM2
         NAT=NATOMS
         M2=0
         GOTO 430
      ENDIF
      DO 5006 I=1,NATOMS
         DO 480 J=1,3
            C(I,J)=CZ(I,J)
480      CONTINUE
5006  CONTINUE
C
C     PRINT COORDINATES
C
      IF (IP.EQ.1) THEN
         WRITE (1,1000)
         DO 300 I=1,NATOMS
            WRITE (1,1010) I,(C(I,J),J=1,3)
300      CONTINUE
C
C*  DISTANCES BETWEEN ATOMS
C
         WRITE (1,5001)
         DO 810 IA=2,NATOMS
            IA1=IA-1
            DO 811 IB=1,IA1
               DIZ=0.0
               DO 473 IAW=1,3
                  FACT1=C(IA,IAW)-C(IB,IAW)
                  DIZ=DIZ+FACT1*FACT1
473            CONTINUE
               DIS(IB)=SQRT(DIZ)
811         CONTINUE
            IFLAG=1
            ILOWER=1
400         IUPPER=ILOWER+14
            IF (IUPPER.GE.IA1) THEN
               IUPPER=IA1
               IFLAG=0
            ENDIF
            WRITE (1,5000) (JA,JA=ILOWER,IUPPER)
            WRITE (1,5003) IA,(DIS(JA),JA=ILOWER,IUPPER)
            IF (IFLAG.NE.0) THEN
               ILOWER=ILOWER+15
               GOTO 400
            ENDIF
810      CONTINUE
      ENDIF
C
C  CONVERT TO ATOMIC UNITS
C
      DO 5007 IA=1,NATOMS
         DO 860 I=1,3
            C(IA,I)=C(IA,I)/0.529
860      CONTINUE
5007  CONTINUE
      RETURN
      END
C=========================================================================
C
C  Numerical Recipe substitute for ZGE (Van HOVE TONG ZGT)
C  Subroutine LUDCMP performs Gaussian elimination as the fist step in the 
C  solution of a system of linear equations. This is used to multiply
C  the inverse of a matrix into a vector, the multiplication being done
C  later by subroutine LUBSKB. For Details see Numerical Recipe.
C  Minor modification by BARBIERI
C
C=========================================================================
      SUBROUTINE LUDCMP(A,N,NP,INDX,D,VV)
      PARAMETER (NMAX=1000)
      COMPLEX A,VV,SUM,DUM,AAMAX
      DIMENSION A(NP,NP),INDX(N),VV(NP)
      TINY=1.e-10
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (CABS(A(I,J)).GT.CABS(AAMAX)) AAMAX=CABS(A(I,J))
11      CONTINUE
cjcm        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*CABS(SUM)
          IF (CABS(DUM).GE.CABS(AAMAX)) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(CABS(A(J,J))-TINY.LT.0.)A(J,J)=TINY*(1.,1.)
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(CABS(A(N,N))-TINY.LT.0.)A(N,N)=TINY*(1.,1.)
      RETURN
      END
C=========================================================================
C
C Numerical recipe substitute for ZSU ( Van Hove-Tong ZST)
C  Subroutine LUBKSB terminates the solution of a system of linear
C  equations initiated by subroutine LUDCMP, by back-substituting the
C  constant vector.
C
C=========================================================================
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      COMPLEX A,B,SUM
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
C=========================================================================
C
C  Subroutine LXGETT generated the needed relationships between
C  different ordering sequences of the (L,M) pairs (L.LE.LMAX, ABS(M)
C  .LE.L). Three orderings are considered, the \NATURAL\ one (N), the
C  \COPLANAR\ one (C) and the \SYMMETRIZED\ one (S). They are tabulated
C  below for the case of LMAX=4.
C
C     I=  0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2
C         1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
C
C  N  L=  0 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4
C     M=  0-1 0 1-2-1 0 1 2-3-2-1 0 1 2 3-4-3-2-1 0 1 2 3 4
C
C  C  L=  0 1 1 2 2 2 3 3 3 3 4 4 4 4 4 1 2 2 3 3 3 4 4 4 4
C     M=  0-1 1-2 0 2-3-1 1 3-4-2 0 2 4 0-1 1-2 0 2-3-1 1 3
C                  L+M=EVEN              *      L+M=ODD
C
C  S  L=  0 2 2 2 4 4 4 4 4 1 1 3 3 3 3 1 3 3 3 2 2 4 4 4 4
C     M=  0-2 0 2-4-2 0 2 4-1 1-3-1 1 3 0-2 0 2-1 1-3-1 1 3
C             L=EVEN          *  L=ODD    * L=ODD *  L=EVEN
C             M=EVEN          *  M=ODD    * M=EVEN*   M=ODD
C
C  To describe the relationships, a particular pair (L,M) in a particul-
C  ar sequence shall be represented here by N(I), C(I) or S(I) (I=1,
C  LMMAX)  e.g. N(10)=(3,-3). The relationships LX,LXI,LT,LXM generated
C  by LXGETT are now defined by
C
C   LX        N(LX(I))        =  S(I).
C   LXI  C(I)            =  S(LXI(I)) IF(I.LE.LEV),
C                           S(LXI(I)+LEV) IF(I.GT.LEV).
C   LT        CSM(N(LT(I)))   =  S(I) Where CSM means change the sign of M. 
C   LXM  N(I)            =  S(LXM(I)).
C
C Output:
C =======
C
C   LX,LXI,LT,LXM        =  Output permutations of (L,M) sequence.
C   LAN     same as LX but gives the L value instead of the N position
C
C Input:
C ======
C
C   LMAX                 =  Largest value of L.
C   LMMAX                =  (LMAX+1)**2.
C
C Modified version of routine LXGENT from the Van Hove/Tong LEED package.
C Modificatoins by WANDER.
C
C=========================================================================
C
      SUBROUTINE LXGENT2(LX,LXI,LT,LXM,LMAX,LMMAX,LAN)
C
      DIMENSION LX(LMMAX),LXM(LMMAX),LXI(LMMAX),LT(LMMAX),LAN(LMMAX)
C
      LEV=(LMAX+1)*(LMAX+2)/2
      LEE=(LMAX/2+1)*(LMAX/2+1)
      LEO=((LMAX+1)/2+1)*((LMAX+1)/2)+LEE
      LOE=((LMAX-1)/2+1)*((LMAX-1)/2+1)+LEO
      LL=0
      L1=0
      LT1=0
      L=-1
1     L=L+1
      M=-L
2     LL=LL+1
      IF (MOD(L+M,2).EQ.1) THEN
         LEV=LEV+1
         LT(LL)=LEV
         IF (MOD(L,2).EQ.0) THEN
            LOE=LOE+1
            LX(LOE)=LL
            LAN(LOE)=L
         ELSE
            LEO=LEO+1
            LX(LEO)=LL
            LAN(LEO)=L
         ENDIF
      ELSE
         LT1=LT1+1
         LT(LL)=LT1
         IF (MOD(L,2).EQ.1) THEN
            LEE=LEE+1
            LX(LEE)=LL
            LAN(LEE)=L
         ELSE
            L1=L1+1
            LX(L1)=LL
            LAN(L1)=L
         ENDIF
      ENDIF
      M=M+1
      IF (L.GE.M) GOTO 2
      IF (L.LT.LMAX) GOTO 1
      DO 12 L=1,LMMAX
         L1=LX(L)
         LT1=LT(L1)
         LXI(LT1)=L
12    CONTINUE
      L1=LEE+1
      DO 13 L=L1,LMMAX
         LXI(L)=LXI(L)-LEE
13    CONTINUE
      L1=LMAX+1
      JLM=1
      DO 17 L=1,L1
         LL=L+L-1
         DO 16 M=1,LL
            JLP=JLM+LL+1-2*M
            DO 14 LT1=1,LMMAX
               IF (JLM.EQ.LX(LT1)) GOTO 15
14          CONTINUE
15          LT(LT1)=JLP
            LXM(JLM)=LT1
            JLM=JLM+1
16       CONTINUE
17    CONTINUE
      RETURN
      END
C===================================================================
C
C Subroutine MSET calculates the information for the truncation of the
C Q matrix.
C
C Input Parameter;
C ================
C
C LSMAX                   =    MAXIMUM L VALUE FOR THE TENSOR Q
C LLCUT                   =    LARGEST VALUE OF L+L' TO BE INCLUDED
C IPCUT                   =    DIMENSION OF MICUT, MJCUT
C Q(MICUT(K),MJCUT(K))    =    THE KTH INDEPENDENT ELEMENT OF Q
C ICUT                    =    NUMBER OF INDEPENDENT ELEMENTS OF Q
C
C =========================================================================
C
      SUBROUTINE MSET(LSMAX,LLCUT,ICUT,MICUT,MJCUT,IPCUT)
C
      DIMENSION MICUT(IPCUT),MJCUT(IPCUT)
C    
100   FORMAT (' NUMBER OF ELEMENTS IN Q IS ',I3)
C
      ICUT=0
      I=0
      DO 101 L=0,LSMAX
         DO 102 M=-L,L
            I=I+1
            IP=0
            DO 103 LP=0,LSMAX
               DO 10 MP=-LP,LP
                  IP=IP+1
                  IF ((L+LP).LE.LLCUT) THEN
                     IF (IP.GE.I) THEN
                        ICUT=ICUT+1
                        MICUT(ICUT)=I
                        MJCUT(ICUT)=IP
                     ENDIF
                  ENDIF
10             CONTINUE
103         CONTINUE
102      CONTINUE
101   CONTINUE
      WRITE (1,100) ICUT
      RETURN
      END
C=======================================================================
C
C  Subroutine PRPGAT performs one propagation through one layer for RFS
C  subroutines like RFSO2 and RFSO3. A layer asymmetrical in +-X can be
C  handled, as well as independent beam sets.
C
C Parameter List;
C ===============
C
C  RA,TA,RB,TB    =  Input diffraction matrices for current layer (R for
C                    reflection, T for transmission, A,B for incidence towards
C                    +-X, resp.).
C  N              =  Total No. of beams at current energy.
C  NM             =  Largest No. of beams in any current beam set (But =N for 
C                    overlayer)
C  AW             =  Working space.
C  I              =  Index of interlayer spacing to which current call to 
C                   PRPGAT will lead.
C  L1,L2          =  Current choice of plane-wave propagators, referring to 
C                   second index of matrix PK. L1 is to describe propagation 
C                   in direction towards interlayer spacing I, L2 in opposite 
C                   direction.
C  CRIT           =  Criterion for penetration convergence (Set in calling 
C                   routine).
C  IR             =  Output flag for penetration convergence.
C  IA             =  +-1 indicates propagation towards +-X.
C  BNORM          =  Output measure of current wavefield amplitude at current 
C                   layer.
C
C Modified version of routine RFSO2 from the VAN HOVE/TONG LEED package.
C Modifications by ROUS and WANDER.
C Modifications by BARBIERI
C
C =========================================================================
C
      SUBROUTINE PRPGAT(RA,TA,RB,TB,N,NM,ANEW,ND,NSL,NL,AW,I,PK,L1,L2,
     & CRIT,IR,IA,BNORM,NROM)
C
      COMPLEX RA(NROM,NM),TA(NROM,NM),ANEW(N,ND),AW(N,2),PK(N,8),CZ
      COMPLEX RB(NROM,NM),TB(NROM,NM)
      DIMENSION NSL(NL)
C
6     FORMAT ('***THIS ORDER TOO DEEP')
C
      CZ=CMPLX(0.0,0.0)
      BNORM=0.0
      NA=0
      DO 4 NN=1,NL
         NB=NSL(NN)
         DO 1 K=1,NB
            KNA=K+NA
C
C  PROPAGATE WAVEFIELD TO CURRENT LAYER FROM THE TWO NEAREST LAYERS
C
            AW(K,1)=PK(KNA,L1)*ANEW(KNA,I-IA)
            AW(K,2)=PK(KNA,L2)*ANEW(KNA,I)
1        CONTINUE
         DO 10 J=1,NB
            JNA=J+NA
            ANEW(JNA,I)=CZ
C
C  SELECT FORMULA ACCORDING TO PROPAGATION DIRECTION
C
            IF (IA.EQ.-1) THEN
               DO 9 K=1,NB
                  ANEW(JNA,I)=ANEW(JNA,I)+TB(JNA,K)*AW(K,1)+RA(JNA,K)
     &             *AW(K,2)
9              CONTINUE
            ELSE
               DO 2 K=1,NB
C
C  TRANSMIT AND REFLECT AND ADD WAVES
C
                  ANEW(JNA,I)=ANEW(JNA,I)+TA(JNA,K)*AW(K,1)+RB(JNA,K)
     &             *AW(K,2)
2              CONTINUE
            ENDIF
C
C  GET MEASURE OF NEW WAVEFIELD MAGNITUDE
C
            FACT1=REAL(ANEW(JNA,I))*REAL(ANEW(JNA,I))
            FACT2=AIMAG(ANEW(JNA,I))*AIMAG(ANEW(JNA,I))
            BNORM=BNORM+FACT1+FACT2
10       CONTINUE
         NA=NA+NB
4     CONTINUE
C
C  CHECK AGAINST PENETRATION LIMIT, EXCEPT WHEN EMERGING (IN WHICH CASE
C  IR=1)
C
      IF (.NOT.((I.LT.ND-1).OR.(IR.EQ.1))) THEN
         WRITE (1,6)
         IR=1
      ENDIF
C
C  CHANGE I TO INDICATE TO WHICH NEW INTERLAYER SPACING THE NEXT CALL
C  TO PRPGAT WILL HAVE TO LEAD
C
      I=I+IA
C
C  CHECK ON PENETRATION CONVERGENCE
C
      IF (BNORM.LE.CRIT) IR=1
      RETURN
      END
C======================================================================
C  SUBROUTINE PSTEMP INCORPORATES THE THERMAL VIBRATION EFFECTS IN THE
C  PHASE SHIFTS, THROUGH A DEBYE-WALLER FACTOR. ISOTROPIC VIBRATION
C  AMPLITUDES ARE ASSUMED.
C   PPP= CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CPPP.
C   N3= NO. OF INPUT PHASE SHIFTS.
C   N2= DESIRED NO. OF OUTPUT TEMPERATURE-DEPENDENT PHASE SHIFTS.
C   N1= N2+N3-1.
C   DR0= 4TH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDES.
C   DR= ISOTROPIC RMS VIBRATION AMPLITUDE AT REFERENCE TEMPERATURE T0.
C   T0= ARBITRARY REFERENCE TEMPERATURE FOR DR.
C   TEMP= ACTUAL TEMPERATURE.
C   E= CURRENT ENERGY (REAL NUMBER).
C   PHS= INPUT PHASE SHIFTS.
C   DEL= OUTPUT (COMPLEX) PHASE SHIFTS.
C   11.01.95 UL: uses BESSEL from TLEED package to clculate spherical 
C   bessel-functions. Yield better convergence for higher vibration amplitudes
C   than original calculation scheme                                     110195
C
C Note:
C =====
C
C Modified version of routine PSTEMP from the VAN HOVE/TONG LEED package.
C Modifications by LOEFFLER.
C
C  DIMENSIONS ARE SET FOR N3.LE.11, N2.LE.11, N1.LE.21
      SUBROUTINE  PSTEMP (PPP, N1, N2, N3, DR0, DR, T0, TEMP, E, PHS,
     1DEL)
      COMPLEX  DEL, SUM, CTAB
      COMPLEX BJ(30)
      COMPLEX  Z, CI, CS, CL
      COMPLEX CEXP,CLOG
      DIMENSION  PPP(N1,N2,N3), PHS(N3), DEL(N2), SUM(20)
      DIMENSION  CTAB(20)
      PI = 3.14159265
      CI = CMPLX(0.0,1.0)
      DO 170 J = 1, N2
         DEL(J) = (0.0,0.0)
 170  CONTINUE
      ALFA = DR * DR * TEMP/T0
      ALFA = 0.166667 * SQRT(ALFA * ALFA + DR0)
      FALFE =  - 4.0 * ALFA * E
      IF (ABS(FALFE)-1.0E-3)  180, 200, 200
 180  CONTINUE
      DO 190 J = 1, N3
         DEL(J) = CMPLX(PHS(J),0.0)
 190  CONTINUE
      GO TO 360
C
COMMENT BJ(N1) IS LOADED WITH SPHERICAL BESSEL FUNCTIONS OF
C       THE FIRST KIND; ARGUMENT Z
  200 Z = FALFE * CI
c
c     use subroutine BESSEL
      call bessel(BJ,Z,N1-1,N1)
c
      CS = CMPLX(1.0,0.0)
      FL = 1.0
      DO 280 I = 1, N1
cc    WRITE(6,*) 'i=',i,' BJ(i)=',BJ(i)
      BJ(I) = EXP(FALFE) * FL * CS * BJ(I)
cc    WRITE(6,*) 'i=',i,' BJ(i)=',BJ(i)
      FL = FL + 2.0
      CS = CS * CI
  280 CONTINUE
C
      FL = 1.0
      DO 290 I = 1, N3
      CTAB(I) = (CEXP(2.0 * PHS(I) * CI) - (1.0,0.0)) * FL
      FL = FL + 2.0
  290 CONTINUE
C
      ITEST = 0
      LLLMAX = N2
      FL = 1.0
      DO 350 LLL = 1, N2
         SUM(LLL) = CMPLX(0.0,0.0)
         DO 305 L = 1, N3
         LLMIN = IABS(L - LLL) + 1
         LLMAX = L + LLL - 1
         DO 300 LL = LLMIN, LLMAX
            SUM(LLL) = SUM(LLL) + PPP(LL,LLL,L) * CTAB(L) * BJ(LL)
 300     CONTINUE
 305  CONTINUE
      DEL(LLL) =  - CI * CLOG(SUM(LLL) + (1.0,0.0))/(2.0,0.0)
      ABSDEL = CABS(DEL(LLL))
      IL = LLL - 1
      IF (ABSDEL-1.0E-2)  320, 310, 310
  310 ITEST = 0
      GO TO 350
  320 IF (ITEST-1)  340, 330, 340
  330 LLLMAX = LLL
      GO TO 360
 340  ITEST = 1
      FL = FL + 2.0
  350 CONTINUE
  360 RETURN
      END
C======================================================================
C
C  PSTEMo is the old version of PSTEMP which is not used any more
C  Subroutine PSTEMP incorporates the thermal vibration effects in the 
C  phase shifts, through a Debye-Waller factor. Isotropic vibration
C  amplitudes are assumed.
C
C Parameter List;
C ==============
C PPP         =    CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CPPP.
C N3          =    NO. OF INPUT PHASE SHIFTS.
C N2          =    DESIRED NO. OF OUTPUT TEMPERATURE-DEPENDENT PHASE SHIFTS.
C N1          =    N2+N3-1.
C DR0         =    4TH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDES.
C DR          =    ISOTROPIC RMS VIBRATION AMPLITUDE AT REFERENCE TEMPERATURE 
C                  T0.
C T0          =    ARBITRARY REFERENCE TEMPERATURE FOR DR.
C TEMP        =    ACTUAL TEMPERATURE.
C E           =    CURRENT ENERGY (REAL NUMBER).
C PHS         =    INPUT PHASE SHIFTS.
C DEL         =    OUTPUT (COMPLEX) PHASE SHIFTS.
C
C Note:
C =====
C
C Modified version of routine PSTEMP from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C=========================================================================
C
      SUBROUTINE PSTEMo(PPP,N1,N2,N3,DR0,DR,T0,TEMP,E,PHS,DEL)
C
      COMPLEX BJ,DEL,SUM,CTAB
      COMPLEX Z,FF,FN,FN1,FN2,CI,CS,CC,CL
      DIMENSION PPP(N1,N2,N3),PHS(N3),BJ(30),DEL(N2),SUM(20)
      DIMENSION CTAB(20)
C
      CI=CMPLX(0.0,1.0)
      DO 170 J=1,N2
         DEL(J)=(0.0,0.0)
170   CONTINUE
      ALFA=DR*DR*TEMP/T0
      ALFA=0.166667*SQRT(ALFA*ALFA+DR0)
      FALFE=-4.0*ALFA*E
      IF (ABS(FALFE).GE.1.0E-3) THEN
C
COMMENT BJ(N1) IS LOADED WITH SPHERICAL BESSEL FUNCTIONS OF
C       THE FIRST KIND; ARGUMENT Z
C
         Z=FALFE*CI
         LLL=N1-1
         CS=CSIN(Z)
         CC=CCOS(Z)
         BJ(1)=CS/Z
         LP=INT(3.5*CABS(Z))
         IF (LLL.GT.LP) THEN
            FF=CMPLX(1.0,0.0)
            CL=FF
            LPP=LP+1
            DO 220 MM=1,LPP
               CL=CL+(2.0,0.0)
               FF=FF*Z/CL
220         CONTINUE
            FL=FLOAT(LP)
            DO 230 J=LPP,LLL
               FL=FL+1.0
               B=1.0/(12.0*FL+42.0)
               FN=1.0-B*Z*Z
               B=1.0/(8.0*FL+20.0)
               FN=1.0-B*Z*Z*FN
               B=1.0/(4.0*FL+6.0)
               FN=1.0-B*Z*Z*FN
               BJ(J+1)=FF*FN
               CL=CL+(2.0,0.0)
               FF=FF*Z/CL
230         CONTINUE
            LLL=LP
            IF (LP.EQ.0) GOTO 270
         ENDIF
         BJ(2)=(BJ(1)-CC)/Z
         FN1=BJ(2)
         FN2=BJ(1)
         AM=3.0
         IF (LP.NE.1) THEN
            ILL=LLL+1
            DO 260 IA=3,ILL
               FN=AM*FN1/Z-FN2
               FN2=FN1
               FN1=FN
               AM=AM+2.0
               BJ(IA)=FN
260         CONTINUE
         ENDIF
270      CS=CMPLX(1.0,0.0)
         FL=1.0
         DO 280 I=1,N1
            BJ(I)=EXP(FALFE)*FL*CS*BJ(I)
            FL=FL+2.0
            CS=CS*CI
280      CONTINUE
         FL=1.0
         DO 290 I=1,N3
            CTAB(I)=(CEXP(2.0*PHS(I)*CI)-(1.0,0.0))*FL
            FL=FL+2.0
290      CONTINUE
         ITEST=0
         FL=1.0
         DO 350 LLL=1,N2
            SUM(LLL)=CMPLX(0.0,0.0)
            DO 361 L=1,N3
               LLMIN=IABS(L-LLL)+1
               LLMAX=L+LLL-1
               DO 300 LL=LLMIN,LLMAX
                  SUM(LLL)=SUM(LLL)+PPP(LL,LLL,L)*CTAB(L)*BJ(LL)
300            CONTINUE
361         CONTINUE
            DEL(LLL)=-CI*CLOG(SUM(LLL)+(1.0,0.0))/(2.0,0.0)
            ABSDEL=CABS(DEL(LLL))
            IL=LLL-1
            IF (ABSDEL.GE.1.0E-2) THEN
               ITEST=0
            ELSEIF (ITEST.EQ.1) THEN
               GOTO 360
            ELSE
               ITEST=1
            ENDIF
            FL=FL+2.0
350      CONTINUE
      ELSE
         DO 190 J=1,N3
            DEL(J)=CMPLX(PHS(J),0.0)
190      CONTINUE
      ENDIF
360   RETURN
      END
C==================================================================== 
C                                                                     
C Routine QGEN calculates the Q matrix for the current layer and Energy.
C 
C Input Parameters;
C =================
C
C GA,GB              = G VECTORS FROM CVEC
C LSMMAX             = (LSMAX+1)**2
C LMMAX              = (LMAX+1)**2
C NLLAY              = INDEX OF CURRENT LAYER
C NLAY               = NUMBER OF SUBPLANES IN COMPOSITE LAYER
C E,VPI              = CURRENT (COMPLEX) ENERGY
C NEXIT              = INDEX OF CURRENT BEAM
C NT0                = TOTAL NUMBER OF EXIT BEAMS
C AK2M,AK3M          = PARALLEL COMPONENTS OF MOMENTUM OF EACH EXIT BEAM
C TV                 = AREA OF UNRECONSTRUCTED UNIT CELL
C Q                  = THE Q MATRIX FOR THE CURRENT LAYER AND EXIT BEAM
C                     DIRECTION AT THE CURRENT ENERGY
C
C =========================================================================
C
      SUBROUTINE QGEN(GA,GB,Q,LMMAX,LSMMAX,NLLAY,NLAY,E,VPI,NEXIT,NT0,
     & AK2M,AK3M,TV)
C
      COMPLEX GA(NLAY,LMMAX,LSMMAX),GB(LMMAX,LSMMAX)
      COMPLEX Q(LSMMAX,LSMMAX)
      COMPLEX CSUM,CI,PREG,CAK,XA
      DIMENSION AK2M(NT0),AK3M(NT0)
C
      CI=CMPLX(0.0,1.0)
      D=0.0
      IF (NEXIT.GT.0) D=AK2M(NEXIT)**2+AK3M(NEXIT)**2
      CAK=CMPLX(2.0*E,-2.0*VPI+0.0000001)
      CAK=CSQRT(CAK)
      XA=CMPLX(2.0*E-D,-2.0*VPI+0.0000001)
      XA=CSQRT(XA)
      PREG=1.0/(2.0*TV*CAK*XA)
      DO 100 I2=1,LSMMAX
         DO 110 I3=1,LSMMAX
            CSUM=CMPLX(0.0,0.0)
            DO 120 I1=1,LMMAX
               CSUM=CSUM+GB(I1,I2)*GA(NLLAY,I1,I3)
120         CONTINUE
            Q(I2,I3)=PREG*CSUM
110      CONTINUE
100   CONTINUE
      RETURN
      END
C=========================================================================
C
C  Subroutine READCT reads in data relevant to the composite layer(s).
C  Input from MAIN program:-
C
C   NLAY    =   Max. No. of subplanes in any composite layer.
C   LMMAX   =   No. of spherical harmonics used =(LMAX+1)**2
C   NST1    =   Number of compoiste layers in input.
C   LAFLAG  =   Number of layers in each composite layer.
C
C IN COMMON BLOCKS
C
C   CPVPOS  =   Atomic coordinates of the atoms in the composite layers
C               before sorting as needed in subroutine LOOKUP
C
C This is a modified version of routine READCL from the Van Hove/Tong LEED
C package. Modifcations by WANDER and BARBIERI.
C 
C
C========================================================================
C
      SUBROUTINE READCT(NLAY,VPOS,FPOS,POSS,CPVPOS,NTAUAW,LPSAW,
     % IPR,LAFLAG,NST1,ASB,VICL,VCL,FRCL,TST,TSTS,ASA,INVECT)
C
      DIMENSION FPOS(NLAY,3),VPOS(NST1,NLAY,3),LPSAW(NST1,NLAY)
      DIMENSION CPVPOS(NST1,NLAY,3),ASB(NST1,3),ASA(10,3)
      DIMENSION LAFLAG(NST1),NTAUAW(NST1),POSS(NLAY,3)
      DIMENSION VICL(NST1),VCL(NST1),FRCL(NST1)
C
      COMMON /ZMAT/IANZ(40),IZ(40,4),BL(40),ALPHA(40),BETA(40),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(40),ALPHAS(40),BETAS
     & (40),PHIR,PHIM1,PHIM2
C
160   FORMAT (3F7.4)
161   FORMAT (3F7.2)
180   FORMAT (/' ASB',1I3,3X,3F9.4)
200   FORMAT (500I3)
215   FORMAT (/' COMPOSITE LAYER VECTOR   ',3F9.4)
216   FORMAT (' NTAU = ',20I3)
217   FORMAT (' PHASE SHIFT ASSIGNMENT IN COMPOSITE LAYER No.',1I2,
     & 3X,200I3)
444   FORMAT (//' COMPOSITE LAYER No. ',1I3)
445   FORMAT (' =================== ')
1020  FORMAT (I3,I4,F11.6,I4,F11.6,I4,F11.6,I4)
C
      IF (NLAY.NE.0) THEN
C
C IFLAG = O USE FPOS
C       = 1 USE GEOMV
C
         READ (5,*) IFLAG
C
C  NTAU= NO. OF CHEMICAL ELEMENTS IN THE COMPOSITE LAYER
C
         READ (5,200) (NTAUAW(I),I=1,NST1)
         IF (IPR.GT.0) WRITE (1,216) (NTAUAW(I),I=1,NST1)
C
C  LPSAW GIVES CHEMICAL IDENTITY FOR EACH SUBPLANE, REFERRING TO ORDER
C  OF INPUT OF PHASE SHIFTS IN SUBROUTINE READT. LPSAW(J,I)=K MEANS I-TH
C  SUBPLANE OF J'TH COMPOSITE LAYER HAS K-TH SET OF PHASE SHIFTS IN INPUT 
C  SEQUENCE
C
         DO 648 K=1,NST1
            READ (5,200) (LPSAW(K,I),I=1,LAFLAG(K))
            IF (IPR.GT.0) WRITE (1,217) K,(LPSAW(K,I),I=1,LAFLAG(K))
648      CONTINUE
C
C USE IFLAG AS POINTER TO INDICATE WETHER INPUT IS VIA Z-MATRIX (IFLAG=1)
C OR VPOS IFLAG=0
C
         IF (IFLAG.EQ.0) THEN
            DO 481 K=1,NST1
               DO 1041 I=1,LAFLAG(K)
C
C  VPOS  3-D VECTORS POINTING FROM AN ARBITRARY REFERENCE POINT TO ONE
C  ATOM IN EACH SUBPLANE (THE REFERENCE POINT WILL BE USED TO PRODUCE
C  A PAIR OF DIFFERENT REFERENCE POINTS IN MPERTI AND MTINV)(ANGSTROM)
C
C                  READ (5,160) (VPOS(K,I,J),J=1,3)
                  READ (5,*) (VPOS(K,I,J),J=1,3)
                  IF (IPR.GT.0) WRITE (1,215) (VPOS(K,I,J),J=1,3)
                  DO 480 J=1,3
                     FPOS(I,J)=VPOS(K,I,J)/0.529
                     VPOS(K,I,J)=VPOS(K,I,J)/0.529
480               CONTINUE
1041           CONTINUE
C               READ (5,160) (ASB(K,I),I=1,3)
               READ (5,*) (ASB(K,I),I=1,3)
               IF (IPR.GT.0) WRITE (1,180) K,(ASB(K,I),I=1,3)
               DO 4666 IK=1,3
                  ASB(K,IK)=ASB(K,IK)/0.529
4666           CONTINUE
C               READ (5,161) FRCL(K),VCL(K),VICL(K)
               READ (5,*) FRCL(K),VCL(K),VICL(K)
               VCL(K)=VCL(K)/27.21
               VICL(K)=VICL(K)/27.21
               NLAY2=LAFLAG(K)
               CALL SORT(FPOS,POSS,NLAY,NLAY2)
               DO 3333 J=1,NLAY2
                  DO 3332 KK=1,3
                     CPVPOS(K,J,KK)=0.529*FPOS(J,KK)
3332               CONTINUE
3333            CONTINUE
481         CONTINUE
         ELSE
C
C  THE FOLLOWING INPUT IS TO BE USED IN CONJUNCTION WITH ROUTINE
C  GEOMV TO DEFINE ATOMIC POSITIONS, ESPECIALLY IN MOLECULES, AND
C  VARIATIONS THEREOF FOR A STRUCTURAL SEARCH. SEE EXPLANATIONS
C  IN GEOMV FOR DETAILS.
C
            DO 1042 II=1,NST1
               READ (5,200) NATOMS
               NZ=NATOMS+3
               DO 30 I=4,NZ
                  READ (5,1020) IANZ(I),IZ(I,1),BLS(I),IZ(I,2),ALPHAS
     &             (I),IZ(I,3),BETAS(I),IZ(I,4)
30             CONTINUE
               READ (5,160) PHIR,PHIM1,PHIM2
C
C  Read the vector connecting the bottom of the II CLayer to the top
C  of the II+1 CLayer (first substrate layer if II=NST1)
C
               WRITE (1,*)
               WRITE (1,444) II
               WRITE (1,445)
               READ (5,160) (ASB(II,I),I=1,3)
               IF (IPR.GT.0) WRITE (1,180) II,(ASB(II,I),I=1,3)
               DO 466 IK=1,3
                  ASB(II,IK)=ASB(II,IK)/0.529
466            CONTINUE
               READ (5,161) FRCL(II),VCL(II),VICL(II)
               VCL(II)=VCL(II)/27.21
               VICL(II)=VICL(II)/27.21
C
C CALCULATE STRUCTURES IN FPOS AND COPY INTO VPOS FOR FUTURE REFERENCE
C
               NLAY2=LAFLAG(II)
               CALL GEOMV(1,FPOS,NLAY2,NLAY)
               DO 1043 J=1,NLAY2
                  DO 348 K=1,3
                     VPOS(II,J,K)=FPOS(J,K)
348               CONTINUE
1043           CONTINUE
C CPVPOS corresponds to sorted coordinates needed in subroutine LOOKUP
               CALL SORT(FPOS,POSS,NLAY,NLAY2)
               DO 333 J=1,NLAY2
                  DO 332 K=1,3
                     CPVPOS(II,J,K)=0.529*FPOS(J,K)
332               CONTINUE
333            CONTINUE
1042        CONTINUE
         ENDIF
      ENDIF
      TSTS=TST
C      AMINA=100.
      AMINA=ASA(1,1)
      DO 248 I=2,INVECT
         AMINA=AMIN1(ASA(I,1),AMINA)
248   CONTINUE
      DO 249 I=1,NST1
         AMINA=AMIN1(ASB(I,1),AMINA)
249   CONTINUE
      TST=ALOG(TST)/AMINA
      TST=TST*TST
      RETURN
      END
C ======================================================================
C 
C Subroutine READPL reads in the information relevant to the Tensor part
C of the LEED calculation.
C
C Parameters from main program;
C =============================
C
C NT0     =  Number of exit beams in calculation
C NLAY    =  Total Number of layers in composite layers
C IPR     =  Print control parameter
C
C Author: WANDER.  Modifications BARBIERI
C
C ======================================================================
C 
      SUBROUTINE READPL(NT0,NSET,PQFEX,NINSET,NDIM,DISP,
     & NLAY,IPR,ALPHA,BETA,GAMMA,ITMAX,
     & FTOL1,FTOL2,MFLAG,LLFLAG,NGRID)
C
      DIMENSION PQFEX(2,NT0),NINSET(20)
      DIMENSION DISP(NLAY,3),LLFLAG(NLAY+1)
      COMMON /RPL/DVOPT
      COMMON /WIV2/PERSH,NIV,NSE1(30),NSE2(30)
cjcm      COMMON /WIV2/PERSH,NIV
C
100   FORMAT (20I3)
110   FORMAT (10F7.4)
115   FORMAT (3F7.6,I3)
116   FORMAT (F7.6,I3)
120   FORMAT (/' Calculating for Exit Beams; ')
130   FORMAT (' =========================== ')
200   FORMAT (/' Initial Coordinates in Search ')
210   FORMAT (' ============================= ')
!220   FORMAT (/' Line Search Information ')
!230   FORMAT (' ======================= ')
!240   FORMAT (/' Moving coordinate ',I3,' with ',I3,' moving atoms ')
!250   FORMAT (/' Moving Atom Identifiers; ')
!260   FORMAT (/' Number of steps in Line Search ',I3,' Step Size',F7.4)
270   FORMAT (/' Search Information; ')
280   FORMAT (/' ALPHA =',F7.4,' BETA =',F7.4,' GAMMA =',F7.4)
290   FORMAT (/' Max. Number of Iterations =',I3)
300   FORMAT (/' Down-hill Simplex Search Requested')
310   FORMAT (/' Powell Direction-set Requested')
320   FORMAT (/' Direction Set with Principal Directions Requested')
!350   FORMAT (/' Line Search Requested')
C
C Loop over the beamsets reading in the beam indices of each desired exit
C beam.
C
      IOFF=0
      DO 1000 I=1,NSET
         READ (4,100) NINSET(I)
         DO 1010 J=1,NINSET(I)
            READ (4,110) (PQFEX(K,IOFF+J),K=1,2)
1010     CONTINUE
         IOFF=IOFF+NINSET(I)
1000  CONTINUE
      IF (IPR.GT.0) THEN
         WRITE (1,120)
         WRITE (1,130)
         DO 1015 I=1,NT0
            WRITE (1,110) (PQFEX(K,I),K=1,2)
1015     CONTINUE
      ENDIF
C
C Read in number of dimensions to be used in the search algorithm, NDIM.
C If NDIM=1, then only the perpendicular component of the displacement
C will be allowed to vary. NDIM=3 will generate a search over all of 
C cartesian space. If NDIM=0, then the program will perform a line
C search over one (or more) coordinates.
C
      READ (4,100) NDIM
C
C Read in initial geometry, initial displacement of the inner potential,
C and whether to include the coordinate in the search.
C
      DO 1020 I=1,NLAY
         READ (4,115) (DISP(I,J),J=1,3),LLFLAG(I)
1020  CONTINUE
         READ (4,116) DVOPT,LLFLAG(NLAY+1)
      IF (IPR.GT.0) THEN
         WRITE (1,200)
         WRITE (1,210)
         DO 1030 I=1,NLAY
            WRITE (1,110) (DISP(I,J),J=1,3)
1030     CONTINUE
      ENDIF
      READ (4,100)MFLAG,NGRID,NIV
      IF(IPR.GT.0)THEN
         IF(MFLAG.EQ.1)THEN
          WRITE(1,300)
         ENDIF
         IF(MFLAG.EQ.2)THEN
           WRITE(1,310)
         ENDIF
         IF(MFLAG.EQ.3)THEN
           WRITE(1,320)
         ENDIF
      ENDIF
      READ (4,100) ITMAX
      READ (4,110) ALPHA,BETA,GAMMA
      READ (4,110) FTOL1,FTOL2
      IF (IPR.GT.0)THEN
          WRITE (1,270)
          WRITE (1,280) ALPHA,BETA,GAMMA
          WRITE (1,290) ITMAX
      ENDIF
      RETURN
      END
C ========================================================================
C
C Subroutine READT is not used any more, and substituted by the
C new subroutines READT1 and READT2 (see below)
C
C Subroutine READT inputs most of the data required by the 'conventional'
C LEED part of TLEED with the exception of data involving the composite 
C layer. It is a modified version of routine READIN from the Van Hove/Tong 
C LEED package. The modifications involve a reordering of the order of data 
C input. Modifications by WANDER.
C
C Input From Main Program;
C ========================
C
C NPSI       =  NUMBER OF ENERGIES AT WHICH PHASE SHIFTS ARE READ IN.
C LMAX       =  MAXIMUM VALUE OF THE ANGULAR MOMENTUM
C NPHASE     =  NUMBER OF PHASE SHIFTS*NUMBER OF ELEMENTS (DIMENSION OF PHSS)
C
C=========================================================================
C
      SUBROUTINE READT(TVA,RAR1,RAR2,ASA,INVECT,TVB,IDEG,NL,V,VL,JJS,
     & TST,THETA,FI,LMMAX,NPSI,ES,PHSS,L1,IPR,NEL)
C
      COMPLEX VL(NL,2)
      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),ASA(10,3),ARB1(2)
      DIMENSION ARB2(2),RBR1(2),RBR2(2)
      DIMENSION CPARB1(2),CPARB2(2)
      DIMENSION V(NL,2),JJS(NL,IDEG),ES(NPSI),PHSS(NPSI,80)
      DIMENSION IT1(5),THDB(5),AM1(5),FPER1(5),FPAR1(5),DR01(5)
      DIMENSION DRPER1(5),DRPAR1(5)
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
      COMMON /MS/LMAX
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /TEMP/IT1,TI,T0,DRPER1,DRPAR1,DR01
      COMMON /LO/CPARB1,CPARB2
C
130   FORMAT (5F9.4)
140   FORMAT (/' PARAMETERS FOR INTERIOR ')
160   FORMAT (3F7.4)
166   FORMAT (3F11.8)
161   FORMAT (3F7.2)
170   FORMAT (/' SURF VECS',2(5X,2F8.4))
171   FORMAT (15X,2F8.4)
172   FORMAT (/' Interlayer Vectors ')
180   FORMAT (/' ASA',1I3,3X,3F7.4)
200   FORMAT (24I3)
210   FORMAT (/' FR = ',F7.4,' ASE = ',F7.4)
280   FORMAT (' ',1I4,F10.3,F7.3,4X,4I3)
281   FORMAT (' TST = ',1F7.4)
285   FORMAT (' THETA  FI = ',2F7.2)
290   FORMAT ('VV',4X,F7.2,'VPIS',4X,F7.2)
305   FORMAT (' IT1 = ',3I3)
325   FORMAT (' THDB = ',1F9.4,' AM = ',1F9.4,' FPER = ',1F9.4,
     & ' FPAR = ',1F9.4,' DR0 = ',1F9.4)
340   FORMAT (16F7.4)
350   FORMAT (/6X,'PHASE SHIFTS')
355   FORMAT (' LMAX = ',1I3)
360   FORMAT ('E = ',1F7.4,'  1ST ELEMENT',3X,16F8.4)
361   FORMAT (' ',11X,' ELEMENT # ',I2,3X,16F8.4)
!362   FORMAT (' ',11X,'  3RD ELEMENT',3X,16F8.4)
C
      PI=3.14159265
      IF (IPR.GT.0) WRITE (1,140)
C
C  FR IS THE FRACTION OF THE TOP LAYER TO SUBSTRATE SPACING THAT IS
C  ASSIGNED TO THE SUBSTRATE (OVER WHICH THE SUBSTRATE MUFFIN-TIN
C  CONSTANT AND DAMPING APPLY).
C  ASE IS THE SPACING BETWEEN THE SURFACE (WHERE MUFFIN-TIN CONSTANT
C  AND DAMPING SET IN) AND THE TOP-LAYER NUCLEI (ANGSTROM).
C
      READ (5,160) FR,ASE
      IF (IPR.GT.0) WRITE (1,210) FR,ASE
      ASE=ASE/0.529
C
C  TST GIVES THE CRITERION FOR THE SELECTION OF BEAMS AT ANY ENERGY
C  TO BE SELECTED A BEAM MAY NOT DECAY TO LESS THAN A FRACTION TST OF
C  ITS INITIAL VALUE ON TRAVELING FROM ANY LAYER TO THE NEXT
C
      READ (5,160) TST
      IF (IPR.GT.0) WRITE (1,281) TST
C
C  (THETA,FI) IS THE DIRECTION OF INCIDENCE. THETA=0 AT NORMAL
C  INCIDENCE. FI=0 FOR INCIDENCE ALONG Y-AXIS IN (YZ) SURFACE PLANE
C  (X-AXIS IS PERPENDICULAR TO SURFACE)(DEGREE)
C
      READ (5,161) THETA,FI
      IF (IPR.GT.0) WRITE (1,285) THETA,FI
      THETA=THETA*PI/180.0
      FI=FI*PI/180.0
C
C  VV IS VACUUM LEVEL WITH RESPECT TO SUBSTRATE MUFFIN-
C  TIN CONSTANT (EV)
C
      READ (5,161) VV,VPIS
      IF (IPR.GT.0) WRITE (1,290) VV,VPIS
      VV=VV/27.21
      VPIS=VPIS/27.21
C
C  LMAX= HIGHEST L-VALUE TO BE CONSIDERED
C
      READ (5,200) LMAX
      IF (IPR.GT.0) WRITE (1,355) LMAX
C
C  LMMAX= NO. OF SPHERICAL HARMONICS TO BE CONSIDERED
C
      LMMAX=(LMAX+1)*(LMAX+1)
C
C  NEL= NO. OF CHEMICAL ELEMENTS IN CALCULATION (.LE.5)
C
      READ (5,200) NEL
      IF (NEL.GT.5) THEN
         WRITE (1,*) ' CURRENT CODE NOT DIMENSIONED FOR NEL.GT.5'
c     jcm replace abort with stop
c     jcm         CALL ABORT()
         STOP
      ENDIF
C
C  IT1(I)=1,(0) MEANS  DO (NOT) INCLUDE THERMAL VIBRATION EFFECTS
C  FOR ATOM TYPE I
C
      READ (5,200) (IT1(I),I=1,NEL)
      IF (IPR.GT.0) WRITE (1,305) (IT1(I),I=1,NEL)
C
C  THDB= DEBYE TEMPERATURE (KELVIN).
C  AM= ATOMIC MASS (IN NUCLEON UNITS).
C  FPER,FPAR ARE ENHANCEMENT FACTORS FOR THE MEAN SQUARE VIBRATION
C  AMPLITUDES PERPENDICULAR AND PARALLEL TO THE SURFACE, RESP.
C  DR0 IS THE RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE (ANGSTROM).
      DO 48 I=1,NEL
         READ (5,130) THDB(I),AM1(I),FPER1(I),FPAR1(I),DR01(I)
         IF (IPR.GT.0) WRITE (1,325) THDB(I),AM1(I),FPER1(I),FPAR1(I),
     &    DR01(I)
48    CONTINUE
C
C  TI IS THE CALCULATION TEMPERATURE
C
      READ (5,130) TI
      T0=TI
      IF (T0.LT.0.001) T0=0.001
C
C  DRPER AND DRPAR ARE RMS VIBRATION AMPLITUDES PERPENDICULAR AND
C  PARALLEL TO THE SURFACE, RESP., EVALUATED AT REFERENCE TEMPERATURE
C  T0=MAX(TI,0.001 K)
C
      DO 49 I=1,NEL
         DR2=1.546E3*T0/(AM1(I)*THDB(I)*THDB(I))
         DRPER1(I)=SQRT(FPER1(I)*DR2)
         DRPAR1(I)=SQRT(FPAR1(I)*DR2)
         DR01(I)=(DR01(I)/0.529)*(DR01(I)/0.529)
         DR01(I)=DR01(I)*DR01(I)
49    CONTINUE
      L1=LMAX+1
      DO 660 I=1,NPSI
C
C  ES= ENERGIES (HARTREES) AT WHICH PHASE SHIFTS ARE INPUT. LINEAR
C  INTERPOLATION OF THE PHASE SHIFTS WILL OCCUR FOR ACTUAL ENERGIES
C  FALLING BETWEEN THE VALUES OF ES (AND LINEAR EXTRAPOLATION ABOVE
C  THE HIGHEST ES)
C
         READ (5,160) ES(I)
C
C  PHSS STORES THE INPUT PHASE SHIFTS (RADIAN)
C
         DO 661 II=1,NEL
           IO=(II-1)*L1
           READ (5,340) (PHSS(I,L),L=1+IO,L1+IO)
661      CONTINUE
660   CONTINUE
      IF (IPR.GT.0) THEN 
        WRITE (1,350)
        DO 670 I=1,NPSI
        DO 671 II=1,NEL
           IO=(II-1)*L1
           IF (II.EQ.1) THEN
              WRITE (1,280)
              WRITE (1,360) ES(I),(PHSS(I,L),L=1+IO,L1+IO)
           ELSE 
              WRITE (1,361)II,(PHSS(I,L),L=1+IO,L1+IO) 
           ENDIF
671   CONTINUE
670   CONTINUE
       ENDIF
C
C  ARA1 AND ARA2 ARE TWO 2-D BASIS VECTORS OF THE SUBSTRATE LAYER
C  LATTICE. THEY SHOULD BE EXPRESSED IN TERMS OF THE PLANAR CARTESIAN
C  Y- AND Z-AXES (X-AXIS IS PERPENDICULAR TO SURFACE)(ANGSTROM)
C
C      READ (5,160) (ARA1(I),I=1,2)
C      READ (5,160) (ARA2(I),I=1,2)
      READ (5,*) (ARA1(I),I=1,2)
      READ (5,*) (ARA2(I),I=1,2)
      IF (IPR.GT.0) THEN
         WRITE (1,170) (ARA1(I),I=1,2)
         WRITE (1,171) (ARA2(I),I=1,2)
      ENDIF
      DO 460 I=1,2
         ARA1(I)=ARA1(I)/0.529
         ARA2(I)=ARA2(I)/0.529
460   CONTINUE
      TVA=ABS(ARA1(1)*ARA2(2)-ARA1(2)*ARA2(1))
      ATV=2.0*PI/TVA
C
C  RAR1 AND RAR2 ARE THE RECIPROCAL-LATTICE BASIS VECTORS CORRESPONDING
C  TO ARA1 AND ARA2
C
      RAR1(1)=ARA2(2)*ATV
      RAR1(2)=-ARA2(1)*ATV
      RAR2(1)=-ARA1(2)*ATV
      RAR2(2)=ARA1(1)*ATV
C
C  ARB1,ARB2,RBR1,RBR2 ARE EQUIVALENT TO ARA1,ARA2,RAR1,
C  RAR2 BUT FOR AN OVERLAYER 
C
C      READ (5,160) (ARB1(I),I=1,2)
C      READ (5,160) (ARB2(I),I=1,2)
      READ (5,*) (ARB1(I),I=1,2)
      READ (5,*) (ARB2(I),I=1,2)
      DO 465 I=1,2
         CPARB1(I)=ARB1(I)
         CPARB2(I)=ARB2(I)
465   CONTINUE
      IF (IPR.GT.0) THEN
         WRITE (1,170) (ARB1(I),I=1,2)
         WRITE (1,171) (ARB2(I),I=1,2)
      ENDIF
      DO 467 I=1,2
         ARB1(I)=ARB1(I)/0.529
         ARB2(I)=ARB2(I)/0.529
467   CONTINUE
      TVB=ABS(ARB1(1)*ARB2(2)-ARB1(2)*ARB2(1))
      ATV=2.0*PI/TVB
      RBR1(1)=ARB2(2)*ATV
      RBR1(2)=-ARB2(1)*ATV
      RBR2(1)=-ARB1(2)*ATV
      RBR2(2)=ARB1(1)*ATV
C
C  INVECT IS NUMBER OF INTERLAYER VECTORS IN INPUT FILE
C
      READ (5,*) INVECT
      IF (IPR.GT.0) WRITE (1,172)
      DO 348 J=1,INVECT
         READ (5,166) (ASA(J,I),I=1,3)
         IF (IPR.GT.0) WRITE (1,180) J,(ASA(J,I),I=1,3)
         DO 466 I=1,3
            ASA(J,I)=ASA(J,I)/0.529
466      CONTINUE
348   CONTINUE
C
C  SLIND COMPARES SUBSTRATE LATTICE AND SUPERLATTICE FOR THE BENEFIT
C  OF FMAT
C
      CALL SLIND(V,VL,JJS,NL,IDEG,2.0E-4)
      RETURN
      END
C ========================================================================
C Subroutines READT1 and READT2 input most of the data required
C by the 'conventional' LEED part of TLEED with the exception of
C data involving the composite layer. It is a modified version
C of routine READIN from the Van Hove/Tong LEED package.
C The modifications involve a reordering of the order of data 
C input. Modifications by WANDER.
C Later modifications by DOELL: Subroutine READT was divided into two parts,
C READT1 and READT2. This allows for the calculation of the superlattice
C characterization NL1 and NL2.
C
C Input From Main Program;
C ========================
C
C NPSI       =  NUMBER OF ENERGIES AT WHICH PHASE SHIFTS ARE READ IN.
C LMAX       =  MAXIMUM VALUE OF THE ANGULAR MOMENTUM
C NPHASE     =  NUMBER OF PHASE SHIFTS*NUMBER OF ELEMENTS (DIMENSION OF PHSS)
C
C=========================================================================
C
      SUBROUTINE READT1(TVA,RAR1,RAR2,INVECT,TVB,IDEG,NL,
     & TST,THETA,FI,LMMAX,NPSI,ES,PHSS,L1,IPR,NEL,
     & IPNL1,IPNL2)
C
      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),ARB1(2)
      DIMENSION ARB2(2),RBR1(2),RBR2(2)
      DIMENSION CPARB1(2),CPARB2(2)
      DIMENSION ES(NPSI),PHSS(NPSI,80)
      DIMENSION IT1(5),THDB(5),AM1(5),FPER1(5),FPAR1(5),DR01(5)
      DIMENSION DRPER1(5),DRPAR1(5)
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
      COMMON /MS/LMAX
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /TEMP/IT1,TI,T0,DRPER1,DRPAR1,DR01
      COMMON /LO/CPARB1,CPARB2
C
130   FORMAT (5F9.4)
140   FORMAT (/' PARAMETERS FOR INTERIOR ')
160   FORMAT (3F7.4)
161   FORMAT (3F7.2)
170   FORMAT (/' SURF VECS',2(5X,2F8.4))
171   FORMAT (15X,2F8.4)
200   FORMAT (24I3)
210   FORMAT (/' FR = ',F7.4,' ASE = ',F7.4)
280   FORMAT (' ',1I4,F10.3,F7.3,4X,4I3)
281   FORMAT (' TST = ',1F7.4)
285   FORMAT (' THETA  FI = ',2F7.2)
290   FORMAT ('VV',4X,F7.2,'VPIS',4X,F7.2)
305   FORMAT (' IT1 = ',3I3)
325   FORMAT (' THDB = ',1F9.4,' AM = ',1F9.4,' FPER = ',1F9.4,
     & ' FPAR = ',1F9.4,' DR0 = ',1F9.4)
340   FORMAT (16F7.4)
350   FORMAT (/6X,'PHASE SHIFTS')
355   FORMAT (' LMAX = ',1I3)
360   FORMAT ('E = ',1F7.4,'  1ST ELEMENT',3X,16F8.4)
361   FORMAT (' ',11X,' ELEMENT # ',I2,3X,16F8.4)
C
      PI=3.14159265
      IF (IPR.GT.0) WRITE (1,140)
C
C  FR IS THE FRACTION OF THE TOP LAYER TO SUBSTRATE SPACING THAT IS
C  ASSIGNED TO THE SUBSTRATE (OVER WHICH THE SUBSTRATE MUFFIN-TIN
C  CONSTANT AND DAMPING APPLY).
C  ASE IS THE SPACING BETWEEN THE SURFACE (WHERE MUFFIN-TIN CONSTANT
C  AND DAMPING SET IN) AND THE TOP-LAYER NUCLEI (ANGSTROM).
C
      READ (5,160) FR,ASE
      IF (IPR.GT.0) WRITE (1,210) FR,ASE
      ASE=ASE/0.529
C
C  TST GIVES THE CRITERION FOR THE SELECTION OF BEAMS AT ANY ENERGY
C  TO BE SELECTED A BEAM MAY NOT DECAY TO LESS THAN A FRACTION TST OF
C  ITS INITIAL VALUE ON TRAVELING FROM ANY LAYER TO THE NEXT
C
      READ (5,160) TST
      IF (IPR.GT.0) WRITE (1,281) TST
C
C  (THETA,FI) IS THE DIRECTION OF INCIDENCE. THETA=0 AT NORMAL
C  INCIDENCE. FI=0 FOR INCIDENCE ALONG Y-AXIS IN (YZ) SURFACE PLANE
C  (X-AXIS IS PERPENDICULAR TO SURFACE)(DEGREE)
C
      READ (5,161) THETA,FI
      IF (IPR.GT.0) WRITE (1,285) THETA,FI
      THETA=THETA*PI/180.0
      FI=FI*PI/180.0
C
C  VV IS VACUUM LEVEL WITH RESPECT TO SUBSTRATE MUFFIN-
C  TIN CONSTANT (EV)
C
      READ (5,161) VV,VPIS
      IF (IPR.GT.0) WRITE (1,290) VV,VPIS
      VV=VV/27.21
      VPIS=VPIS/27.21
C
C  LMAX= HIGHEST L-VALUE TO BE CONSIDERED
C
      READ (5,200) LMAX
      IF (IPR.GT.0) WRITE (1,355) LMAX
C
C  LMMAX= NO. OF SPHERICAL HARMONICS TO BE CONSIDERED
C
      LMMAX=(LMAX+1)*(LMAX+1)
C
C  NEL= NO. OF CHEMICAL ELEMENTS IN CALCULATION (.LE.5)
C
      READ (5,200) NEL
      IF (NEL.GT.5) THEN
         WRITE (1,*) ' CURRENT CODE NOT DIMENSIONED FOR NEL.GT.5'
c     jcm replace abort with stop CALL ABORT()
         STOP
      ENDIF
C
C  IT1(I)=1,(0) MEANS  DO (NOT) INCLUDE THERMAL VIBRATION EFFECTS
C  FOR ATOM TYPE I
C
      READ (5,200) (IT1(I),I=1,NEL)
      IF (IPR.GT.0) WRITE (1,305) (IT1(I),I=1,NEL)
C
C  THDB= DEBYE TEMPERATURE (KELVIN).
C  AM= ATOMIC MASS (IN NUCLEON UNITS).
C  FPER,FPAR ARE ENHANCEMENT FACTORS FOR THE MEAN SQUARE VIBRATION
C  AMPLITUDES PERPENDICULAR AND PARALLEL TO THE SURFACE, RESP.
C  DR0 IS THE RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE (ANGSTROM).
      DO 48 I=1,NEL
         READ (5,130) THDB(I),AM1(I),FPER1(I),FPAR1(I),DR01(I)
         IF (IPR.GT.0) WRITE (1,325) THDB(I),AM1(I),FPER1(I),FPAR1(I),
     &    DR01(I)
48    CONTINUE
C
C  TI IS THE CALCULATION TEMPERATURE
C
      READ (5,130) TI
      T0=TI
      IF (T0.LT.0.001) T0=0.001
C
C  DRPER AND DRPAR ARE RMS VIBRATION AMPLITUDES PERPENDICULAR AND
C  PARALLEL TO THE SURFACE, RESP., EVALUATED AT REFERENCE TEMPERATURE
C  T0=MAX(TI,0.001 K)
C
      DO 49 I=1,NEL
         DR2=1.546E3*T0/(AM1(I)*THDB(I)*THDB(I))
         DRPER1(I)=SQRT(FPER1(I)*DR2)
         DRPAR1(I)=SQRT(FPAR1(I)*DR2)
         DR01(I)=(DR01(I)/0.529)*(DR01(I)/0.529)
         DR01(I)=DR01(I)*DR01(I)
49    CONTINUE
      L1=LMAX+1
      DO 660 I=1,NPSI
C
C  ES= ENERGIES (HARTREES) AT WHICH PHASE SHIFTS ARE INPUT. LINEAR
C  INTERPOLATION OF THE PHASE SHIFTS WILL OCCUR FOR ACTUAL ENERGIES
C  FALLING BETWEEN THE VALUES OF ES (AND LINEAR EXTRAPOLATION ABOVE
C  THE HIGHEST ES)
C
         READ (5,160) ES(I)
C
C  PHSS STORES THE INPUT PHASE SHIFTS (RADIAN)
C
         DO 661 II=1,NEL
           IO=(II-1)*L1
           READ (5,340) (PHSS(I,L),L=1+IO,L1+IO)
661      CONTINUE
660   CONTINUE
      IF (IPR.GT.0) THEN 
        WRITE (1,350)
        DO 670 I=1,NPSI
        DO 671 II=1,NEL
           IO=(II-1)*L1
           IF (II.EQ.1) THEN
              WRITE (1,280)
              WRITE (1,360) ES(I),(PHSS(I,L),L=1+IO,L1+IO)
           ELSE 
              WRITE (1,361)II,(PHSS(I,L),L=1+IO,L1+IO) 
           ENDIF
671   CONTINUE
670   CONTINUE
       ENDIF
C
C  ARA1 AND ARA2 ARE TWO 2-D BASIS VECTORS OF THE SUBSTRATE LAYER
C  LATTICE. THEY SHOULD BE EXPRESSED IN TERMS OF THE PLANAR CARTESIAN
C  Y- AND Z-AXES (X-AXIS IS PERPENDICULAR TO SURFACE)(ANGSTROM)
C
C      READ (5,160) (ARA1(I),I=1,2)
C      READ (5,160) (ARA2(I),I=1,2)
      READ (5,*) (ARA1(I),I=1,2)
      READ (5,*) (ARA2(I),I=1,2)
      IF (IPR.GT.0) THEN
         WRITE (1,170) (ARA1(I),I=1,2)
         WRITE (1,171) (ARA2(I),I=1,2)
      ENDIF
      DO 460 I=1,2
         ARA1(I)=ARA1(I)/0.529
         ARA2(I)=ARA2(I)/0.529
460   CONTINUE
      TVA=ABS(ARA1(1)*ARA2(2)-ARA1(2)*ARA2(1))
      ATV=2.0*PI/TVA
C
C  RAR1 AND RAR2 ARE THE RECIPROCAL-LATTICE BASIS VECTORS CORRESPONDING
C  TO ARA1 AND ARA2
C
      RAR1(1)=ARA2(2)*ATV
      RAR1(2)=-ARA2(1)*ATV
      RAR2(1)=-ARA1(2)*ATV
      RAR2(2)=ARA1(1)*ATV
C
C  ARB1,ARB2,RBR1,RBR2 ARE EQUIVALENT TO ARA1,ARA2,RAR1,
C  RAR2 BUT FOR AN OVERLAYER 
C
C      READ (5,160) (ARB1(I),I=1,2)
C      READ (5,160) (ARB2(I),I=1,2)
      READ (5,*) (ARB1(I),I=1,2)
      READ (5,*) (ARB2(I),I=1,2)
      DO 465 I=1,2
         CPARB1(I)=ARB1(I)
         CPARB2(I)=ARB2(I)
465   CONTINUE
      IF (IPR.GT.0) THEN
         WRITE (1,170) (ARB1(I),I=1,2)
         WRITE (1,171) (ARB2(I),I=1,2)
      ENDIF
      DO 467 I=1,2
         ARB1(I)=ARB1(I)/0.529
         ARB2(I)=ARB2(I)/0.529
467   CONTINUE
      TVB=ABS(ARB1(1)*ARB2(2)-ARB1(2)*ARB2(1))
      ATV=2.0*PI/TVB
      RBR1(1)=ARB2(2)*ATV
      RBR1(2)=-ARB2(1)*ATV
      RBR2(1)=-ARB1(2)*ATV
      RBR2(2)=ARB1(1)*ATV
C
C calculate NL1 and NL2 from above informations
      CALL NL1NL2(IPR,IPNL1,IPNL2,NL)
C
      RETURN
      END
C ========================================================================
C
C Subroutines READT1 and READT2 input most of the data required
C by the 'conventional' LEED part of TLEED with the exception of
C data involving the composite layer. It is a modified version
C of routine READIN from the Van Hove/Tong LEED package.
C The modifications involve a reordering of the order of data 
C input. Modifications by WANDER.
C Later modifications by DOELL: Subroutine READT was divided into two parts,
C READT1 and READT2. This allows for the calculation of the superlattice
C characterization NL1 and NL2.
C
C
C=========================================================================
C
      SUBROUTINE READT2(ASA,INVECT,IDEG,NL,V,VL,JJS,IPR)
C
      COMPLEX VL(NL,2)
      DIMENSION ARA1(2),ARA2(2),ASA(10,3),ARB1(2)
      DIMENSION ARB2(2),RBR1(2),RBR2(2)
      DIMENSION V(NL,2),JJS(NL,IDEG)
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
166   FORMAT (3F11.8)
172   FORMAT (/' Interlayer Vectors ')
180   FORMAT (/' ASA',1I3,3X,3F7.4)
C
C
C  INVECT IS NUMBER OF INTERLAYER VECTORS IN INPUT FILE
C
      READ (5,*) INVECT
      IF (IPR.GT.0) WRITE (1,172)
      DO 348 J=1,INVECT
         READ (5,166) (ASA(J,I),I=1,3)
         IF (IPR.GT.0) WRITE (1,180) J,(ASA(J,I),I=1,3)
         DO 466 I=1,3
            ASA(J,I)=ASA(J,I)/0.529
466      CONTINUE
348   CONTINUE
C
C  SLIND COMPARES SUBSTRATE LATTICE AND SUPERLATTICE FOR THE BENEFIT
C  OF FMAT
C
      CALL SLIND(V,VL,JJS,NL,IDEG,2.0E-4)
      RETURN
      END
C =========================================================================
C
c subroutine nl1nl2 determines NL, NL1 and NL2 from the vectors ARA1,
c ARA2, ARB1 and ARB2.
c author: R. Doell, 2/15/96
C
C =========================================================================

      SUBROUTINE NL1NL2(IPR,IPNL1,IPNL2,NL)

      REAL ARA1(2),ARA2(2),ARB1(2),ARB2(2),VECT(2),VECTS(2)
      REAL RBR1(2),RBR2(2)
      REAL DETA,DETB,DET1,DET2,RHELP
      INTEGER IPR,IPNL1,IPNL2,NL1,NL2,NL,I1,I2,INL1,INL2
      INTEGER IC1,IC2,IC3,IC4,COEFF1,COEFF2
      INTEGER IMAT(-50:50,-50:50)

      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2

c-----------------------------------------------------------------------

c>> calculate DETB = det(ARB1,ARB2)
      DETB = ARB1(1)*ARB2(2) - ARB1(2)*ARB2(1)
         IF (DETB .EQ. 0.) THEN
             WRITE (1,'(A30,A20)') ' ERROR: ARB1 and ARB2 are not ',
     &                             'linearly independent'
             STOP
         ENDIF

cccc       calculate DETA = det(ARA1,ARA2)
           DETA = ARA1(1)*ARA2(2) - ARA1(2)*ARA2(1)
           IF (DETA .EQ. 0.) THEN
             WRITE (1,'(A30,A20)') ' ERROR: ARA1 and ARA2 are not ',
     &                             'linearly independent'
             STOP
         ENDIF

c-----------------------------------------------------------------------

c>> calculate NL
         TVA = ABS (ARA1(1)*ARA2(2) - ARA1(2)*ARA2(1))
         TVB = ABS (ARB1(1)*ARB2(2) - ARB1(2)*ARB2(1))
      NL = NINT (TVB / TVA)
         IF (ABS(TVB/TVA-FLOAT(NL)) .GT. 0.00001) THEN
             WRITE (1,'(A31,A24)') ' ERROR: size of superstructure ',
     &                             'does not match substrate'
             STOP
         ENDIF
         IF (NL .GT. 100) THEN
             WRITE (1,'(A28,A35)') ' ERROR: correct dimension of',
     &             ' I1PT and I2PT in SUBROUTINE NL1NL2'
             STOP
         ENDIF

c-----------------------------------------------------------------------

c>> check NL1-NL2 combinations consistent with NL
      DO 100 INL1=1,NL
         RHELP = FLOAT(NL)/FLOAT(INL1) - FLOAT(NL/INL1)

        IF (RHELP .LT. 0.0001) THEN
            INL2 = NL/INL1

cccc          preset matrix IMAT
              DO 102 I1=-50,50
              DO 101 I2=-50,50
                 IMAT(I1,I2)=0
101           CONTINUE
 102          CONTINUE


c>>   check all points in the area spanned by +- 2*ARB1 +- 2*ARB2
c>>   around the origin

          DO 211 IC1=-2,2
             DO 201 IC2=-2,2
               VECTS(1) = FLOAT(IC1)*ARB1(1) + FLOAT(IC2)*ARB2(1)
               VECTS(2) = FLOAT(IC1)*ARB1(2) + FLOAT(IC2)*ARB2(2)

               DO 212 IC3=0,INL1-1
                  DO 202 IC4=0,INL2-1
                     VECT(1) = VECTS(1)+FLOAT(IC3)*ARA1(1)+
     &                    FLOAT(IC4)*ARA2(1)
                     VECT(2) = VECTS(2)+FLOAT(IC3)*ARA1(2)+
     &                    FLOAT(IC4)*ARA2(2)


c>>    solve the equation system
cccc        VECT(1) = COEFF1*ARA1(1) + COEFF2*ARA2(1)
cccc        VECT(2) = COEFF1*ARA1(2) + COEFF2*ARA2(2)

cccc        calculate determinants of (ARA2,VECT) and (ARA1,VECT)
                     DET1 = ARA2(1)*VECT(2) - ARA2(2)*VECT(1)
                     DET2 = ARA1(1)*VECT(2) - ARA1(2)*VECT(1)

                     COEFF1 = NINT(-DET1 / DETA)
                     COEFF2 = NINT(DET2 / DETA)

                     IF ((COEFF1.GT.50).OR.(COEFF2.GT.50).OR.
     &                    (COEFF1.LT.-50).OR.(COEFF2.LT.-50)) THEN
                        WRITE (1,'(A28,A25)') ' ERROR: check ',
     &                   'dimensions of IMAT in SUBROUTINE NL1NL2'
                        STOP
                     ENDIF


c>>   check whether point was pointed at before, if so:
c>>   current INL1 and INL2 don't work...

                     IF (IMAT(COEFF1,COEFF2) .GT. 0) THEN
                        GOTO 100
                     ELSE
                        IMAT(COEFF1,COEFF2) = 1
                     ENDIF
 202              CONTINUE
 212           CONTINUE
 201        CONTINUE
 211     CONTINUE

c-----------------------------------------------------------------------

c>> set NL1 and NL2 to INL1 and INL2, resp.
      GOTO 300
      ENDIF

100   CONTINUE
300   CONTINUE

      NL1 = INL1
      NL2 = INL2

c-----------------------------------------------------------------------

c>> perform some checks
      IF (NL .NE. (NL1*NL2)) THEN
             WRITE (1,'(A30,A25)') ' ERROR: check lattice vectors ',
     &                             'ARA1, ARA2, ARB1 and ARB2'
             STOP
      ENDIF
      IF (NL1 .GT. IPNL1) THEN
             WRITE (1,'(A28,A17,I3)') ' ERROR: the parameter IPNL1 ',
     &         'must be at least ', NL1
             STOP
      ENDIF
c
      IF (NL2 .GT. IPNL2) THEN
             WRITE (1,'(A28,A17,I3)') ' ERROR: the parameter IPNL2 ',
     &         'must be at least ', NL2
             STOP
      ENDIF

c-----------------------------------------------------------------------

c>> write values for NL1 and NL2, if desired
      IF (IPR .GT. 0) THEN
         WRITE (1,'(A6,I3,A9,I3)') ' NL1 =',NL1,'    NL2 =',NL2
      ENDIF

c-----------------------------------------------------------------------

      RETURN
      END
C =========================================================================
C 
C  Subroutine SB computes spherical Bessel functions.
C 
C Modified version of routine SB from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C If the modulus of X is small one should use Miller's device to compute
C accurate real part of the Hankel function. Otherwise the upward recursion
C relation works better.
C Modifications by Barbieri
C
C =========================================================================
C
      SUBROUTINE SB(X,HH,N3)
C
      COMPLEX A,B,C,HH,X,ZJ0
      DIMENSION HH(N3)
      COMPLEX ZJ(20),ZY(20),ZJS(20)
C
330   FORMAT ('SB INFINITE')
C  
C LSTART determines where to start the downward recursion relation for
C ZJ using Miller's device. (see BESSEL)
      LSTART=13
      A=(0.0,1.0)
      C=X
      B=A*C
      F=CABS(X)
      IF (F.LE.1.0E-38) THEN
         WRITE (1,330)
      ELSE
C This has been tested to work for up to L=10
       IF(F.LE.2.) THEN
C Upward recurrence for ZY (recurrence relation is stable for ZY)
         ZY(1)=-CCOS(C)/C
         ZY(2)=-CCOS(C)/C**2 -CSIN(C)/C
         DO 1530 J=3,N3
            ZY(J)=(2.0*(J-2)+1.0)/C*ZY(J-1)-ZY(J-2)
1530     CONTINUE
C Donward recurrence for ZJ (see BESSEL)
         ZJS(LSTART)=(0.,0.)
         ZJS(LSTART-1)=(1.0,0.)
C
C GENERATE BESSEL FUNCTIONS FOR L=0,LSTART-2 BY BACKWARD
C RECURRENCE.
C
         DO 100 IL=0,LSTART-3
            L=LSTART-3-IL
            ZJS(L+1)=FLOAT(2*(L+1)+1)*ZJS(L+2)/X-ZJS(L+3)
100      CONTINUE
C
C EVALUATE BESSEL FN FOR L=0 EXPLICITLY
C
         ZJ0=CSIN(X)/X
C
C NORMALISE PREVIOUSLY CALCULATED BESSEL FNS
C
         DO 110 L=0,N3-1
            ZJ(L+1)=ZJS(L+1)*ZJ0/ZJS(1)
            HH(L+1)=ZJ(L+1)+A*ZY(L+1)
110      CONTINUE
       ELSE
         HH(1)=CEXP(B)*(-A)/C
         HH(2)=CEXP(B)*(1.0/(-C)-A/C**2)
         DO 1520 J=3,N3
            HH(J)=(2.0*(J-2)+1.0)/C*HH(J-1)-HH(J-2)
1520     CONTINUE
       ENDIF
      ENDIF
      RETURN
      END
C =========================================================================
C
C  Subroutine SETPOS is too trivial to explain 
C
C =========================================================================
      SUBROUTINE SETPOS(FPOS,VPOS,LPSAW,LL,I,NCL,NLMAX,NST1)
      DIMENSION FPOS(NLMAX,3),VPOS(NST1,NLMAX,3),LPSAW(NST1,NLMAX)
         DO 549 J=1,3
            FPOS(I,J)=VPOS(NCL,I,J)
549      CONTINUE
         LL=LPSAW(NCL,I)      
      RETURN
      END
C =========================================================================
C
C  Subroutine SH computes spherical harmonics.
C
C Modified version of routine SH from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE SH(NHARM,Z,FI,Y)
C
      COMPLEX A,AA,ANOR,BB,BB1,Q1A,Y,YSTAR,YY,Z,ZNW
      DIMENSION Y(NHARM,NHARM)
C
      AM=1.0
      A=(0.0,1.0)*FI
      RZ=REAL(Z)
      ZNW=CSQRT((1.0,0.0)-Z*Z)
      ANORA=0.2820948
      YY=(1.0,0.0)
      DO 1510 L=1,NHARM
         BM=FLOAT(L)-1.0
         Q1A=CEXP(BM*A)
         RAWF=1.0
         DO 48 IAW=1,L+1
            RAWF=RAWF*(-1.0)
48       CONTINUE
         ANOR=ANORA*Q1A*RAWF
         ANORA=ANORA*SQRT(1.0+0.5/(BM+1.0))/(2.0*BM+1.0)
         IF (.NOT.((ABS(ABS(RZ)-1.0).LT.1.0E-10).AND.L.EQ.1)) THEN
            YY=1.0
            DO 49 IAW=1,INT(BM+0.1)
               YY=YY*ZNW
49          CONTINUE
         ENDIF
         AA=AM*YY
         AM=AM*(2.0*BM+1.0)
         BB=AM*Z*YY
         DO 1500 LL=L,NHARM
            BN=FLOAT(LL)
            Y(LL,L)=AA*ANOR
            IF (L.NE.1) THEN
               YSTAR=Y(LL,L)/(Q1A*Q1A)
               Y(L-1,LL)=YSTAR
            ENDIF
            BB1=BB
            BB=BB1*Z+(BB1*Z-AA)*(BN+BM)/(BN-BM+1.0)
            AA=BB1
            ANOR=ANOR*SQRT((1.0+1.0/(BN-0.5))*((BN-BM)/(BN+BM)))
1500     CONTINUE
1510  CONTINUE
      RETURN
      END
C =========================================================================
C
C  Subroutine SH computes spherical harmonics.
C
C Modified version of routine SH from the VAN HOVE/TONG LEED package.
C Modifications by Barbieri.
C
C =========================================================================
C
      SUBROUTINE SHR(NHARM,Z,FI,Y)
C
      COMPLEX A,Q1A,Q2A,Y,YSTAR
      DIMENSION Y(NHARM,NHARM)
C
      AM=1.0
      A=(0.0,1.0)*FI
      RZ=Z
      ZNW=SQRT(1.0-Z*Z)
      ANORA=0.2820948
      YY=1.0
      RAWF=-1.0
      DO 1510 L=1,NHARM
         BM=FLOAT(L)-1.0
         Q1A=CEXP(BM*A)
         Q2A=(Q1A*Q1A)
         RAWF=RAWF*(-1.0)
         ANOR=ANORA*RAWF
         ANORA=ANORA*SQRT(1.0+0.5/(BM+1.0))/(2.0*BM+1.0)
         IF (.NOT.((ABS(ABS(RZ)-1.0).LT.1.0E-10).AND.L.EQ.1)) THEN
            YY=1.0
            DO 49 IAW=1,INT(BM+0.1)
               YY=YY*ZNW
49          CONTINUE
         ENDIF
         AA=AM*YY
         AM=AM*(2.0*BM+1.0)
         BB=AM*Z*YY
         DO 1500 LL=L,NHARM
             BN=FLOAT(LL)
            Y(LL,L)=AA*ANOR*Q1A
            IF (L.NE.1) THEN
               YSTAR=Y(LL,L)/Q2A
               Y(L-1,LL)=YSTAR
            ENDIF
            BB1=BB
            BB=BB1*Z+(BB1*Z-AA)*(BN+BM)/(BN-BM+1.0)
            AA=BB1
            ANOR=ANOR*SQRT((1.0+1.0/(BN-0.5))*((BN-BM)/(BN+BM)))
1500     CONTINUE
1510  CONTINUE
      RETURN
      END
C===================================================================== 
C 
C Subroutine SHORT dumps the truncated Q matrix to the transfer file.
C
C Input Parameters;
C =================
C 
C Q             =  Q MTRIX FO CURRENT BEAM, LAYER AND ENERGY
C LSMMAX        =  (LSMAX+1)**2
C IPR           =  PRINT CONTROL PARAMETER
C E,VV          =  CURRENT ENERGY
C AK2M,AK3M     =  PARALLEL COMPONENTS OF MOMENTUM FOR EACH EXIT BEAM
C NT0           =  NUMBER OF EXIT BEAMS IN CALCULATION
C NEXIT         =  INDEX OF CURRENT EXIT BEAM
C ICUT          =  NUMBER OF INDEPENDANT ELEMENTS IN Q
C MICUT,MJCUT   =  LOCATION OF INDEPENDANT ELEMENTS IN Q MATRIX
C
C =========================================================================
C
C      SUBROUTINE SHORT(Q,LSMMAX,NFILE,E,VV,AK2M,AK3M,NT0,NEXIT, changed by ZZ on 09/07/04
      SUBROUTINE SHORT_own(Q,LSMMAX,NFILE,E,VV,AK2M,AK3M,NT0,NEXIT,
     & ICUT,QCUT,MICUT,MJCUT)
C
      COMPLEX Q(LSMMAX,LSMMAX),QCUT(ICUT)
      DIMENSION MICUT(ICUT),MJCUT(ICUT)
      DIMENSION AK2M(NT0),AK3M(NT0)
C
!1000  FORMAT (500(2E13.5,/))
C
      A=2.0*E-2.0*VV-AK2M(NEXIT)**2-AK3M(NEXIT)**2
      DO 100 IC=1,ICUT
         I=MICUT(IC)
         IP=MJCUT(IC)
         QCUT(IC)=Q(I,IP)
         IF (IP.NE.I) QCUT(IC)=Q(I,IP)+Q(IP,I)
100   CONTINUE
C      IF (A.GT.0) WRITE (NFILE,1000) (QCUT(I),I=1,ICUT)
      IF (A.GT.0) WRITE (NFILE) (QCUT(I),I=1,ICUT)
      RETURN
      END
C=========================================================================
C
C Slind sets up a matrix JJS(JS,J) containing details of
C how the sublattices JS are transformed into one another by
C rotations through J*2.0*PI/IDEG. V(JS,2) contains the adding      
C vectors defining the sublattices JS in polar form.
C
C  AUTHOR  PENDRY.
C
C Parameter List;
C ===============
C
C V                 =    Vectors pointing to different sublattices (In polar
C                       coordinates).
C VL                =    Working space.
C JJS               =    Output relationship of sublattices under rotation.
C NL                =    No. of sublattices.
C IDEG              =    Degree of symmetry of lattice  IDEG-fold rotation axis.
C EPSD              =    Small parameter.
C
C IN COMMON BLOCKS
C
C BR1,BR2           =    Basis vectors of substrate layer lattice.
C AR1,AR2,RAR1,RAR2 =    Basis vectors of superlattice in direct and
C                        reciprocal space.
C NL1,NL2           =    Superlattice characterization codes (See TLEED1).
C
C=========================================================================
C       
      SUBROUTINE SLIND(V,VL,JJS,NL,IDEG,EPSD)
C
      COMPLEX VL,VLA,VLB,CI
      COMPLEX CEXP
      DIMENSION V(NL,2),JJS(NL,IDEG),VL(NL,2),AR1(2),AR2(2)
      DIMENSION BR1(2),BR2(2),RAR1(2),RAR2(2)
C
      COMMON /SL/BR1,BR2,AR1,AR2,RAR1,RAR2,NL1,NL2
C
      CI=CMPLX(0.0,1.0)
      PI=3.14159265
C
C  SET UP VECTORS V DEFINING SUBLATTICES AND QUANTITIES VL FOR LATER
C  REFERENCE.
C
      I=1
      S1=0.0
      DO 560 J=1,NL1
         S2=0.0
         DO 550 K=1,NL2
            ADR1=S1*BR1(1)+S2*BR2(1)
            ADR2=S1*BR1(2)+S2*BR2(2)
            V(I,1)=SQRT(ADR1*ADR1+ADR2*ADR2)
            V(I,2)=0.0
            IF (V(I,1).GT.0) V(I,2)=ATAN2(ADR2,ADR1)
            VL(I,1)=CEXP(CI*(ADR1*RAR1(1)+ADR2*RAR1(2)))
            VL(I,2)=CEXP(CI*(ADR1*RAR2(1)+ADR2*RAR2(2)))
            I=I+1
            S2=S2+1.0
550      CONTINUE
         S1=S1+1.0
560   CONTINUE
C
C  ROTATE EACH VECTOR V AND FIND TO WHICH V IT BECOMES EQUIVALENT IN
C  TERMS OF THE QUANTITIES VL. THIS EQUIVALENCE MEANS BELONGING TO THE
C  SAME SUBLATTICE.
C
      AINC=2.0*PI/FLOAT(IDEG)
      DO 591 I=1,NL
         ADR=V(I,1)
         ANG=V(I,2)
         DO 590 K=1,IDEG
            ANG=ANG+AINC
            CANG=COS(ANG)
            SANG=SIN(ANG)
            A=RAR1(1)*CANG+RAR1(2)*SANG
            B=RAR2(1)*CANG+RAR2(2)*SANG
            VLA=CEXP(CI*ADR*A)
            VLB=CEXP(CI*ADR*B)
            DO 570 J=1,NL
               TEST=CABS(VLA-VL(J,1))+CABS(VLB-VL(J,2))
               IF (TEST-5.0*EPSD.LE.0) GOTO 580
570         CONTINUE
580         JJS(I,K)=J
590      CONTINUE
591   CONTINUE
      RETURN
      END
C===========================================================================
C
C Subroutine SORT reorganizes the coordinates of the compiste layer in the     
C same manner as SRTLAY does.
C
C Modified version of SRTLAY from the Van Hove/Tong LEED package. Modifications
C by WANDER.
C
C===========================================================================
C
      SUBROUTINE SORT(POS,POSS,NLMX,NLAY)
C
      DIMENSION POS(NLMX,3),POSS(NLMX,3),POSA(3)
C
1000  FORMAT (/,10X,'COORDINATES AFTER SORTING',/,13X,'X',14X,'Y',14X,
     & 'Z')
1010  FORMAT (7X,F10.4,2(5X,F10.4))
C
      DO 1112 I=1,NLAY
         DO 1 J=1,3
            POSS(I,J)=POS(I,J)
1        CONTINUE
1112  CONTINUE
C
C  ANALYSE ORDER OF SUBPLANE POSITIONS ALONG X-AXIS AND REORDER IN
C  ASCENDING POSITION ALONG +X-AXIS (EQUALLY POSITIONED SUBPLANES ARE
C  NOT PERMUTED)
C
      NLAY1=NLAY-1
      DO 7 I=1,NLAY1
         II=I+1
         KM=I
         PM=POSS(I,1)
         DO 2 K=II,NLAY
            IF (POSS(K,1).LT.PM) THEN
               PM=POSS(K,1)
               KM=K
            ENDIF
2        CONTINUE
         IF (KM.NE.I) THEN
            DO 3 J=1,3
               POSA(J)=POSS(KM,J)
3           CONTINUE
            DO 5 KK=II,KM
               K=KM+II-KK
               DO 4 J=1,3
                  POSS(K,J)=POSS(K-1,J)
4              CONTINUE
5           CONTINUE
            DO 6 J=1,3
               POSS(I,J)=POSA(J)
6           CONTINUE
         ENDIF
7     CONTINUE
      WRITE (1,1000)
      DO 1111 I=1,NLAY
         WRITE (1,1010) (POSS(I,J)*0.529,J=1,3)
         DO 1113 II=1,3
            POS(I,II)=POSS(I,II)
1113     CONTINUE
1111  CONTINUE
      RETURN
      END
C =========================================================================
C
C  Subroutine SPHRM computes spherical harmonics.
C
C Parameter List;
C ===============
C
C LMAX          =   LARGEST VALUE OF L
C YLM           =   OUTPUT COMPLEX SPHERICAL HARMONIC
C LMMAX         =   (LMAX+1)**2
C CT            =   COS(THETA)  (COMPLEX)
C ST            =   SIN(THETA)  (COMPLEX)
C CF            =   CEXP(I*FI)
C
C Modified version of PENDRY'S routine SPHRM. Modifications by WANDER,BARBIERI
C
C =========================================================================
C
      SUBROUTINE SPHRM(LMAX,YLM,LMMAX,CT,ST,CF)
C
      COMPLEX YLM
      COMPLEX CT,ST,CF,SF,SA
      DIMENSION FAC1(20),FAC3(20),FAC2(300),YLM(LMMAX)
C
      PI=4.0*ATAN(1.0)
      LM=0
      CL=0.0
      A=1.0
      B=1.0
      ASG=1.0
      LL=LMAX+1
      DO 550 L=1,LL
         FAC1(L)=ASG*SQRT((2.0*CL+1.0)*A/(4.0*PI*B*B))
         FAC3(L)=SQRT(2.0*CL)
         CM=-CL
         LN=L+L-1
         DO 540 M=1,LN
            LO=LM+M
            FAC2(LO)=SQRT((CL+1.0+CM)*(CL+1.0-CM)/((2.0*CL+3.0)*(2.0*CL
     &       +1.0)))
            CM=CM+1.0
540      CONTINUE
         CL=CL+1.0
         A=A*2.0*CL*(2.0*CL-1.0)/4.0
         B=B*CL
         ASG=-ASG
         LM=LM+LN
550   CONTINUE
      LM=1
      CL=1.0
      ASG=-1.0
      SF=CF
      SA=CMPLX(1.0,0.0)
      YLM(1)=CMPLX(FAC1(1),0.0)
      DO 560 L=1,LMAX
         LN=LM+L+L+1
         YLM(LN)=FAC1(L+1)*SA*SF*ST
         YLM(LM+1)=ASG*FAC1(L+1)*SA*ST/SF
         YLM(LN-1)=-FAC3(L+1)*FAC1(L+1)*SA*SF*CT/CF
         YLM(LM+2)=ASG*FAC3(L+1)*FAC1(L+1)*SA*CT*CF/SF
         SA=ST*SA
         SF=SF*CF
         CL=CL+1.0
         ASG=-ASG
         LM=LN
560   CONTINUE
      LM=1
      LL=LMAX-1
      DO 580 L=1,LL
         LN=L+L-1
         LM2=LM+LN+4
         LM3=LM-LN
         DO 570 M=1,LN
            LO=LM2+M
            LP=LM3+M
            LQ=LM+M+1
            YLM(LO)=-(FAC2(LP)*YLM(LP)-CT*YLM(LQ))/FAC2(LQ)
570      CONTINUE
         LM=LM+L+L+1
580   CONTINUE
      RETURN
      END
C  file LEEDSATL.SB4  Feb. 29, 1996
C
C**************************************************************************
C  Symmetrized Automated Tensor LEED (SATLEED):  subroutines, part 4
C  Version 4.1 of Automated Tensor LEED
C
C======================================================================  
C Subroutine SPLINE compute the second derivative Y2 of the tabulated phase shifts
C contained in the array Y corresponding to the energy values contained in X.
C 
C This subroutine is lifted from  Numerical Recipes 
C
C References 'NUMERICAL RECIPES' W.H.Press,B.P.Flannery, et al. Cambridge
C             University Press
C
C==============================================================================
      SUBROUTINE SPLINE(x,y,n,yp1,ypn,y2)
      parameter (nmax=500)
      dimension x(n),y(n),y2(n),u(nmax)
      if (yp1.gt..99e10) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e10) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
C ============================================================================
C Subroutine SPLINT performs the cubic spline interpolation giving the phase
C shift y at energy x, from the tabulated arrays ya and xa.
C This subroutine is lifted from  Numerical Recipes but KLO and KHI are
C now input parameters and the routine is accordingly modified
C
C References 'NUMERICAL RECIPES' W.H.Press,B.P.Flannery, et al. Cambridge
C             University Press
C
C Modifications by BARBIERI
C==============================================================================
      SUBROUTINE SPLINT(xa,ya,y2a,n,x,y,yp,klo,khi)
      dimension xa(n),ya(n),y2a(n)
      h=xa(khi)-xa(klo)
cjcm      if (h.eq.0.) pause 'bad xa input. in SPLINT'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      yp=(ya(khi)-ya(klo))/h -(3.*a**2 -1.)*h*y2a(klo)/6. +
     *  (3.*b**2 -1.)*h*y2a(khi)/6. 
      return
      end
C =========================================================================
C
C  Subroutine TAUMAT2 computes the matrix TAU (Intra-subplane multiple
C  scattering in (L,M)-representation) for a single Bravais-lattice
C  layer for each of several chamical elements. The matrix elements are
C  ordered thus: The (L,M) sequence is the 'symmetrized' one (CF. book).
C  TAU is split into block-diagonal parts corresponding to even L+M
C  and odd L+M of dimensions LEV and LOD, resp. In the matrix TAU the
C  part for L+M= odd is left-justified and the matrix reduced in width
C  for storage efficiency. The TAU matrices for different chemical
C  elements are arranged under each other in columnar fashion;at the top
C  is the TAU based on the atomic T-matrix elements TSF(1,L), followed by
C  the TAU's based on TSF(2,L) etc.
C
C  It also collects the matrix inverse of TAU needed 
C  for the tensor 
C 
C Parameter List;
C ===============
C
C TAU          =  OUPUT MATRIX (SEE ABOVE)(symmetrized order)
C TAUINV       =   t^-1*(1-tG)=tau^-1=TAUINV)(symmetrized order)
C                  we have to collect it for all the different elements
C                  in all composite layers
C LMT          =  NTAU*LMMAX
C NTAU         =  NUMBER OF CHEMICAL ELEMENTS TO BE USED
C X            =  WORKING SPACE
C TSF          =  ATOMIC T-MATRIX ELEMENTS
C FLMS         =  LATTICE SUMS FROM SUBROUTINE FMAT
C NL           =  NUMBER OF SUBLATTICES TO BE CONSIDERED IN FMAT (ONLY THE
C                FIRST LATTICE SUM IS USED BY TAUMAT)
C CLM              =  CLEBSCH-GORDAN COEFFICIENTS FROM SUBROUTINE CELMG
C NLM          =  DIMENSION OF CLM
C LXI          =  PERMUTATION OF (L,M) SEQUENCE FROM SUBROUTINE LXGENT
C NT           =  NUMBER OF BEAMS IN CALCULATION AT CURRENT ENERGY
C PQ           =  BEAM LIS
C NA           =  OFFSET IN CURRENT CALL TO MTINVT
C NLL          =  NUMBER OF SUBLATTICES INCLUDED IN CURRENT CALL TO MTINVT
C FLM              =  WORKING SPACE
C
C In Common Blocks;
C =================
C
C E,VPI        =  CURRENT COMPLEX ENERGY
C Modified version of routine TAUMAT from the VAN HOVE/TONG LEED package.
C
C =========================================================================
C
      SUBROUTINE TAUMAT2(TAU,LMT,NTAU,X,LEV,LEV2,LOD,TSF,LMMAX,LMAX,
     & FLMS,NL,KLM,LM,CLM,NLM,LXI,NT,PQ,NA,NLL,FLM,NTAUSH,TAUINV,
     & NLTU,LAN)
C
      DIMENSION CLM(NLM),LXI(LMMAX),PQ(2,NT),LAN(LMMAX)
      DIMENSION BR1(2),BR2(2),AR1(2),AR2(2),RAR1(2),RAR2(2)
      COMPLEX AK,CZ,TAU(LMT,LEV),X(LEV,LEV2),TSF(6,16),FLMS(NL,KLM)
      COMPLEX DET,FLM(KLM),CI,XA
      COMPLEX TAUINV(NLTU,LEV)
C      COMPLEX TCH(50,50)
C
      COMMON E,AK2,AK3,VPI
      COMMON /SL/BR1,BR2,AR1,AR2,RAR1,RAR2,NL1,NL2
C
      CZ=(0.0,0.0)
      CI=(0.0,1.0)
      AK=-0.5/CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      DO 100 K=1,KLM
         FLM(K)=CZ
100   CONTINUE
      BK2=PQ(1,1+NA)
      BK3=PQ(2,1+NA)
      JS=1
      S1=0.
      DO 130 J=1,NL1
         S2=0.
         DO 120 K=1,NL2
            ADR1=S1*BR1(1)+S2*BR2(1)
            ADR2=S1*BR1(2)+S2*BR2(2)
            ABR=ADR1*BK2+ADR2*BK3
            XA=CEXP(ABR*CI)
            DO 110 I=1,KLM
               FLM(I)=FLM(I)+FLMS(JS,I)*XA
110         CONTINUE
            IF (NLL.EQ.1) GOTO 144
            JS=JS+1
            S2=S2+1.
120      CONTINUE
         S1=S1+1.
130   CONTINUE
144   CONTINUE
      DO 141 IT=1,NTAU
         IT2=IT+NTAUSH
         DO 13 IL=1,2
            DO 142 I=1,LEV
               DO 1 J=1,LEV2
                  X(I,J)=CZ
1              CONTINUE
142         CONTINUE
            LL=LOD
            IF (IL.EQ.2) LL=LEV
C
C  GENERATE MATRIX 1-X FOR L+M= ODD (IL=1), LATER FOR L+M= EVEN (IL=2)
C  The output X is equal to 1-tG of Van Hove-Tong eq. 45
C
            CALL XMT(IL,FLM,X,LEV,LL,TSF,IT2,LM,LXI,LMMAX,KLM,
     &       CLM,NLM,NST)
            LD2=(IT2-1)*LMMAX
            IF (IL.EQ.1) LD2=LD2+LEV
            DO 1422 I=1,LL
               IF(IL.EQ.1) LLL=LAN(I+LEV)+1
               IF(IL.EQ.2) LLL=LAN(I)+1
               DO 1100 J=1,LL
                     TAUINV(I+LD2,J)=X(I,J)
C                  IF (CABS(TSF(IT2,LLL)).GE.1.0E-06) THEN
C                     TAUINV(I+LD2,J)=X(I,J)/(AK*TSF(IT2,LLL))
C                  IF (CABS(TSF(IT2,LLL)).GE.1.0E-06.AND
C     &               .I.EQ.J) THEN
C                     TAUINV(I+LD2,J)=1./(AK*TSF(IT2,LLL))
C                  ELSE
C                     TAUINV(I+LD2,J)=CZ
C                  ENDIF
1100           CONTINUE
1422         CONTINUE
C
C  PREPARE QUANTITIES INTO WHICH INVERSE OF 1-X WILL BE MULTIPLIED
C
            IF (IL.LT.2) THEN
               IS=0
               LD1=0
               L=1
7              LD=LD1+1
               LD1=LD+L-1
               DO 8 I=LD,LD1
                  X(I,I+LOD)=AK*TSF(IT2,L+1)
8              CONTINUE
               L=L+2
               IF (L.LE.LMAX) GOTO 7
               IS=IS+1
               L=2
               IF (IS.LE.1) GOTO 7
            ELSE
               IS=0
               LD1=0
               L=0
3              LD=LD1+1
               LD1=LD+L
               DO 4 I=LD,LD1
                  X(I,I+LEV)=AK*TSF(IT2,L+1)
4              CONTINUE
               L=L+2
               IF (L.LE.LMAX) GOTO 3
               IS=IS+1
               L=1
               IF (IS.LE.1) GOTO 3
            ENDIF
            LL2=LL+LL
C
C  PERFORM INVERSION AND MULTIPLICATION
C
            CALL CXMTXT(X,LEV,LL,LL,LL2,MARK,DET,-1)
            LD=(IT-1)*LMMAX
            IF (IL.EQ.1) LD=LD+LEV
            DO 143 I=1,LL
               DO 11 J=1,LL
C
C  PUT RESULT IN TAU 
C
                  TAU(LD+I,J)=X(I,J+LL)
11             CONTINUE
143         CONTINUE
13       CONTINUE
141   CONTINUE
      RETURN
      END
C =========================================================================
C
C  Subroutine TAUT2 produces the quantites YLM(G) in the proper
C  symmetrized order
C
C Parameter List;
C ===============
C
C TAUG,TAUGM    =  OUTPUT RESULTING VECTORS FOR K(G+) AND K(G-),RESP.
C LTAUG         =  (NTAU+NINV)*LMMAX (NINV TO BE 0 IN CALL FROM MTINV).
C LT            =  PERMUTATION OF (LM) ORDER, FROM LXGENT.
C LMNI          =  NINV*LMMAX, USED ONLY IF LINV=1.
C LMNI2         =  2*LMNI, USED ONLY IF LINV=1.
C JGP           =  CURRENT INCIDENT BEAM.
C
C
C Modified version of routine TAUY from the VAN HOVE/TONG LEED package.
C
C =========================================================================
C
      SUBROUTINE TAUT2(TAUG,TAUGM,LTAUG,LMT,LEV,CYLM,NT,LMMAX,LT,
     & NTAU,LOD,LEE,LOE,JGP,TSF,LAN,NTAUSH)
C
      COMPLEX CZ,CF,AK
      COMPLEX CYLM(NT,LMMAX)
      COMPLEX TAUG(LTAUG),TAUGM(LTAUG),TSF(6,16)
      DIMENSION LT(LMMAX),LAN(LMMAX)
      COMMON E,AK2,AK3,VPI
C
      CZ=(0.0,0.0)
      AK=-0.5/CSQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
C
C  PERFORM MATRIX PRODUCT TAU*YLM(G+-) FOR EACH CHEMICAL ELEMENT
C  for this modified version use TAU=t instead 
C
      DO 275 I=1,NTAU
         IT2=I+NTAUSH
         IS=(I-1)*LMMAX
         DO 250 JLM=1,LEV
               LLL=LAN(JLM)+1
               KLP=LT(JLM)
               CF=CYLM(JGP,KLP)*(AK*TSF(IT2,LLL))
               IF (JLM.GT.LEE) THEN
                  TAUG(IS+JLM)=-CF
                  TAUGM(IS+JLM)=-CF
C                  TAUG(IS+JLM)=-CYLM(JGP,KLP)
C                  TAUGM(IS+JLM)=-CYLM(JGP,KLP)
               ELSE
                  TAUG(IS+JLM)=CF
                  TAUGM(IS+JLM)=CF
               ENDIF
250      CONTINUE
         IS=IS+LEV
         DO 270 JLM=1,LOD
               LLL=LAN(JLM+LEV)+1
               KLP=LT(JLM+LEV)
               CF=CYLM(JGP,KLP)*(AK*TSF(IT2,LLL))
               IF (JLM.GT.LOE) THEN
                  TAUG(IS+JLM)=-CF
                  TAUGM(IS+JLM)=CF
               ELSE
                  TAUG(IS+JLM)=CF
                  TAUGM(IS+JLM)=-CF
               ENDIF
270      CONTINUE
275   CONTINUE
      RETURN
      END
C =========================================================================
C
C  Subroutine TFOLD produces, for subroutine MSMF, the individual
C  reflection and transmission matrix elements for a Bravais-Lattice
C  layer. A Debye-Waller factor is included explicitly, if required 
C  (IT.NE.0).
C
C Parameter List;
C ===============
C
C JG                    =   INDEX OF CURRENT SCATTERED BEAM.
C NA                    =   OFFSET OF PRESENT BEAM SET IN LIST PQ.
C YLM                   =   RESULT OF (1-X)**(-1)*Y*T, FROM MSMF.
C CYLM                  =   SET OF SPHERICAL HARMONICS, FROM MSMF.
C N                     =   NOT USED.
C RA,TA                 =   OUTPUT MATRIX ELEMENTS.
C
C In Common Blocks;
C =================
C
C GP                    =   CURRENT INCIDENT BEAM.
C LM                    =   LMAX+1.
C IT                    =   0  NO EXPLICIT DEBYE-WALLER FACTOR TO BE USED 
C                         (COMPLEX PHASE SHIFTS INSTEAD).
C IT                    =   1  USE EXPLICIT DEBYE-WALLER FACTORS.
C TEMP,T0,DRX,DRY,DR0   =   THERMAL VIBRATION DATA, SEE MAIN PROGRAM.
C
C This is a modified version of routine MFOLD from the VAN HOVE/TONG
C LEED package. Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE TFOLD(JG,NA,E,VPI,AK2,AK3,YLM,CYLM,LMMAX,PQ,NT,RA,
     & TA)
C
      COMPLEX CYLM,YLM,CZ,RU,CI,RA,TA,ST,SM,SL,CM,CY,CTR,CTT
      COMPLEX R,T
      COMPLEX EDW
      DIMENSION GP(2),PQ(2,NT),EDW(2)
      DIMENSION CYLM(NT,LMMAX),YLM(LMMAX)
      DIMENSION G(2)
C
      COMMON /MFB/GP,LM,LAY
      COMMON /BT/IT,TEMP,T0,DRX,DRY,DR0
C
      CZ=(0.0,0.0)
      RU=(1.0,0.0)
      CI=(0.0,1.0)
      IGC=0
      JGA=JG+NA
      G(1)=PQ(1,JGA)
      G(2)=PQ(2,JGA)
      IF (IT.EQ.0) THEN
         DO 150 I=1,2
            EDW(I)=RU
150      CONTINUE
      ELSE
C
C  COMPUTE DEBYE-WALLER FACTORS
C
         CALL DEBWAT(G,GP,E,VPI,AK2,AK3,TEMP,T0,DRX,DRY,DR0,EDW)
      ENDIF
      R=CZ
      T=CZ
C
C  START SUMMATION OVER (L,M)
C
      ST=RU
      SL=RU
      CM=RU
      JLM=1
      L1=1
      DO 240 L=1,LM
         SM=SL
         LL=L+L-1
         L1=L1+LL
         DO 230 M=1,LL
            MM=L1-M
            CY=CYLM(JG,MM)
            MM=L1-LL-1+M
            CTT=YLM(JLM)*ST
            CTR=CTT*SL
            CTT=CTT*SM
            R=R+CTR*CY
            T=T+CTT*CY
            JLM=JLM+1
            CM=CM*CI
            SM=-SM
230      CONTINUE
         CM=CONJG(CM)
         SL=-SL
         ST=ST*CI
240   CONTINUE
      RA=R*EDW(1)
      TA=T*EDW(2)
      IF (IGC.GT.0) THEN
         GP(2)=-GP(2)
         AK3=-AK3
      ENDIF
      RETURN
      END
C =========================================================================
C
C  Subroutine TRINT computes the reflected beam intensities from the
C  (complex) reflected amplitudes. Included are angular prefactors (energy
C dependant). TRINT prints out non-zero intensities only.
C
C Parameter List;
C ===============
C
C   N         =      NO. OF BEAMS AT CURRENT ENERGY.
C   WV        =      INPUT REFLECTED AMPLITUDES.
C   AT        =      OUTPUT REFLECTED INTENSITIES.
C   PQ        =      LIST OF RECIPROCAL-LATTICE VECTORS G (BEAMS).
C   PQF       =      SAME AS PQ, BUT IN UNITS OF THE RECIPROCAL-LATTICE 
C                    CONSTANTS.
C   VV        =      VACUUM LEVEL ABOVE SUBSTRATE MUFFIN-TIN CONSTANT.
C   THETA,FI  =      POLAR AND AZIMUTHAL ANGLES OF INCIDENT BEAM.
C   EEV       =      CURRENT ENERGY IN EV ABOVE VACUUM LEVEL.
C   A         =      STRUCTURAL PARAMETER OR OTHER IDENTIFIER TO BE PUNCHED 
C                    ON CARDS.
C
C In Common Blocks;
C =================
C
C   E         =      CURRENT ENERGY IN HARTREES ABOVE SUBSTRATE MUFFIN-TIN 
C                    CONSTANT.
C   VPI       =      IMAGINARY PART OF CURRENT ENERGY.
C   CK2,CK3   =      PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C   AS             =      NOT USED.
C
C Modified version of routine RINT from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE TRINT(N,WV,PQ,PQF,VV,THETA,FI,IPR)
C
      DIMENSION WV(N),PQ(2,N),PQF(2,N)
      COMPLEX WV
C
      COMMON /X4/E,VPI,CK2,CK3
C
10    FORMAT (' E=',1F10.4,' THETA=',1F10.4,' FI=',1F10.4,/)
20    FORMAT ('    PQ1    PQ2      REF INT ')
30    FORMAT (2F7.3,1E13.5,I7)
C
      AK=SQRT(AMAX1(2.0*E-2.0*VV,0.0))
      BK2=AK*SIN(THETA)*COS(FI)
      BK3=AK*SIN(THETA)*SIN(FI)
C
C  C IS K(PERP) IN VACUUM FOR INCIDENT BEAM
C
      C=AK*COS(THETA)
      C=AK*AK-BK2*BK2-BK3*BK3
      C=SQRT(C)
      TH1=THETA*180.0/3.14159265
      FI1=FI*180.0/3.14159265
      IF (IPR.GT.0) THEN
         WRITE (1,10) E,TH1,FI1
         WRITE (1,20)
      ENDIF
      DO 120 J=1,N
         AK2=PQ(1,J)+BK2
         AK3=PQ(2,J)+BK3
         A=2.0*E-2.0*VV-AK2*AK2-AK3*AK3
C         AT(J)=0.0
         AT=0.0
C
C  SKIP PRINT OUTPUT FOR NON-EMERGING BEAMS
C
         IF (A.GT.0) THEN
C
C  A IS K(PERP) IN VACUUM FOR SCATTERED BEAMS
C
            A=SQRT(A)
            WR=REAL(WV(J))
            WI=AIMAG(WV(J))
C
C  AT IS REFLECTED INTENSITY (FOR UNIT INCIDENT CURRENT)
C
C            AT(J)=(WR*WR+WI*WI)*A/C
            AT=(WR*WR+WI*WI)*A/C
            IF (IPR.GT.0) WRITE (1,30) PQF(1,J),PQF(2,J),AT
         ENDIF
120   CONTINUE
      RETURN
      END
C=======================================================================
C
C  Subroutine TSCATF interpolates tabulated phase shifts and produces
C  the atomic T-matrix elements (output in AF and TSF0). These are also
C  corrected for thermal vibrations (output in CAF and TSF). AF and CAF
C  are meant to be used bu subroutine MSMFT, TSF0 and TSF by MTINVT.
C  If NFLAGINT=1 the energy at which the phase shift is needed is outside
C  the range of the input tabulated phase shifts
C  Dimensions are such that LMAX=9, NEL max =3, max # of input energies=90
C Parameter List;
C ===============
C
C IEL         = CHEMICAL ELEMENT TO BE TREATED NOW, IDENTIFIED BY THE INPUT
C               SEQUENCE ORDER OF THE PHASE SHIFTS (IEL=1,2 OR 3).
C L1          = LMAX+1.
C ES          = LIST OF ENERGIES AT WHICH PHASE SHIFTS ARE TABULATED.
C PHSS        = TABULATED PHASE SHIFTS.
C NPSI        = NO. OF ENERGIES AT WHICH PHASE SHIFTS ARE GIVEN.
C IT(i)          = specifies whether temperature effects are included
C               for element i
C EB-V        = CURRENT ENERGY (V CAN BE USED TO DESCRIBE LOCAL VARIATIONS
C               OF THE MUFFIN-TIN CONSTANT).
C PPP         = CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CPPP.
C NN1         = NN2+NN3-1.
C NN2         = NO. OF OUTPUT TEMPERATURE-CORRECTED PHASE SHIFTS DESIRED.
C NN3         = NO. OF INPUT PHASE SHIFTS.
C DR0         = FOURTH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE.
C DRPER       = RMS VIBRATION AMPLITUDE PERPENDICULAR TO SURFACE.
C DRPAR       = RMS VIBRATION AMPLITUDE PARALLEL TO SURFACE.
C T0          = TEMPERATURE AT WHICH DRPER AND DRPAR HAVE BEEN COMPUTED.
C T           = CURRENT TEMPERATURE.
C PHSS2       = second derivative of the tabulated phase shifts(as a function
C               of the energy) needed for cubic spline interpolation
C
C Modfied version of routine TSCATF from the VAN HOVE/TONG package.
C Modifications by WANDER.
C Modifications of the interpolation scheme by BARBIERI
C
C=========================================================================
C
      SUBROUTINE TSCATF(IEL,L1,ES,PHSS,PHSS2,NPSI,IT,EB,V,PPP,
     & NN1,NN2,NN3,DR0,DRPER,DRPAR,T0,T,TSF0,TSF,AF,CAF,NFLAGINT)
C
      COMPLEX CI,DEL(16),CA,AF(L1),CAF(L1),TSF0(6,16),TSF(6,16)
      DIMENSION PHSS(NPSI,80),PHSS2(NPSI,80),ES(NPSI),PPP(NN1,NN2,NN3)
      DIMENSION PHS(16),PHSL(90),PHSL2(90)
      DIMENSION DR0(5),DRPER(5),DRPAR(5),IT(5)
C
C
700   FORMAT (' TOO LOW ENERGY FOR AVAILABLE PHASE SHIFTS')
C     
      NFLAGINT=0
      CI=(0.0,1.0)
      E=EB-V
      IF (E.LT.ES(1)) THEN
         WRITE (1,700)
         STOP
      ELSE
C
C  FIND SET OF PHASE SHIFTS APPROPRIATE TO DESIRED CHEMICAL ELEMENT
C  AND INTERPOLATE TO CURRENT ENERGY (OR EXTRAPOLATE TO ENERGIES
C  ABOVE THE RANGE GIVEN FOR THE PHASE SHIFTS)
C
         IO=(IEL-1)*L1
         I=1
720      IF ((E-ES(I))*(E-ES(I+1)).LE.0) GOTO 750
         I=I+1
         IF (I.LT.NPSI) GOTO 720
         NFLAGINT=1
         I=I-1
750      KLO=I
         KHI=I+1
         DO 760 L=1,L1
            DO 360 II=1,NPSI
               PHSL(II)=PHSS(II,L+IO)
               PHSL2(II)=PHSS2(II,L+IO)
360         CONTINUE
            CALL SPLINT(ES,PHSL,PHSL2,NPSI,E,PHSOUT,PHP
     &      ,KLO,KHI)
            PHS(L)=PHSOUT
760      CONTINUE
C
C
C  PRODUCE AND STORE TEMPERATURE-INDEPENDENT T-MATRIX ELEMENTS
C  ONLY IF IT(IEL)=0
C
         IF(IT(IEL).EQ.0) THEN
         DO 790 L=1,L1
            A=PHS(L)
            AF(L)=A*CI
            AF(L)=CEXP(AF(L))
            A=SIN(A)
            AF(L)=A*AF(L)
            TSF0(IEL,L)=AF(L)
            TSF(IEL,L)=AF(L)
790      CONTINUE
         ENDIF
C
C  AVERAGE ANY ANISOTROPY OF RMS VIBRATION AMPLITUDES
C
         DR=SQRT((DRPER(IEL)*DRPER(IEL)+2.0*DRPAR(IEL)*DRPAR(IEL))/3.0)
C
C  COMPUTE TEMPERATURE-DEPENDENT PHASE SHIFTS (DEL)
C
         DR01=DR0(IEL)
C
C  PRODUCE AND STORE TEMPERATURE-DEPENDENT T-MATRIX ELEMENTS
C  ONLY IF IT=1
C
         IF(IT(IEL).EQ.1) THEN
         CALL PSTEMP(PPP,NN1,NN2,NN3,DR01,DR,T0,T,E,PHS,DEL)
         DO 840 L=1,L1
            CA=DEL(L)
            CAF(L)=CA*CI
            CAF(L)=CEXP(CAF(L))
            CAF(L)=-CI*(CAF(L)*CAF(L)-1.0)/2.0
            TSF(IEL,L)=CAF(L)
C            c1=.5
C            c2=1.-c1
C            IF(IEL.EQ.2) TSF(1,L)=c1*TSF(1,L)+c2*TSF(2,L)
840      CONTINUE
         ENDIF
         RETURN
      ENDIF
      END
C======================================================================
C
C WAVE2 sets the parallel wavevector of the incident and time reversed
C exit beams for each energy.
C                                                                      
C Modified version of subroutine WAVE2 (written by ROUS). Modifications by
C WANDER.
C
C=========================================================================
C
      SUBROUTINE WAVE2(AK2,AK3,THETA,FI,E,VV,AK21,AK31,AK2M,AK3M,NT0,
     & RAR1,RAR2,PQFEX,PSQ,NEXIT,SPQF,NT,NBIN)
C
      DIMENSION AK2M(NT0),AK3M(NT0),PQFEX(2,NT0),SPQF(2,NT)
      DIMENSION PSQ(2,NT0),RAR1(2),RAR2(2)
C
33    FORMAT (' WAVE2: ERROR. EXIT BEAM NO',I4,
     & ' IS NOT IN INPUT    LIST')
C
      IF (NEXIT.EQ.0) THEN
         AK=SQRT(2.0*(E-VV))*SIN(THETA)
         AK2=AK*COS(FI)
         AK3=AK*SIN(FI)
         NBIN=1
         DO 100 IG=1,NT0
            DO 110 I=1,2
               PSQ(I,IG)=PQFEX(1,IG)*RAR1(I)+PQFEX(2,IG)*RAR2(I)
110         CONTINUE
            AK2M(IG)=-AK2-PSQ(1,IG)
            AK3M(IG)=-AK3-PSQ(2,IG)
100      CONTINUE
      ELSE
C
C NEXIT>0 SO LOCATE THIS BEAM IN MAIN BEAM LIST.
C
C FIRST CHECK FOR EMERGENCE OF THIS BEAM
C
         A=2.0*(E-VV)-AK2M(NEXIT)*AK2M(NEXIT)-AK3M(NEXIT)*AK3M(NEXIT)
         IF (A.GT.0) THEN
C
C BEAM EMERGENCES. WHICH BEAM DOES THIS CORRESPOND TO
C
            DO 120 I=1,NT
               DIFF1=ABS(PQFEX(1,NEXIT)+SPQF(1,I))
               DIFF2=ABS(PQFEX(2,NEXIT)+SPQF(2,I))
               IF ((DIFF1.LT.1.0E-02).AND.(DIFF2.LT.1.0E-02)) GOTO 131
120         CONTINUE
C
C CAN'T FIND THIS BEAM
C
            WRITE (1,33) NEXIT
            STOP
C
C SET BEAM LABEL FOR THIS BEAM
C
131         NBIN=I
         ELSE
C
C SET NBIN TO NEGATIVE NUMBER TO INDICATE NON-EMERGENCE
C
            NBIN=-1
         ENDIF
      ENDIF
      AK21=AK2
      AK31=AK3
      RETURN
      END
C =========================================================================
C
C  Subroutine XMA produces the intra-layer multiple-scattering matrix X
C  for a single Bravais-Lattice layer. The (L,M) sequence is the
C  \SYMMETRIZED\ one and X is split into two parts corresponding to
C  L+M= Even (XEV) and L+M= Odd (XOD), as a result of block-diagonalization.
C
C Parameter List;
C ===============
C
C FLM       = LATTICE SUM FROM SUBROUTINES FMAT AND MSMF.
C XEV,XOD   = OUTPUT MATRIX X IN TWO PARTS.
C LEV       = (LMAX+1)*(LMAX+2)/2.
C LOD       = LMAX*(LMAX+1)/2.
C AF           = NOT USED.
C CAF       = TEMPERATURE-DEPENDENT ATOMIC SCATTERING T-MATRIX ELEMENTS.
C LM        = LMAX+1.
C LX,LXI    = PERMUTATIONS OF (L,M) SEQUENCE, FROM SUBROUTINE LXGENT.
C LMMAX     = (LMAX+1)**2.
C KLM       = (2*LMAX+1)*(2*LMAX+2)/2.
C CLM       = CLEBSCH-GORDON COEFFICIENTS, FROM SUBROUTINE CELMG.
C NLM       = DIMENSION OF CLM (SEE TLEED1).
C
C Modified version of routine XM from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE XM(FLM,XEV,XOD,LEV,LOD,CAF,LM,LXI,LMMAX,KLM,CLM,
     & NLM)
C
      INTEGER LXI
      COMPLEX XEV,XOD,FLM,CZERO,ACC,CAF
      DIMENSION CLM(NLM),XEV(LEV,LEV),XOD(LOD,LOD),FLM(KLM)
      DIMENSION CAF(LM),LXI(LMMAX)
C
      LMAX=LM-1
      CZERO=CMPLX(0.0,0.0)
      DO 611 I=1,LEV
         DO 440 J=1,LEV
            XEV(I,J)=CZERO
440      CONTINUE
611   CONTINUE
      DO 612 I=1,LOD
         DO 450 J=1,LOD
            XOD(I,J)=CZERO
450      CONTINUE
612   CONTINUE
      L2MAX=LMAX+LMAX
      NODD=LOD
      JSET=1
      MM=LEV
      N=1
C
C  FIRST XOD IS CREATED
C
460   J=1
      L=JSET
470   M=-L+JSET
      JL=L+1
480   K=1
      LPP=JSET
490   MPP=-LPP+JSET
500   MPA=IABS(MPP-M)
      LPA=IABS(LPP-L)
      IF (LPA.GT.MPA) MPA=LPA
      MP1=MPP-M+L2MAX+1
      LP1=L+LPP+1
      ACC=CZERO
530   JLM=(LP1*LP1+MP1-L2MAX)/2
      ACC=ACC+CLM(N)*FLM(JLM)
      N=N+1
      LP1=LP1-2
      IF (LP1-1-MPA.GE.0) GOTO 530
      JX=LXI(J+MM)
      KX=LXI(K+MM)
      XEV(JX,KX)=ACC*CAF(JL)
      K=K+1
      MPP=MPP+2
      IF (LPP.GE.MPP) GOTO 500
      LPP=LPP+1
      IF (LMAX.GE.LPP) GOTO 490
      J=J+1
      M=M+2
      IF (L.GE.M) GOTO 480
      L=L+1
      IF (LMAX.GE.L) GOTO 470
      IF (JSET.GT.0) THEN
         DO 613 J=1,NODD
            DO 600 K=1,NODD
               XOD(J,K)=XEV(J,K)
600         CONTINUE
613      CONTINUE
         JSET=0
         MM=0
C
C  NOW RETURN TO CREATE XEV
C
         GOTO 460
      ENDIF
      RETURN
      END
C =========================================================================
C
C  Subroutine XMT produces the intra-layer multiple scattering matrix
C  1-X for subroutine TAUMAT, output in the TONG convention using
C  input in the PENDRY convention. XMT must be called twice, first to
C  produce the values for L+M= odd, then for L+M= even (in that order).
C
C Parameter List;
C ===============
C
C IL        =  1 FOR L+M= ODD.
C           =  2 FOR L+M= EVEN.
C X         =  OUTPUT MATRIX 1-X, BLOCK-DIAGONALIZED, 2ND BLOCK LEFT-ADJUSTED.
C LL        =  INPUT EITHER LOD OR LEV.
C TSF       =  ATOMIC T-MATRIX ELEMENTS.
C NTAU      =  NO. OF CHEMICAL ELEMENTS CONSIDERED.
C IT        =  INDEX OF CURRENT CHEMICAL ELEMENT (.LE.NTAU).
C N         =  RUNNING INDEX OF CLM  MAY NOT BE RESET BETWEEN THE TWO CALLS TO
C              XMT IN TAUMAT.
C
C Modified version of routine XMT from the VAN HOVE/TONG LEED package.
C Modifications by WANDER.
C
C =========================================================================
C
      SUBROUTINE XMT(IL,FLM,X,LEV,LL,TSF,IT,LM,LXI,LMMAX,KLM,
     & CLM,NLM,N)
C
      COMPLEX TSF(6,16)
      complex acc
      COMPLEX X,FLM,CZERO,CI,ST,SU,RU
      DIMENSION CLM(NLM),X(LEV,LL),FLM(KLM),LXI(LMMAX)
C
      LMAX=LM-1
      RU=(1.0,0.0)
      CZERO=CMPLX(0.0,0.0)
      CI=CMPLX(0.0,-1.0)
      L2MAX=LMAX+LMAX
C
C  IF IL=1, CONSIDER L+M= ODD ONLY
C  IF IL=2, CONSIDER L+M= EVEN ONLY
C
      IF (IL.LT.2) THEN
         JSET=1
         MM=LEV
         N=1
      ELSE
         JSET=0
         MM=0
      ENDIF
      J=1
      L=JSET
470   M=-L+JSET
      MEXP=MOD(L,4)
      ST=CI**MEXP
480   K=1
      LPP=JSET
490   MPP=-LPP+JSET
      JLPP=LPP+1
      MEXP=MOD(LPP,4)
      SU=CI**MEXP/ST
500   MPA=IABS(MPP-M)
      LPA=IABS(LPP-L)
      IF (LPA.GT.MPA) MPA=LPA
      MP1=MPP-M+L2MAX+1
      LP1=L+LPP+1
      ACC=CZERO
530   JLM=(LP1*LP1+MP1-L2MAX)/2
      ACC=ACC+CLM(N)*FLM(JLM)
      N=N+1
      LP1=LP1-2
      IF (LP1-1-MPA.GE.0) GOTO 530
      JX=LXI(J+MM)
      KX=LXI(K+MM)
      X(KX,JX)=-ACC*TSF(IT,JLPP)*SU
      IF (J.EQ.K) X(KX,JX)=X(KX,JX)+RU
      K=K+1
      MPP=MPP+2
      IF (LPP.GE.MPP) GOTO 500
      LPP=LPP+1
      IF (LMAX.GE.LPP) GOTO 490
      J=J+1
      M=M+2
      IF (L.GE.M) GOTO 480
      L=L+1
      IF (LMAX.GE.L) GOTO 470
      RETURN
      END
C =========================================================================
C
C Subroutine AVEINT computes the average intensity of the peaks in each
C beam and returns it as EEAVE 
C
C =========================================================================
C
      SUBROUTINE AVEINT(A,AP,NBD,NE1,NE2,NB,IEERG,EEAVE,NFL)
C
      DIMENSION A(NBD,IEERG)
      DIMENSION NE1(NBD),NE2(NBD),AP(NBD,IEERG),EEAVE(30)
C    
C we go through the derivative of each beam and record the intensity
C when it becomes zero and it corresponds to a peak
C
       EEAVE(NB+1)=.0
       DO 10 I=1,NB
         ZMAX=0.
         EAV=0.
         ICOUNT=0
         J1=NE1(I)
         IF(NFL.EQ.0) J1=1
C
C if exp. and theoret. beams do not have energy range in common skip
C the calculation (of course this should not be the norm!)
C
         IF(NFL.EQ.1.AND.NE1(I).GE.NE2(I)) GOTO 10
         DO 100 J=J1,NE2(I)-1
           IF(AP(I,J)*AP(I,J+1).LT.0..AND.AP(I,J).GT.0.) THEN
c
c try to eliminate false peaks in the noise region
c
             IF(A(I,J).GT.ZMAX)ZMAX=A(I,J)
             IF (A(I,J).GT.ZMAX/50.) THEN
               EAV=EAV+A(I,J)
               ICOUNT=ICOUNT+1
             ENDIF
           ENDIF
100      CONTINUE            
         FL=FLOAT(ICOUNT)
         IF (NFL.EQ.1.AND.FL.EQ.0) THEN 
             EEAVE(I)=.0001
         ELSE
             EEAVE(I)=EAV/FL
         ENDIF
         EEAVE(NB+1)=EEAVE(NB+1)+EEAVE(I)
10     CONTINUE
       EEAVE(NB+1)=EEAVE(NB+1)/FLOAT(NB)
       RETURN
       END
C=========================================================================
C
C Subroutine BESSEL generates the spherical Bessel functions for a complex
C argument Z, L=0 to LMAX using Miller's device.
C
C See Abramowitz and Stegun P.452.
C
C Input Parameters;
C =================
C
C Z          =    Complex argument of the Bessel function
C LMAX       =    Maximum vlue of the angular momentum
C LMAX1      =    LMAX+1 
C BJ         =    Output Bessel function of order L=0 to LMAX
C
C
C Note
C ====
C
C When this routine is run on a micro with limited accuracy, care is needed
C to ensure that numerical overflow does not occur. This is done by setting
C LSTART (the starting L value for backwards recurrence). LSTART must be 
C sufficiently large to ensure that the Bessel functions are accurate, but
C small enough so that overflow does not occur. 
C
C Author: ROUS. Modifications by WANDER.
C
C=====================================================================
C
      SUBROUTINE BESSEL(BJ,Z,LMAX,LMAX1)
C
      PARAMETER (LSTART=11)
c      COMPLEX BJ(LMAX1),Z,BJ0,BJS(20),CI
      COMPLEX BJ(LMAX1),BJ0,BJS(20),CI
      COMPLEX Z
C
C SET CONSTANTS
C
      CI=CMPLX(0.0,1.0)
C
C ASSUME BESSEL FN FOR L=LSTART+1 IS ZERO & FOR L=LSTART =1.0
C
      BJS(LSTART)=CMPLX(0.0,0.0)
      BJS(LSTART-1)=CMPLX(1.0,0.0)
C
C GENERATE BESSEL FUNCTIONS FOR L=0,LSTART-2 BY BACKWARD
C RECURRENCE.
C
      DO 100 IL=0,LSTART-3
         L=LSTART-3-IL
         BJS(L+1)=FLOAT(2*(L+1)+1)*BJS(L+2)/Z-BJS(L+3)
100   CONTINUE
C
C EVALUATE BESSEL FN FOR L=0 EXPLICITLY
C
      BJ0=CSIN(Z)/Z
C
C NORMALISE PREVIOUSLY CALCULATED BESSEL FNS
C
      DO 110 L=0,LMAX
         BJ(L+1)=BJS(L+1)*BJ0/BJS(1)
110   CONTINUE
      RETURN
      END
C=====================================================================
C
C  SUBROUTINE BINSRX FINDS A REQUIRED INTERPOLATION INTERVAL
C  BY BINARY SEARCH (SUCCESSIVE HALVING OF INITIAL INTERVAL)
C
C=====================================================================
C
      SUBROUTINE BINSRX(IL, IH, X, WORX, LENGTH)
C
      DIMENSION WORX(LENGTH)
C
      I = (LENGTH + 1)/2
      IHI =  (I + 1)/2
      IF(X.LT. WORX(1) .OR. X .GE. WORX(LENGTH)) GO TO 100
20    CONTINUE
      IF(I .LT. 1) GO TO 40
      IF(I .GE. LENGTH) GO TO 30
      IF(X .GE.WORX(I) .AND. X .LE. WORX(I + 1)) GO TO 50
      IF(X .LT. WORX(I)) GO TO 30
40    I = IHI + I
      IHI = (IHI + 1) / 2
      GO TO 20
30    I = I - IHI
      IHI =(IHI + 1) / 2
      GO TO 20
50    IL = I
      IH = I + 1
      GO TO 110
100   IF(X .LT. WORX(1)) GO TO 105
      IL = LENGTH - 1
      IH = LENGTH
      GO TO 110
105   IL = 1
      IH = 2
110   RETURN
      END
C ============================================================================
C
C FUNCTION BRENT takes the triplet AX,BX,CX which brackets the minimum
C of the R factor computed from the point PCOM along the direction XICOM
C and evaluate this minimum. This is returned as XMIN, with the value of the
C R factor being returned as BRENT.TOL is the fractional accuracy as passed
C from LINMIN
C
C References 'NUMERICAL RECIPES' W.H.Press,B.P.Flannery, et al. Cambridge
C             University Press
C 
C Modifications by BARBIERI
C
C Input Parameters;
C =================
C
C NLAY                   = NUMBER OF LAYERS IN COMPOSITE LAYER
C NDIM                   = DIMENSIONALITY OF SEARCH (1=X AXIS ONLY
C                                                    3=X,Y,Z AXES)
C PCOM                      = ARRAY CONTAINING THE INITIAL STARTING POINT
C                          IN THE SEARCH
C XICOM                     = ARRAY STORING THE INITIAL DIRECTION SET
C DISP                   = COORDINATES INPUT BY USER
C NNDIM                  = TOTAL NUMBER OF DIMENSIONS IN SEARCH (=NLAY*NDIM)
C ADISP                  = GEOMETRY OF CURRENT POINT IN SAME FORMAT AS DISP
C                          (USED AS INPUT TO FUNCV)
C DVOPT                  = SHIFT IN THE INNER POTENTIAL FOR THE STARTING 
C                          CONFIGURATION (in COMMON /RPL )
C LLFLAG                  = INDICATES WHETHER THE LAYER (OR NON STRUCTURAL
C                          PARAMETER) COORDINATES HAVE TO BE VARIED 
C                          IN THE SEARCH
C VOPT                   = INNER POTENTIAL OF CURRENT INPUT (used as input
C                          in FUNCV), in COMMON /NSTR
C ITMAX                  = MAXIMUM NUMBER OF ITERATIONS TO BE PERFORMED BY
C                          POWELL
C ISTART                 = 0 USE PARAMETERS FROM TLEED4.I
C                          1 RESTART SEARCH FROM COORDINATES IN RESTART.D
C IPR                    = PRINT CONTROL PARAMETER
C ILOOK                  = LOOKUP TABLE FOR DOMAIN AVERAGING
C ACOORD                 = STORAGE FOR SYMMETRY EQUIVALENT COORDINATES
C MICUT,MJCUT,PSQ,JYLM,  = DUMMY ARRAYS TO BE PASSED TO ROUTINE FUNCV
C BJ,YLM,QS,AT,XISTS,
C XIST,INBED,IEERG,AE,EE,
C NEE,NBEA,BENAME,IPR,XPL,
C YPL,NNN,AP,APP,YE,TSE,
C TSE2,TSEP,TSEP2,TSEPP,
C TSEPP2,TSEY2,WR,WB,IBP,
C ETH
C in Common:
C IFUNC                  =Total number of function evaluations
C
C P and XI are returned as the best point found, and the latest directiob set
C respectively
C
C
C =======================================================================
C
      FUNCTION BRENT(ax,bx,cx,fax,fbx,fcx,tol,xmin,pcom,xicom,NLAY,
     & NDIM,DISP,NNDIM,ADISP,ILOOK,ACOORD,MICUT,
     & MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     & NDIML,DIREC,NLTIN,LPOINT)
C
      PARAMETER (ibmax=100,cgold=.3819660,zeps=1.0e-3)
      PARAMETER (MAXC=100,accfun=3.e-3)
C
C ACCFUN corresponds to the square root of the accuracy with
C which the function to be minimized can be computed. A limited accuracy 
C can be due to machine precision or numerical noise inherent in the
C way the function is computed 
C
      DIMENSION PCOM(MAXC),XICOM(MAXC),LLFLAG(NLAY+NNST)
      DIMENSION DIREC(NLAY,2,2),NDIML(NLAY),LSFLAG(NLAY)
      DIMENSION DISP(NLAY,3),ADISP(NLAY,3),XT(MAXC)
      DIMENSION ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1)
      COMPLEX JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /RPL/DVOPT
      COMMON /NSTR/VOPT,NNST,NNSTEF
      a=ax
      b=cx
      fa=fax
      fb=fcx
      if (ax.gt.cx) then
         a=cx
         fa=fcx
         b=ax
         fb=fax
      endif
      w=a
      v=b
      if (fa.gt.fb) then
         w=b
         v=a
      endif
      x=bx
      fx=fbx
      fv=max(fa,fb)
      fw=min(fa,fb)
      xm=0.5*(a+b)
      if(x.ge.xm) then
         e=a-x
      else
         e=b-x
      endif
      d=cgold*e
      do 11 iter=1,ibmax
        xm=0.5*(a+b)
C        tol1=tol+zeps
        tolf=accfun*abs(x)
        tolf2=2.*tolf
        tol1=tol*abs(x)+zeps
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tolf) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.
     *        p.ge.q*(b-x)) goto 1
          d=p/q
          u=x+d
C          if(u-a.lt.tolf2 .or. b-u.lt.tolf2) d=sign(tolf,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=cgold*e
2       if(abs(d).ge.tolf) then
          u=x+d
        else
          u=x+sign(tolf,d)
        endif
        do 4 j=1,NNDIM
          xt(j)=pcom(j)+u*xicom(j)
4       continue
        CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
        IFUNC=IFUNC+1
        fu=fval
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
cjcm      pause 'brent exceed maximum iterations.'
3     xmin=x
      brent=fx
      return
      end
C ==========================================================================
C
C Subroutine COMNEI finds the common energy interval between theory and
C experiment.
C
C From the VAN HOVE R-factor program
C
C ==========================================================================
C
      SUBROUTINE COMNEI(EE,NBED,NEE,ET,NBTD,NET,IBE,IBT,V0,EINCR,NE1,
     & NE2,NT1,NT2,EET,IEERG)
C
      DIMENSION EE(NBED,IEERG),NEE(NBED),ET(NBTD,IEERG),NET(NBTD)
C
      NE=NEE(IBE)
      NT=NET(IBT)
      DE1=ET(IBT,1)+V0-EE(IBE,1)
      DE2=ET(IBT,NT)+V0-EE(IBE,NE)
      IF (DE1.LT.0.) THEN
         NE1=1
         NT1=INT((-DE1/EINCR)+0.0001)+1
      ELSE
         NE1=INT((DE1/EINCR)+0.0001)+1
         NT1=1
      ENDIF
      IF (DE2.LT.0.) THEN
         NE2=NE-INT((-DE2/EINCR)+0.0001)
         NT2=NT
      ELSE
         NE2=NE
         NT2=NT-INT((DE2/EINCR)+0.0001)
      ENDIF
      IF (NE2.GT.NE1) EET=EE(IBE,NE2)-EE(IBE,NE1)
      RETURN
      END
C =========================================================================
C
C Subroutine DELXGEN generates the change in plane wave amplitude produced
C by the displacements held in ACOORD.
C
C Input Parameters;
C =================
C
C QS              =  Tensor
C ACOORD          =  Coordinates of each symmetrically equivalent domain
C MICUT,MJCUT     =  Coordinates of tensor element before truncation. Used
C                    to reconstruct tensor.
C NLAY            =  Number of layers in calculation
C NDOM            =  Number of symmetry equivalent domains
C YLM,BJ,JYLM     =  Work Space
C DELXI           =  Changes in plane wave amplitude for each beam
C E,VPIS          =  Current (complex) energy
C IOFF            =  Current offset in truncated tensor
C NEXIT           =  Index of current exit beam
C 
C In Common Blocks;
C =================
C
C Block: TLVAL
C ------------
C
C LSMAX        =    Maximum size of L used in single centre expansion
C LSMMAX       =    (LSMAX+1)**2
C ICUT         =    Size of Q for each energy point/layer/beam
C LSMAX1       =    LSMAX+1
C NT0          =    Number of exit beams
C IQSIZ        =    Total size of Q. (=ICUT*NLAY*NT0*NERG)
C
C Modified version of routine MATNEW by ROUS and WANDER. Modifications
C by WANDER,BARBIERI.
C
C =========================================================================
C
      SUBROUTINE DELX2(QS,ACOORD,MICUT,MJCUT,NLAY,YLM,BJ,
     & JYLM,E,VPIS,IOFF,J,NLTIN,LPOINT,XIST,PRE2,
     & NERG,IE)
C
      DIMENSION ACOORD(12,NLAY,3)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),C(3)
      DIMENSION LPOINT(NLTIN)
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2),RBR1(2),RBR2(2)
cc      COMPLEX QS(IQSIZ),YLM(LSMMAX),BJ(LSMAX1),JYLM(LSMMAX)
      COMPLEX QS(IQSIZ),YLM(LSMMAX),JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX DELXI,CAPPA,Z,CI,ST,PRE
      COMPLEX XIST(NT0,NERG),AMAT,TMINUS
      DIMENSION PRE2(50)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
C Set constants
C
      NL=NL1*NL2
      PI=4.0*ATAN(1.0)
      PI4=16.0*ATAN(1.0)
      CI=CMPLX(0.0,1.0)
      DO 1002 JJ=1,NT0
           XIST(JJ,IE)=(0.,0.)
1002  CONTINUE
C
C Begin loop over subplanes of composite layer
C
         DO 100 IJK2=1,NLTIN
            IJK=LPOINT(IJK2)
C
C Pull out current coordinate into temporary storage, and generate
C total vector displacement.
C
           CL=0.0
           DO 120 K=1,3
              C(K)=ACOORD(J,IJK,K)/0.529
               CL=CL+C(K)*C(K)
120         CONTINUE 
C
C Set prefactors. We require bessel functions for Z=CAPPA*MOD(C).
C
            CAPPA=CMPLX(2.0*E,-2.0*VPIS+0.000001)
            Z=CSQRT(CAPPA)*SQRT(CL)
C
C Is this atom moving?
C
            IF (CL.LT.0.0000001) THEN
C
C If not, clear the functions BJ and YLM. Set only the first element of
C each array.
C
              DO 130 L=0,LSMAX
                 BJ(L+1)=CMPLX(0.0,0.0)
130            CONTINUE
              BJ(1)=CMPLX(1.0,0.0)
               DO 140 LM=1,LSMMAX
                YLM(LM)=CMPLX(0.0,0.0)
140            CONTINUE
              YLM(1)=1.0/SQRT(4.0*PI)
           ELSE
C
C The atom is moving, so set the bessel function (BJ) and spherical harmonic
C (YLM)
C
               CALL BESSEL(BJ,Z,LSMAX,LSMAX1)
              CALL HARMON(YLM,C,LSMAX,LSMMAX)
            ENDIF
C
C Generate the structure factor JYLM(l,m)=BJ(l+1)*CI**(-l)*YLM(l,-m)
C The JYLM are independent on the beams
C
C
           I=0
            ST=CI
           DO 150 L=0,LSMAX
              ST=ST/CI
              PRE=BJ(L+1)*ST
               DO 160 M=-L,L
                 I=I+1
                 IM=I-2*M
                 JYLM(I)=PRE*YLM(IM)
160            CONTINUE
150         CONTINUE
C
C Generate the matrix element by multipling Q on both the right and left
C by JYLM
C
C                 AMAT=CMPLX(0.0,0.0)
C          IIOF=IOFF
          DO 110 JJ=1,NT0
             AMAT=CMPLX(0.0,0.0)
             IF(PRE2(JJ).NE.0) THEN
                TMINUS=QS(IOFF+1)/(PI4)
                DO 170 IC=1,ICUT
                   I2=MICUT(IC)
                   I3=MJCUT(IC)
C                   AMAT(JJ)=AMAT(JJ)+JYLM(I2)*QS(IOFF+IC)*JYLM(I3)
                   AMAT=AMAT+JYLM(I2)*QS(IOFF+IC)*JYLM(I3)
170             CONTINUE
                IOFF=IOFF+ICUT
C           XIST(JJ,IE)=XIST(JJ,IE)+(AMAT(JJ)-TMINUS(JJ))/FLOAT(NL)
           XIST(JJ,IE)=XIST(JJ,IE)+(AMAT-TMINUS)/FLOAT(NL)
             ENDIF
110       CONTINUE
100      CONTINUE
      RETURN
      END
C ========================================================================
C
C Subroutine EPSZJ calculates the Zanazzi-Jona EPS (= MAX[ABS{derivative}])
C
C From the VAN HOVE R-factor program.
C
C ========================================================================
C
      SUBROUTINE EPSZJ(AEP,NBED,IB,NE1,NE2,EPS,IEERG)
C
      DIMENSION AEP(NBED,IEERG)
C
      EPS=0.0
      DO 10 IE=NE1,NE2
         A=ABS(AEP(IB,IE))
         IF (A.GT.EPS) EPS=A
10    CONTINUE
      RETURN
      END
C ============================================================================
C
C Subroutine EXPAN reads in the experimental data, and performs a preliminary
C analysis of this data.
C
C Note:
C EEAVE  = average energy of each beam after averaging. The dimension
C          20 is the maximum number of experimental beams after averaging
C
C
C Author: Wander. Based on Van Hove R-Factor program.
C         modifications: Barbieri
C
C ============================================================================
C
      SUBROUTINE EXPAN(INBED,IEERG,AE,EE,NEE,NBEA,BENAME,IPR,XPL,YPL,
     & NNN,AEP,AEPP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,VI)
C
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION AEP(INBED,IEERG),AEPP(INBED,IEERG),YE(INBED,IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10)
      DIMENSION TSEPP2(INBED),TSEY2(INBED)
C
C      CHARACTER*4 TITLE(20)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA,FI
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /REXP/EEINCR
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
      COMMON /TIT/TITLE(20)
C
100   FORMAT(20A4)
110   FORMAT(20I3)
120   FORMAT(/' ERROR IN READE:MORE EXPERIMENTAL BEAMS THAN HAVE BEEN 
     & DIMENSIONED')
C
C  READ AND PRINT DESCRIPTION OF EXPERIMENT
C
      READ (11,100)(TITLE(I),I=1,20)
      WRITE (1,100)(TITLE(I),I=1,20)
      READ (11,110)NBED
      IF(NBED.GT.INBED)THEN
         WRITE(1,120)
         STOP
      ENDIF
C
C READ IN EXPERIMENTAL DATA
C
      CALL READE(AE,EE,NBED,NEE,NBEA,BENAME,IPR,IEERG)
cjcm      write(*,*) 'READE:', NEE
c      write(*,*) EE
C
C  IF AN IRREGULAR INPUT ENERGY GRID WAS USED, FIRST INTERPOLATE ALL
C  EXP. DATA TO THE A GRID specified by the input EEINCR
C
      IF (IRGEXP.EQ.1) THEN
cjcm         write(*,*) 'Calling INTPOL: IRGEXP .eq. 1'
cjcm         write(*,*) 'INTPOL:', NEE
         CALL INTPOL(AE,AEP,AEPP,NBED,NEE,NBED,EE,
     %        EEINCR,IPR,XPL,YPL,IEERG)
      ENDIF
C
C  AVERAGE DATA FROM DIFFERENT EXPERIMENTS AND ORDER BY INCREASING ENERGY
C
      CALL EXPAV(AE,EE,NBED,NEE,BENAME,NBEA,NBE,IPR,XPL,NNN,IEERG)
cjcm      write(*,*) 'EXPAV:', NEE
C
C  ESTIMATE THE AVERAGE ENERGY OF EACH EXPERIMENTAL BEAM
C 
C      EGRID=EEINCR
C      IF (IRGEXP.EQ.1) EGRID=EINCR
C      DO 11 I=1,NBE
C          CALL VARSUM(AE,AE,AE,AE,NBED,1,I,1,1,NEE(I),0,
C     &       EGRID,0.,0.,1,EEAVE(I),YPL,IEERG)
C          DE=EE(I,NEE(I))-EE(I,1)
C          EEAVE(I)=EEAVE(I)/DE
C
C NBEA WAS RETURNED eEQUAL TO ZERO BY EXPAV. NOW WE RESET IT SO THAT
C IT NUMBERS THE NON EQUIVALENT BEAMS AFTER AVERAGING
C
C          NBEA(I)=I
C11    CONTINUE
C
C  SMOOTH EXP. DATA THE NUMBER OF TIMES DESIRED
C
         IF (ISMOTH.NE.0) THEN
            DO 10 I=1,ISMOTH
               CALL SMOOTH(AE,EE,NBED,NBE,NEE,IPR,IEERG)
10          CONTINUE
cjcm            write(*,*) 'SMOOTH:', NEE
         ENDIF
C
C  INTERPOLATE EXP. DATA TO WORKING GRID (MULTIPLES OF EINCR EV), UNLESS
C  DONE BEFORE
C
         IF (IRGEXP.EQ.0) THEN
cjcm            write(*,*) 'Calling INTPOL: IRGEXP .eq. 0'
cjcm            write(*,*) 'INTPOL:', NEE
            CALL INTPOL(AE,AEP,AEPP,NBED,NEE,NBE,EE,EINCR,
     &           IPR,XPL,YPL,IEERG)
         ENDIF
C
C  RENORMALIZE EXP. BY DIVIDING BY FITTED A*EXP(-ALPH*(E+10))
C
         IF (IREN.EQ.1) THEN
            CALL RENORM(AE,NBED,NBE,NEE,EE,IPR,YPL,IEERG)
         ENDIF
C
C  ESTIMATE THE AVERAGE INTENSITY OF THE PEAKS IN EACH EXPERIMENTAL BEAM
C 
         CALL AVEINT(AE,AEP,NBED,NEE,NEE,NBE,IEERG,EEAVE,0)
C
C  PRODUCE (1ST AND) 2ND DERIVATIVES OF EXP. DATA
C
C         CALL DERL(AE,NEE,NBED,NBE,AEP,EINCR,IEERG)
C         CALL DERL(AEP,NEE,NBED,NBE,AEPP,EINCR,IEERG)
C
C  PRODUCE PENDRY Y FUNCTION FOR EXP. DATA
C
         CALL YPEND(AE,AEP,NBED,NBE,NEE,EE,YE,
     %    VI,IPR,IEERG,EEAVE)
C
C  PRODUCE SOME INTEGRALS OVER EXP. DATA
C
C      DO 1112 J=1,NBE
C      DO 1111 I=1,NEE(J)
C         WRITE (*,*) EE(J,I),YE(J,I)
C         WRITE (*,*) EE(J,I),AE(J,I)
C1111  CONTINUE
C      WRITE (*,*)
C1112  CONTINUE
         DO 20 IB=1,NBE
            IE2=NEE(IB)
            IF (IE2.NE.0) THEN
               TSE(IB)=0.
               TSE2(IB)=0.
               TSEP(IB)=0.
               TSEP2(IB)=0.
               TSEPP(IB)=0.
               TSEPP2(IB)=0.
               TSEY2(IB)=0.
               CALL VARSUM(AE,AE,AE,AE,NBED,1,IB,1,1,IE2,0,
     &          EINCR,0.,0.,1,TSE(IB),YPL,IEERG)
               IF (WR(3).GE.1.E-6) CALL VARSUM(AE,AE,AE,AE,NBED,1,
     &          IB,1,1,IE2,0,EINCR,0.,0.,2,TSE2(IB),YPL,IEERG)
               IF (WR(4).GE.1.E-6) CALL VARSUM(AEP,AE,AE,AE,NBED,
     &          1,IB,1,1,IE2,0,EINCR,0.,0.,3,TSEP(IB),YPL,IEERG)
               IF (WR(5).GE.1.E-6) CALL VARSUM(AEP,AE,AE,AE,NBED,
     &          1,IB,1,1,IE2,0,EINCR,0.,0.,2,TSEP2(IB),YPL,IEERG)
               IF (.NOT.(WR(6).LT.1.E-6.AND.WR(9).LT.1.E-6)) CALL 
     &          VARSUM(AEPP,AE,AE,AE,NBED,1,IB,1,1,IE2,0,
     &          EINCR,0.,0.,3,TSEPP(IB),YPL,IEERG)
               IF (WR(7).GE.1.E-6) CALL VARSUM(AEPP,AE,AE,AE,NBED,
     &          1,IB,1,1,IE2,0,EINCR,0.,0.,2,TSEPP2(IB),YPL,IEERG)
               IF (WR(10).GE.1.E-6) CALL VARSUM(YE,AE,AE,AE,NBED,
     &          1,IB,1,1,IE2,0,EINCR,0.,0.,2,TSEY2(IB),YPL,IEERG)
            ENDIF
20       CONTINUE
         RETURN
         END
C ============================================================================
C
C Subroutine EXPAV averages together symmetrically equivalent experimental
C beams as specified by NBEA
C
C ============================================================================
C
      SUBROUTINE EXPAV(AE,EE,NBED,NEE,BENAME,NBEA,NBE,IPR,A,NV,IEERG)
C
      DIMENSION AE(NBED,IEERG),EE(NBED,IEERG),NEE(NBED),NBEA(NBED)
      DIMENSION A(IEERG),NV(IEERG),BENAME(5,NBED)
C
110   FORMAT (' EXP. ENERG. AND INTENS. AFTER AVERAGING IN BEAM',1I4,
     & /,50(5(1F7.2,1E13.4,3X),/))
C
      NBE=1
      DO 90 IB=1,NBED
         IF (NBEA(IB).NE.0) THEN
            EMIN=1.E6
            EMAX=0.
C
C  LOOP OVER BEAMS TO BE AVERAGED TOGETHER WITH BEAM IB TO OBTAIN
C  MINIMUM AND MAXIMUM ENERGIES
C
            DO 45 IB1=IB,NBED
               IF (NBEA(IB1).EQ.NBEA(IB)) THEN
                  N=NEE(IB1)
                  DO 43 IE=1,N
                     IF (EMIN.GT.EE(IB1,IE)) THEN
                        EMIN=EE(IB1,IE)
                        IBMIN=IB1
                     ENDIF
                     IF (EMAX.LT.EE(IB1,IE)) EMAX=EE(IB1,IE)
43                CONTINUE
               ENDIF
45          CONTINUE
C
C  FIND ENERGY INCREMENT
C
            EMIP=1.E6
            N=NEE(IBMIN)
            DO 48 IE=1,N
               EP=EE(IBMIN,IE)
               IF (EP.LT.EMIP.AND.EP.GT.(EMIN+.01)) EMIP=EP
48          CONTINUE
            DE=EMIP-EMIN
            NEMAX=INT((EMAX-EMIN)/DE+0.0001)+1
            DO 50 IE=1,NEMAX
               NV(IE)=0
               A(IE)=0.
50          CONTINUE
            IEMAX=0
            NBEAT=NBEA(IB)
C
C  LOOP OVER SAME BEAMS AGAIN, AVERAGING OVER THESE BEAMS AND
C  REORDERING ENERGIES
C
            DO 70 IB1=IB,NBED
               IF (NBEA(IB1).EQ.NBEAT) THEN
                  NEMAX=NEE(IB1)
                  DO 60 IE=1,NEMAX
                     IEN=INT((EE(IB1,IE)-EMIN)/DE+0.0001)+1
                     IF (IEN.GT.IEMAX) IEMAX=IEN
                     A(IEN)=A(IEN)+AE(IB1,IE)
                     NV(IEN)=NV(IEN)+1
60                CONTINUE
                  NBEA(IB1)=0
               ENDIF
70          CONTINUE
            DO 80 IE=1,IEMAX
               EE(NBE,IE)=EMIN+FLOAT(IE-1)*DE
               AE(NBE,IE)=A(IE)/FLOAT(NV(IE))
80          CONTINUE
            NEE(NBE)=IEMAX
C
C  KEEP NAME OF FIRST BEAM ENCOUNTERED IN SET OF BEAMS TO BE AVERAGED
C
            DO 85 I=1,5
               BENAME(I,NBE)=BENAME(I,IB)
85          CONTINUE
            NBE=NBE+1
         ENDIF
90    CONTINUE
      NBE=NBE-1
      IF (IPR.GE.2) THEN
         DO 100 IB=1,NBE
            N=NEE(IB)
            WRITE (1,110) IB,(EE(IB,IE),AE(IB,IE),IE=1,N)
100      CONTINUE
      ENDIF
      RETURN
      END
C =======================================================================
C 
C Subroutine FUNCV returns the value of the R-factor at a point in parameter
C space specified by COORD and DISP.
C
C Input Parameters;
C =================
C
C FVAL           =  Value of the R-factor at the point specified by PR and
C                   DISP
C NLAY           =  Number of subplanes in composite layer
C ADISP          =  Current coordinates of required structure assembled from
C                   PR and DISP (Output from routine SETCOR)
C VOPT           =  Optimal value of the inner potential
C ILOOK          =  Look up table for domain averaging
C ACOORD         =  Storage for symmetry equivalent corrdinates
C MICUT,MJCUT    =  Location of tensor element in tensor before truncation
C                   (used to reconstruct tensor)
C PSQ            =  Labels of exit beams in reciprocal space
C JYLM           =  Work space
C BJ             =  Work space
C YLM            =  Work Space
C XISTS          =  Reference structure plane wave amplitudes
C XIST           =  Total plane wave amplitudes for the current structure
C AT             =  Intensities generated from amplitudes
C INBED          =  Number of experimental beams
C IEERG          =  Total energy range after interpolation
C AE             =  Experimental intensities
C EE             =  Experimental energies
C NEE            =  Number of data points in each experimental beam
C NBEA           =  Beam averaging information
C BENAME         =  Identifier for each experimental beam
C IPR            =  Print control parameter
C XPL,YPL,NNN    =  Work space
C AP             =  First derivative of experimental intensities
C APP            =  Second derivative of experimental intensities
C YE             =  Experimental Y-function
C TSE,TSE2.TSEP  =  Integrals over experimental data
C TSEP2, TSEPP,
C TSEPP2,TSEY2
C WR             =  R-factor weights for 10 R-factor average
C WB             =  Beam weights within each R-factor
C IBP            =  Theoretical beam averaging information
C
C COMMON BLOCKS
C =============
C
C VMIN,VMAXM,DV  =  Range of search over inner potential (block VINY)
C EINCR          =  Grid step to be used after interpolation (block VINY)
C EI,EF,DE       =  Energy range of calculation (block ENY)
C NERG           =  Number of energy points (block ENY)
C NSYM           =  Symmetry code of surface (block ENY)
C NDOM           =  Number of symmetry equivalent domains (block ENY)
C VV             =  Real part of inner potential (block ENY)
C LSMAX          =  Max L value to be used in single centre expansion 
C                   (block TLVAL)
C LSMMAX         =  (LSMAX+1)**2 (block TLVAL)
C ICUT           =  Number of tensor elements for each energy, each layer 
C                   and each beam (block TLVAL)
C LSMAX1         =  LSMAX+1 (block TLVAL)
C NT0            =  Number of exit beams (block TLVAL)
C IQSIZ          =  Total size of tensor after truncation 
C                   =(ICUT*NLAY*NERG*NT0) (block TLVAL)
C
C
C ==========================================================================
C
      SUBROUTINE FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
C
      DIMENSION ADISP(NLAY,3),ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0),AT(NT0,NERG)
      DIMENSION ETH(NT0,IEERG),IBP(NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),YT(NT0,IEERG)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),PRE2(50)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1),DELXI,AAK,AAJ,PRE,CI
      COMPLEX JYLM(LSMMAX),DELXI,AAK,AAJ,PRE,CI
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA,FI
      COMMON /ENY/EI,EF,DE,NERG,NSYM,NDOM,VV,VPIS
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
C
C Set constants
C
      CI=CMPLX(0.0,1.0)
C
C Find maximum X displacement to correct for perpendicular displacements
C of the surface barrier (program will place surface barrier a distance
C ASE above the top atom). The Following is not correct 
C when there is more than one atom per unit cell in the outermost layer
C
      AMAXD=100.0
      IOFF=0
      DO 101 I=1,1
        IF(ADISP(I,1).LE.AMAXD)AMAXD=ADISP(I,1)
101   CONTINUE
C
C First create symmetry equivalent displacements for all specified
C domains
C
      CALL DISDOM(ADISP,ACOORD,ILOOK,NLAY,0)
C
C Begin loop over energies.
C
C      I2=MCLOCK()
      DO 100 I=1,NERG
C
C EEV is energy in eV, E is shifted energy in Hartrees
C
         EEV=EI+FLOAT(I-1)*DE
        E=EEV/27.21+VV
C
C Set parallel and perpendicular components of momentum for incident 
C direction
C
         AKZ=SQRT(2.0*(E-VV))*COS(THETA)
         AK=SQRT(2.0*(E-VV))*SIN(THETA)
         AK2=AK*COS(FI)
         AK3=AK*SIN(FI)
C
C Compute scaling needed to correct for shift in surface barrier
C (This is for incident beam)
C
        AAK=CMPLX(2.0*(E-VV)-AK2*AK2-AK3*AK3,-2.0*VPIS)
        AAK=CSQRT(AAK)
C
C Set up indices of exit beams in units of reciprocal lattice vectors
C Also compute beam dependent prefactors and emergence condition for all
C beams
        DO 110 J=1,NT0
C
C Set parallel and perpendicular componenets of momentum for each 
C exit beam
C
            ETH(J,I)=EEV
            AK2M=-AK2-PSQ(1,J)
            AK3M=-AK3-PSQ(2,J)
C
C Compute scaling needed to correct for shift in surface barrier
C (This is for exit direction)
C
           AAJ=CMPLX(2.0*(E-VV)-AK2M*AK2M-AK3M*AK3M,-2.0*VPIS)
           AAJ=CSQRT(AAK)
C
C Does this beam emerge?
C
            AAA=2.0*(E-VV)-AK2M*AK2M-AK3M*AK3M
C
C Work out scaling needed to account for shift in surface barrier
C (PRE=(inward contribution + outward contribution for this beam)*shift)
C
            PRE=(AAK+AAJ)*AMAXD
            PRE=CEXP(-CI*PRE)
            IF(AAA.GT.0)THEN
               AAB=SQRT(AAA)
               PRE2(J)=AAB/AKZ/FLOAT(NDOM)*CABS(PRE)*CABS(PRE)
            ELSE
               PRE2(J)=0.0
            ENDIF
            AT(J,I)=0.0
110      CONTINUE
C
C Work out the change in this beam (DELXI) for the current
C structure, including all atoms and domains.
C
         DO 111 K=1,NDOM
            IOFFOL=IOFF
            CALL DELX2(QS,ACOORD,MICUT,MJCUT,NLAY,YLM,BJ,
     &        JYLM,E,VPIS,IOFFOL,K,NLTIN,LPOINT,
     &        XIST,PRE2,NERG,I)
C
C Set up the new plane wave amplitude XIST as the sum of the reference
C structure plane wave amplitude (XISTS) and the change created by moving
C to the current structure (XIST from DELX).
C
            DO 1000 J=1,NT0
               XIST(J,I)=XISTS(J,I)+XIST(J,I)
               AT(J,I)=AT(J,I)+CABS(XIST(J,I))*CABS(XIST(J,I))
1000        CONTINUE
111      CONTINUE  
         IOFF=IOFFOL 
C
C Finally include prefactors 
C
         DO 1001 J=1,NT0
            AT(J,I)=AT(J,I)*PRE2(J)
1001     CONTINUE
C              AT(J,I)=AT(J,I)
C               AT(J,I)=AT(J,I)
100   CONTINUE
C      I3=MCLOCK()
C      FTIME=FLOAT(I3-I2)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU,'funcv'
      CALL RFAC(AT,ETH,INBED,IEERG,AE,EE,NEE,IPR,XPL,
     & YPL,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,
     & IBP,NT0,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,
     & RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,FVAL,VOPT,YT)
C      I4=MCLOCK()
C      FTIME=FLOAT(I4-I3)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU,'rfac'
      RETURN
      END
C =======================================================================
C
C Subroutine GRID evaluates the RFACTOR along the principal directions 
C at 2*MSTEP+1 points centered on P (spacing = ZSTEP) and plot it together
C with the predicted parabola
C
C =======================================================================
C
Cga made into a function 
Cga      SUBROUTINE GRID(ZSTEP,MSTEP,NLAY,NDIM,DISP,XI,
      function GRID(ZSTEP,MSTEP,NLAY,NDIM,DISP,XI,
     & NNDIM,ADISP,IPR,PT,XIT,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,NGRID,
     & LSFLAG,NDIML,DIREC,NLTIN,LPOINT)
C
C MAXC is the maximum dimension of the parameter space where we perform
C optimazation of the R factor
C
      PARAMETER (MAXC=100 )
C
      DIMENSION DISP(NLAY,3),PT(NNDIM),XI(NNDIM,NNDIM)
      DIMENSION XIT(NNDIM),LLFLAG(NLAY+NNST),EIGEN(MAXC)
      DIMENSION DIREC(NLAY,2,2),NDIML(NLAY),LSFLAG(NLAY)
      DIMENSION XX1(MAXC)
      DIMENSION ADISP(NLAY,3)
      DIMENSION ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1)
      COMPLEX JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /RPL/DVOPT
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
C
!501   FORMAT (' ======================= ')
!505   FORMAT (/' CONVERGENCE TOLERANCE ACHIEVED =',F7.6)
!506   FORMAT (/' COORDINATES AT MINIMUM;'/)
!510   FORMAT (3F9.4)
!514   FORMAT (1F7.4)
!515   FORMAT (' OPTIMUM VALUE OF INNER POTENTIAL =',1F7.4)
C
C define the vector P corresponding to the minimum configuration,
C the set of directions along which to evaluate 
C the Rfactor, and the eigenvalues of the hessian at the minimum
C
      pi=atan(1.0)*4.
      CALL SRETV(FMIN,XI,EIGEN,PT,LLFLAG,LSFLAG,NNDIM,
     & JJDIM,NLAY,ADISP,NDIML,DIREC)
      READ(10,*) (EEAVT(J),J=1,NBE+1)
      IF(NGRID.EQ.2) THEN
        DO 999 I=1,NNDIM
           READ (4,*) (XI(K,I),K=1,NNDIM)
999     CONTINUE
      ENDIF
C
C Looop over the different principal directions, NGRID=0 corresponds
C to cartesian coordinates, NGRID=1 to principal directions,NGRID=2
C to user specified directions,NGRID=3 to a user modified program.
C
      IFUNC=1
      IF(NGRID.EQ.0) THEN
          DO 100 I=1,NNDIM
          DO 110 J=1,NNDIM
             IF (J.EQ.I) THEN
                XI(I,J)=1.0
             ELSE
                XI(I,J)=.0
             ENDIF
110       CONTINUE
100       CONTINUE
      ENDIF
C
C This is for the user to play with.  
C The example provided is for  a tilting-angle error bar. 
C
      IF(NGRID.EQ.3) THEN
          blold=1.4814
          BL=sqrt((PT(1)-PT(3)-blold)**2 + (PT(2)-Pt(4))**2)
          bld=blold-bl
          PT(1)=PT(3)
          PT(2)=PT(4)
          DO 88 I=1,Mstep
            th=float(i)/180.*pi/4.
            DO 881 K=1,NNDIM
              XX1(K)=PT(K)
881         CONTINUE
            XX1(1)=PT(1)+bld+bl*(1.-cos(th))
            XX1(2)=PT(2)+bl*sin(th)
          CALL SETCOR2(LSFLAG,NDIML,DIREC,XX1,DISP,ADISP,
     &     NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IFUNC=IFUNC+1
          VOPT=VOPT+DVOPT+VV*27.21
          WRITE (15,*)float(i)/4.,FVAL
88        continue   
          goto 55
      ENDIF
C
C end of the user play.
C
      DO 10 II=1,NNDIM
C
C define the displacement vector
C
        DO 120 I=1,NNDIM
           XIT(I)=XI(I,II)*ZSTEP
120     CONTINUE
C
C evaluate Rfactor
C
        DO 130 I=-MSTEP,MSTEP
          FI=FLOAT(I)
          DO 140 K=1,NNDIM
            XX1(K)=PT(K)+FI*XIT(K)
140       CONTINUE
          CALL SETCOR2(LSFLAG,NDIML,DIREC,XX1,DISP,ADISP,
     &     NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IFUNC=IFUNC+1
          VOPT=VOPT+DVOPT+VV*27.21
          WRITE (15,*) ZSTEP*FI,FVAL
C         WRITE (16,514) FVAL
C         WRITE (2,515) VOPT
130      continue
         WRITE (15,*) 
         IF(NGRID.EQ.1) THEN
           DO 20 KK=-25,25
             FI=FLOAT(KK)
             X=ZSTEP*FI
             Y=EIGEN(II)*X*X/2.
             WRITE (15,*) X,Y+FMIN
20         CONTINUE
           WRITE (15,*) 
         ENDIF 
10    CONTINUE
Cga
Cga55    RETURN
 55   continue
      grid = vopt
      return
      END
C===================================================================== 
C
C Subroutine HARMON generated the spehrical harmonics YLM for a real
C vector C.
C
C Input Parameters;
C =================
C 
C C(i), i=1,3         =   Displacement vector
C LMAX                =   Maximum value of the angular momentum
C LMMAX               =   (LMAX+1)*(LMAX+1)
C YLM                 =   Ouput spherical harmonic
C
C Author: ROUS. Modifications by WANDER.
C
C=====================================================================
C
      SUBROUTINE HARMON(YLM,C,LMAX,LMMAX)
C
      COMPLEX CT,ST,CF,YLM(LMMAX)
      DIMENSION C(3)
C
      CD=C(2)*C(2)+C(3)*C(3)
      YA=SQRT(CD+C(1)*C(1))
      B=0.0
      CF=CMPLX(1.0,0.0)
      IF (CD.GT.1.0E-7) THEN
         B=SQRT(CD)
         CF=CMPLX(C(2)/B,-C(3)/B)
      ENDIF
      CT=C(1)/YA
      ST=B/YA
      CALL SPHRM(LMAX,YLM,LMMAX,CT,ST,CF)
      RETURN
      END
C =========================================================================
C
C Subroutine INTPOL interpolates LEED intensities onto a working grid with 
C step size EINCR (eV)
C
C =========================================================================
C
      SUBROUTINE INTPOL(A,AP,APP,NBD,NE,NB,E,EINCR,IPR,X,WORYT,IEERG)
C
      DIMENSION WORYT(IEERG),X(IEERG),A(NBD,IEERG),E(NBD,IEERG)
      DIMENSION NE(NBD),WORYT2(2000),WORYT4(2000),AP(NBD,IEERG)
      DIMENSION APP(NBD,IEERG),APPP(20,2000)
C
C The dimension of APPP limits the number of inequivalent beams to 20
C 2000 is the maximum number of energy points in a single beam after
C interpolation
C
      COMMON /ENY/EI,EF,DE,NERG,NSYM,NDOM,VV,VPIS
6     FORMAT (
     & '** IN PRESENT BEAM TOO FEW ENERGY VALUES FOR INTERPOLATION')
80    FORMAT ('INTENSITIES AFTER INTERPOLATION IN BEAM',1I4,/,50(5
     & (1F7.2,1E13.4,3X),/))
C
      ITIL=0
      ITIH=0
      E2=EINCR/2.
C      E2=1.
      DO 60 IB=1,NB
         NEM=NE(IB)
cjcm         write(*,*) 'INTPOL: NE(IB)  = ',IB, NEM 
C
C  FIND FIRST NON-ZERO INTENSITY (FOR THEORY, WHERE NON-EMERGENCE OF
C  CURRENT BEAM CAN OCCUR)
C
         DO 30 IE=1,NEM
            IMIN=IE
            IF (A(IB,IE).GT.1.E-6) GOTO 91
30       CONTINUE
cjcm         write(*,*) 'INTPOL(30): IMIN = ',IMIN 
cjcm         write(*,*) 'INTPOL(30): NEM  = ',NEM 
91       CONTINUE
         IF (IMIN.NE.NEM) THEN
            LMIN=INT((E(IB,IMIN)-E2)/EINCR)+1
            LMIN=MAX0(LMIN,0)
C XMIN (XMAX)is the minimum (maximum) energy on the working grid 
C for which the computed intensities are available
            XMIN=FLOAT(LMIN)*EINCR
cjcm            write(*,*) 'INTPOL: E(IB,NEM)  = ', E(IB,NEM)
cjcm            write(*,*) 'INTPOL: E2, EINCR  = ', E2, EINCR
            LMAX=INT((E(IB,NEM)+E2)/EINCR)
            XMAX=FLOAT(LMAX)*EINCR
C
C  NPTS IS NO. OF POINTS USED ON THE INTERPOLATION GRID
C TO GET CONTINUES RFACTOR WHEN VARYING VOPT ONE SHOULD HAVE
C THEORETICAL DATA UP TO E=EMAXexp +|VOPT| where VOPT is the maximum
C expected shift in the inner potential. THE SAME IS TRUE FOR
C THE MINIMUM THEORETICAL Emin=EMINexp -|VOPT| 
C
            NPTS=LMAX-LMIN+1
cjcm debug
c     
cjcm            write(*,*) 'LMAX, LMIN, NPTS, IMIN', LMAX, LMIN, NPTS, IMIN
            DO 5 I=IMIN,NEM
               X(I-IMIN+1)=E(IB,I)
               WORYT(I-IMIN+1)=A(IB,I)
5           CONTINUE
            NEM=NEM-IMIN+1
            IF (NEM.GE.1) THEN
CTEST
            YP1=1.E10
            YPN=1.E10
C XVAL can be < or > of X(1) because of the shifts in X due to a varying
C VOPT. For purpose of interpolation however we want it to be > than X(1)
            XMIN=XMIN+EINCR
            XVAL=XMIN
cjcm            IF(XVAL.LT.X(1)) PAUSE 'AGAIN!!'
            KLO=1
            KHI=2
C this shift allows all points in the grid to be included in the energy 
C interval for which we have data
            NPTS=NPTS-2
            NE(IB)=NPTS
c            write(*,*) 'INTPOL: NE(IB)  = ',IB, NE(IB) 
            CALL SPLINE(X,WORYT,NEM,YP1,YPN,WORYT2)
            CALL SPLINE(X,WORYT2,NEM,YP1,YPN,WORYT4)
               DO 10 I=1,NPTS
28          IF(XVAL.GE.X(KLO).AND.XVAL.LT.X(KHI))GOTO 27 
                  KLO=KLO+1
                  KHI=KHI+1
cjcm                  IF(KHI.GT.NEM)PAUSE 'PROBLEM'
                  GOTO 28
27             CALL SPLINT(X,WORYT,WORYT2,NEM,XVAL,A(IB,I),
     %            AP(IB,I),KLO,KHI)
               CALL SPLINT(X,WORYT2,WORYT4,NEM,XVAL,APP(IB,I),
     %            APPP(IB,I),KLO,KHI)
C
C  INTERPOLATE (AND SET NEGATIVE INTENSITIES TO ZERO)
C
C              CALL YVAL(A(IB,I),AP(IB,I),XVAL,WORYT,X,NEM,ITIL,ITIH)
C                  IF (A(IB,I).LT.0.0) THEN
C                    A(IB,I)=0.0
C                    AP(IB,I)=0.0
C                  ENDIF
                  XVAL=XVAL+EINCR
10             CONTINUE
               E(IB,1)=XMIN
               DO 50 IE=2,NPTS
                  E(IB,IE)=E(IB,IE-1)+EINCR
50             CONTINUE
               GOTO 60
            ENDIF
         ENDIF
C
C  TOO FEW ENERGY VALUES FOR INTERPOLATION. THIS BEAM WILL BE SKIPPED
C  FROM NOW ON
C
         NE(IB)=0
         IF (IPR.GE.2) THEN
            WRITE (1,6)
         ENDIF
60    CONTINUE
      IF (IPR.GE.2) THEN
         DO 70 IB=1,NB
            N=NE(IB)
            IF (N.NE.0) WRITE (1,80) 
     %       IB,(E(IB,IE),A(IB,IE),AP(IB,IE),IE=1,N)
70       CONTINUE
      ENDIF
      RETURN
      END
C =======================================================================
C
C Subroutine JACOBI diagonalizes the symmetric matrix A and returns
C eigenvalues (in D) and eigenvectors (in V)
C From Numerical Recipes
C
C =======================================================================
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.
11      CONTINUE
        V(IP,IP)=1.
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
cjcm      PAUSE '50 iterations should never happen'
      RETURN
      END
C =======================================================================
C
C Subroutine LINMIN prepares for the one dimensional minimization. 
C Given an NNDIMensional point P and an NNDIMemsional direction XIT
C moves and reset P to where the R-factor(P) takes on aminimum
C along the direction XIT from P, and replaces XIT by the actual vector
C displacement that P was moved. Returns FRET as the value of the
C Rfactor at the returned location P.
C
C References 'NUMERICAL RECIPES' W.H.Press,B.P.Flannery, et al. Cambridge
C             University Press
C 
C Modifications by BARBIERI
C
C Input Parameters;
C =================
C
C AX                     = First bracketing point candidate 
C XX                     = Second bracketing point candidate
C BX                     = Third bracketing point candidate
C NLAY                   = NUMBER OF LAYERS IN COMPOSITE LAYER
C NDIM                   = DIMENSIONALITY OF SEARCH (1=X AXIS ONLY
C                                                    3=X,Y,Z AXES)
C P                      = ARRAY CONTAINING THE INITIAL STARTING POINT
C                          IN THE SEARCH
C DISP                   = COORDINATES INPUT BY USER
C NNDIM                  = TOTAL NUMBER OF DIMENSIONS IN SEARCH (=NLAY*NDIM)
C ADISP                  = GEOMETRY OF CURRENT POINT IN SAME FORMAT AS DISP
C                          (USED AS INPUT TO FUNCV)
C DVOPT                  = SHIFT IN THE INNER POTENTIAL FOR THE STARTING 
C                          CONFIGURATION (in COMMON /RPL )
C LLFLAG                  = INDICATES WHETHER THE LAYER (OR NON STRUCTURAL
C                          PARAMETER) COORDINATES HAVE TO BE VARIED 
C                          IN THE SEARCH
C LSFLAG                 =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C                          NSYM) IN THE COMPOSITE LAYER. LLFLAG(i)=LLFLAG(j)
C                          INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C                          AS EQUIVALENT IN THE SEARCH
C NDIML                  =  ARRAY GIVING THE EFFECTIVE DIMENSIONALITY OF ATOM
C                           j (ACTUALLY TO BE USED IN CONJUNCTION WITH LLFLAG)
C DIREC                  =  SET OF DIRECTIONS
C VOPT                   = INNER POTENTIAL OF CURRENT INPUT (used as input
C                          in FUNCV), in COMMON /NSTR
C FTOL2                  = CONVERGENCE CRITERIA FOR POWELL SEARCH
C ITMAX                  = MAXIMUM NUMBER OF ITERATIONS TO BE PERFORMED BY
C                          POWELL
C ISTART                 = 0 USE PARAMETERS FROM TLEED4.I
C                          1 RESTART SEARCH FROM COORDINATES IN RESTART.D
C IPR                    = PRINT CONTROL PARAMETER
C ILOOK                  = LOOKUP TABLE FOR DOMAIN AVERAGING
C ACOORD                 = STORAGE FOR SYMMETRY EQUIVALENT COORDINATES
C MICUT,MJCUT,PSQ,JYLM,  = DUMMY ARRAYS TO BE PASSED TO ROUTINE FUNCV
C BJ,YLM,QS,AT,XISTS,
C XIST,INBED,IEERG,AE,EE,
C NEE,NBEA,BENAME,IPR,XPL,
C YPL,NNN,AP,APP,YE,TSE,
C TSE2,TSEP,TSEP2,TSEPP,
C TSEPP2,TSEY2,WR,WB,IBP,
C ETH
C in Common:
C IFUNC                  =Total number of function evaluations
C
C ============================================================================
      SUBROUTINE LINMIN(XMIN,AX,XX,BX,FRET,P,XIT,NLAY,NDIM,NNDIM,ADISP,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,DISP,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     & NDIML,DIREC,NLTIN,LPOINT)
C
C MAXC is the maximum dimension of parameter space, TOL is the
C tolerance in parameter space ( or fractional tolerance: see BRENT ) 
C with which we accept this one dimensional minimization
C
      PARAMETER (MAXC=100 ,tol=1.e-4)
C
      DIMENSION PCOM(MAXC),XICOM(MAXC)
      DIMENSION DIREC(NLAY,2,2),NDIML(NLAY),LSFLAG(NLAY)
      DIMENSION XIT(NNDIM),DISP(NLAY,3),P(NNDIM)
      DIMENSION ADISP(NLAY,3),LLFLAG(NLAY+NNST)
      DIMENSION ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1)
      COMPLEX JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /RPL/DVOPT
      COMMON /NSTR/VOPT,NNST,NNSTEF
C
C PCOM and XICOM are needed to compute the R-factor from P along the direction
C XIT
C
C      I2=MCLOCK()
      do 11 j=1,NNDIM
        pcom(j)=p(j)
        xicom(j)=xit(j)
11    continue
C
C bracket the minimum
C
C the value at x=0. is XMIN. Use -xx and xx as first attempt to bracket
      ax=-xx
      CALL MNBRAK(XMIN,ax,xx,bx,fa,fx,fb,pcom,xicom,NLAY,NDIM,
     & DISP,NNDIM,ADISP,ILOOK,ACOORD,MICUT,
     & MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     & NDIML,DIREC,NLTIN,LPOINT)
C
C perform the minimization from PCOM alon XICOM
C
C however, if we are already close enough to the minimum (as specified
C by typicall error bars for LEED analysis) we will stop here
C 
C       if(abs(ax-bx).le..02) then
C          FRET=fx
C          XMIN=xx
C          goto 99
C       endif
      FRET=BRENT(ax,xx,bx,fa,fx,fb,tol,xmin,pcom,xicom,NLAY,NDIM,
     & DISP,NNDIM,ADISP,ILOOK,ACOORD,MICUT,
     & MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     & NDIML,DIREC,NLTIN,LPOINT)
      do 12 j=1,NNDIM
        xit(j)=xmin*xit(j)
        p(j)=p(j)+xit(j)
12    continue
C      I3=MCLOCK()
C      FTIME=FLOAT(I3-I2)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU,'total linmin time'
      return
      end
C ============================================================================
C
C Subroutine MNBRACK takes the initial points AX,XX,BX and returns
C a new triplet of points which bracket a minimum of the R-factor
C from the point PCOM along the direction XICOM
C
C References 'NUMERICAL RECIPES' W.H.Press,B.P.Flannery, et al. Cambridge
C             University Press
C 
C Modifications by BARBIERI
C
C Input Parameters;
C =================
C
C NLAY                   = NUMBER OF LAYERS IN COMPOSITE LAYER
C NDIM                   = DIMENSIONALITY OF SEARCH (1=X AXIS ONLY
C                                                    3=X,Y,Z AXES)
C PCOM                   = ARRAY CONTAINING THE POINT IN PARAMETER SPACE
C XICOM                  = DIRECTION ALONG WHICH TO EVALUATE THE Rfactor
C DISP                   = COORDINATES INPUT BY USER
C NNDIM                  = TOTAL NUMBER OF DIMENSIONS IN SEARCH (=NLAY*NDIM)
C ADISP                  = GEOMETRY OF CURRENT POINT IN SAME FORMAT AS DISP
C                          (USED AS INPUT TO FUNCV)
C DVOPT                  = SHIFT IN THE INNER POTENTIAL FOR THE STARTING 
C                          CONFIGURATION (in COMMON /RPL )
C LLFLAG                  = INDICATES WHETHER THE LAYER (OR NON STRUCTURAL
C                          PARAMETER) COORDINATES HAVE TO BE VARIED 
C                          IN THE SEARCH
C LSFLAG                 =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C                          NSYM) IN THE COMPOSITE LAYER. LLFLAG(i)=LLFLAG(j)
C                          INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C                          AS EQUIVALENT IN THE SEARCH
C NDIML                  =  ARRAY GIVING THE EFFECTIVE DIMENSIONALITY OF ATOM
C                           j (ACTUALLY TO BE USED IN CONJUNCTION WITH LLFLAG)
C DIREC                  =  SET OF DIRECTIONS
C VOPT                   = INNER POTENTIAL OF CURRENT INPUT (used as input
C                          in FUNCV), in COMMON /NSTR
C ITMAX                  = MAXIMUM NUMBER OF ITERATIONS TO BE PERFORMED BY
C                          POWELL
C ISTART                 = 0 USE PARAMETERS FROM TLEED4.I
C                          1 RESTART SEARCH FROM COORDINATES IN RESTART.D
C IPR                    = PRINT CONTROL PARAMETER
C ILOOK                  = LOOKUP TABLE FOR DOMAIN AVERAGING
C ACOORD                 = STORAGE FOR SYMMETRY EQUIVALENT COORDINATES
C MICUT,MJCUT,PSQ,JYLM,  = DUMMY ARRAYS TO BE PASSED TO ROUTINE FUNCV
C BJ,YLM,QS,AT,XISTS,
C XIST,INBED,IEERG,AE,EE,
C NEE,NBEA,BENAME,IPR,XPL,
C YPL,NNN,AP,APP,YE,TSE,
C TSE2,TSEP,TSEP2,TSEPP,
C TSEPP2,TSEY2,WR,WB,IBP,
C ETH
C in Common:
C IFUNC                  =Total number of function evaluations
C
C
C  ax,bx,cx AND the function values fa,fb,fc are returned such that
C  they bracket a minimum
C
C =======================================================================
C
      SUBROUTINE MNBRAK(XMIN,ax,bx,cx,fa,fb,fc,pcom,xicom,NLAY,NDIM,
     & DISP,NNDIM,ADISP,ILOOK,ACOORD,MICUT,
     & MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     & NDIML,DIREC,NLTIN,LPOINT)
C
      PARAMETER (gold=1.618034, glimit=100., tiny=1.e-8)
      PARAMETER (MAXC=100 )
C
      DIMENSION PCOM(MAXC),XICOM(MAXC),LLFLAG(NLAY+NNST)
      DIMENSION DIREC(NLAY,2,2),NDIML(NLAY),LSFLAG(NLAY)
      DIMENSION ADISP(NLAY,3),DISP(NLAY,3),XT(MAXC)
      DIMENSION ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1)
      COMPLEX JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /RPL/DVOPT
      COMMON /NSTR/VOPT,NNST,NNSTEF
      do 11 j=1,NNDIM
        xt(j)=pcom(j)+ax*xicom(j)
11    continue
C
C Set coordinate matrix ADISP from temporary store P
C
      CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      IFUNC=IFUNC+1
      fa=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
      do 2 j=1,NNDIM
        xt(j)=pcom(j)+bx*xicom(j)
2     continue
      CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      IFUNC=IFUNC+1
      fb=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
      if(fb.gt.XMIN.AND.fa.gt.XMIN) then
        cx=bx
        fc=fb
        bx=0.
        fb=XMIN 
        return
      endif
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+gold*(bx-ax)
      do 3 j=1,NNDIM
        xt(j)=pcom(j)+cx*xicom(j)
3     continue
      CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      IFUNC=IFUNC+1
      fc=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
        ulim=bx+glimit*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          do 4 j=1,NNDIM
            xt(j)=pcom(j)+u*xicom(j)
4         continue
          CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IFUNC=IFUNC+1
          fu=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            go to 1
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            go to 1
          endif
          u=cx+gold*(cx-bx)
          do 5 j=1,NNDIM
            xt(j)=pcom(j)+u*xicom(j)
5         continue
          CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IFUNC=IFUNC+1
          fu=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
        else if((cx-u)*(u-ulim).gt.0.)then
          do 6 j=1,NNDIM
            xt(j)=pcom(j)+u*xicom(j)
6         continue
          CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IFUNC=IFUNC+1
          fu=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+gold*(cx-bx)
            fb=fc
            fc=fu
            do 7 j=1,NNDIM
              xt(j)=pcom(j)+u*xicom(j)
7           continue
            CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
            IFUNC=IFUNC+1
            fu=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          do 8 j=1,NNDIM
            xt(j)=pcom(j)+u*xicom(j)
8         continue
          CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IFUNC=IFUNC+1
          fu=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
        else
          u=cx+gold*(cx-bx)
          do 9 j=1,NNDIM
            xt(j)=pcom(j)+u*xicom(j)
9         continue
          CALL SETCOR2(LSFLAG,NDIML,DIREC,XT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IFUNC=IFUNC+1
          fu=fval
cjcm      write(*,*) 'BRENT: ', IFUNC, fval
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        go to 1
      endif
      return
      end
C  file LEEDSATL.SB5  Feb. 29, 1996
C
C**************************************************************************
C  Symmetrized Automated Tensor LEED (SATLEED):  subroutines, part 5
C  Version 4.1 of Automated Tensor LEED
C
C =========================================================================
C
C Subroutine OPSIGN determines the fractional energy range where the
C theoretcial and experimental data sets have slopes of opposite sign.
C
C From the VAN HOVE R-factor program.
C
C ========================================================================
C
      SUBROUTINE OPSIGN(AEP,NBED,IBE,ATP,NBTD,IBT,EINCR,IE1,IE2,NV,ROS,
     &  IEERG)
C
      DIMENSION AEP(NBED,IEERG),ATP(NBTD,IEERG)
C
      EINCR2=0.5*EINCR
      ROS=0.
      DO 10 IE=IE1,IE2
         IES=IE+NV
         IF ((AEP(IBE,IE)*ATP(IBT,IES)).LE.0.) THEN
            ROS=ROS+EINCR
            IF (IE.EQ.IE1.OR.IE.EQ.IE2) ROS=ROS-EINCR2
         ENDIF
10    CONTINUE
      RETURN
      END
C =======================================================================
C
C Subroutine POWELL performs the direction set search for the R-factor
C minimum. At exit,the vector P is set to the best point found, XI is the 
C the current direction set, FRET is the returned function value at P, and
C ITER is the number of iteration taken. The routine LINMIN is used 
C
C References 'NUMERICAL RECIPES' W.H.Press,B.P.Flannery, et al. Cambridge
C             University Press
C 
C Author BARBIERI
C
C Input Parameters;
C =================
C
C NLAY                   = NUMBER OF LAYERS IN COMPOSITE LAYER
C NDIM                   = DIMENSIONALITY OF SEARCH (1=X AXIS ONLY
C                                                    3=X,Y,Z AXES)
C P                      = ARRAY CONTAINING THE INITIAL STARTING POINT
C                          IN THE SEARCH
C XI                     = ARRAY STORING THE INITIAL DIRECTION SET
C DISP                   = COORDINATES INPUT BY USER
C NNDIM                  = TOTAL NUMBER OF DIMENSIONS IN SEARCH (=NLAY*NDIM)
C ADISP                  = GEOMETRY OF CURRENT POINT IN SAME FORMAT AS DISP
C                          (USED AS INPUT TO FUNCV)
C DVOPT                  = SHIFT IN THE INNER POTENTIAL FOR THE STARTING 
C                          CONFIGURATION (in COMMON /RPL )
C LLFLAG                  = INDICATES WHETHER THE LAYER (OR NON STRUCTURAL
C                          PARAMETER) COORDINATES HAVE TO BE VARIED 
C                          IN THE SEARCH
C LSFLAG                 =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C                          NSYM) IN THE COMPOSITE LAYER. LLFLAG(i)=LLFLAG(j)
C                          INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C                          AS EQUIVALENT IN THE SEARCH
C LRFLAG                 =  0 OR 1, CORRESPONDING TO TENSOR SEARCH AND
C                           WITHOUT TENSOR SEARCH, SEPERATELY.
C NDIML                  =  ARRAY GIVING THE EFFECTIVE DIMENSIONALITY OF ATOM
C                           j (ACTUALLY TO BE USED IN CONJUNCTION WITH LLFLAG)
C DIREC                  =  SET OF DIRECTIONS
C VOPT                   = INNER POTENTIAL OF CURRENT INPUT (used as input
C                          in FUNCV), in COMMON /NSTR
C FTOL2                  = CONVERGENCE CRITERIA FOR POWELL SEARCH
C ITMAX                  = MAXIMUM NUMBER OF ITERATIONS TO BE PERFORMED BY
C                          POWELL
C ISTART                 = 0 USE PARAMETERS FROM TLEED4.I
C                          1 RESTART SEARCH FROM COORDINATES IN RESTART.D
C IPR                    = PRINT CONTROL PARAMETER
C ILOOK                  = LOOKUP TABLE FOR DOMAIN AVERAGING
C ACOORD                 = STORAGE FOR SYMMETRY EQUIVALENT COORDINATES
C MICUT,MJCUT,PSQ,JYLM,  = DUMMY ARRAYS TO BE PASSED TO ROUTINE FUNCV
C BJ,YLM,QS,AT,XISTS,
C XIST,INBED,IEERG,AE,EE,
C NEE,NBEA,BENAME,IPR,XPL,
C YPL,NNN,AP,APP,YE,TSE,
C TSE2,TSEP,TSEP2,TSEPP,
C TSEPP2,TSEY2,WR,WB,IBP,
C ETH
C in Common:
C IFUNC                  =Total number of function evaluations
C
C P and XI are returned as the best point found, and the latest directiob set
C respectively
C
C
C =======================================================================
C
Cga      SUBROUTINE POWELL(P,XI,NLAY,NDIM,DISP,NNDIM,
      function POWELL(P,XI,NLAY,NDIM,DISP,NNDIM,
     & ADISP,FTOL2,ASTEP,VSTEP,ITMAX,ISTART,IPR,PTT,PT,XIT,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,
     & NDIML,DIREC,NLTIN,LPOINT,LRFLAG)
C
C MAXC is the maximum dimension of the parameter space where we perform
C optimazation of the R factor
C
      PARAMETER (MAXC=100)
C
      DIMENSION XI(NNDIM,NNDIM),DISP(NLAY,3),PT(NNDIM)
      DIMENSION P(NNDIM),XIT(NNDIM),LLFLAG(NLAY+NNST)
      DIMENSION DIREC(NLAY,2,2),NDIML(NLAY)
      DIMENSION D(MAXC),LSFLAG(NLAY)
      DIMENSION ADISP(NLAY,3),PTT(NNDIM)
      DIMENSION ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1)
      COMPLEX JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /RPL/DVOPT
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
C
500   FORMAT (/' DIRECTION SET SEARCH COMPLETE ')
501   FORMAT (' ======================= ')
502   FORMAT (' DIMENSIONALITY OF SEARCH=',I5)
503   FORMAT (' ISMOTH=',I5)
505   FORMAT (/' CONVERGENCE TOLERANCE ACHIEVED =',F7.6)
506   FORMAT (/' COORDINATES AT MINIMUM;  LSFLAG'/)
510   FORMAT (3F12.8,I4)
513   FORMAT (' NUMBER OF FUNCTIONAL EVALUATIONS =',I5)
514   FORMAT (/' OPTIMUM R-FACTOR =',1F7.4)
515   FORMAT (' OPTIMUM VALUE OF INNER POTENTIAL =',1F7.4)
516   FORMAT (/' NUMBER OF ITERATIONS =',1I3)
520   FORMAT (' SEARCH HAS EXCEEDED MAXIMUM NUMBER OF ITERATIONS')
!521   FORMAT (' OTHER VERTICES HAVE COORDINATES; ')
!522   FORMAT (' R-FACTOR =',1F7.4)
!523   FORMAT (I4,F7.4,70F7.4)

cjcm  del2 unitialized?
      del2 = 0.0
cjcm      write(*,*) 'powell: ******************************************'
      
      IF(LRFLAG.EQ.0) THEN
C?
      SCAL=VSTEP/ASTEP
      IF (ISTART.EQ.0) THEN
C
C First initialize the set of directions in XI
C
        DO 100 I=1,NNDIM
          DO 110 J=1,NNDIM
             IF (J.EQ.I) THEN
                XI(I,J)=1.0
             ELSE
                XI(I,J)=.0
             ENDIF
110       CONTINUE
100     CONTINUE 
C 
C Define initial point in configuration space. N1 takes care of coordinates
C which have to be moved together 
C
         JJDIM=0
         KK=0
         IF (NNST.GE.1) THEN
            DO 103 I=1,NNST
              IF (LLFLAG(I+NLAY).NE.0) THEN
                 KK=KK+1
                 P(KK)=0.
                 D(KK)=1.
              ENDIF
103         CONTINUE
         ENDIF
         N1=0
         DO 105 I=1,NLAY
C
C Is this layer to be included in the search?
C
            IF (LSFLAG(I).GT.N1) THEN
C
C yes. Then record its coordinates according to the effective
C dimensionality of the layer.(Notice that if the dimensionality
C is 2 then the Y coordinate read from the input actually
C corresponds to the coefficient of the unit vector specified
C by DIREC; but this will be important in SETCOR)
C
              N1=N1+1
              DO 115 J=1,NDIML(N1)
                 KK=KK+1
                 P(KK)=DISP(I,J)
C
C D() is not used in POWELL but in POWELL2.  
C
                 D(KK)=1.
115           CONTINUE
            ENDIF
105      CONTINUE
C
C set nonstructural parameters
C
      ELSE
C
C If restarting, upload coordinates from restart file and initialize
C bracketing triple
C
         CALL SRETV(FMIN,XI,D,P,LLFLAG,LSFLAG,
     &    NNDIM,JJDIM,NLAY,ADISP,NDIML,DIREC)
      ENDIF
      IFUNC=0
      iter=0
C
C Set coordinate matrix ADISP from temporary store P
C
C      write(0,*)'before first setcor2'
C      write(0,*) 'P=',P
C      write(0,*) 'DISP=',DISP
C      write(0,*) 'ADISP=',ADISP
      CALL SETCOR2(LSFLAG,NDIML,DIREC,P,DISP,ADISP,
     & NLAY,NNDIM)
C      write(0,*)'after first setcor2'
C      write(0,*) 'P=',P
C      write(0,*) 'DISP=',DISP
C      write(0,*) 'ADISP=',ADISP
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      FRET=FVAL

      IFUNC=IFUNC+1
      do 11 j=1,NNDIM   
        pt(j)=p(j)
11    continue
C
C reinitialize after NNDIM iterations the set of directions in XI
C This is quite important because it prevents that directions  become
C linearly dependent.  
C
1     iter=iter+1
      NN=MOD(ITER,NNDIM)
      IF(NN.EQ.0) THEN
        DO 102 I=1,NNDIM
          DO 112 J=1,NNDIM
             IF (J.EQ.I) THEN
                XI(I,J)=1.0
             ELSE
                XI(I,J)=.0
             ENDIF
112       CONTINUE
102     CONTINUE
      ENDIF
      fp=fret
      ibig=0
      del=0.
C
C loop over the different directions..... 
C
      do 13 i=1,NNDIM   
        do 12 j=1,NNDIM  
          xit(j)=xi(j,i)
12      continue
        AX=.0
        XX=ASTEP
C        BX=.02
C
C compute the minimum of the R factor along  direction I
C
       FPTT=FRET
       FMN=FRET
C      I2=MCLOCK()
      CALL LINMIN(FMN,AX,XX,BX,FRET,P,XIT,NLAY,NDIM,NNDIM,ADISP,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,DISP,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     &      NDIML,DIREC,NLTIN,LPOINT)
cjcm      write(*,*) 'powell: after linmin, i, fret', i, fret
C      I3=MCLOCK()
C      FTIME=FLOAT(I3-I2)
C      CPU=FTIME/100.0
C      WRITE(*,*) CPU,'linmin'
C
C and record the direction of maximum decrease of the function
C
      if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          del2=abs(fp-fret)
C          write(*,*) 'ibig=',i
          ibig=i
      endif

13    continue
      FRAC=2.*abs(fp-fret)/(abs(fp)+abs(fret))
C
C are we done?
C
cjcm
cjcm      write(*,*) 'powell: ', iter, fret, frac
      if(FRAC.LE.FTOL2)GOTO 1001
      if(iter.eq.itmax)GOTO 1002  
C
C check to see whether to replace the direction of maximum decrease in the
C direction set
C
      anorm=0.
      do 14 j=1,NNDIM
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
        anorm=anorm+xit(j)**2
14    continue
      anorm=sqrt(anorm)
      CALL SETCOR2(LSFLAG,NDIML,DIREC,PTT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      FPTT=FVAL
      IFUNC=IFUNC+1
C
C one reason not to replace it
C
      if(fptt.ge.fp)go to 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del2)**2-del2*(fp-fptt)**2
C
C another reason not to replace it
C
      if(t.ge.0.)go to 1
C
C the new direction is good and can be used, but first scale the vector
C
      do 156 j=1,NNDIM
        xit(j)=xit(j)/anorm
        xi(j,ibig)=xit(j)
156   continue
        AX=.0
        XX=ASTEP
      CALL LINMIN(FVAL,AX,XX,BX,FRET,P,XIT,NLAY,NDIM,NNDIM,ADISP,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,DISP,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     & NDIML,DIREC,NLTIN,LPOINT)
C      do 15 j=1,NNDIM
C        xi(j,ibig)=xit(j)
C15    continue
      go to 1
1002  CONTINUE
C
C Iteration count exceeded, so dump coordinates of the minimum
C to retsart file
C
      WRITE (2,520)
      CALL SETCOR2(LSFLAG,NDIML,DIREC,P,DISP,ADISP,
     & NLAY,NNDIM)
      CALL SDUMP(XI,D,FRET,ADISP,NNDIM,NDIM,NLAY)
      RETURN
C
C Convergence achieved in R factor value. 
C
 1001 CONTINUE
cjcm      
cjcm      write(*,*) 'powell: converged ', iter, fret, frac
cjcm      write(*,*) 'powell: ******************************************'
      CALL SETCOR2(LSFLAG,NDIML,DIREC,P,DISP,ADISP,
     & NLAY,NNDIM)

C  Recompute intensity for the minimum configuration and write the
C  IV curves for the best structure
C
      IFUNC1=IFUNC
      IFUNC=0
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
C      WRITE(*,*) 'rfacter for minimum ',FVAL
      CALL SDUMP(XI,D,FVAL,ADISP,NNDIM,NDIM,NLAY)
      WRITE(10,*) (EEAVT(J),J=1,NBE+1)
      WRITE(10,*) (EEAVE(J),J=1,NBE+1)
      WRITE (2,500)
      WRITE (2,501)
      CALL WRIV(AT,ETH,AE,EE,IEERG,NT0,NBED,IBK,WR,
     & ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,BENAME,
     & NBE,NST1,NST2)
      VOPT=VOPT+DVOPT+VV*27.21
      WRITE (2,505) FRAC
      WRITE (2,502) NNDIM
      WRITE (2,503) ISMOTH
      WRITE (2,506)
      DO 150 I=1,NLAY
         WRITE (2,510) (ADISP(I,J),J=1,3),LSFLAG(I)
150   CONTINUE
      WRITE (2,516) ITER
      WRITE (2,513) IFUNC1
      WRITE (2,514) FRET
      WRITE (2,515) VOPT
C?    if lrflag .ne. o
             ELSE
C?
C
C  Compute intensity for the original configuration and write the
C  IV curves for the original structure without any minimization.

C
      IFUNC1=IFUNC
      IFUNC=0
      CALL FUN2(FVAL,NLAY,DISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      CALL SDUMP(XI,D,FVAL,DISP,NNDIM,NDIM,NLAY)
      WRITE(10,*) (EEAVT(J),J=1,NBE+1)
      WRITE(10,*) (EEAVE(J),J=1,NBE+1)
      WRITE (2,500)
      WRITE (2,501)
      CALL WRIV(AT,ETH,AE,EE,IEERG,NT0,NBED,IBK,WR,
     & ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,BENAME,
     & NBE,NST1,NST2)
      VOPT=VOPT+DVOPT+VV*27.21
      WRITE (2,505) FRAC
      WRITE (2,502) NNDIM
      WRITE (2,503) ISMOTH
      WRITE (2,506)
      DO 151 I=1,NLAY
         WRITE (2,510) (DISP(I,J),J=1,3),LSFLAG(I)
151   CONTINUE
      WRITE (2,516) ITER
      WRITE (2,513) IFUNC1
      WRITE (2,514) FVAL
c     WRITE (2,514) FRET
      WRITE (2,515) VOPT

C      WRITE(0,*) 'rfacter = ',FVAL
      fret = fval 
             ENDIF
C?
Cga

      powell = fret

      RETURN
      END
C =======================================================================
C
C Subroutine POWELL2 searchs for the R-factor minimum by performing
C one dimensional minimization along a set of orthogonal direction which
C are rotated to match the principal directions at the minimum as the 
C search proceeds. 
C At exit,the vector P is set to the best point found, XI is the 
C the current direction set, D contains the eigenvalues giving the curvatures
C at the minimum,FRET is the returned function value at P, and
C ITER is the number of iteration taken. The routine LINMIN is used 
C
C References 'NUMERICAL RECIPES' W.H.Press,B.P.Flannery, et al. Cambridge
C             University Press
C 
C Modifications by BARBIERI
C
C Input Parameters;
C =================
C
C NLAY                   = NUMBER OF LAYERS IN COMPOSITE LAYER
C NDIM                   = DIMENSIONALITY OF SEARCH (1=X AXIS ONLY
C                                                    3=X,Y,Z AXES)
C P                      = ARRAY CONTAINING THE INITIAL STARTING POINT
C                          IN THE SEARCH
C XI                     = ARRAY STORING THE INITIAL DIRECTION SET
C DISP                   = COORDINATES INPUT BY USER
C NNDIM                  = TOTAL NUMBER OF DIMENSIONS IN SEARCH (=NLAY*NDIM)
C ADISP                  = GEOMETRY OF CURRENT POINT IN SAME FORMAT AS DISP
C                          (USED AS INPUT TO FUNCV)
C DVOPT                  = SHIFT IN THE INNER POTENTIAL FOR THE STARTING 
C                          CONFIGURATION (in COMMON /RPL )
C LLFLAG                  = INDICATES WHETHER THE LAYER (OR NON STRUCTURAL
C                          PARAMETER) COORDINATES HAVE TO BE VARIED 
C                          IN THE SEARCH
C LSFLAG                 =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C                          NSYM) IN THE COMPOSITE LAYER. LLFLAG(i)=LLFLAG(j)
C                          INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C                          AS EQUIVALENT IN THE SEARCH
C NDIML                  =  ARRAY GIVING THE EFFECTIVE DIMENSIONALITY OF ATOM
C                           j (ACTUALLY TO BE USED IN CONJUNCTION WITH LLFLAG)
C DIREC                  =  SET OF DIRECTIONS
C VOPT                   = INNER POTENTIAL OF CURRENT INPUT (used as input
C                          in FUNCV), in COMMON /NSTR
C FTOL2                  = CONVERGENCE CRITERIA FOR POWELL SEARCH
C ITMAX                  = MAXIMUM NUMBER OF ITERATIONS TO BE PERFORMED BY
C                          POWELL
C ISTART                 = 0 USE PARAMETERS FROM TLEED4.I
C                          1 RESTART SEARCH FROM COORDINATES IN RESTART.D
C IPR                    = PRINT CONTROL PARAMETER
C ILOOK                  = LOOKUP TABLE FOR DOMAIN AVERAGING
C ACOORD                 = STORAGE FOR SYMMETRY EQUIVALENT COORDINATES
C MICUT,MJCUT,PSQ,JYLM,  = DUMMY ARRAYS TO BE PASSED TO ROUTINE FUNCV
C BJ,YLM,QS,AT,XISTS,
C XIST,INBED,IEERG,AE,EE,
C NEE,NBEA,BENAME,IPR,XPL,
C YPL,NNN,AP,APP,YE,TSE,
C TSE2,TSEP,TSEP2,TSEPP,
C TSEPP2,TSEY2,WR,WB,IBP,
C ETH
C in Common:
C IFUNC                  =Total number of function evaluations
C
C P and XI are returned as the best point found, and the latest directiob set
C respectively
C
C
C =======================================================================
C
Cga changed into function
Cga      SUBROUTINE POWELL2(P,XI,NLAY,NDIM,DISP,NNDIM,
      function POWELL2(P,XI,NLAY,NDIM,DISP,NNDIM,
     & ADISP,FTOL2,ASTEP,VSTEP,ITMAX,ISTART,IPR,PTT,XIT,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,
     & NDIML,DIREC,NLTIN,LPOINT)
C
C MAXC is the maximum dimension of the parameter space where we perform
C optimazation of the R factor
C
      PARAMETER (MAXC=100 )
C
      DIMENSION A(MAXC,MAXC),XTEM(MAXC,MAXC)
      DIMENSION D(MAXC),XINEW(MAXC,MAXC)
      DIMENSION DIREC(NLAY,2,2),NDIML(NLAY),LSFLAG(NLAY)
      DIMENSION XI(NNDIM,NNDIM),DISP(NLAY,3)
      DIMENSION P(NNDIM),XIT(NNDIM),LLFLAG(NLAY+NNST)
      DIMENSION XX1(MAXC),XX2(2,MAXC),XX3(2,MAXC)
      DIMENSION ADISP(NLAY,3),PTT(NNDIM)
      DIMENSION ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1)
      COMPLEX JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /RPL/DVOPT
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
C
500   FORMAT (/' DIRECTION SET SEARCH COMPLETE ')
501   FORMAT (' ======================= ')
502   FORMAT (' DIMENSIONALITY OF SEARCH=',I5)
503   FORMAT (' ISMOTH=',I5)
505   FORMAT (/' CONVERGENCE TOLERANCE ACHIEVED =',F7.6)
506   FORMAT (/' COORDINATES AT MINIMUM;  LSFLAG'/)
507   FORMAT (/' PRICIPAL DIRECTIONS AT MINIMUM;'/)
508   FORMAT (/' EIGENVALUES AT MINIMUM;'/)
510   FORMAT (3F9.4,I4)
511   FORMAT (10F9.4)
513   FORMAT (' NUMBER OF FUNCTIONAL EVALUATIONS =',I5)
514   FORMAT (/' OPTIMUM R-FACTOR =',1F7.4)
515   FORMAT (' OPTIMUM VALUE OF INNER POTENTIAL =',1F7.4)
516   FORMAT (/' NUMBER OF ITERATIONS =',1I3)
520   FORMAT (' SEARCH HAS EXCEEDED MAXIMUM NUMBER OF ITERATIONS')
!521   FORMAT (' OTHER VERTICES HAVE COORDINATES; ')
!522   FORMAT (' R-FACTOR =',1F7.4)
!523   FORMAT (I4,F7.4,70F7.4)
C
      SCAL=VSTEP/ASTEP
      IF (ISTART.EQ.0) THEN
C
C First initialize the set of directions in XI 
C
        DO 100 I=1,NNDIM
            DO 110 J=1,NNDIM
               IF (J.EQ.I) THEN
                  XI(I,J)=1.0
               ELSE
                  XI(I,J)=.0
               ENDIF
110         CONTINUE
100     CONTINUE 
C 
C Define initial point in configuration space 
C
         JJDIM=0
         KK=0
C
C set nonstructural parameters
C
         IF (NNST.GE.1) THEN
           DO 103 I=NLAY+1,NLAY+NNST
             IF (LLFLAG(I).NE.0) THEN
               KK=KK+1
               P(KK)=0.
               D(KK)=1.
             ENDIF
103        CONTINUE
         ENDIF
         N1=0
         DO 105 I=1,NLAY
C
C Is this layer to be included in the search?
C
            IF (LSFLAG(I).GT.N1) THEN
C
C yes. 
C
              N1=N1+1
              DO 115 J=1,NDIML(N1)
                 KK=KK+1
                 P(KK)=DISP(I,J)
C
C  initialize eigenvalues D() 
C
                 D(KK)=1.
115           CONTINUE
            ENDIF
105      CONTINUE
      ELSE
C
C If restarting, upload coordinates,directions and tentative eigenvalues
C from restart file 
C
         CALL SRETV(FMIN,XI,D,P,LLFLAG,LSFLAG,
     &    NNDIM,JJDIM,NLAY,ADISP,NDIML,DIREC)
      ENDIF
      IFUNC=0
      iter=0
C
C Set coordinate matrix ADISP from temporary store P
C
      CALL SETCOR2(LSFLAG,NDIML,DIREC,P,DISP,ADISP,
     & NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      fp=FVAL
      FRET=FVAL
      IFUNC=IFUNC+1
C
C Now set the scale for computing gradient and Hessian of the Rfactor
C at the minimum. If we consider the variation of R along the principal
C direction I corresponding to an eigenvalue LAMBDAI, we will probe tha 
C function R a distance X away from the minimum P where  
C LAMBDAI*X^2=.01*R(P)
C
      eps=sqrt(.01*fp)
1     iter=iter+1
      IF (iter.gt.1) THEN
C
C All this work to compute A, now we use also b to extrapolate the
C position of the minimum  and accept it if is lower than current
C minimum. f=f(p) -bx +xAx/2 gives a minimum at x=A^(-1)b but one has
C to be careful about the coordinate system. Now we define the
C extrapolated minimum
C
         do 19 j=1,NNDIM
             ptt(j)=.0
             do 191 m=1,NNDIM
             do 192 n=1,NNDIM
                 ptt(j)=ptt(j)+xi(j,m)*XINEW(n,m)*XX1(n)/D(m)
192          continue
191          continue 
             ptt(j)=p(j)+ptt(j)
19       continue
C
C and check if it gives a lower value for function
C 
         CALL SETCOR2(LSFLAG,NDIML,DIREC,PTT,DISP,ADISP,
     & NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
          IF (FVAL.lt.fpl) THEN
             do 193 i=1,NNDIM
                p(i)=ptt(i)
193          continue 
             fp=fval
             FRET=fp
          ELSE
             fp=fpl
          ENDIF 
      ENDIF
C
C now go on with the search and
C loop over the different directions..... 
C
      do 13 i=1,NNDIM   
        do 12 j=1,NNDIM  
          xit(j)=xi(j,i)
12      continue
        AX=.0
        XX=.01
        BX=.02
C
C compute the minimum of the R factor along  direction I
C
C      I2=MCLOCK()
      FPTT=FRET
      FMN=FRET
      CALL LINMIN(FMN,AX,XX,BX,FRET,P,XIT,NLAY,NDIM,NNDIM,ADISP,
     & ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,NERG,
     & AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,DISP,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,IPR,
     & NDIML,DIREC,NLTIN,LPOINT)
C
13    continue
      fpl=FRET
C
C Fp is the value of the function at p at the beginning of the loop,
C fpl is the last computed value of the function at p
C
      FRAC=2.*abs(fp-fpl)/(abs(fp)+abs(fpl))
C
C Now estimate the Hessian A at P
C The second derivatives are estimated by probing the function 
C eps away from P, but taking into account that different directions
C have different scales associated with them (typically the scale of
C x displacements are smaller than that of y,z displacements, the different
C scale for the inner potential is  partially taken care of by SCAL)
C Both the Hessian and the gradient are computed to order O(eps^4)
      do 14 i=1,NNDIM
         do 141 kk=1,2
            eps=-eps
            epst=eps
            do 142 kkk=1,2
              do 145 ii=1,NNDIM
                 ptt(ii)=p(ii)+epst*xi(ii,i)/sqrt(ABS(D(i)))
145           continue
      CALL SETCOR2(LSFLAG,NDIML,DIREC,PTT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
              IFUNC=IFUNC+1
              IF (kk.EQ.1) THEN
                XX2(kkk,i)=FVAL
              ELSE
                XX3(kkk,i)=FVAL
              ENDIF
              epst=eps*float(kkk+1)
142        continue
141      continue
C
C compute b to order O(eps^4)
C
         fac1=2./3.
         fac2=1./12.
         XX1(i)=sqrt(ABS(D(i)))/eps*(fac1*(XX2(1,i)-XX3(1,i))-
     &               fac2*(XX2(2,i)-XX3(2,i)))
14    continue
      do 15 i=1,NNDIM
         do 16 j=1,NNDIM
            do 155 ii=1,NNDIM
               ptt(ii)=p(ii)+eps*(xi(ii,i)/sqrt(ABS(D(i)))
     &                 -xi(ii,j)/sqrt(ABS(D(j))))
155         continue
      CALL SETCOR2(LSFLAG,NDIML,DIREC,PTT,DISP,ADISP,
     & NLAY,NNDIM)
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
            A(i,j)=FVAL-XX3(1,i)-XX2(1,j)+fpl
            IFUNC=IFUNC+1
16       continue
15    continue
C
C evaluate the Hessian A to order O(eps^4)
C      
      do 17 i=1,NNDIM
         do 18 j=i,NNDIM
            A(i,j)=-(A(i,j)+A(j,i))*sqrt(ABS(D(i)*D(j)))/(2.*eps**2)
            A(j,i)=A(i,j)
18       continue
17    continue
C
C Now diagonalize A.  The rotation matrix needed to reset the direction set 
C to the principal directions of A is XINEW (the columns are the normalized 
C eigenvectors of A)
C
      CALL JACOBI(A,NNDIM,MAXC,D,XINEW,NROT)
C and fix possible inconsistencies
      do 194 i=1,NNDIM
         if (D(i).lt.0.) D(i)=1.
194   continue
C
C rotate old directions to approximate principal directions
C (the Hessian matrix was expressed in the basis
C specified by the old direction set, hence XINEW is also giving the
C coordinates of the new direction set in that basis)
C
      do 111 i=1,NNDIM
         do 112 j=1,NNDIM
            XTEM(i,j)=0.
            do 113 k=1,NNDIM
              XTEM(i,j)=XTEM(i,j)+XINEW(k,j)*XI(i,k)
113         continue
112      continue
111   continue
      do 114 i=1,NNDIM
         do 116 j=1,NNDIM
            XI(i,j)=XTEM(i,j)
116      continue
114   continue
C
C are we done?
C
      if(FRAC.LE.FTOL2)GOTO 1001
      if(iter.eq.itmax)GOTO 1002  
      go to 1
1002  CONTINUE
C
C Iteration count exceeded, so dump coordinates of the minimum
C to retsart file
C
      WRITE (2,520)
      CALL SETCOR2(LSFLAG,NDIML,DIREC,P,DISP,ADISP,
     & NLAY,NNDIM)
      CALL SDUMP(XI,D,FPL,ADISP,NNDIM,NDIM,NLAY)
      RETURN
C
C Convergence achieved in R factor value, dump minimum configuration,
C principal direction, eigenvalues etc to RESTART.D file which can then 
C be used by GRID. 
C
1001  CONTINUE
      CALL SETCOR2(LSFLAG,NDIML,DIREC,P,DISP,ADISP,
     & NLAY,NNDIM)
C
C  Recompute intensity for the minimum configuration and write the
C  IV curves for the best structure
C
      IFUNC1=IFUNC
      IFUNC=0
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      CALL SDUMP(XI,D,FVAL,ADISP,NNDIM,NDIM,NLAY)
      WRITE(10,*) (EEAVT(J),J=1,NBE+1)
      WRITE(10,*) (EEAVE(J),J=1,NBE+1)
      WRITE (2,500)
      WRITE (2,501)
      CALL WRIV(AT,ETH,AE,EE,IEERG,NT0,NBED,IBK,WR,
     & ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,BENAME,
     & NBE,NST1,NST2)
      VOPT=VOPT+DVOPT+VV*27.21
      WRITE (2,505) FRAC
      WRITE (2,502) NNDIM
      WRITE (2,503) ISMOTH
      WRITE (2,506)
      DO 150 I=1,NLAY
         WRITE (2,510) (ADISP(I,J),J=1,3),LSFLAG(I)
150   CONTINUE
      WRITE (2,507)
      DO 1513 I=1,NNDIM
         WRITE (2,511) (XI(J,I),J=1,NNDIM)
1513   CONTINUE
      WRITE (2,508)
      WRITE (2,511) (D(J),J=1,NNDIM)
      WRITE (2,516) ITER
      WRITE (2,513) IFUNC1
      WRITE (2,514) FRET
      WRITE (2,515) VOPT
Cga
      powell2 = fret
      RETURN
      END
C ======================================================================
C
C Subroutine QUADEQ fits a quadratic equation to three evenly spaced
C points in parameter space, seperated by a distance STEP.
C
C Input Parameters;
C =================
C
C VAL      =  FUNCTIONALLY VALUE AT THE THREE POINTS TO BE FITTED
C A,B,C    =  QUADRATIC COEFFICIENTS
C STEP     =  STEP SIZE BETWEEN POINTS
C
C Author: WANDER
C
C ======================================================================
C
      SUBROUTINE QUADEQ(VAL,A,B,C,STEP,XVAL)
C
      DIMENSION VAL(3)
C
C Calculate Values of A and B in Y=A*X^2+B*X+C
C
      FAC=STEP*STEP
      A=(VAL(1)-VAL(2))/FAC
      A=A-(VAL(1)-VAL(3))/(2*FAC)
      B=(VAL(2)-VAL(1))/STEP-A*(2*XVAL-STEP)
      C=VAL(2)-A*XVAL*XVAL-B*XVAL
      RETURN
      END
C ============================================================================
C
C Subroutine READE loads the experimental data from the file on UNIT=11.
C
C Input Parameters;
C -----------------
C
C AE        =  Experimental intensities
C EE        =  Energy range
C NBED      =  Number of experimental beams
C NEE       =  Number of energy points in each beam
C NBEA      =  Beam averaging information
C BENAME    =  Label of beam in reciprocal space
C IPR       =  Print control parameter
C NBMAX     =  Number of symmetry inequivalent beams
C
C =============================================================================
C
      SUBROUTINE READE(AE,EE,NBED,NEE,NBEA,BENAME,IPR,IEERG)
C
      DIMENSION AE(NBED,IEERG),EE(NBED,IEERG),NEE(NBED),NBEA(NBED)
      DIMENSION BENAME(5,NBED)
      CHARACTER(LEN=20) FMT
      COMMON /REXP/EEINCR
C
30    FORMAT (5(25I3,/))
40    FORMAT (I4,' EXP. BEAMS TO BE READ IN')
50    FORMAT ('AVERAGING SCHEME OF EXP. BEAMS',5(25I3,/))
60    FORMAT (20A4)
10    FORMAT (20A4)
!70    FORMAT ('EXP. BEAM ',1I3,' (',5A4,')')
35    FORMAT (1I3,1E13.4)
36    FORMAT (5F7.4)
80    FORMAT ('  EXP. ENERGIES AND INTENSITIES',/,50(5(1F7.2,1E13.4,
     & 3X),/))
C
C  READ IN EXP. BEAM AVERAGING INFORMATION. IF NBEA(I)=NBEA(J), EXP.
C  BEAMS I AND J WILL BE AVERAGED TOGETHER. THE RELATION NBEA(J).GT.
C  NBEA(I) IF J.GT.I MUST HOLD. IF NBEA(I)=0, EXP. BEAM I WILL BE SKIPPED
C  LATER ON
C
      READ (11,30) (NBEA(I),I=1,NBED)
      IF(IPR.GT.0)THEN
         WRITE(1,40)NBED
         WRITE(1,50)(NBEA(I),I=1,NBED)
      ENDIF
C
C  READ INPUT FORMAT OF EXP. INTENSITIES
C
C      READ (11,60) (FMT(I),I=1,20)
      READ (11,60) FMT
      READ (11,36) EEINCR
C
      DO 90 IB=1,NBED
         READ (11,10) (BENAME(I,IB),I=1,5)
         IF(IPR.GT.0)THEN
            WRITE(1,10)(BENAME(I,IB),I=1,5)
         ENDIF
C
C  READ IN NO. OF DATA POINTS TO BE INPUT FOR THIS BEAM AND CONSTANT
C  CORRECTION FACTOR FOR INTENSITIES. THIS FACTOR IS MEANT TO ALLOW
C  NORMALIZATION TO INTENSITIES OF THE ORDER OF 1 (NOT NECESSARY, BUT
C  SAFE), AND TO MATCH UP CURVES TO BE AVERAGED TOGETHER WHEN THEIR
C  ENERGY RANGES DIFFER (TO AVOID DISCONTINUITIES AT ENERGIES WHERE
C  THE NUMBER OF CURVES AVERAGED TOGETHER CHANGES)
C
         READ (11,35) NEE(IB),FAC
         N=NEE(IB)
C
C  READ (AND MAYBE PRINT) EXP. INTENSITIES
C

cjcm replace FMT with 99 and 2F12.3 which is what is in the exp.d file 9/29/2016
         READ (11,99) (EE(IB,IE),AE(IB,IE),IE=1,N)
 99      FORMAT(2F12.3)
         IF (IPR.GE.0) THEN
            WRITE (1,80) (EE(IB,IE),AE(IB,IE),IE=1,N)
         ENDIF
         DO 86 IE=1,N
            AE(IB,IE)=AE(IB,IE)*FAC
86       CONTINUE
90    CONTINUE
      RETURN
      END
C =========================================================================
C
C Subroutine RENORM fits the function A*EXP(-ALPH*(E+10)) to given
C data, by making it pass through the centres of gravity of the first    
C and second halves of the data, and divides the data by that function.
C
C =========================================================================
C
      SUBROUTINE RENORM(A,NBD,NB,NE,E,IPR,Y,IEERG)
C
      DIMENSION A(NBD,IEERG),NE(NBD),E(NBD,IEERG),Y(IEERG)
C
15    FORMAT (' RENORM. IN BEAM',1I4,' - AA =',1E15.4,', ALPH =',
     & 1E15.4)
60    FORMAT ('ENERG. AND INTENS. AFTER RENORMALIZATION IN BEAM',1I3,
     & /,50(5(1F7.2,1E13.4,3X),/))
C
      DO 30 IB=1,NB
         NT=NE(IB)
         IF (NT.GT.1) THEN
            N1=NT/2
            N2=N1+1
            IF (2*N1.NE.NT) N2=N2+1
            E1=0.5*(E(IB,N1)+E(IB,1))
            E2=0.5*(E(IB,NT)+E(IB,N2))
C
C  DO NOT RENORMALIZE IF ONLY SHORT ENERGY RANGE AVAILABLE
C
            IF (2.*(E2-E1).GE.100.) THEN
               DO 12 IE=1,NT
                  Y(IE)=A(IB,IE)
12             CONTINUE
               CALL SUM(Y,1,1,1.,1,N1,S1,IEERG)
               CALL SUM(Y,1,1,1.,N2,NT,S2,IEERG)
               S1=S1/FLOAT(N1-1)
               S2=S2/FLOAT(N1-1)
               ALPH=(ALOG(S1/S2))/(E2-E1)
C
C  PREVENT NEGATIVE ALPH
C
               IF (ALPH.LT.0.) ALPH=0.
               AA=S1*EXP(ALPH*(E1+10.))
               IF (IPR.GE.2) WRITE (1,15) IB,AA,ALPH
               DO 20 IE=1,NT
                  A(IB,IE)=A(IB,IE)/AA*EXP(ALPH*(E(IB,IE)+10.))
20             CONTINUE
            ENDIF
         ENDIF
30    CONTINUE
      IF (IPR.GE.2) THEN
         DO 50 IB=1,NB
            N=NE(IB)
            IF (N.NE.0) WRITE (1,60) IB,(E(IB,IE),A(IB,IE),IE=1,N)
50       CONTINUE
      ENDIF
      RETURN
      END
C ========================================================================
C
C Subroutine RETRV uploads the plane wave amplitudes and tensor elements
C from the transfer file.
C
C Input Parameters
C ================
C
C NCUT        =  DIMENSION OF QS
C ICUT        =  NUMBER OF ELEMENTS IN Q
C NLAY        =  NUMBER OF LAYERS IN COMPOSITE LAYER
C NT0         =  NUMBER OF EXIT BEAMS
C XISTS       =  PLANE WAVE AMPLITUDES
C QS          =  TENSOR
C THETA,FI    =  ANGLE OF INCIDENCE
C PSQ,PQFEX   =  INDICES OF EXIT BEAMS IN RECIPROCAL SPACE
C RAR1,RAR2   =  RECIPROCAL LATTICE VECTORS
C
C COMMON BLOCK
C ============
C
C EI,EF,DE    =  Energy range (block ENY)
C NERG        =  Number of energy points (block ENY)
C NSYM        =  Symmetry code of surface (block ENY)
C NDOM        =  Number of domains on surface (block ENY)
C VV          =  Real part of inner potential (block ENY)
C
C AUTHOR: WANDER
C
C ========================================================================
C
      SUBROUTINE RETR2(NCUT,ICUT,NLAY,NT0,XISTS,QS,THETA,FI,PSQ,PQFEX,
     & RAR1,RAR2)
C
      DIMENSION PSQ(2,NT0),PQFEX(2,NT0),RAR1(2),RAR2(2)
      DIMENSION NEM(50)
      COMPLEX QS(NCUT),XISTS(NT0,NERG)
C
      COMMON /ENY/EI,EF,DE,NERG,NSYM,NDOM,VV,VPIS
C
!1000  FORMAT (3F7.4)
!1010  FORMAT (500(2E13.5,/))
1020  FORMAT (' **** ERROR: ENERGIES FROM TLEED5 AND SHORT.T DO NOT 
     & AGREE ')
C
C Set offset to zero
C
      IOFF=0
      REWIND (22)
C
C Start loop over energies
      write (1,9010)
9010  format (' read energies:')
C
      DO 100 I=1,NERG
         ENERG=(EI+FLOAT(I-1)*DE)/27.21+VV
C         READ (22,1000) E
         READ (22) E
cjcm
         write (1,9011) I,E, ENERG   
9011     format (' i,e,energy=',i5,f10.3,f10.3)
C
C Check that current enrgy point agrees with the transfer file energy point
C
         IF (ABS(E-ENERG).GT.0.0001) GOTO 1021
C
C Retrieve plane wave amplitudes for each exit beam
C
C         READ (22,1010) (XISTS(II,I),II=1,NT0)
         READ (22) (XISTS(II,I),II=1,NT0)
C
C Set parallel momentum components
C
         AK=SQRT(2.0*(E-VV))*SIN(THETA)
         AK2=AK*COS(FI)
         AK3=AK*SIN(FI)
C
C Set up indices of exit beams in units of reciprocal lattice vectors
C and check how many emerging beams at this energy
C
         NBEM=0
         DO 110 J=1,NT0
            DO 120 K=1,2
               PSQ(K,J)=PQFEX(1,J)*RAR1(K)+PQFEX(2,J)*RAR2(K)
120         CONTINUE
            AK2M=-AK2-PSQ(1,J)
            AK3M=-AK3-PSQ(2,J)
            A=2.0*(E-VV)-AK2M*AK2M-AK3M*AK3M
            NEM(J)=0
C
C Does this beam emerge?
C
            IF (A.GT.0) THEN
               NBEM=NBEM+1
               NEM(J)=1
            ENDIF
110      CONTINUE
        
C
C retrieve the tensor for each emerging beam and for each layer 
C
         IOF2=0
         DO 1100 J=1,NT0
            IF(NEM(J).EQ.1) THEN
               IIOF=IOFF
               DO 130 K=1,NLAY
C                  READ (22,1010) (QS(IOFF+IJ),IJ=1,ICUT)
C                  READ (22) (QS(IOFF+IJ),IJ=1,ICUT)
C                  READ (22,1010) (QS(IIOF+IJ),IJ=1,ICUT)
                  READ (22) (QS(IIOF+IJ),IJ=1,ICUT)
C
C Increment offset
C
                  IIOF=IIOF+NBEM*ICUT
                  IOF2=IOF2+ICUT
130            CONTINUE
               IOFF=IOFF+ICUT
            ENDIF
1100     CONTINUE
         IOFF=IOFF+IOF2-NBEM*ICUT
C
C Next energy point
C
100   CONTINUE
      RETURN
1021  WRITE (1,1020)
      STOP
      END
C ==========================================================================
C 
C Subroutine RFAC generates the R-factor for the current structure. 
C
C Author: Modified version of VAN HOVE'S R-factor program. Modifications
C by WANDER and BARBIERI.
C
C R1,R2, etc are the Rfactors corresponding to each beam.
C AR(11) are the averaged R factors.
C ===========================================================================
C
      SUBROUTINE RFAC(AT,ETH,INBED,IEERG,AE,EE,NEE,IPR,XPL,
     & YPL,AEP,AEPP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,
     & IBP,NT0,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,
     & RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,BARAV,BV0,YT)
C
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG),AE(INBED,IEERG)
      DIMENSION EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION XPL(IEERG),AEP(INBED,IEERG)
      DIMENSION AEPP(INBED,IEERG),YE(INBED,IEERG),TSE(INBED),TSE2(INBED)
      DIMENSION TSEP(INBED),TSEP2(INBED),TSEPP(INBED),TSEPP2(INBED)
      DIMENSION WR(10),WB(NT0),TSEY2(INBED)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
C++++
      DIMENSION IBP(NT0)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA,FI
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /ENY/EI,EF,DE,NERG,NSYM,NDOM,VINER,VIMAG
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /RPL/DVOPT
      COMMON /REXP/EEINCR
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /WIV2/PERSH,NIV,NSE1(30),NSE2(30)
C
C SET AVERAGE R-FACTOR TO A LARGE NUMBER
C
      BARAV=100.0
C
C PERFORM DOMAIN-AVERAGING OF THEORETCIAL DATA
C
      CALL RINTAV(AT,NT0,IBP,IEERG,NERG,NBMAX,IBK)
      DO 100 I=1,NBMAX
         NET(I)=NERG
100   CONTINUE
C
C SHIFT THE ENERGY ACCORDING TO THE SHIFT OF THE INNER POTENTIAL
C DETERMINED BY THE SEARCH ALGORITHM
C
cjcm      write(*,*) 'RFAC: shifting energy, check for <0 ??'
      DO 101 I=1,NBMAX
        DO 102 J=1,NERG
           ETH(I,J)=ETH(I,J)-BV0-DVOPT
cjcm check for negative energy and write out warning
           if(ETH(I,J) .LT. 0.0) THEN
              write(*,*) 'RFAC: ***** WARNING: ETH < 0 *****'
              write(*,*) 'RFAC: i,j, eth(i,j), bv0, dvopt:',
     &         I, J, ETH(I,J), BV0, DVOPT
           endif
102     CONTINUE
101   CONTINUE
C
C Now interpolate on the same grid used when smoothing the experimental 
C data (EEINCR)
C and smooth as done with the experimental data
C
         IF (ISMOTH.NE.0) THEN
cjcm         write(*,*) 'Calling INTPOL: ISMOTH= ', ISMOTH
cjcm         write(*,*) 'INTPOL: before call', NET
            CALL INTPOL(AT,ATP,ATPP,NT0,NET,NBMAX,ETH,EEINCR,
     &      IPR,XPL,YPL,IEERG)
cjcm         write(*,*) 'INTPOL: after call, before smooth', NET
            DO 10 I=1,ISMOTH
               CALL SMOOTH(AT,ETH,NT0,NBMAX,NET,IPR,IEERG)
10          CONTINUE
cjcm         write(*,*) 'INTPOL: after smooth', NET
         ENDIF
C
C INTERPOLATE THEORY ONTO WORKING GRID
C
c         write(*,*) 'Calling INTPOL: INTP THEORY ONTO WORKING GRID '
c         write(*,*) 'INTPOL:', NET
      CALL INTPOL(AT,ATP,ATPP,NT0,NET,NBMAX,ETH,EINCR,
     & IPR,XPL,YPL,IEERG)
C
C  RENORMALIZE TH. BY DIVIDING BY FITTED A*EXP(-ALPH*(E+10))
C
      IF (IREN.EQ.1) THEN
         CALL RENORM(AT,NT0,NBMAX,NET,ETH,IPR,YPL,IEERG)
      ENDIF
C
C  DETERMINE ENERGY INTERVAL COMMON TO EXP. AND THEORY. THIS INTERVAL,
C  OF LENGTH EET, IS BOUNDED BY THE GRID POINTS (NE1,NE2) AND (NT1,NT2)
C  FOR EXP. AND THEORY, RESP.
C
      DO 129 IBE=1,NBE
            IBT=IBK(IBE)
            NST1(IBE)=0
            NST2(IBE)=0
            NSE1(IBE)=0
            NSE2(IBE)=0
            CALL COMNEI(EE,NBED,NEE,ETH,NT0,NET,IBE,IBT,VO,EINCR,
     &      NE1,NE2,NT1,NT2,EET(IBE),IEERG)
            NST1(IBE)=NT1
            NST2(IBE)=NT2
            NSE1(IBE)=NE1
            NSE2(IBE)=NE2
129   CONTINUE
C
C  ESTIMATE THE AVERAGE INTENSITY OF THE PEAKS IN EACH TH. BEAM
C  SO TO SHIFT INTENSITIES WHEN COMPUTING THE PENDRY R FACTOR
C  THIS IS DONE ONLY AT THE BEGINNING OF THE CALCULATION, AND
C  AT THE END
C
         IF(IFUNC.EQ.0) THEN
C           CALL AVEINT(AT,ATP,NT0,NET,NBMAX,ETH,IEERG,EEAVT)
           CALL AVEINT(AT,ATP,NT0,NST1,NST2,NBMAX,
     &       IEERG,EEAVT,1)
         ENDIF
C
C  1ST AND 2ND DERIVATIVES  are now computed by INTPOL
C
C  PRODUCE PENDRY Y FUNCTION FOR TH. DATA
C
      CALL YPEND(AT,ATP,NT0,NBMAX,NET,ETH,YT,
     %  VPIS,IPR,IEERG,EEAVT)
C
C  PRODUCE SOME INTEGRALS OVER THEOR. DATA
C
      DO 20 IB=1,NBMAX
         IE2=NET(IB)
         IF (IE2.NE.0) THEN
            CALL VARSUM(AT,AT,AT,AT,NT0,1,IB,1,1,IE2,0,EINCR,0.,0.,
     &       1,TST(IB),YPL,IEERG)
            CALL VARSUM(YT,AT,AT,AT,NT0,1,IB,1,1,IE2,0,EINCR,0.,0.,
     &       2,TSTY2(IB),YPL,IEERG)
         ENDIF
20    CONTINUE
C
C      VO=BV0
         AROS=0.0
         AR1=0.0
         AR2=0.0
         ARP1=0.0
         ARP2=0.0
         ARPP1=0.0
         ARPP2=0.0
         ARRZJ=0.0
         ARMZJ=0.0
         ARPE=0.0
         ARAV=0.0
         ERANG=0.0
C         DO 62 IB=1,NT0
C62       CONTINUE
C
C START LOOP OVER EXPERIMENTAL BEAMS
C 
         DO 130 IBE=1,NBE
            RAV(IBE)=-1.0
C
C IBT INDICATES THE THEORETICAL BEAM CORRESPONDING TO THIS EXPERIMENTAL BEAM
C
            IBT=IBK(IBE)
            IF (IBT.NE.0.AND.NEE(IBE).NE.0.AND.NET(IBT).NE.0) THEN
               ROS(IBE)=0.0
               R1(IBE)=0.0
               R2(IBE)=0.0
               RP1(IBE)=0.0
               RP2(IBE)=0.0
               RPP1(IBE)=0.0
               RPP2(IBE)=0.0
               RRZJ(IBE)=0.0
               RMZJ(IBE)=0.0
               RPE(IBE)=0.0
               RAV(IBE)=0.0
C
C  FIX ENERGY INTERVAL COMMON TO EXP. AND THEORY. THIS INTERVAL,
C  OF LENGTH EET, IS BOUNDED BY THE GRID POINTS (NE1,NE2) AND (NT1,NT2)
C  FOR EXP. AND THEORY, RESP.
C
               VO=0.
               NT1=NST1(IBE)
               NT2=NST2(IBE)
               NE1=NSE1(IBE)
               NE2=NSE2(IBE)
C               CALL COMNEI(EE,NBED,NEE,ETH,NT0,NET,IBE,IBT,VO,EINCR,
C    &          NE1,NE2,NT1,NT2,EET(IBE),IEERG)
               IF (NT2.GT.NT1) THEN
C
C  FIND EPS = MAX(ABS(DERIV. OF EXP.)) AND SAME FOR THEORY, FOR ZANAZZI-JONA 
C  R-FACTOR
C
                  CALL EPSZJ(AEP,NBED,IBE,NE1,NE2,EPSE,IEERG)
                  CALL EPSZJ(ATP,NT0,IBT,NT1,NT2,EPST,IEERG)
C
C  IF ANY IV-CURVE IS TRUNCATED, THE INTEGRALS PERFORMED BEFORE SHOULD
C  BE REDUCED ACCORDINGLY
C
                  NE=NEE(IBE)
                  NT=NET(IBT)
                  SS=0.0
                  SS2=0.0
                  SSP=0.0
                  SSP2=0.0
                  SSPP=0.0
                  SSPP2=0.0
                  SSY2=0.0
                  SU=0.0
                  SU2=0.0
                  SUP=0.0
                  SUP2=0.0
                  SUPP=0.0
                  SUPP2=0.0
                  SUY2=0.0
                  IF (NE1.NE.1) THEN
                     CALL VARSUM(AE,AE,AE,AE,NBED,1,IBE,1,1,NE1,0,
     &                EINCR,0.,0.,1,SS,YPL,IEERG)
                     IF (WR(3).GE.1.E-6) CALL VARSUM(AE,AE,AE,AE,NBED,
     &                1,IBE,1,1,NE1,0,EINCR,0.,0.,2,SS2,YPL,IEERG)
                     IF (WR(4).GE.1.E-6) CALL VARSUM(AEP,AE,AE,AE,NBED,
     &                1,IBE,1,1,NE1,0,EINCR,0.,0.,3,SSP,YPL,IEERG)
                     IF (WR(5).GE.1.E-6) CALL VARSUM(AEP,AE,AE,AE,NBED,
     &                1,IBE,1,1,NE1,0,EINCR,0.,0.,2,SSP2,YPL,IEERG)
                     IF (.NOT.(WR(6).LT.1.E-6.AND.WR(9).LT.1.E-6)) CALL 
     &                VARSUM(AEPP,AE,AE,AE,NBED,1,IBE,1,1,NE1,0,EINCR,
     &                0.,0.,3,SSPP,YPL,IEERG)
                     IF (WR(7).GE.1.E-6) CALL VARSUM(AEPP,AE,AE,AE,NBED,
     &                1,IBE,1,1,NE1,0,EINCR,0.,0.,2,SSPP2,YPL,IEERG)
                     IF (WR(10).GE.1.E-6) CALL VARSUM(YE,AE,AE,AE,NBED,
     &                1,IBE,1,1,NE1,0,EINCR,0.,0.,2,SSY2,YPL,IEERG)
                  ENDIF
                  IF (NE2.NE.NE) THEN
                     CALL VARSUM(AE,AE,AE,AE,NBED,1,IBE,1,NE2,NE,0,
     &                EINCR,0.,0.,1,SU,YPL,IEERG)
                     IF (WR(3).GE.1.E-6) CALL VARSUM(AE,AE,AE,AE,NBED,
     &                1,IBE,1,NE2,NE,0,EINCR,0.,0.,2,SU2,YPL,IEERG)
                     IF (WR(4).GE.1.E-6) CALL VARSUM(AEP,AE,AE,AE,NBED,
     &                1,IBE,1,NE2,NE,0,EINCR,0.,0.,3,SUP,YPL,IEERG)
                     IF (WR(5).GE.1.E-6) CALL VARSUM(AEP,AE,AE,AE,NBED,
     &                1,IBE,1,NE2,NE,0,EINCR,0.,0.,2,SUP2,YPL,IEERG)
                     IF (.NOT.(WR(6).LT.1.E-6.AND.WR(9).LT.1.E-6)) CALL 
     &                VARSUM(AEPP,AE,AE,AE,NBED,1,IBE,1,NE2,NE,0,EINCR,
     &                0.,0.,3,SUPP,YPL,IEERG)
                     IF (WR(7).GE.1.E-6) CALL VARSUM(AEPP,AE,AE,AE,NBED,
     &                1,IBE,1,NE2,NE,0,EINCR,0.,0.,2,SUPP2,YPL,IEERG)
                     IF (WR(10).GE.1.E-6) CALL VARSUM(YE,AE,AE,AE,NBED,
     &                1,IBE,1,NE2,NE,0,EINCR,0.,0.,2,SUY2,YPL,IEERG)
                  ENDIF
                  SE=TSE(IBE)-SS-SU
                  SE2=TSE2(IBE)-SS2-SU2
                  SEP=TSEP(IBE)-SSP-SUP
                  SEP2=TSEP2(IBE)-SSP2-SUP2
                  SEPP=TSEPP(IBE)-SSPP-SUPP
                  SEPP2=TSEPP2(IBE)-SSPP2-SUPP2
                  SEY2=TSEY2(IBE)-SSY2-SUY2
                  SS=0.0
                  SSY2=0.0
                  SU=0.0
                  SUY2=0.0
                  IF (NT1.NE.1) THEN
                     CALL VARSUM(AT,AT,AT,AT,NT0,1,IBT,1,1,NT1,0,EINCR,
     &                0.,0.,1,SS,YPL,IEERG)
                      IF (WR(10).GE.1.E-6) CALL VARSUM(YT,AT,AT,AT,NT0,
     &                 1,IBT,1,1,NT1,0,EINCR,0.,0.,2,SSY2,YPL,IEERG)
                  ENDIF
                  IF (NT2.NE.NT) THEN
                     CALL VARSUM(AT,AT,AT,AT,NT0,1,IBT,1,NT2,NT,0,EINCR,
     &                0.,0.,1,SU,YPL,IEERG)
                      IF (WR(10).GE.1.E-6) CALL VARSUM(YT,AT,AT,AT,NT0,
     &                 1,IBT,1,NT2,NT,0,EINCR,0.,0.,2,SUY2,YPL,IEERG)
                  ENDIF
                  ST=TST(IBT)-SS-SU
                  STY2=TSTY2(IBT)-SSY2-SUY2
C
C  NORMALIZATION FACTOR EXP.-THEORY
C
                  C=SE/ST
                  CEPST=C*EPST
C
C  PRODUCE INTEGRALS INVOLVING BOTH EXP. AND THEORY
C
                  NV=NT1-NE1
                  IF (WR(2).GE.1.E-6) CALL VARSUM(AE,AT,AE,AE,NBED,NT0,
     &             IBE,IBT,NE1,NE2,NV,EINCR,0.,C,4,S,YPL,IEERG)
                  IF (WR(3).GE.1.E-6) CALL VARSUM(AE,AT,AE,AE,NBED,NT0,
     &             IBE,IBT,NE1,NE2,NV,EINCR,0.,C,5,S2,YPL,IEERG)
                  IF (WR(4).GE.1.E-6) CALL VARSUM(AEP,ATP,AE,AE,NBED,
     &             NT0,IBE,IBT,NE1,NE2,NV,EINCR,0.,C,4,SP,YPL,IEERG)
                  IF (WR(5).GE.1.E-6) CALL VARSUM(AEP,ATP,AE,AE,NBED,
     &             NT0,IBE,IBT,NE1,NE2,NV,EINCR,0.,C,5,SP2,YPL,IEERG)
                  IF (WR(6).GE.1.E-6) CALL VARSUM(AEPP,ATPP,AE,AE,NBED,
     &             NT0,IBE,IBT,NE1,NE2,NV,EINCR,0.,C,4,SPP,YPL,IEERG)
                  IF (WR(7).GE.1.E-6) CALL VARSUM(AEPP,ATPP,AE,AE,NBED,
     &             NT0,IBE,IBT,NE1,NE2,NV,EINCR,0.,C,5,SPP2,YPL,IEERG)
                  IF (WR(8).GE.1.E-6) CALL VARSUM(AEP,ATP,AEPP,ATPP,
     &             NBED,NT0,IBE,IBT,NE1,NE2,NV,EINCR,EPSE,C,6,SRZJ,YPL,
     &             IEERG)
                  IF (WR(9).GE.1.E-6) CALL VARSUM(AEP,ATP,AEPP,ATPP,
     &             NBED,NT0,IBE,IBT,NE1,NE2,NV,EINCR,CEPST,C,6,SMZJ,
     &             YPL,IEERG)
                  IF (WR(10).GE.1.E-6) CALL VARSUM(YE,YT,YE,YE,NBED,
     &             NT0,IBE,IBT,NE1,NE2,NV,EINCR,0.,1.,5,SY2,YPL,IEERG)
C
C  PRODUCE R-FACTORS (ALL ARE NORMALIZED TO ABOUT 1 FOR ANTICORRELATED
C  CURVES,I.E. FOR (SIN(E))**2 COMPARED WITH (COS(E))**2 OVER ONE PERIOD)
C
                  IF (WR(1).GE.1.E-6) THEN
C
C  R-FACTOR BASED ON FRACTION OF ENERGY INTERVAL WITH EXP. AND THEOR.
C  SLOPES OF OPPOSITE SIGNS
C
                     CALL OPSIGN(AEP,NBED,IBE,ATP,NT0,IBT,EINCR,NE1,NE2,
     &                NV,ROS(IBE),IEERG)
                     ROS(IBE)=ROS(IBE)/EET(IBE)
                     AROS=AROS+WB(IBE)*EET(IBE)*ROS(IBE)
                  ENDIF
                  IF (WR(2).GE.1.E-6) THEN
C
C  R-FACTOR BASED ON INTEGRAL OF ABS(EXP-C*TH)
C
                     R1(IBE)=S/SE
                     AR1=AR1+WB(IBE)*EET(IBE)*R1(IBE)
                  ENDIF
                  IF (WR(3).GE.1.E-6) THEN
C
C  R-FACTOR BASED ON INTEGRAL OF (EXP-C*TH)**2
C
                     R2(IBE)=S2/SE2
                     AR2=AR2+WB(IBE)*EET(IBE)*R2(IBE)
                  ENDIF
                  IF (WR(4).GE.1.E-6) THEN
C
C  R-FACTOR AS R1 BUT USING 1ST DERIVATIVES
C
                     RP1(IBE)=SP/SEP
                     ARP1=ARP1+WB(IBE)*EET(IBE)*RP1(IBE)
                  ENDIF
                  IF (WR(5).GE.1.E-6) THEN
C
C  R-FACTOR AS R2 BUT USING 1ST DERIVATIVES
C
                     RP2(IBE)=SP2/SEP2
                     ARP2=ARP2+WB(IBE)*EET(IBE)*RP2(IBE)
                  ENDIF
                  IF (WR(6).GE.1.E-6) THEN
C
C  R-FACTOR AS R1 BUT USING 2ND DERIVATIVES
C
                     RPP1(IBE)=SPP/SEPP
                     ARPP1=ARPP1+WB(IBE)*EET(IBE)*RPP1(IBE)
                  ENDIF
                  IF (WR(7).GE.1.E-6) THEN
C
C  R-FACTOR AS R2 BUT USING 2ND DERIVATIVES
C
                     RPP2(IBE)=SPP2/SEPP2
                     ARPP2=ARPP2+WB(IBE)*EET(IBE)*RPP2(IBE)
                  ENDIF
                  IF (WR(8).GE.1.E-6) THEN
C
C  REDUCED R-FACTOR ACCORDING TO ZANAZZI-JONA 
C
                     RRZJ(IBE)=SRZJ/(0.027*SE)
                     ARRZJ=ARRZJ+WB(IBE)*EET(IBE)*RRZJ(IBE)
                  ENDIF
                  IF (WR(9).GE.1.E-6) THEN
C
C  R-FACTOR ACCORDING TO MODIFIED ZANAZZI-JONA FORMULA .
C  EPS BASED ON EXP. INT. IS REPLACED BY EPS BASED ON THEOR. INT.  THE
C  NORMALIZ. INTEGRAL SE BASED ON INT. IS REPLACED BY AN EQUIVALENT
C  INTEGRAL BASED ON THE 2ND DERIVATIVE OF THE INT., THEREBY REMOVING
C  THE NEED FOR A REDUCED R-FACTOR (NO 0.027)
C
                     RMZJ(IBE)=SMZJ/SEPP
                     ARMZJ=ARMZJ+WB(IBE)*EET(IBE)*RMZJ(IBE)
                  ENDIF
                  IF (WR(10).GE.1.E-6) THEN
C
C  R-FACTOR ACCORDING TO PENDRY 
C
                     RPE(IBE)=SY2/(SEY2+STY2)
                     ARPE=ARPE+WB(IBE)*EET(IBE)*RPE(IBE)
                  ENDIF
C
C  AVERAGE OVER ABOVE R-FACTORS
C  FIRST CALCULATE TOTAL WEIGHT FOR NORMALIZATION
C
                  WS=0.
                  DO 220 I=1,10
                     WS=WS+WR(I)
220               CONTINUE
                  RAV(IBE)=(WR(1)*ROS(IBE)+WR(2)*R1(IBE)+WR(3)*R2
     &             (IBE)+WR(4)*RP1(IBE)+WR(5)*RP2(IBE)+WR(6)*RPP1
     &             (IBE)+WR(7)*RPP2(IBE)+WR(8)*RRZJ(IBE)+WR(9)
     &             *RMZJ(IBE)+WR(10)*RPE(IBE))/WS
                  ARAV=ARAV+WB(IBE)*EET(IBE)*RAV(IBE)
                  ERANG=ERANG+WB(IBE)*EET(IBE)
               ENDIF
C
C EITHER IBE OR IBT IS 0. THEREFOE, DON'T BOTHER WITH THIS BEAM
C
            ELSE
               NST1(IBE)=0
               NST2(IBE)=0
            ENDIF
C
C END OF LOOP OVER BEAMS
C
130      CONTINUE
         AR(1)=AROS/ERANG
         AR(2)=AR1/ERANG
         AR(3)=AR2/ERANG
         AR(4)=ARP1/ERANG
         AR(5)=ARP2/ERANG
         AR(6)=ARPP1/ERANG
         AR(7)=ARPP2/ERANG
         AR(8)=ARRZJ/ERANG
         AR(9)=ARMZJ/ERANG
         AR(10)=ARPE/ERANG
         AR(11)=ARAV/ERANG
         IF (AR(11).LT.BARAV) THEN
            BARAV=AR(11)
C            BV0=VO
         ENDIF
!jcm         write(*,*) 'RFAC: BARAV(FVAL)', BARAV
      RETURN
      END
C ==========================================================================
C Subroutine RFIN reads in the information required by the R-factor code.
C
C Input Parameters
C ----------------
C
C IBP     =  THEORETICAL BEAM AVERAGING INFORMATION
C NT0     =  NUMBER OF EXIT BEAMS
C WB      =  WEIGHTING INFORMATION FOR BEAMS IN EACH R-FACTOR
C WR      =  WEIGHTING INFORMATION FOR R-FACTORS IN 10 R-FACTOR AVERAGE
C
C AUTHOR:WANDER
C
C ==========================================================================
C
      SUBROUTINE RFIN(IBP,NT0,WB,WR,IPR)
C
      DIMENSION IBP(NT0),WR(10),WB(NT0)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA,FI
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
C
100   FORMAT (40I3)
!110   FORMAT (40F7.4)
140   FORMAT (/' Beam Averaging Information; ')
150   FORMAT (' =========================== ')
160   FORMAT (/' R-Factor Information; ')
170   FORMAT (' ===================== ')
180   FORMAT (/' Interpolation Step Size;',1F7.4)
!190   FORMAT (' Inner Potential Range;',2F7.4,' Step Size ',1F7.4)
C
C Read in the information needed to link together the experimental and
C theoretical beams. If IBP(i)=IBP(j) then theoretical beams i and j will
C be avveraged and compared with experimental beam IBP(i).
C
      MAXB=0
C++++ READ (12,100) (IBP(J),J=1,NT0)
      READ (12,*) (IBP(J),J=1,NT0)
      IF (IPR.GT.0) THEN
         WRITE (1,140)
         WRITE (1,150)
         WRITE (1,100) (IBP(J),J=1,NT0)
      ENDIF
      DO 148 I=1,NT0
         IF(IBP(I).GT.MAXB)MAXB=IBP(I)
148   CONTINUE
C
C Read in energy grid needed in R-factor calculation
C
C++++ READ (12,110) EINCR
      READ (12,*) EINCR
      IF (IPR.GT.0) THEN
         WRITE (1,160)
         WRITE (1,170)
         WRITE (1,180) EINCR
      ENDIF
C
C Read in weighting information for R-factors.
C WB=Beam weighting in each R-factor
C WR=R-factor weighting in 10 R-factor average
C
C++++ READ (12,110)(WB(I),I=1,MAXB)
C++++ READ (12,110)(WR(I),I=1,10)
      READ (12,*)(WB(I),I=1,MAXB)
      READ (12,*)(WR(I),I=1,10)
C
C Read in remaining parameters;
C ISMOTH=Desired number of experimental smoothing 
C IREN=1(0) Renormalize (or not) the experimental and theoretical curves.
C IRGEXP=0 if experimental data is on a uniform grid,
C        1 if data is on a non-uniform grid. In this case, an initial 
C          interpolation will be done. This is not recommended, as it may 
C          enhance experimental noise.
C
C++++ READ (12,100)ISMOTH,IREN,IRGEXP 
      READ (12,*)ISMOTH,IREN,IRGEXP 
      RETURN
      END
C ===========================================================================
C
C Subroutine RINTAV averages LEED intensities over different domains.
C
C Author: WANDER
C
C ===========================================================================
C
      SUBROUTINE RINTAV(ATH,NT0,IBP,IEERG,NERG,NBMAX,IBK)
C
      DIMENSION ATH(NT0,IEERG),IBP(NT0),TEMP(30,221),ITEMP(30),IBK(NT0)
C
C SET NBMAX TO THE MAXIMUM NUMBER OF INEQUIVALENT TH. BEAMS (AFTER DOMAIN
C AVERAGING
C
      NBMAX=0
      DO 100 I=1,NT0
         IF (IBP(I).GT.NBMAX) NBMAX=IBP(I)
100   CONTINUE
C
C TEMP is the temporary storage for the averaged th. beams
C
      DO 110 I=1,NBMAX
         DO 120 K=1,NERG
            TEMP(I,K)=0
120      CONTINUE
         ICOUNT=0
         DO 220 J=1,NT0
            IF (IBP(J).EQ.I) THEN
               DO 230 K=1,NERG
                    TEMP(I,K)=TEMP(I,K)+ATH(J,K)
230            CONTINUE
               ICOUNT=ICOUNT+1             
            ENDIF
220      CONTINUE
         IF (ICOUNT.GT.0) THEN
            DO 240 K=1,NERG
              TEMP(I,K)=TEMP(I,K)/FLOAT(ICOUNT)
240         CONTINUE
            ITEMP(I)=I
         ELSE
            ITEMP(I)=0
         ENDIF            
110   CONTINUE
      DO 300 I=1,NBMAX
         IBK(ITEMP(I))=I
         DO 310 J=1,NERG
            ATH(I,J)=TEMP(I,J)
310      CONTINUE
300   CONTINUE
      RETURN
      END
C =======================================================================
C
C Subroutine SDUMP dumps the current search parameters corresponding to
C the minimum to a restart file,
C from which they can be retrieved and used again.
C
C Input Parameters;
C =================
C 
C COORD        =   COORDINATE OF THE VECTOR TO BE DUMPED
C NNDIM        =   SECOND DIMENSION OF COORD
C
C AUTHOR: BARBIERI
C
C =======================================================================
C
      SUBROUTINE SDUMP(XI,EIGEN,FMIN,ADISP,NNDIM,NDIM,NLAY)
C
      PARAMETER (MAXC=100 )
C
      DIMENSION ADISP(NLAY,3),XI(NNDIM,NNDIM)
      DIMENSION EIGEN(MAXC)
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /POW/IFUNC,MFLAG,SCAL
C
C Rewind file to over write current contents (saves storage space)
C
      REWIND (UNIT=10)
C
C Dump dimensionality of search and search flag
C
      WRITE (10,*) NNDIM,MFLAG,NDIM
      DO 25 I=1,NLAY 
           WRITE (10,*) (ADISP(I,J),J=1,3)
25    CONTINUE
         WRITE (10,*) VOPT
         WRITE (10,*) FMIN
         IF(MFLAG.GE.2) THEN
           DO 16 I=1,NNDIM
             WRITE (10,*) (XI(J,I),J=1,NNDIM)
16         CONTINUE
             WRITE (10,*) (EIGEN(J),J=1,NNDIM)
         ENDIF
      RETURN
      END
C ==========================================================================
C                 
C Subroutine SETCOR generates a set of coordinates for input to FUNCV.
C
C Input Parameters:
C =================
C
C PR       = CURRENT VERTEX OF SIMPLEX
C DISP     = USER DEFINED COOERDINATES (FROM TLLED4.I)
C ADISP    = COORDINATES OF CURRENT VERTEX IN SAME FORMAT AS DISP
C NLAY     = NUMBER OF LAYERS UNDER CONSIDERATION
C NDIM     = CURRENT DIMENSIONALITY OF SEARCH (1=X AXIS ONLY
C                                          3=X,Y,Z AXES)
C NNDIM    = TOTAL DIMENSIONALITY OF SEARCH (=NLAY*NDIM)
C LSFLAG   =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C             NSYM) IN THE COMPOSITE LAYER. LLFLAG(i)=LLFLAG(j)
C             INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C             AS EQUIVALENT IN THE SEARCH
C LLFLAG   =   ARRAYS SPECIFYING WHETHER THE ATOM HAS TO BE INCLUDED
C              (LLFLAG(i)=1) IN THE SEARCH OR NOT (LLFLAG(i)=0)
C NDIML    =  ARRAY GIVING THE EFFECTIVE DIMENSIONALITY OF ATOM
C             j (ACTUALLY TO BE USED IN CONJUNCTION WITH LLFLAG)
C DIREC    =  SET OF DIRECTIONS SPECIFYING A BASIS OF UNIT VECTORS
C             DEFINED ON EACH ATOM AND RESPECTING THE SYMMETRY  
C             IMPOSED ON THE SEARCH
C
C AUTHOR: BARBIERI
C ==========================================================================
C
      SUBROUTINE SETCOR2(LSFLAG,NDIML,DIREC,PR,DISP,ADISP,
     % NLAY,NNDIM)
C
      DIMENSION PR(NNDIM),DISP(NLAY,3),ADISP(NLAY,3)
      DIMENSION NDIML(NLAY),DIREC(NLAY,2,2)
      DIMENSION LSFLAG(NLAY),N2(20)
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /POW/IFUNC,MFLAG,SCAL
C
      K=NNSTEF
      N1=0
      DO 100 I=1,NLAY
C
C Is this layer included in the search?
C
         IF (LSFLAG(I).NE.0) THEN
C
C Yes. But is it equivalent to some of the previous ones ?
C
            IF (LSFLAG(I).GT.N1) THEN
C
C No. Then its coordinates are in PR . Get them and take care of 
C equivalent atoms also
C
              N1=N1+1
              KK=0
              DO 120 II=I+1,NLAY
                IF (LSFLAG(II).EQ.N1) THEN
                   KK=KK+1
                   N2(KK)=II
                ENDIF
120           CONTINUE
C
C Now KK is the number of atoms different from I but equivalent to it.
C N2 points to their layer index.
C
C vertical coordinate
C
              K=K+1
              ADISP(I,1)=PR(K)
              DO 125 II=1,KK
                   ADISP(N2(II),1)=PR(K)
125           CONTINUE
C
C in plane coordinates
C
              IF(NDIML(N1).GT.1) THEN
                ADISP(I,2)=.0
                ADISP(I,3)=.0
                   DO 114 II=1,KK
                      ADISP(N2(II),2)=.0
                      ADISP(N2(II),3)=.0
114                CONTINUE
                DO 110 J=1,NDIML(N1)-1
                  K=K+1
                  ADISP(I,2)=ADISP(I,2)+PR(K)*DIREC(I,J,1)
                  ADISP(I,3)=ADISP(I,3)+PR(K)*DIREC(I,J,2)
                  DO 115 II=1,KK
                     ADISP(N2(II),2)=ADISP(N2(II),2)+PR(K)
     &                *DIREC(N2(II),J,1)
                     ADISP(N2(II),3)=ADISP(N2(II),3)+PR(K)
     &                *DIREC(N2(II),J,2)
115               CONTINUE
110             CONTINUE
              ENDIF
            ENDIF
         ELSE
C
C This layer is not included in the saerch. The coordinate is determined
C by DISP
C
            DO 135 KK=1,3
               ADISP(I,KK)=DISP(I,KK)
135         CONTINUE
         ENDIF
100   CONTINUE
C
C set non structural parameters
C
      IF (NNSTEF.GE.1) THEN
        VOPT=PR(1)
        IF(MFLAG.GE.2) VOPT=VOPT*SCAL
      ELSE
        VOPT=0.
      ENDIF
      RETURN
      END
C =========================================================================
C
C Subroutine SIMPL performs the down-hill simplex search for the R-factor
C minimum. At exit, all corners of the simplex will be within FTOL of the
C minimum value.
C
C Input Parameters;
C =================
C 
C NLAY                   = NUMBER OF LAYERS IN COMPOSITE LAYER
C NDIM                   = DIMENSIONALITY OF SEARCH (1=X AXIS ONLY
C                                                    3=X,Y,Z AXES)
C COORD                  = ARRAY STORING THE COORDINATES OF THE VERTICES 
C                         OF THE SIMPLEX
C DISP                   = COORDINATES INPUT BY USER
C DVOPT                  = SHIFT IN THE INNER POTENTIAL FOR THE STARTING 
C                          CONFIGURATION (in COMMON /RPL )
C LLFLAG                  = INDICATES WHETHER THE LAYER (OR NON STRUCTURAL
C                          PARAMETER) COORDINATES HAVE TO BE VARIED 
C                          IN THE SEARCH
C NNDIM                  = TOTAL NUMBER OF DIMENSIONS IN SEARCH 
C LSFLAG                 =  ARRAY SPECIFYING EQUIVALENT ATOMS (ACCORDING TO
C                          NSYM) IN THE COMPOSITE LAYER. LLFLAG(i)=LLFLAG(j)
C                          INDICATES THAT i and j ATOMS HAVE TO BE CONSIDERED
C                          AS EQUIVALENT IN THE SEARCH
C NDIML                  =  ARRAY GIVING THE EFFECTIVE DIMENSIONALITY OF ATOM
C                           j (ACTUALLY TO BE USED IN CONJUNCTION WITH LLFLAG)
C DIREC                  =  SET OF DIRECTIONS
C                          
C NNDIM1                 = NUMBER OF VERICES OF SIMPLEX (=NNDIM+1)
C ASTEP                  = INITIAL LENGTH OF SIDES OF SIMPLEX
C VSTEP                  = INITIAL LENGTH OF SIDE OF SIMPL. IN THE INNER
C                          POTENTIAL DIRECTION   
C PR                     = TEMPORARY STORAGE FOR NEW SIMPEX VERTEX
C ADISP                  = GEOMETRY OF CURRENT POINT IN SAME FORMAT AS DISP
C                              (USED AS INPUT TO FUNCV)
C VOPT                   = INNER POTENTIAL OF CURRENT INPUT (used as input
C                          in FUNCV), in COMMON /NSTR
C FTOL1                  = CONVERGENCE CRITERIA FOR SIMPLEX SEARCH
C ITMAX                  = MAXIMUM NUMBER OF ITERATIONS TO BE PERFORMED BY 
C                          SIMPL
C PBAR                   = VECTOR AVERAGE OF VERTICES OF SIMPLEX. USED TO CHOOSE
C                         NEXT DISTORTION DIRECTION OF SIMPLEX
C VAL                    = VALUE OF THE R-FACTOR AT THE POINT ADISP RETURNED BY
C                        ROUTINE FUNCV
C AALPHA,BBETA,GGAMMA    = PARAMETERS GOVERNING DISTORTIONS OF SIMPLEX
C PRR                    = TEMPORARY STORAGE FOR NEW SIMPEX VERTEX
C ISTART                 = 0 USE PARAMETERS FROM TLEED4.I
C                          1 RESTART SEARCH FROM COORDINATES IN RESTART.D
C IPR                    = PRINT CONTROL PARAMETER
C ILOOK                  = LOOKUP TABLE FOR DOMAIN AVERAGING
C ACOORD                 = STORAGE FOR SYMMETRY EQUIVALENT COORDINATES
C
C MICUT,MJCUT,PSQ,JYLM,  = DUMMY ARRAYS TO BE PASSED TO ROUTINE FUNCV
C BJ,YLM,QS,AT,XISTS,
C XIST,INBED,IEERG,AE,EE,
C NEE,NBEA,BENAME,IPR,XPL,
C YPL,NNN,AP,APP,YE,TSE,
C TSE2,TSEP,TSEP2,TSEPP,
C TSEPP2,TSEY2,WR,WB,IBP,
C ETH
C
C AUTHOR: WANDER,BARBIERI   based on NUMER. RECIP.
C =======================================================================
C
Cga made into function call
Cga      SUBROUTINE SIMPL(NLAY,NDIM,COORD,DISP,NNDIM,NNDIM1,ASTEP,VSTEP,
      function SIMPL(NLAY,NDIM,COORD,DISP,NNDIM,NNDIM1,ASTEP,VSTEP,
     & PR,ADISP,FTOL1,ITMAX,PBAR,VAL,AALPHA,BBETA,GGAMMA,PRR,ISTART,
     & IPR,ILO,ILOOK,ACOORD,MICUT,MJCUT,PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,
     & NERG,AT,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,XPL,YPL,NNN,AP,APP,YE,
     & TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,ETH,
     & ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,R2,RP1,RP2,RPP1,
     & RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,YT,LLFLAG,LSFLAG,
     & NDIML,DIREC,NLTIN,LPOINT)
C
      PARAMETER (MAXC=100 )
C
C     MAXC is the maximum number of coordinates of the parameter space
C          where we perform the search
C
      
      DIMENSION COORD(NNDIM1,NNDIM),DISP(NLAY,3)
      DIMENSION DIREC(NLAY,2,2),NDIML(NLAY)
      DIMENSION PR(NNDIM),ZMIN(MAXC),LLFLAG(NLAY+NNST),LSFLAG(NLAY)
      DIMENSION ADISP(NLAY,3),VAL(NNDIM1),PBAR(NNDIM),PRR(NNDIM)
      DIMENSION ILOOK(12,NLAY),ACOORD(12,NLAY,3)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG)
      DIMENSION MICUT(ICUT),MJCUT(ICUT),PSQ(2,NT0)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPOINT(NLTIN)
C++++ 
      DIMENSION AP(INBED,IEERG),APP(INBED,IEERG),YE(INBED,IEERG)
      DIMENSION XI(NNDIM,NNDIM)
C++++
C     COMPLEX JYLM(LSMMAX),BJ(LSMAX1)
      COMPLEX JYLM(LSMMAX)
      COMPLEX BJ(LSMAX1)
      COMPLEX YLM(LSMMAX),QS(IQSIZ),XISTS(NT0,NERG),XIST(NT0,NERG)
C
      COMMON /TLVAL/LSMAX,LSMMAX,ICUT,LSMAX1,NT0,IQSIZ
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /RPL/DVOPT
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
C
500   FORMAT (/' SIMPLEX SEARCH COMPLETE ')
501   FORMAT (' ======================= ')
502   FORMAT (' DIMENSIONALITY OF SEARCH=',I5)
503   FORMAT (' ISMOTH=',I5)
505   FORMAT (/' CONVERGENCE TOLERANCE ACHIEVED =',F7.6)
506   FORMAT (/' COORDINATES AT MINIMUM;  LSFLAG '/)
510   FORMAT (3F9.4,I4)
513   FORMAT (' NUMBER OF FUNCTIONAL EVALUATIONS =',I4)
514   FORMAT (/' OPTIMUM R-FACTOR =',1F7.4)
515   FORMAT (' OPTIMUM VALUE OF INNER POTENTIAL =',1F7.4)
516   FORMAT (/' NUMBER OF ITERATIONS =',1I3)
517   FORMAT (/' NUMBER OF REPEATED SEARCHES =',1I3)
520   FORMAT (' SEARCH HAS EXCEEDED MAXIMUM NUMBER OF ITERATIONS')
521   FORMAT (' OTHER VERTICES HAVE COORDINATES; ')
522   FORMAT (' R-FACTOR =',1F7.4)
523   FORMAT ('INNER POTENTIAL =',1F7.4)
C 
C The search reaches its first step when a certain convergence (FTOL1) 
C is reached in R-factor value.
C To control somewhat also the convergence in parameter space we are going 
C to proceed to a next step in which we repeat the search again 
C starting from the presumed minimum and compare it to the new point obtained.
C The search ends when on each coordinate we get a convergence of order
C sqrt(xtols/NNDIM).(IPAR counts the number of iter. of this second step)
C This is however a delicate point because the search is usually
C much more sensitive to some coordinates than to others and in this way
C we are considering all coordinates to be equivalent  
C 
      FNNDIM=FLOAT(NNDIM)
      XTOL=1.E-4*FNNDIM 
      VTOL=1.E-3
C
C IF NNDIM=0 return information about the input structure
C
      IF(NNDIM.EQ.0) GOTO 133
C
C First initialize vertices of simplex. Set origin from DISP.
C The last parameter is the inner potential
C 
      IF (ISTART.EQ.0) THEN
         KK=0
C
C set nonstructural parameters
C
         IF (NNST.GE.1) THEN
           DO 103 I=NLAY+1,NLAY+NNST
             IF (LLFLAG(I).NE.0) THEN
               KK=KK+1
               COORD(1,KK)=0.
             ENDIF
103        CONTINUE
         ENDIF
         N1=0
         DO 100 I=1,NLAY
C
C Is this layer to be included in the search?
C
            IF (LSFLAG(I).GT.N1) THEN
C
C yes. Then record its coordinates according to the effective 
C dimensionality of the layer.(Notice that if the dimensionality
C is 2 then the Y coordinate read from the input actually
C corrisponds to the coefficient of the unit vector specified
C by DIREC; but this will be important in SETCOR)
C
              N1=N1+1
              DO 110 J=1,NDIML(N1)
                 KK=KK+1
                 COORD(1,KK)=DISP(I,J)
110           CONTINUE
            ENDIF
100      CONTINUE
      ELSE
C
C If restarting, upload coordinates from restart file
C
C++++
C I don't know whether this change is correct (M. Gierer)
C        CALL SRETV(FMIN,XI,ZMIN,P,LLFLAG,LSFLAG,
C    &    NNDIM,JJDIM,NLAY,ADISP,NDIML,DIREC)
         CALL SRETV(FMIN,XI,ZMIN,PR,LLFLAG,LSFLAG,
     &    NNDIM,JJDIM,NLAY,ADISP,NDIML,DIREC)
         DO 101 K=1,NNDIM
             COORD(1,K)=PR(K)
101      CONTINUE         
      ENDIF
C
C store the initial minimum coordinates
C
      DO 102 K=1,NNDIM
          ZMIN(K)=COORD(1,K)
102   CONTINUE         
C
C Set remaining vertices by addition of vector ASTEP to leading diagonal
C Add VSTEP for the inner potential
      IFUNC=0
      ICOUNT=0
      IPAR=0
1002  CONTINUE
      IPAR=IPAR+1
      DO 115 J=1,NNDIM
         DO 120 K=2,NNDIM1
             COORD(K,J)=COORD(1,J)
             IF (K-1.EQ.J.AND.J.LE.NNDIM-NNSTEF) 
     &       COORD(K,J)=COORD(K,J)+ASTEP
             IF (NNSTEF.GE.1.AND.K-1.EQ.J.AND.J.EQ.NNDIM) 
     &       COORD(K,J)=COORD(K,J)+VSTEP
120      CONTINUE
115   CONTINUE
C
C Now set VAL to be value of simplex at each vertex. Begin by copying
C coordinates of each vertex into temporary storage
C
      DO 130 I=1,NNDIM1
         DO 135 J=1,NNDIM
            PR(J)=COORD(I,J)
135      CONTINUE
C
C Set coordinate matrix ADISP from temporary store PR
C
         CALL SETCOR2(LSFLAG,NDIML,DIREC,PR,DISP,ADISP,
     %         NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
         IFUNC=IFUNC+1
         VAL(I)=FVAL
130   CONTINUE
C
C Begin search. First isolate the highest and lowest points of the simplex
C (IHI and ILO respect.) and the second highest point INHI.
C
1000  ILO=1
      IHI=1
      INHI=1
      VHI=VAL(1)
      VNHI=VAL(1)
      VLO=VAL(1)
      DO 140 I=2,NNDIM1
         IF (VAL(I).GT.VHI) THEN
            VHI=VAL(I)
            IHI=I
         ELSEIF (VAL(I).LT.VLO) THEN
            VLO=VAL(I)
            ILO=I
         ENDIF
         IF (VAL(I).GT.VNHI.AND.VAL(I).LT.VHI) THEN
            VNHI=VAL(I)
            INHI=I
         ENDIF
140   CONTINUE
C
C Compute range of R-factor spread as fractional value
C
      FRAC=ABS(VHI-VLO)/ABS(VHI+VLO)
      IF (FRAC.LE.FTOL1) GOTO 1001
      IF (ICOUNT.LE.ITMAX) THEN
         ICOUNT=ICOUNT+1
         DO 160 I=1,NNDIM
            PBAR(I)=0.0
160      CONTINUE
         DO 170 I=1,NNDIM1
            IF (I.NE.IHI) THEN
C
C Average all points except the highest
C
               DO 180 J=1,NNDIM
                  PBAR(J)=PBAR(J)+COORD(I,J)
180            CONTINUE
            ENDIF
170      CONTINUE
C
C Reflect the simplex from the highest point
C
         DO 190 J=1,NNDIM
            PBAR(J)=PBAR(J)/NNDIM
            PR(J)=(1.0+AALPHA)*PBAR(J)-AALPHA*COORD(IHI,J)
190      CONTINUE
C
C Evaluate the R-factor at the new point
C
C Set coordinate matrix ADISP
C
         CALL SETCOR2(LSFLAG,NDIML,DIREC,PR,DISP,ADISP,
     %         NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
         IFUNC=IFUNC+1
C
C Is the new point lower than the best?
C
         YPR=FVAL
         IF (YPR.LE.VAL(ILO)) THEN
C
C If so, try a further enhancement by a factor GGAMMA
C
            DO 200 J=1,NNDIM
               PRR(J)=GGAMMA*PR(J)+(1.0-GGAMMA)*PBAR(J)
200         CONTINUE
C
C Evaluate at the new point
C Set coordinate matrix ADISP
C
            CALL SETCOR2(LSFLAG,NDIML,DIREC,PRR,DISP,ADISP,
     %         NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
            IFUNC=IFUNC+1
C
C Did this extra move work?
C
            YPRR=FVAL 
            IF (YPRR.LT.VAL(ILO)) THEN
C
C If so, use this new point to replace the high point
C
               DO 210 J=1,NNDIM
                  COORD(IHI,J)=PRR(J)
210            CONTINUE
               VAL(IHI)=YPRR
            ELSE
C
C If the extra extrapolation failed, just use the reflacted point
C
               DO 220 J=1,NNDIM
                  COORD(IHI,J)=PR(J)
220            CONTINUE
               VAL(IHI)=YPR
            ENDIF
         ELSEIF (YPR.GE.VAL(INHI)) THEN
            IF (YPR.LT.VAL(IHI)) THEN
C
C The new point is the second worst found so far
C
               DO 230 J=1,NNDIM
                  COORD(IHI,J)=PR(J)
230            CONTINUE
               VAL(IHI)=YPR
            ENDIF
            DO 240 J=1,NNDIM
C
C Look for an intermediate lower point by contracting the simplex
C
               PRR(J)=BBETA*COORD(IHI,J)+(1.0-BBETA)*PBAR(J)
240         CONTINUE
C
C Set coordinate matrix ADISP
C
            CALL SETCOR2(LSFLAG,NDIML,DIREC,PRR,DISP,ADISP,
     %         NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
            IFUNC=IFUNC+1
            YPRR=FVAL
            IF (YPRR.LT.VAL(IHI)) THEN
C
C The contraction gives an improvement so accept it
C
               DO 250 J=1,NDIM
                  COORD(IHI,J)=PRR(J)
250            CONTINUE
               VAL(IHI)=YPRR
            ELSE
C
C Highest point seems to unremovable by reflection. Instead, contract
C around the lowest point
C
               DO 260 I=1,NNDIM1
                  IF (I.NE.ILO) THEN
                     DO 270 J=1,NNDIM
                        PR(J)=0.5*(COORD(I,J)+COORD(ILO,J))
                        COORD(I,J)=PR(J)
270                  CONTINUE
C
C Compute new value at this point
C Set coordinate matrix ADISP
C
                CALL SETCOR2(LSFLAG,NDIML,DIREC,PR,DISP,ADISP,
     %         NLAY,NNDIM)
C
C Generate value of R-Factor at point ADISP
C
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
                     IFUNC=IFUNC+1
                     VAL(I)=FVAL
                  ENDIF
260            CONTINUE
            ENDIF
         ELSE
C
C Nothing special about the current point. Use it to replace the high
C point and carry on as before
C
            DO 280 J=1,NNDIM
               COORD(IHI,J)=PR(J)
280         CONTINUE
            VAL(IHI)=YPR
         ENDIF
         GOTO 1000
      ENDIF
C
C Iteration count exceeded, so dump coordinates of the minimum
C to retsart file
C
      WRITE (2,520)
      DO 146 I=1,NNDIM
         PR(I)=COORD(ILO,I)
146   CONTINUE
      CALL SETCOR2(LSFLAG,NDIML,DIREC,PR,DISP,ADISP,
     & NLAY,NNDIM)
      CALL SDUMP(COORD,ZMIN,VAL(ILO),ADISP,NNDIM,NDIM,NLAY)
      RETURN
C
C Convergence achieved in R factor value. Is this true also in parameter
C space value or has the simplex contracted so much in volume to produce
C a fake convergence?
C
1001  CONTINUE
      DO 145 I=1,NNDIM
         PR(I)=COORD(ILO,I)
145   CONTINUE
C
C     check the distance between the previous "minimum" and the new
C     "converged" minimum. Convergence is checked on structural parameters
C     only.
C
      DIST=.0
      DO 147 I=1,NNDIM
         IF (I.NE.NNDIM) THEN
            DIST=DIST+(PR(I)-ZMIN(I))**2
         ELSE
            DISTV=(PR(I)-ZMIN(I))**2
         ENDIF
         ZMIN(I)=PR(I)
         COORD(1,I)=PR(I)
147   CONTINUE
C
C     If we don't have an acceptable minimum let's start again
C
      IF (DIST.GE.XTOL.OR.DISTV.GE.VTOL) GOTO 1002
C
C  Dump relevant information in RESTART.D
C
133   CALL SETCOR2(LSFLAG,NDIML,DIREC,PR,DISP,ADISP,
     %         NLAY,NNDIM)
C
C  Recompute intensity for the minimum configuration and write the 
C  IV curves for the best structure
C
      IFUNC1=IFUNC
      IFUNC=0
      CALL FUN2(FVAL,NLAY,ADISP,VOPT,ILOOK,ACOORD,MICUT,MJCUT,
     & PSQ,JYLM,BJ,YLM,QS,XISTS,XIST,AT,INBED,IEERG,AE,EE,NEE,NBEA,
     & BENAME,IPR,XPL,YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,
     & TSEPP2,TSEY2,WR,WB,IBP,ETH,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,
     & IBK,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,
     & YT,NLTIN,LPOINT)
      CALL SDUMP(COORD,ZMIN,FVAL,ADISP,NNDIM,NDIM,NLAY)
      WRITE(10,*) (EEAVT(J),J=1,NBE+1)
      WRITE(10,*) (EEAVE(J),J=1,NBE+1)
      WRITE (2,500)
      WRITE (2,501)
      CALL WRIV(AT,ETH,AE,EE,IEERG,NT0,NBED,IBK,WR,
     & ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,BENAME,
     & NBE,NST1,NST2)
      VOPT=VOPT+DVOPT+VV*27.21
      IF (NNDIM.EQ.0) THEN
        WRITE (2,514) FVAL
        GOTO 301
      ENDIF
      WRITE (2,505) FRAC
      WRITE (2,502) NNDIM
      WRITE (2,503) ISMOTH
      WRITE (2,506)
      DO 150 I=1,NLAY
         WRITE (2,510) (ADISP(I,J),J=1,3),LSFLAG(I)
150   CONTINUE
      WRITE (2,516) ICOUNT
      WRITE (2,517) IPAR
      WRITE (2,513) IFUNC1
      WRITE (2,514) VAL(ILO)
      WRITE (2,515) VOPT
      IF (IPR.GT.0) THEN
         WRITE (2,521)
         DO 300 I=1,NNDIM1
            IF (I.NE.ILO) THEN
               DO 310 J=1,NNDIM
                  PR(J)=COORD(I,J)
310            CONTINUE
               CALL SETCOR2(LSFLAG,NDIML,DIREC,PR,DISP,ADISP,
     %         NLAY,NNDIM)
               VOPT=VOPT+DVOPT+VV*27.21
               DO 1514 K=1,NLAY
                  WRITE (2,510) (ADISP(K,J),J=1,3),LSFLAG(K)
1514            CONTINUE
               WRITE (2,522) VAL(I)
               WRITE (2,523) VOPT
            ENDIF
300      CONTINUE
      ENDIF
Cga
Cga301   RETURN
301   continue
      simpl = val(ilo)
      RETURN
      END
C =========================================================================
C
C Subroutine SMOOTH smoothes a given set of data by weighted three-point
C smoothing.
C
C =========================================================================
C
      SUBROUTINE SMOOTH(AE,EE,NBED,NBE,NEE,IPR,IEERG)
C
      DIMENSION AE(NBED,IEERG),EE(NBED,IEERG),NEE(NBE)
C 
60    FORMAT ('EXP. ENERG. AND INTENS. AFTER SMOOTHING IN BEAM',1I3,
     & /,50(5(1F7.2,1E13.4,3X),/))
      DO 40 IB=1,NBE
         N=NEE(IB)-1
C
C  IF TOO FEW POINTS, NO SMOOTHING POSSIBLE
C
         IF (N.GE.3) THEN
            AM=AE(IB,1)
            DO 30 IE=2,N
               AF=AE(IB,IE)
               E21=EE(IB,IE)-EE(IB,IE-1)
               E32=EE(IB,IE+1)-EE(IB,IE)
               IF (ABS(E21-E32).LT.0.0001) THEN
                  AE(IB,IE)=0.25*(2.*AF+AM+AE(IB,IE+1))
               ELSE
                  AE(IB,IE)=0.5*(AF+(E32*AM+E21*AE(IB,IE+1))/(E21+E32))
               ENDIF
               AM=AF
30          CONTINUE
         ENDIF
40    CONTINUE
      IF (IPR.GE.2) THEN
         DO 50 IB=1,NBE
            N=NEE(IB)
            IF (N.NE.0) WRITE (1,60) IB,(EE(IB,IE),AE(IB,IE),IE=1,N)
50       CONTINUE
      ENDIF
      RETURN
      END
C =======================================================================
C
C Subroutine SRETV retrieves the current search parameters from the 
C restart file, to which they were written by SDUMP.
C NLAY and NNST should be unchanged. IT is required that 
C NDIMO<NDIM and DVOPT is not changed from the previous run
C
C
C COORD       =   ARRAY CONTAINING COORDINATES OF THE STARTING CONFIGURATION
C                 READ FROM restart.d
C
C Author : Barbieri
C =======================================================================
C
      SUBROUTINE SRETV(FMIN,XI,EIGEN,P,LLFLAG,LSFLAG,
     & NNDIM,JJDIM,NLAY,ADISP,NDIML,DIREC)
C
      PARAMETER (MAXC=100 )
C
      DIMENSION P(NNDIM),LLFLAG(NLAY+NNST),ADISP(NLAY,3),NDIML(NLAY)
      DIMENSION LSFLAG(NLAY),DIREC(NLAY,2,2)
      DIMENSION XI(NNDIM,NNDIM),EIGEN(MAXC)
      COMMON /NSTR/VOPT,NNST,NNSTEF
      COMMON /POW/IFUNC,MFLAG,SCAL
C
1200  FORMAT (' CURRENT DIMENSION OF PARAMETER SPACE; ',I3)
1300  FORMAT (' DIMENSION FROM FILE RESTART.D; ',I3)
C
C Rewind file to start
C
      REWIND (UNIT=10)
C
C Upload total dimensionality , MFLAG, and NDIM used in previous run 
C
      READ (10,*) JJDIM,MFLAGO,NDIMO
      DO 25 I=1,NLAY
           READ (10,*) (ADISP(I,J),J=1,3)
25    CONTINUE
C
C
      READ (10,*) VOPT
      READ (10,*) FMIN
      WRITE (2,1200) NNDIM
      WRITE (2,1300) JJDIM
C
C generate COORD 
C
      KK=NNSTEF
      N1=0
      DO 105 I=1,NLAY
C
C Is this layer to be included in the search?
C
         IF (LSFLAG(I).GT.N1) THEN
C
C yes
C
              N1=N1+1
              DO 115 J=1,NDIML(N1)
                 KK=KK+1
                 P(KK)=ADISP(I,J)
C
C if ndiml=2 project the displacement onto the basis vector used in
C SETCOR
C
                 IF (NDIML(N1).EQ.2.AND.J.EQ.2) THEN
                   ZMOD=(DIREC(I,1,1)**2+DIREC(I,1,2)**2)
                   P(KK)=(ADISP(I,2)*DIREC(I,1,1)+ADISP(I,3)
     &              *DIREC(I,1,2))/ZMOD
                 ENDIF
115           CONTINUE
          ENDIF
105   CONTINUE
C
C set nonstructural parameters
C
      KK=0
      IF (NNST.GE.1) THEN
           DO 103 I=NLAY+1,NLAY+NNST
             IF (LLFLAG(I).NE.0) THEN
               KK=KK+1
               P(KK)=VOPT
               IF(MFLAG.GE.2) P(KK)=P(KK)/SCAL
             ENDIF
103        CONTINUE
      ENDIF
      IF(MFLAG.GE.2) THEN
           DO 120 I=1,NNDIM
              DO 130 J=1,NNDIM
                 IF (J.EQ.I) THEN
                    XI(I,J)=1.0
                 ELSE
                    XI(I,J)=.0
                 ENDIF
130           CONTINUE
              EIGEN(I)=1.
120      CONTINUE
      ENDIF
C
C now retrieve directions and eigenvalues if needed
C
      IF(NNDIM.EQ.JJDIM) THEN
        IF(MFLAG.EQ.MFLAGO.OR.MFLAG.EQ.4) THEN
             DO 16 I=1,NNDIM
               READ (10,*) (XI(J,I),J=1,NNDIM)
16           CONTINUE
             READ (10,*) (EIGEN(J),J=1,NNDIM)
        ENDIF
      ENDIF
      RETURN
      END
C====================================================================
C
C  SUBROUTINE STFPTS FINDS, GIVEN THE INTERPOLATION INTERVAL, THE
C  FOUR NEAREST GRID POINTS AND THE CORRESPONDING ORDINATE VALUES
C
C====================================================================
C
      SUBROUTINE STFPTS(IL,IH,WORX,WORY,LENGTH)
C
      DIMENSION WORX(LENGTH),WORY(LENGTH),TEMP(2, 4)
C
      COMMON / DATBLK / X0, X1, X2, X3, Y0, Y1, Y2, Y3
      I = IL - 1
      IF(IL .LE. 1) I = I + 1
      IF(IH .GE. LENGTH) I = I - 1
      DO 10 K = 1, 4
      N = K + I - 1
      TEMP(1,K) = WORX(N)
      TEMP(2, K) = WORY(N)
10    CONTINUE
      X0 = TEMP(1, 1)
      X1 = TEMP(1, 2)
      X2 = TEMP(1, 3)
      X3 = TEMP(1, 4)
      Y0 = TEMP(2, 1)
      Y1 = TEMP(2, 2)
      Y2 = TEMP(2, 3)
      Y3 = TEMP(2, 4)
      RETURN
      END
C =========================================================================
C
C Subroutine SUM integrates by the simple trapezoid rule (after Zanazzi-Jona)
C
C =========================================================================
C
      SUBROUTINE SUM(Y,NBD,IB,H,I1,I2,S,IEERG)
C
      DIMENSION Y(NBD,IEERG)
C
      A=0.
      S=0.
      DO 10 J=I1,I2
         A=A+Y(IB,J)
10    CONTINUE
      S=A-0.5*(Y(IB,I1)+Y(IB,I2))
      S=S*H
      RETURN
      END
C =========================================================================
C
C Subroutine VARSUM integrates over various combinations of the input 
C functions (tabulated) A1,A2,B1,B2, depending on the value of NF.
C NV is a relative shift of the X-axis between functions. IE1,IE2 are
C the integration limits. With the integrand of the Zanazzi-Jona R-factor
C a 10 times denser grid is used and first interpolated to;
C
C   NF     INTEGRAND
C
C    1       A1
C    2       A1**2
C    3       ABS(A1)
C    4       ABS(A1-C*A2)
C    5       (A1-C*A2)**2
C    6       ABS(B1-C*B2)*ABS(A1-C*A2)/(ABS(A1)+EPS)
C
C =========================================================================
C
      SUBROUTINE VARSUM(A1,A2,B1,B2,NBD1,NBD2,IB1,IB2,
     & IE1,IE2,NV,EINCR,EPS,C,NF,S,Y,IEERG)
C
      DIMENSION A1(NBD1,IEERG),A2(NBD2,IEERG),
     & B1(NBD1,IEERG),B2(NBD2,IEERG),Y(IEERG),Y1(221),
     & Y2(221),Y3(221),Y4(221),YY(2210)
C
      N=0
C
C  FOR ZANAZZI-JONA R-FACTOR INTERPOLATION ONTO 10-FOLD DENSER GRID
C  IS MADE (This is not needed any more (Barbieri))
C
      IF (NF.EQ.7) THEN
         DO 110 IE=IE1,IE2
            N=N+1
            IES=IE+NV
            Y(N)=FLOAT(N-1)*EINCR
            Y1(N)=A1(IB1,IE)
            Y2(N)=A2(IB2,IES)
            Y3(N)=B1(IB1,IE)
            Y4(N)=B2(IB2,IES)
110      CONTINUE
         DE=0.1*EINCR
         NN=10*(N-1)+1
         DO 120 IE=1,NN
            X=FLOAT(IE-1)*DE
            ITIL=0
            ITIH=0
            CALL YVAL(AA1,AA1P,X,Y1(1),Y,N,ITIL,ITIH)
            ITIL=0
            ITIH=0
            CALL YVAL(AA2,AA2P,X,Y2(1),Y,N,ITIL,ITIH)
            ITIL=0
            ITIH=0
            CALL YVAL(AB1,AB1P,X,Y3(1),Y,N,ITIL,ITIH)
            ITIL=0
            ITIH=0
            CALL YVAL(AB2,AB2P,X,Y4(1),Y,N,ITIL,ITIH)
            YY(IE)=ABS(AB1-C*AB2)*ABS(AA1-C*AA2)/(ABS(AA1)+EPS)
120      CONTINUE
         CALL SUM(YY,1,1,DE,1,NN,S,IEERG)
      ELSE
         DO 80 IE=IE1,IE2
            N=N+1
            IES=IE+NV
            IF (NF.EQ.2) THEN
               Y(N)=A1(IB1,IE)**2
            ELSEIF (NF.EQ.3) THEN
               Y(N)=ABS(A1(IB1,IE))
            ELSEIF (NF.EQ.4) THEN
               Y(N)=ABS(A1(IB1,IE)-C*A2(IB2,IES))
            ELSEIF (NF.EQ.5) THEN
               Y(N)=(A1(IB1,IE)-C*A2(IB2,IES))**2
            ELSEIF (NF.EQ.6) THEN
               Y(N)=ABS(B1(IB1,IE)-C*B2(IB2,IES))*ABS(A1(
     &          IB1,IE)-C*A2(IB2,IES))/(ABS(A1(IB1,IE))+EPS)
            ELSE
               Y(N)=A1(IB1,IE)
            ENDIF
80       CONTINUE
         CALL SUM(Y,1,1,EINCR,1,N,S,IEERG)
      ENDIF
      RETURN
      END
C=======================================================================
C
C  SUBROUTINE XNTERP PERFORMS 3RD-ORDER POLYNOMIAL INTERPOLATION
C  AND RETURNS THE FUNCTION (Y) AND ITS DERIVATIVE (YP) AT THE INPUT
C  POINT X
C
C=======================================================================
C
      SUBROUTINE XNTERP(Y,YP,X)
C
      COMMON / DATBLK / X0, X1, X2, X3, Y0, Y1, Y2, Y3
C
      TERM = Y0
      FACT1 = X - X0
      F1P = X - X1
      F2P = X - X2
      F12P=F1P*F2P
      F02P=FACT1*F2P
      D2=FACT1+F1P
      FACT2 = (Y1 - Y0) / (X1 - X0)
      TERM = TERM + FACT1 * FACT2
      TERMP = FACT2
      FACT1 = FACT1 * (X - X1)
      D3=FACT1+F02P+F12P
      FACT2 = ((Y2 - Y1)/(X2 - X1) - FACT2) / (X2 - X0)
      TERM = TERM + FACT1 * FACT2
      TERMP = TERMP + D2*FACT2
      FACT1 = FACT1 * (X - X2)
      TEMP = ((Y3 - Y2)/(X3 - X2) - (Y2 - Y1)/(X2 - X1))/(X3 - X1)
      FACT2 = (TEMP - FACT2) / (X3 - X0)
      Y = TERM + FACT1 * FACT2
      YP = TERMP + D3*FACT2
      RETURN
      END
C=======================================================================
C
C  SUBROUTINE XNTRP2 PERFORMS 2ND OR 1ST ORDER POLYNOMIAL INTERPOLATION
C
C=======================================================================
C
      SUBROUTINE XNTRP2(X,Y,YP,XS,YS,N)
C
      DIMENSION XS(N),YS(N)
C
      IF (N.GT.2) GO TO 10
      Y=(YS(2)-YS(1))*(X-XS(1))/(XS(2)-XS(1))+YS(1)
      YP=(YS(2)-YS(1))/(XS(2)-XS(1))
      RETURN
10    A=(YS(1)-YS(2))/(XS(1)-XS(2))/(XS(2)-XS(3))-
     &  (YS(1)-YS(3))/(XS(1)-XS(3))/(XS(2)-XS(3))
      B=(YS(1)-YS(2))/(XS(1)-XS(2))-A*(XS(1)+XS(2))
      C=YS(1)-(A*XS(1)+B)*XS(1)
      Y=C+X*(B+A*X)
      YP=B+2.*A*X
      RETURN
      END
C =======================================================================
C 
C Subroutine WRIV dumps the IV curves for the best structure found
C
C Input Parameters;
C =================
C
C VOPT           =  Optimal value of the inner potential
C AT             =  Theoretical Intensities of the best configuration
C AE             =  Experimental intensities
C ETH            =  Theoretical energies
C EE             =  Experimental energies
C NBED           =  Number of experimental beams
C NT0            =  Number of exit beams 
C IEERG          =  Maximum Number of energy points after interpolation
C NET            =  Number of data points in each THEORETICAL beam
C NEE            =  Number of data points in each EXPERIMENTAL beam
C NBEA           =  Beam averaging information
C BENAME         =  Identifier for each experimental beam
C IBK            =  IBK(I) is the Theoretical beam corresponding to
C                   the experimental beam I
C
C COMMON BLOCKS
C =============
C
C NBMAX          =  Number of beams after averaging  (block WIV)
C EEAVE          =  Average intensity of the peaks in the experimental beams
C EEAVT          =  Average intensity of the peaks in the theoretical beams
C
C Author: Barbieri
C
C ==========================================================================
C
      SUBROUTINE WRIV(AT,ETH,AE,EE,IEERG,NT0,NBED,IBK,
     & WR,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,BENAME,
     & NBE,NST1,NST2)
C
      DIMENSION ETH(NT0,IEERG),AT(NT0,IEERG)
      DIMENSION AE(NBED,IEERG),EE(NBED,IEERG)
      DIMENSION IBK(NT0),NST1(NT0),NST2(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),RFAC(30)
      DIMENSION BENAME(5,NBED),WR(10)
      CHARACTER NC(9)
      CHARACTER(LEN=2) IV,NC2(21)
CGPS      CHARACTER*4 IVNAME
      CHARACTER(LEN=100) IVNAME,tleed4i
C
      common /pdir/tleed4i
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
      COMMON /WIV2/PERSH,NIV,NSE1(30),NSE2(30)
      COMMON /TIT/TITLE(20)
C
      DATA NC(1)/'1'/NC(2)/'2'/NC(3)/'3'/NC(4)/'4'/NC(5)/'5'/
      DATA NC(6)/'6'/NC(7)/'7'/NC(8)/'8'/NC(9)/'9'/
      DATA NC2(1)/'10'/NC2(2)/'11'/NC2(3)/'12'/NC2(4)/'13'/
      DATA NC2(5)/'14'/NC2(6)/'15'/NC2(7)/'16'/NC2(8)/'17'/
      DATA NC2(9)/'18'/NC2(10)/'19'/NC2(11)/'20'/
C
110   FORMAT ('TitleText: ',5A4,' Beam ',5A4,' Rfac',I2,
     & '=',F6.4)
111   FORMAT (' RFAC',I2)
500   FORMAT (/' THE NUMBER OF THEORETICAL AND EXPERIMENTAL BEAMS 
     & DISAGREE ')
!100   FORMAT(20A4)
      IF (NBMAX.NE.NBE) THEN
         WRITE (48,500)
         RETURN
      ENDIF
C
C Find out which Rfactor has been computed
C
      NCOUNT=0
      DO 8 I=1,10
         IF (WR(I).GT.1.e-6) THEN
           NCOUNT=NCOUNT+1 
           NRFAC=I 
         ENDIF
         IF(NCOUNT.GT.1) NRFAC=11
8     CONTINUE
      WRITE(2,111) NRFAC
C
C and record that Rfactor
C
      DO 9 I=1,NBE
         RFAC(I)=0.
         IF (WR(1).GT.1.e-6) RFAC(I)=ROS(I)+RFAC(I)
         IF (WR(2).GT.1.e-6) RFAC(I)=R1(I)+RFAC(I)
         IF (WR(3).GT.1.e-6) RFAC(I)=R2(I)+RFAC(I)
         IF (WR(4).GT.1.e-6) RFAC(I)=RP1(I)+RFAC(I)
         IF (WR(5).GT.1.e-6) RFAC(I)=RP2(I)+RFAC(I)
         IF (WR(6).GT.1.e-6) RFAC(I)=RPP1(I)+RFAC(I)
         IF (WR(7).GT.1.e-6) RFAC(I)=RPP2(I)+RFAC(I)
         IF (WR(8).GT.1.e-6) RFAC(I)=RRZJ(I)+RFAC(I)
         IF (WR(9).GT.1.e-6) RFAC(I)=RMZJ(I)+RFAC(I)
         IF (WR(10).GT.1.e-6) RFAC(I)=RPE(I)+RFAC(I)
         RFAC(I)=RFAC(I)/FLOAT(NCOUNT)
9     CONTINUE
C
C write down the experimental beams(remembering the shift in intensity
C for the Pendry Rfactor)
C
CGPS+
      I=1
888      if (tleed4i(I:(I+6)).NE.'tleed4i') then
             I=I+1
             goto 888
         else
             ilength=I-1
         end if
         
      IV='iv'
      NBEAMS=NBE
      NOUT=48
      NBB=NOUT+NBEAMS-1
      DO 1 I=1,NBE
          IF (NOUT.LE.NBB) THEN
              IF (I.LT.10) IVNAME=tleed4i(1:ilength)//IV//NC(I)
              IF (I.GE.10) IVNAME=tleed4i(1:ilength)//IV//NC2(I-9)
c              IF (I.LT.10) IVNAME='/scratch/zz217/run1/'//IV//NC(I)
c              IF (I.GE.10) IVNAME='/scratch/zz217/run1/'//IV//NC2(I-9)
              OPEN(UNIT=NOUT,FILE=IVNAME,STATUS='UNKNOWN')
Cga              OPEN(UNIT=NOUT,FILE=IVNAME,STATUS='UNKNOWN')
CGSS              OPEN(UNIT=NOUT,FILE='/dev/null',STATUS='UNKNOWN')
C              OPEN(UNIT=NOUT,FILE='/dev/null',STATUS='UNKNOWN')
              WRITE (NOUT,110) (TITLE(II),II=1,5),
     &        (BENAME(II,I),II=1,5),NRFAC,RFAC(I) 
              WRITE (NOUT,*) '"IV exp'
          ENDIF  
C         DO 2 K=1,NEE(I)
          DO 2 K=NSE1(I),NSE2(I)
                E=EE(I,K) 
                ZI=AE(I,K)
C
C shift according to an average intensities for all beams
C
C                IF(NRFAC.EQ.10) ZI=ZI+EEAVE(NBE+1)*PERSH
C
C or to the average intensity of each beam. It should match with the
C choice in YPEND
C
                IF(NRFAC.EQ.10) ZI=ZI+EEAVE(I)*PERSH
            IF (NOUT.LE.NBB) WRITE (NOUT,*) E,ZI
2         CONTINUE
          IF (NOUT.LE.NBB)WRITE (NOUT,*)
          NOUT=NOUT+1
1     CONTINUE
C
C and the theoretical beams after rescaling 
C
      NOUT=48
      DO 3 I=1,NBE
        IBT=IBK(I)
C
C SCAL=(I)/(IBT) if we scale each beam separately (NIV=0)
C SCAL=(NBE+1)/(NBE+1)(NIV=1) for only one scale factor
C
        IF (NIV.EQ.0.AND.EEAVT(IBT).LE.1.E-8) GOTO 3
        SCAL=EEAVE(I)/EEAVT(IBT)
C        SCAL=10000.
        IF (NIV.EQ.1) SCAL=EEAVE(NBE+1)/EEAVT(NBE+1)
        IF (NOUT.LE.NBB) WRITE (NOUT,*) '"IV theory'
C          DO 4 K=1,NET(IBT)
          DO 4 K=NST1(I),NST2(I)
            E=ETH(IBT,K) 
            ZI=AT(IBT,K)*SCAL
C
C Same as above for the shift
C
C            IF(NRFAC.EQ.10) ZI=((ZI/SCAL)+
C     &         EEAVT(NBE+1)*PERSH)*SCAL
            IF(NRFAC.EQ.10) ZI=((ZI/SCAL)+
     &         EEAVT(IBT)*PERSH)*SCAL
            IF (NOUT.LE.NBB) WRITE (NOUT,*) E,ZI
4         CONTINUE
          IF (NOUT.LE.NBB)WRITE (NOUT,*)
          NOUT=NOUT+1
3     CONTINUE
      RETURN
      END
C =========================================================================
C
C Subroutine YPEND calculates the Pendry Y function Y=(A/AP)/((A/AP)**2+VI)
C where AP/A is the logarithmic derivative of the tabulated function A.
C
C =========================================================================
C
      SUBROUTINE YPEND(A,AP,NBD,NB,NE,E,Y,VI,IPR,IEERG,EEAVE)
C
      DIMENSION A(NBD,IEERG),AP(NBD,IEERG),NE(NBD)
      DIMENSION E(NBD,IEERG),Y(NBD,IEERG)
      DIMENSION EEAVE(30)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA,FI
      COMMON /WIV2/PERSH,NIV,NSE1(30),NSE2(30)
C
40    FORMAT (' PENDRY Y FUNCTION IN BEAM',1I3,/,50(5(1F7.2,1E13.4,3X)
     & ,/))
C 
C THEORETICAL AND EXPERIMENTAL ENERGIES IN E() ARE IN EV's
C
      VIEFF=VI*27.21
      PERSH=.05
C      YSH=EEAVE(NB+1)*PERSH
      DO 25 IB=1,NB
C
C WE SHIFT THE INTENSITY OF EACH BEAM BY A PERCENTAGE SHIFT PERSH
C
         YSH=EEAVE(IB)*PERSH
         N=NE(IB)
         IF (N.NE.0) THEN
            DO 20 IE=1,N
               AF=A(IB,IE)+YSH
               IF (AF.LT.1.E-5) THEN
                  APF=AP(IB,IE)
                  IF (ABS(APF).GT.1.E-6) THEN
                     AF=(AF)/APF
                       Y(IB,IE)=1./(AF+VIEFF**2/AF)
                  ELSE
                     Y(IB,IE)=0.
                  ENDIF
               ELSE
                  AF=AP(IB,IE)/(AF)
                  Y(IB,IE)=AF/(1.+AF**2*VIEFF**2)
               ENDIF
20          CONTINUE
         ENDIF
25    CONTINUE
      IF (IPR.GE.2) THEN
         DO 30 IB=1,NB
            N=NE(IB)
            IF (N.NE.0) WRITE (1,40) IB,(E(IB,IE),Y(IB,IE),IE=1,N)
30       CONTINUE
      ENDIF
      RETURN
      END
C=======================================================================
C
C  SUBROUTINE YVAL INTERPOLATES TO GET THE FUNCTION AND ITS FIRST
C  DERIVATIVE (Y and YP resp.)
C
C=======================================================================
C
      SUBROUTINE YVAL(Y,YP,X, WORY, WORX, LENGTH,ITIL,ITIH)
C
      DIMENSION WORY(LENGTH), WORX(LENGTH)
C
      COMMON /DATBLK/ X0,X1,X2,X3,Y0,Y1,Y2,Y3
C
C  IF FEWER THAN FOUR GRID POINTS AVAILABLE, USE 2ND OR 1ST ORDER
C  POLYNOMIAL INTERPOLATION
C
      IF (LENGTH.LT.4) GO TO 10
C
C  FIND REQUIRED INTERPOLATION INTERVAL
C
      CALL BINSRX(IL, IH, X, WORX, LENGTH)
C
C  SKIP NEXT STEP IF SAME INTERVAL IS FOUND AS LAST TIME
C
      IF(IL .EQ. ITIL .AND. IH .EQ. ITIH) GO TO 5
C
C  FIND FOUR NEAREST GRID POINTS AND CORRESPONDING INTENSITIES
C  FOR 3RD-ORDER POLYNOMIAL INTERPOLATION
C
      CALL STFPTS(IL, IH, WORX, WORY, LENGTH)
C
C  DO ACTUAL 3RD-ORDER POLYNOMIAL INTERPOLATION
C
5     CALL XNTERP(Y,YP,X)
      ITIH = IH
      ITIL = IL
      GO TO 20
10    CALL XNTRP2(X,Y,YP,WORX,WORY,LENGTH)
20    RETURN
      END
