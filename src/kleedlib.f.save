
        SUBROUTINE ANGLE_2D ( p1, p2, p3, angle_rad_2d )

c*******************************************************************************
c
c! ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
c
c  Discussion:
c
c    Except for the zero angle case, it should be true that
c
c      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
c
c        P1
c        /
c       /    
c      /     
c     /  
c    P2--------->P3
c
c  Modified:
c
c    15 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), define the rays
c    P1 - P2 and P3 - P2 which define the angle.
c
c    Output, real ( kind = 8 ) ANGLE_RAD_2D, the angle swept out by the rays,
c    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
c    length, then ANGLE_RAD_2D is set to 0.
c
c  implicit none

        PARAMETER (PI=3.141592653589793D+00)

        REAL      P,P1,P2,P3
        DIMENSION P(2),P1(2),P2(2),P3(2)  




        p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) + 
     &  ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )


        p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) - 
     &         ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )


        IF ((p(1).EQ.0).AND.(P(2).EQ.0)) THEN
           angle_rad_2d = 0.0D+00
         ELSE
        angle_rad_2d = atan2 ( p(2), p(1) )
        ENDIF

        IF (ANGLE_RAD_2D.LT.0.0D+00) THEN
            angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
        ELSE
        ENDIF

        return
        end
C =========================================================================
C
C Subroutine AVEINT computes the average intensity of the peaks in each
C beam and returns it as EEAVE 
C
C =========================================================================
C
      SUBROUTINE AVEINT(A,AP,NBD,NE1,NE2,NB,E,IEERG,EEAVE,NFL)
C
      DIMENSION A(NBD,IEERG),E(NBD,IEERG)
      DIMENSION NE1(NBD),NE2(NBD),AP(NBD,IEERG),EEAVE(30)
C    
C we go through the derivative of each beam and record the intensity
C when it becomes zero and it corresponds to a peak
C
       EEAVE(NB+1)=.0
       DO 10 I=1,NB
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
               EAV=EAV+A(I,J)
               ICOUNT=ICOUNT+1
           ENDIF
100      CONTINUE 
         FL=FLOAT(ICOUNT)
         IF (NFL.EQ.1.AND.FL.EQ.0) THEN 
c             EEAVE(I)=.0001
             EEAVE(I)=1.e-11
         ELSE
             EEAVE(I)=EAV/FL
         ENDIF
         EEAVE(NB+1)=EEAVE(NB+1)+EEAVE(I)
10     CONTINUE
       EEAVE(NB+1)=EEAVE(NB+1)/FLOAT(NB)
       RETURN
       END

C============================================================================
C
C  Subroutine BEAMT selects those beams from the input list that are
C  needed at the current energy, based on the parameter TST which limits
C  the decay of plane waves from one layer to the next (the interlayer
C  spacing must already be incorporated in TST, which is done by subroutine
C  READT)
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
C  NP            =   LARGEST NO. OF BEAMS SELECTED IN ANY ONE BEAM SET
C  IPR           =   PRINT CONTROL PARAMETER
C
C Modified version of routine BEAMS from the VAN HOVE/TONG LEED package.
C Modificatiosn by WANDER.
C
C============================================================================
C
      SUBROUTINE BEAMT(KNBS,KNB,SPQ,SPQF,KNT,AK2,AK3,E,TST,NB,PQ,PQF,
     & NT,NP,IPR)
C
      DIMENSION KNB(KNBS),SPQ(2,KNT),SPQF(2,KNT),NB(KNBS)
      DIMENSION PQ(2,KNT),PQF(2,KNT)
C
705   FORMAT (1H0,1I3,13H BEAMS USED  ,8I4)
706   FORMAT (1H ,10(9(2X,2F6.3),/))
C
      KNBJ=0
      NT=0
      NP=0
      DO 704 J=1,KNBS
         N=KNB(J)
         NB(J)=0
         DO 703 K=1,N
            KK=K+KNBJ
            FACT1=(AK2+SPQ(1,KK))*(AK2+SPQ(1,KK))
            FACT2=(AK3+SPQ(2,KK))*(AK3+SPQ(2,KK))
            IF ((2.0*E-FACT1-FACT2)+TST.GE.0) THEN
               NB(J)=NB(J)+1
               NT=NT+1
               DO 702 I=1,2
                  PQ(I,NT)=SPQ(I,KK)
                  PQF(I,NT)=SPQF(I,KK)
702            CONTINUE
            ENDIF
703      CONTINUE
         KNBJ=KNBJ+KNB(J)
         NP=MAX0(NP,NB(J))
704   CONTINUE
      IF (IPR.GT.0) THEN
         WRITE (1,705) NT,(NB(J),J=1,KNBS)
         WRITE (1,706) ((PQF(I,K),I=1,2),K=1,NT)
        do 24 k=1,nt
         WRITE (25,*) E,K,PQF(1,K),PQF(2,K)
24        continue
      ENDIF
      RETURN
      END

C======================================================================
C
C  Subroutine BEMGEN generates the list of beams for the TENSOR LEED
C  program. The beams are ordered by increasing ABS(G) and grouped
C  by subsets induced by the superlattice.
C 
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
C AUTHOR: VAN HOVE/Modified by WANDER
C
C======================================================================
C
      SUBROUTINE BEMGEN(TST,EF,SPQF,SPQ,KNBS,KNB,RAR1,RAR2,KNT,IPR,TVA,
     & DFLAG,KNBMAX,G)
C
      INTEGER LATMAT(2,2),KNB(20),DFLAG
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2),SPQF(2,KNBMAX)
      DIMENSION SPQ(2,KNBMAX),RAR1(2),RAR2(2),RBR1(2),RBR2(2)
      DIMENSION ALMR(2,2),G(2,KNBMAX),ST(4),DG(2)
C
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
1001  FORMAT (' BEAMGENERATION CODE')
1002  FORMAT (' ===================')
165   FORMAT (28H LATMAT HAS ZERO DETERMINANT)
270   FORMAT (23H MORE BEAMS THAN KNBMAX)
C
      IF (IPR.GT.0) THEN
         WRITE (1,1001)
         WRITE (1,1002)
      ENDIF
C
C GENERATE LATICE MATRIX
C
      IF (DFLAG.EQ.0) THEN
         CALL MATGEN(ARA1,ARA2,ARB1,ARB2,LATMAT)
      ELSE
         CALL MATGEN(ARA1,ARA2,ARA1,ARA2,LATMAT)
      ENDIF


C
C SET UP CONSTANTS 
C
      PI=3.1415926535
C
C  ALMR IS MATRIX RELATING SUBSTRATE TO OVERLAYER LATTICES
C

      DET=LATMAT(1,1)*LATMAT(2,2)-LATMAT(1,2)*LATMAT(2,1)
cjcm      write(*,*) 'BEMGEN: DET = ', DET        
      IF (ABS(DET).GE.1.E-5) THEN
         ALMR(1,1)=FLOAT(LATMAT(2,2))/DET
         ALMR(2,2)=FLOAT(LATMAT(1,1))/DET
         ALMR(1,2)=-FLOAT(LATMAT(2,1))/DET
         ALMR(2,1)=-FLOAT(LATMAT(1,2))/DET
         GMAX2=2.*EF/27.21+TST
C
C  GENERATE ALL BEAMS WITHIN A BEAM CIRCLE OF RADIUS SQRT(GMAX2),
C  LIMITED TO KNBMAX BEAMS
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
               IF (KNT.GE.KNBMAX) GOTO 260
               NOP=1
               IF (NI1.EQ.0) NIT=NI2+4
               G(1,KNT)=GT1
               G(2,KNT)=GT2
               SPQF(1,KNT)=ALMR(1,1)*II1+ALMR(2,1)*II2
               SPQF(2,KNT)=ALMR(1,2)*II1+ALMR(2,2)*II2
            ENDIF
280      CONTINUE
 290         IF (NI2.LE.NIT) GOTO 190
         IF (NOP.EQ.1) GOTO 180
C
C  ORDER BEAMS BY INCREASING ABS(G)
C
cjcm         write(*,*) 'order beams by abs(g)'
         KNT1=KNT-1
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
C  ORDER BEAMS BY BEAM SET
C
cjcm         write(*,*) 'BEMGEN: Order beams by beam set'
         
         TWPI=2.*3.1415926535
         I=1
         KNBS=1
330      CONTINUE
         KNB(KNBS)=1
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
                        I2=I+2
                        DO 340 KK=I2,J
                           K=J+I2-KK
                           SPQF(1,K)=SPQF(1,K-1)
                           SPQF(2,K)=SPQF(2,K-1)
                           G(1,K)=G(1,K-1)
                           G(2,K)=G(2,K-1)
340                     CONTINUE
                        G(1,I+1)=ST(1)
                        G(2,I+1)=ST(2)
                        SPQF(1,I+1)=ST(3)
                        SPQF(2,I+1)=ST(4)
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
         DO 10 J=1,KNT
            DO 20 I=1,2
               SPQ(I,J)=SPQF(1,J)*RAR1(I)+SPQF(2,J)*RAR2(I)
20          CONTINUE
10       CONTINUE
         RETURN
         
 260         CONTINUE
cjcm         write(*,*) 'BEMGEN: at 260'
         WRITE (1,270)
      ELSE
         WRITE (1,165)
      ENDIF
cjcm      write(*,*) 'BEMGEN: Should not be here'
        STOP
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
102   IL = LENGTH - 1
      IH = LENGTH
      GO TO 110
105   IL = 1
      IH = 2
110   RETURN
      END

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
      DOUBLEPRECISION FAC(100),a,b,bn,c,cn,blmt
C
      COMMON /F/FAC
C
40    FORMAT (28H INVALID ARGUMENTS FOR BLMT ,6(I3,1H,))
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
                              A=((-1.0)**IA1)/FAC(IA3+1)*FAC(IA2+1)/FAC
     &                         (IA6+1)*FAC(IA4+1)
                              A=A/FAC(IA7+1)*FAC(IA5+1)/FAC(IA8+1)*FAC
     &                         (IS+1)/FAC(IA9+1)
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
                                 BN=SIGN/FAC(IT)*FAC(IB1+1)/FAC(IB3+1)
     &                            *FAC(IB2+1)
                                 BN=BN/FAC(IB4+1)/FAC(IB5+1)
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
                              C=CN/FAC(IC2+1)*FAC(IC1+1)/FAC(IC4+1)*FAC
     &                         (IC3+1)
                              C=C/FAC(IC6+1)*FAC(IC5+1)
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

C=========================================================================
C
C  Subroutine CAAA computes Clebsch-Gordan coefficients for use by
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
      SUBROUTINE CAAA(CAA,NCAA,LMMAX)
C
      DIMENSION CAA(NCAA)
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
            IF (I.EQ.J) THEN
               IF (M1.GT.0) GOTO 50
            ELSEIF (L1.LT.L2.OR.(L1.EQ.L2.AND.(IABS(M1).LT.IABS(M2)))) 
     &       THEN
               GOTO 50
            ENDIF
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
      DOUBLEPRECISION yint,del,pref,prefa,prefb,preft,sum
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
      DOUBLEPRECISION FAC(100),blmt,a,b
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
        INTEGER IBT
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
      DOUBLEPRECISION F(100),SUM,A
C
      COMMON /F/F
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
                        SUM=SUM-SIGN*F(IA1)*F(IA4)/(F(IA2)*F(IA3)*F
     &                   (IA5)*F(IT))
450                  CONTINUE
                     IA1=2+IS-IM1
                     IA2=2+IS-IM2
                     IA3=2+IS-IM3
                     IA4=3+2*(IS-IM3)
                     A=-(-1.0)**(IS-IM2)*F(IA4)*F(IS+1)*F(IM3)*SUM/(F
     &                (IA1)*F(IA2)*F(IA3)*F(2*IS+2))
                  ENDIF
               ENDIF
               PPP(I1,I2,I3)=A
460         CONTINUE
462      CONTINUE
461   CONTINUE
      RETURN
      END


C Subroutine DELXGENTFOL generates the composite layer contribution to the
c outcoming plane wave amplitude in the KINEMATIC LIMIT
C
C Input Parameters:
C =================
C
C E,VPIS          =  Current (complex) energy
C NLAY            =  Number of composite layer (CL) subplanes that are
C                   displaced from reference position
C AK2,AK3,AAJ     = Parallel and perpendicular components of momentum for
C                    incident direction
C AK2M,AK3M,AAK   = Parallel and perpendicular components of momentum for
C                   exit beams 
C ILEMOLTF          = Element type of atoms in NLAY composite layers
C WPOSTF          = Atomic positions of all subplanes in composite layers
C TVA                    = Area of substrate unit cell
C ADISP            =  displacements of overlayer atoms from reference positions 
C F                  = Atomic form-factors
C
C Output Parameters:
C ==================

C XISTF           = Plane wave amplitude of exit beam

C Author: Garcia-Lekue
C =========================================================================
C
      SUBROUTINE DELXGENTFOL_TEMP(NLAY,E,VPIS,AK2,AK3,AK2M,AK3M,
     & NL1,NL2,IELEMOLTF,AAK,AAJ,WPOSTF,TVA,ADISP,XISTF,F)
C
      DIMENSION C(3),ADISP(NLAY,3)
      DIMENSION IELEMOLTF(NLAY)
      DIMENSION POS(3)
      DIMENSION WPOSTF(100,3)
      COMPLEX CI
      COMPLEX AAK,AAJ
      COMPLEX ATOYF,XISTF
      COMPLEX CDELTAKZ,CDELTAKPOS
      COMPLEX F(3)

      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
cjcm        write(*,*) 'entering delxgentfol'

C
C Set constants
C
      NL=NL1*NL2
      PI=4.0*ATAN(1.0)
      CI=CMPLX(0.0,1.0)

C Calculate difference between incident and outgoing beams
        CDELTAKZ=AAK-AAJ
        DELTAK2=AK2-AK2M
        DELTAK3=AK3-AK3M
C
C Initialize plane wave amplitude
         XISTF=CMPLX(0.0,0.0)
C
C Begin loop over subplanes of composite layer
C
            IELEM=1
         DO 100 IJK=1,NLAY

cjcm            write(*,*) 'inside do loop, ijk = ', ijk, ' nlay = ', nlay

C Assign element type to each sublayer 
            IEL=IELEMOLTF(IELEM)
cjcm            write(*,*) 'iel = ', iel

C WPOSTF gives the reference structure, i.e., the input atomic
C positions in unit cell 
C ADISP gives the displacements of overlayer atoms from
C initial positions

cjcm        write(*,*) '120 do loop'
        DO 120 K=1,3
           C(K)=ADISP(IJK,K)/0.529
120        CONTINUE

        POS(1)=WPOSTF(IJK,1) + C(1) 
        POS(2)=WPOSTF(IJK,2) + C(2) 
        POS(3)=WPOSTF(IJK,3) + C(3) 

cjcm        write(*,*) 'compute cdeltakpos'
        CDELTAKPOS=(CDELTAKZ*POS(1))+(DELTAK2*POS(2))+
     1                  (DELTAK3*POS(3))


cjcm        write(*,*) 'compute atoyf'
        ATOYF=F(IEL)*CEXP(CI*CDELTAKPOS)
        ATOYF=ATOYF*(-CI/((AAJ)*TVA*FLOAT(NL)))

        XISTF=XISTF+ATOYF

         IELEM=IELEM+1
 100        CONTINUE

      RETURN
      END

C
C Subroutine DELXGENTFSL generates the substrate layer contribution to the
c outcoming plane wave amplitude in the KINEMATIC LIMIT
C
C Input Parameters:
C =================
C
C E,VPIS          =  Current (complex) energy
C AK2,AK3,AAJ     = Parallel and perpendicular components of momentum for
C                   incident direction
C AK2M,AK3M,AAK   = Parallel and perpendicular components of momentum for
C                   exit beams
C NTAU            = Total number of chemical elements
C POSTF           = Positions of substrate atoms in each substrate layer
C TVA             = Area of substrate unit cell
C FS              = Atomic form-factor for substrate atoms
C
C Output Parameters:
C ==================
C
C XISTF           = Plane wave amplitude of exit beam
C
C Author: Garcia-Lekue
C =========================================================================
C
      SUBROUTINE DELXGENTFSL_TEMP(E,VPIS,NLAY,
     & AK2,AK3,AK2M,AK3M,
     & NL1,NL2,AAK,AAJ,TVA,POSTF,XISTF,FS)
C
      DIMENSION C(3)
      DIMENSION POS(3),POSTF(30,3)
      COMPLEX CAPPA,CI
      COMPLEX FS
      COMPLEX AAK,AAJ
      COMPLEX ATOYF,XISTF
      COMPLEX CDELTAKZ,CDELTAKPOS

      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)


C
C Set constants
C
      NL=NL1*NL2
      PI=4.0*ATAN(1.0)
      CI=CMPLX(0.0,1.0)


C Calculate difference between incident and outgoing beams
        CDELTAKZ=AAK-AAJ
        DELTAK2=AK2-AK2M
        DELTAK3=AK3-AK3M
C Initialize plane wave amplitude
         XISTF=CMPLX(0.0,0.0)
         ATOYF=CMPLX(0.0,0.0)
C
C Begin loop over planes of substrate
C
         DO 100 IJK=1,NL
C
C Input atomic positions in unit cell (one atom per subplane)
C
        POS(1)=POSTF(IJK,1)
        POS(2)=POSTF(IJK,2)
        POS(3)=POSTF(IJK,3)

        CDELTAKPOS=((CDELTAKZ*POS(1))+(DELTAK2*POS(2))+
     1                  (DELTAK3*POS(3)))
        ATOYF=FS*CEXP(CI*CDELTAKPOS)
        ATOYF=ATOYF*(-CI/((AAJ)*TVA*FLOAT(NL)))

        XISTF=XISTF+ATOYF
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
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /REXP/EEINCR
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
      COMMON /TIT/TITLE(20)
C
100   FORMAT(20A4)
110   FORMAT(20I3)
120   FORMAT(/' ERROR IN READE:MORE EXPERIMENTAL BEAMS THAN HAVE BEEN 
     & DIMENSIONED')
CAGL hurrengo lerroa
50    FORMAT (30HAVERAGING SCHEME OF EXP. BEAMS,5(25I3,/))

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
      CALL READE(AE,EE,NBED,NEE,NBEA,BENAME,IPR,IEERG,NBMAX)
C
C  IF AN IRREGULAR INPUT ENERGY GRID WAS USED, FIRST INTERPOLATE ALL
C  EXP. DATA TO THE A GRID specified by the input EEINCR
C
      IF (IRGEXP.EQ.1) THEN
        CALL INTPOL(AE,AEP,AEPP,NBED,NEE,NBED,EE,
     %    EEINCR,IPR,XPL,YPL,IEERG)
      ENDIF
C
C  AVERAGE DATA FROM DIFFERENT EXPERIMENTS AND ORDER BY INCREASING ENERGY
C
      CALL EXPAV(AE,EE,NBED,NEE,BENAME,NBEA,NBE,IPR,XPL,NNN,IEERG)
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
         ENDIF
C
C  INTERPOLATE EXP. DATA TO WORKING GRID (MULTIPLES OF EINCR EV), UNLESS
C  DONE BEFORE
C
         IF (IRGEXP.EQ.0) THEN
          CALL INTPOL(AE,AEP,AEPP,NBED,NEE,NBE,EE,EINCR,
     &    IPR,XPL,YPL,IEERG)
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
         CALL AVEINT(AE,AEP,NBED,NEE,NEE,NBE,EE,IEERG,EEAVE,0)
C
C  PRODUCE (1ST AND) 2ND DERIVATIVES OF EXP. DATA
C
C         CALL DERL(AE,NEE,NBED,NBE,AEP,EINCR,IEERG)
C         CALL DERL(AEP,NEE,NBED,NBE,AEPP,EINCR,IEERG)
C
C  PRODUCE PENDRY Y FUNCTION FOR EXP. DATA
C
         CALL YPEND(AE,AEP,AEPP,NBED,NBE,NEE,EE,YE,
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
110   FORMAT (48H EXP. ENERG. AND INTENS. AFTER AVERAGING IN BEAM,1I4,
     & /,50(5(1F7.2,1E13.4,3X),/))
30    FORMAT (5(25I3,/))

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

C=========================================================================
C
C  Function FACT computes the factorial using an asymptotic expansion.
C  Modified version of function FACT from the Van Hove/Tong LEED package.
C  Modifications by WANDER.
C
C=========================================================================
C
      FUNCTION FACT(L)
C
      DOUBLEPRECISION DFACT,X
C
      IF (L.GT.4) THEN
         X=L+1
         DFACT=DEXP(-X)*(10.0D0**((X-0.5D0)*DLOG10(X)-(X-1.0D0)))*
     &    (DSQRT(6.283185307179586D0))*(1.0+(1.0/(12.0*X))+(1.0/
     &    (288.0D0*(X**2)))-(139.0D0/(51840.0D0*(X**3)))-(571.0D0/
     &    (2488320.0D0*(X**4))))
         FACT=SNGL(DFACT)
      ELSE
         IF (L.EQ.0) FACT=1.0
         IF (L.EQ.1) FACT=0.1
         IF (L.EQ.2) FACT=0.02
         IF (L.EQ.3) FACT=6.0*0.001
         IF (L.EQ.4) FACT=24.0*0.0001
      ENDIF
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
      SUBROUTINE GEOMV(IP,VPOS,NLAY)
C
      DIMENSION VPOS(60,3)
C
      COMMON /ZMAT/IANZ(40),IZ(40,4),BL(40),ALPHA(40),BETA(40),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(40),ALPHAS(40),BETAS
     & (40),PHIR,PHIM1,PHIM2
C
1000  FORMAT (//,
     & 57H PARAMETER OUT OF RANGE OF DEFINED Z MATRIX  $$$ STOP $$$)
1001  FORMAT (9H0NATOMS(=,I4,12H) .NE.NLAY(=,I4,
     & 21H) IN GEOMV  *********)
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
      CALL GEOMXY(IP,VPOS,NLAY)
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
      DIMENSION C(60,3),R(3),U(3)
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
      SUBROUTINE GEOMXY(IP,C,NLAY)
C
      DIMENSION ALPHA(60),BETA(60)
      DIMENSION A(60),B(60),CZ(60,3),D(60),U1(3),U2(3),U3(3),U4(3)
      DIMENSION VJ(3),VP(3),V3(3)
      DIMENSION C(60,3),DIS(60)
C
      COMMON /ZMAT/IANZ(60),IZ(60,4),BL(60),ALPH(60),BET(60),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(60),ALPHAS(60),BETAS
     & (60),PHIR,PHIM1,PHIM2
C
      DATA ZERO/0.0/,ONE/1.0/,TWO/2.0/
      DATA TENM5/1.0E-5/,TENM6/1.0E-6/
      DATA TORAD/1.74532925199E-02/
      DATA PI/3.141592654/
C
1000  FORMAT (//,15X,11HCOORDINATES,/,13X,1HX,14X,1HY,14X,1HZ)
1010  FORMAT (1X,1I4,2X,F10.5,2(5X,F10.5))
1020  FORMAT (//,30X,8HZ MATRIX//1X,9HIN IN1 AN,3X,2HZ1,4X,2HBL,19X,
     & 2HZ2,2X,5HALPHA,17X,2HZ3,3X,4HBETA,17X,2HZ4,/)
1030  FORMAT (3I3)
1040  FORMAT (3I3,I5,G14.7,1H(,I3,1H))
1050  FORMAT (3I3,I5,G14.7,1H(,I3,4H)   ,I5,G14.7,1H(,I3,1H))
1060  FORMAT (3I3,I5,G14.7,1H(,I3,4H)   ,I5,G14.7,1H(,I3,4H)   ,I5,
     & G14.7,1H(,I3,4H)   ,I5)
5000  FORMAT (11X,15I7)
5003  FORMAT (5X,I3,5X,15F7.3)
5001  FORMAT (//5X,15HDISTANCE MATRIX,/)
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

C=============================================================================
C
C Subroutine INTENTF generates the plane wave amplitudes and the 
C intensities of exit beams in the KINEMATIC LIMIT
C Theoretical IV spectra are printed out in ivth# files
C Experimental IV spectra are printed out in ivexp# files                         
C
C Input parameters: 
C =================
C
C NLAY           =  Number of composite layer (CL) subplanes that are
C                   displaced from reference position 
C DISP                  =  Displacement of CL atoms from their reference position
C PSQ, PQFEX     =  Labels of exit beams in reciprocal space
C NTAU                 =  Total number of chemical elements  
C NT0                  =  Number of exit beams
C PHSSEL        =   Complex phase shift
C EI,EF,DE         =  Initial energy point, final energy point and energy   
C                   step within the calculation
C NL1, NL2       =  Superlattice characterization codes
C ILEMOLTF       = Element type of atoms in NLAY composite layers
C WPOSTF           =  Atomic positions of all subplanes in composite layers
C TVA                 =  Area of substrate unit cell
C SPOSTF1        =  Atomic positions of atoms on first substrate layer 
C ASA                 =  Substrate interlayer vectors
C INVECT         =  Number of substrate interlayer vectors
C ITEMP           = Parameter controlling temperature effect
C 
C Output parameters:
C ==================
C
C XIST           =  Total plane wave amplitudes for the current structure
C AT             =  Intensities generated from amplitudes
C COMMON BLOCKS
C =============
C
C VMIN,VMAXM,DV  =  Range of search over inner potential (block VINY)
C EINCR          =  Grid step to be used after interpolation (block VINY)
C VV, VPIS       =  Real and imaginary parts of inner potential (block ENY)
C
C Author: Garcia-Lekue  
C ==========================================================================
C
CAGL      SUBROUTINE INTENTF(NLAY,ADISP,PSQ,NTAU,NT0,
       FUNCTION VINTENTF(NLAY,ADISP,PSQ,NTAU,NT0,
     & PHSSEL,EI,EF,DE,
     & NL1,NL2,IELEMOLTF,WPOSTF,TVA,SPOSTF1,PQFEX,ASA,INVECT,
     & INBED,IEERG,AE,EE,NEE,NBEA,BENAME,IPR,AP,APP,YE,
     & SE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,WR,WB,IBP,NERG,L1,ITEMP)
C
      DIMENSION ADISP(NLAY,3)
      DIMENSION PSQ(2,NT0),PQFEX(2,NT0)
      DIMENSION AT(NT0,NERG)  
      DIMENSION AT2(NT0,IEERG)  
      DIMENSION ETH(NT0,IEERG)
	DIMENSION ASA(10,3)
cjcm need to add dimension for missing parameters to rfactf
	DIMENSION AP(INBED,IEERG), APP(INBED,IEERG), YE(INBED,IEERG)
      DIMENSION AE(INBED,IEERG),EE(INBED,IEERG),NEE(INBED),YPL(IEERG)
      DIMENSION NBEA(INBED),BENAME(5,INBED),XPL(IEERG),NNN(IEERG)
      DIMENSION TSE(INBED),TSE2(INBED),TSEP(INBED),TSEP2(INBED)
      DIMENSION TSEPP(INBED),WR(10),WB(NT0),TSEPP2(INBED)
      DIMENSION TSEY(INBED), TSEY2(INBED),IBP(NT0)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
      DIMENSION LPS2(NLAY),IELEMOLTF(NLAY)
      DIMENSION POSTF(30,3), WPOSTF(100,3)
      DIMENSION SPOSTF1(30,3)
      DIMENSION SPOSTF(20,30,3)
      DIMENSION PLGND(L1)
      COMPLEX F(3),FS,FJ
      COMPLEX AAK,AAJ,PRE,CI
      COMPLEX XIST(NT0,NERG)
      COMPLEX PHSSEL(NERG,NTAU,L1),PHSS(NTAU,L1)
      COMPLEX XISTF,XISTFOL,XISTFS,XISTFSTOT,XISTFS1
      COMPLEX CAPPA
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
      COMMON /ADS/ASL,FR,ASE,VPIS,VPIO,VO,VV
      COMMON /RPL/DVOPT
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /REXP/EEINCR
      COMMON /WIV/NBMAX,EEAVE(30),EEAVT(30)
      COMMON /TIT/TITLE(20)
      COMMON /POW/IFUNC,MFLAG,SCAL
      COMMON /WIV2/PERSH,NIV,NSE1(30),NSE2(30)
      COMMON /NSTR/VOPT,NNST,NNSTEF
C
1000  FORMAT (/,10X,25HCOORDINATES AFTER SORTING,/,13X,1HX,14X,1HY,14X,
     & 1HZ)
1010  FORMAT (7X,F10.7,2(5X,F10.7))
500   FORMAT (/' DIRECTION SET SEARCH COMPLETE ')
501   FORMAT (' ======================= ')
502   FORMAT (' DIMENSIONALITY OF SEARCH=',I5)
503   FORMAT (' ISMOTH=',I5)
505   FORMAT (/' CONVERGENCE TOLERANCE ACHIEVED =',F7.6)
506   FORMAT (/' COORDINATES AT MINIMUM;  LSFLAG'/)
513   FORMAT (' NUMBER OF FUNCTIONAL EVALUATIONS =',I5)
514   FORMAT (/' OPTIMUM R-FACTOR =',1F7.4)
515   FORMAT (' OPTIMUM VALUE OF INNER POTENTIAL =',1F7.4)
516   FORMAT (/' NUMBER OF ITERATIONS =',1I3)
510   FORMAT (3F12.8,I4)


cjcm        write(*,*) 'entering  vintentf'
C
C Set constants
C
        PI=4.0*ATAN(1.0)
        CI=CMPLX(0.0,1.0)

        ANORM1=1.0E-6
c        ANORM1=1.0E-16


C
C Find maximum X displacement to correct for perpendicular displacements
C of the surface barrier (program will place surface barrier a distance
C ASE above the top atom). The Following is not correct 
C when there is more than one atom per unit cell in the outermost layer
C
      AMAXD=100.0
      DO 101 I=1,1
        IF(ADISP(I,1).LE.AMAXD)AMAXD=ADISP(I,1)
101   CONTINUE

C
C Begin loop over energies.
C
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
         AKZ=SQRT(2.0*(E-VV))*COS(THETA(1))
         AK=SQRT(2.0*(E-VV))*SIN(THETA(1))
         AK2=AK*COS(FI(1))
         AK3=AK*SIN(FI(1))
C
C Compute scaling needed to correct for shift in surface barrier
C (This is for incident beam)
C
c 
        AAK=CMPLX(2.0*(E)-AK2*AK2-AK3*AK3,-2.0*VPIS+0.000001)
        AAK=CSQRT(AAK)
C
C Set up indices of exit beams in units of reciprocal lattice vectors
C
        NDI=1
        NSH=NT0
        DO 110 J=1,NT0
           XIST(I,J)=CMPLX(0.0,0.0)
           IF(J.GT.NSH) THEN
              NDI=NDI+1
              NSH=NSH+NT0
              AKZ=SQRT(2.0*(E-VV))*COS(THETA(NDI))
              AK=SQRT(2.0*(E-VV))*SIN(THETA(NDI))
              AK2=AK*COS(FI(NDI))
              AK3=AK*SIN(FI(NDI))
              AAK=CMPLX(2.0*(E-VV)-AK2*AK2-AK3*AK3,-2.0*VPIS+0.000001)
              AAK=CSQRT(AAK)
           ENDIF
C
C Set parallel and perpendicular componenets of momentum for each 
C exit beam
C
            ETH(J,I)=EEV

            AK2M=AK2+PSQ(1,J)
            AK3M=AK3+PSQ(2,J)
C
C Compute scaling needed to correct for shift in surface barrier
C (This is for exit direction)
C
           AAJ=CMPLX(2.0*(E)-AK2M*AK2M-AK3M*AK3M,-2.0*VPIS+0.000001)
           AAJ=-CSQRT(AAJ)
C
C Does this beam emerge?
C
cTF            AAA=2.0*(E)-AK2M*AK2M-AK3M*AK3M
            AAA=2.0*(E-VV)-AK2M*AK2M-AK3M*AK3M
           IF(AAA.GT.0)THEN
             AAB=SQRT(AAA)

C Phase shift for each energy value
C
        DO 33 II=1,NTAU
           DO 34 LL=1,L1
                 PHSS(II,LL)=PHSSEL(I,II,LL)
34        CONTINUE
33      CONTINUE

C Calculate modulus of  wave vector of incident beam
        CAPPA=CMPLX(2.0*E,-2.0*VPIS+0.000001)
        CAPPA=CSQRT(CAPPA)
c Calculate the scattering angle between incoming and outgoing beams
c NOTE THAT WE TAKE REAL IN&OUT MOMENTA!!!

        DRAAK=REAL(AAK)
        DRAAJ=REAL(AAJ)
        PRODKKP=(DRAAK*DRAAJ)+(AK2*AK2M)+(AK3*AK3M)
        CMODK=SQRT(DRAAK**2+AK2**2+AK3**2)
        CMODKP=SQRT(DRAAJ**2+AK2M**2+AK3M**2)

        COSKKP=PRODKKP/(CMODK*CMODKP)
        ARCOS=ACOS(COSKKP)

C Calculate Legendre Polynomials

        XL=COSKKP
        LMAX=L1-1
        CALL LGND(LMAX,XL,PLGND)

C Calculate zero order form factor for all atom types
        DO 25 IEL=1,NTAU
            F(IEL)=CMPLX(0.0,0.0)
        IF (ITEMP.EQ.0) THEN
           F(IEL)=(-2.0*PI/CABS(CAPPA))*CSIN(PHSS(IEL,1))*
     1          CEXP(CI*PHSS(IEL,1))
        ELSE
        DO 26 LL=1,L1
            FJ=(-2.0*PI/CABS(CAPPA))*CSIN(PHSS(IEL,LL))*
     1          CEXP(CI*PHSS(IEL,LL))*
     1          ((2.0*LL)-1.0)*PLGND(LL)
            F(IEL)=F(IEL)+FJ
26      CONTINUE
        ENDIF
25      CONTINUE
C The substrate form factor is given by IEL=NTAU 
        FS=CMPLX(0.0,0.0)
        FS=F(NTAU)
       
C
C Work out scaling needed to account for shift in surface barrier
C (PRE=(inward contribution + outward contribution for this beam)*shift)
C
              PRE=(AAK-AAJ)*AMAXD
              PRE=CEXP(-CI*PRE)
C
C Initialize intensity of exit beam (J) at energy E (I)
C
             AT(J,I)=0.

c OVERLAYER (OL)
c --------------
c Initialize OL contribution to diffracted amplitude
              XISTFOL=CMPLX(0.0,0.0)

cjcm              write(*,*) 'before delxgentfol'

             CALL DELXGENTFOL_TEMP(NLAY,E,VPIS,
     &                AK2,AK3,AK2M,AK3M,
     &               NL1,NL2,IELEMOLTF,AAK,AAJ,WPOSTF,TVA,ADISP,
     &               XISTFOL,F)
cjcm              write(*,*) 'after delxgentfol'


c SUBSTRATE (SL)
c --------------

C Initialize SL contribution to diffracted amplitude
            XISTFSTOT=CMPLX(0.0,0.0)
            DELXISLTOT=CMPLX(0.0,0.0)

C Include only integer beams
        AMOD1=AMOD(2*PQFEX(1,J),2.0)        
        AMOD2=AMOD(2*PQFEX(2,J),2.0)
        IF ((AMOD1.NE.0.).OR.(AMOD2.NE.0.)) THEN         
        GOTO 340
        ELSE

c Start loop over substrate layers
c The 1st SL is always included
             ANORM2=0.0
            XISTFS1=CMPLX(0.0,0.0)
c In each SL there are NL substrate atoms, and their positions
c are obtained by translating the substrate atoms in the 1st layer
c (SPOSTF1) along +X axis, using the interlayer vector ASA
             NL=NL1*NL2
        IS2=1
        DO 730        INL=1,NL
         DO 731 IJ=1,3
            SPOSTF(IS2,INL,IJ)=SPOSTF1(INL,IJ)
731        CONTINUE
730        CONTINUE

        DO 756 INL=1,NL
                DO 757 IJ=1,3
                  POSTF(INL,IJ)=SPOSTF(IS2,INL,IJ)
757      CONTINUE
756      CONTINUE


cjcm        write(*,*) 'before delxgentfsl'
        CALL DELXGENTFSL_TEMP(E,VPIS,NLAY,AK2,AK3,AK2M,AK3M,
     &      NL1,NL2,AAK,AAJ,TVA,POSTF,XISTFS1,FS)

         XISTFSTOT=XISTFS1

c Now start summing up subsequent substrate layers

380        CONTINUE

        IS2=IS2+1
         XISTFS=CMPLX(0.0,0.0) 
c Single interlayer vector
        IF (INVECT.EQ.1) THEN
         DO 745 INL=1,NL
            DO 748 IJ=1,3
             SPOSTF(IS2,INL,IJ)=ASA(1,IJ)+SPOSTF(IS2-1,INL,IJ)
748     CONTINUE
745     CONTINUE
c Two interlayer vectors
         ELSEIF (INVECT.EQ.2) THEN
            IF(MOD(IS2,2).EQ.0.) THEN
                 DO 744 INL=1,NL
                  DO 743 IJ=1,3
              SPOSTF(IS2,INL,IJ)=SPOSTF(IS2-1,INL,IJ)+ASA(1,IJ)        
743     CONTINUE
744     CONTINUE
            ELSE
                 DO 742 INL=1,NL
                  DO 741 IJ=1,3
               SPOSTF(IS2,INL,IJ)=SPOSTF(IS2-1,INL,IJ)+ASA(2,IJ)        
741              CONTINUE
742              CONTINUE
            ENDIF
         ENDIF                 
 
        
         DO 56 INL=1,NL
                DO 57 IJ=1,3
                  POSTF(INL,IJ)=SPOSTF(IS2,INL,IJ)
57        CONTINUE
56        CONTINUE

         CALL DELXGENTFSL_TEMP(E,VPIS,NLAY,AK2,AK3,AK2M,AK3M,
     &        NL1,NL2,AAK,AAJ,TVA,POSTF,XISTFS,FS)

        XISTFSTOT=XISTFSTOT+XISTFS

        ANORM2=CABS(XISTFS)*CABS(XISTFS)
        ANORM1=ANORM1+ANORM2
c Has convergence with respect to the number of SL-s been achieved?
        IF (ANORM2/ANORM1-0.000001.LE.0) GOTO 340 
c        IF (ANORM2/ANORM1-0.0000000000001.LE.0) GOTO 340 
c If not, sum the contribution from next 'deeper' substrate layer
        GOTO 380

        ENDIF

c Convergence with respect to SL number achieved 
340        CONTINUE
C Set up the new plane wave amplitude XIST as the sum of the 
C plane wave amplitudes coming from overlayer atoms (XISTOL) and the
C plane wave amplitudes coming from substrate atoms (XISTFSTOT) 
C
                XIST(J,I)=XISTFOL+XISTFSTOT
C
C Finally, convert plane wave amplitudes to intensities, including
C appropriate prefactors. (AAB/AKZ accounts of relative beam cross-sections
C while PRE accounts for shift in the surface barrier)
C
                  AT(J,I)=CABS(XIST(J,I))
     &                *CABS(XIST(J,I))
               AT(J,I)=AT(J,I)*(AAB)/AKZ
              AT(J,I)=AT(J,I)*CABS(PRE)*CABS(PRE)
            ELSE
               AT(J,I)=0.0
           ENDIF
110      CONTINUE
100   CONTINUE

c Generate the R-factor for the current structure.

        DO 252 J=1,NT0
         DO 251 I=1,IEERG
        IF (I.LE.NERG) THEN
        AT2(J,I)=AT(J,I)
        ELSE
        AT2(J,I)=0.
        ENDIF
251     CONTINUE
252     CONTINUE

cjcm        write(*,*) 'before rfactf'
      CALL RFACTF(AT2,ETH,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,
     & IPR,XPL,
     & YPL,NNN,AP,APP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,
     & TSEY2,WR,WB,
     & IBP,NT0,TSEY,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,
     & ROS,R1,R2,RP1,
     & RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,FVAL,VOPT,YT,NERG)



c Dump the IV curves for the current structure.

c      WRITE (2,1000)
c      DO 1111 I=1,NLAY
c         WRITE (2,1010) (WPOSTF(I,J)*0.529,J=1,3)
c1111  CONTINUE

        WRITE (2,500)
      WRITE (2,501)
cjcm        write(*,*) 'before wrivtf'
        CALL WRIVTF(AT2,ETH,AE,EE,NET,NEE,IEERG,NT0,NBED,VOPT,
     & IBK,WR,WB,
     & ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,BENAME,
     & NBE,NST1,NST2)
c      VOPT=VOPT+DVOPT+VV*27.21
      WRITE (2,505) FRAC
      WRITE (2,502) NNDIM
      WRITE (2,503) ISMOTH
      WRITE (2,506)
      DO 150 I=1,NLAY
         WRITE (2,510) (ADISP(I,J),J=1,3), 1
150   CONTINUE
      WRITE (2,516) ITER
      WRITE (2,513) IFUNC1
      WRITE (2,514) FVAL
      WRITE (2,515) VOPT


c     WRITE (*,*) 'optimum R-factor', FVAL

      WRITE (2,514) FVAL
        vintentf = fval
c        write(*,*) 'vintentf',vintentf

      RETURN
      END


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
C DIMENSION LIMITATIONS:
C The dimension of APPP limits the number of inequivalent beams to 20
C 2000 is the maximum number of energy points in a single beam after
C interpolation
C
      COMMON /ENY/EI,EF,DE,NERG,NSYM,NDOM,VV,VPIS
6     FORMAT (
     & 59H0** IN PRESENT BEAM TOO FEW ENERGY VALUES FOR INTERPOLATION)
80    FORMAT (40H0INTENSITIES AFTER INTERPOLATION IN BEAM,1I4,/,50(5
     & (1F7.2,1E13.4,3X),/))
C
      ITIL=0
      ITIH=0
      E2=EINCR/2.
C      E2=1.

      DO 60 IB=1,NB
         NEM=NE(IB)
C
C  FIND FIRST NON-ZERO INTENSITY (FOR THEORY, WHERE NON-EMERGENCE OF
C  CURRENT BEAM CAN OCCUR)
C
         DO 30 IE=1,NEM
            IMIN=IE
cONA            IF (A(IB,IE).GT.1.E-6) GOTO 91
             IF (A(IB,IE).GT.1.E-16) GOTO 91
30       CONTINUE
91       IF (IMIN.NE.NEM) THEN
            LMIN=INT((E(IB,IMIN)-E2)/EINCR)+1
            LMIN=MAX0(LMIN,0)
C XMIN (XMAX)is the minimum (maximum) energy on the working grid 
C for which the computed intensities are available
            XMIN=FLOAT(LMIN)*EINCR
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
            KLO=1
            KHI=2
C this shift allows all points in the grid to be included in the energy 
C interval for which we have data
            NPTS=NPTS-2
            NE(IB)=NPTS
            CALL SPLINE(X,WORYT,NEM,YP1,YPN,WORYT2)
            CALL SPLINE(X,WORYT2,NEM,YP1,YPN,WORYT4)
               DO 10 I=1,NPTS
28          IF(XVAL.GE.X(KLO).AND.XVAL.LT.X(KHI))GOTO 27 
                  KLO=KLO+1
                  KHI=KHI+1
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

         SUBROUTINE LGND(LMAX,X,P)
C
C Subroutine to generate Legendre polynomials P_L(X)
C for L = 0,1,...,LMAX with given X.
C
         DIMENSION P(LMAX+1)

         P(1) = 1.
         P(2) = X
         DO 100 L = 1, LMAX-1
            P(L+2) = ((2.0*L+1)*X*P(L+1)-L*P(L))/(L+1)
100         CONTINUE
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
      SUBROUTINE LXGENT(LX,LXI,LT,LXM,LMAX,LMMAX)
C
      DIMENSION LX(LMMAX),LXM(LMMAX),LXI(LMMAX),LT(LMMAX)
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
         ELSE
            LEO=LEO+1
            LX(LEO)=LL
         ENDIF
      ELSE
         LT1=LT1+1
         LT(LL)=LT1
         IF (MOD(L,2).EQ.1) THEN
            LEE=LEE+1
            LX(LEE)=LL
         ELSE
            L1=L1+1
            LX(L1)=LL
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
C
C======================================================================
C
      SUBROUTINE MATGEN(ARA1,ARA2,ARB1,ARB2,LATMAT)
C
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2)
      DIMENSION JAW(4),A(4)
      INTEGER LATMAT(2,2)
      COMPLEX WORK(4,4),CONT(4)
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
      CALL ZGE(WORK,JAW,4,4,EMACH)
      CALL ZSU(WORK,JAW,CONT,4,4,EMACH)
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

C ========================================================================
C
C Subroutine OPSIGN determines the fractional energy range where the
C theoretcial and experimental data sets have slopes of opposite sign.
C
C From the VAN HOVE R-factor program.
C
C ========================================================================
C
      SUBROUTINE OPSIGN(AEP,NBED,IBE,ATP,NBTD,IBT,EINCR,IE1,IE2,NV,ROS)
C
      DIMENSION AEP(NBED,1),ATP(NBTD,1)
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

C======================================================================
C
C  Subroutine PSTEMT incorporates the thermal vibration effects in the 
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
      SUBROUTINE PSTEMP(PPP,N1,N2,N3,DR0,DR,T0,TEMP,E,PHS,DEL)
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

c AGL On the boundaries, we need to consider only non-equivalent atoms
        SUBROUTINE QUAD_2D_BOUND (Q1,Q2,Q3,Q4,PB,NBOUND,
     &                        PBBNE,NBNE)
c*******************************************************************************
c
cQUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
c
c  Modified:
c
c    15 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real ( kind = 8 ) Q(2,4), the vertices of the quadrilateral.
c
c    Input, real ( kind = 8 ) P(2), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if the point is in the quadrilateral.
 
c        implicit none


        REAL DPB1,DPB2
        REAL P, Q1,Q2,Q3,Q4
        REAL RX,RY,RXM11,RXP11,RYP12,RXM12
        REAL RXM21,RXP21,RYP22,RXM22
        DIMENSION P(2), Q1(2),Q2(2),Q3(2),Q4(2)
        DIMENSION PB(200,2),AB1(2),AB2(2),PBC(4,2)
        DIMENSION PBB(200,2)
        DIMENSION PBBNE(200,2)

        
c
c  This will only handle convex quadrilaterals.
c


        Q12MOD=SQRT((Q1(1)-Q2(1))**2+(Q1(2)-Q2(2))**2)
        Q13MOD=SQRT((Q1(1)-Q3(1))**2+(Q1(2)-Q3(2))**2)
        Q14MOD=SQRT((Q1(1)-Q4(1))**2+(Q1(2)-Q4(2))**2)
        Q23MOD=SQRT((Q2(1)-Q3(1))**2+(Q2(2)-Q3(2))**2)
        Q24MOD=SQRT((Q2(1)-Q4(1))**2+(Q2(2)-Q4(2))**2)
        Q34MOD=SQRT((Q3(1)-Q4(1))**2+(Q3(2)-Q4(2))**2)


c AGL The vectors defining the lattice are
        AB1(1)=Q2(1)-Q1(1)
        AB1(2)=Q2(2)-Q1(2)
        AB2(1)=Q4(1)-Q1(1)
        AB2(2)=Q4(2)-Q1(2)



c AGL On the boundaries, we need to consider only non-equivalent atoms

        IBC=1
        IBB=1
        DO 101 IB1=1,NBOUND
c AGL First, get corner points
c ----------------------------
           PBQ11=PB(IB1,1)-Q1(1)
           PBQ12=PB(IB1,2)-Q1(2)
           PBQ21=PB(IB1,1)-Q2(1)
           PBQ22=PB(IB1,2)-Q2(2)
           PBQ31=PB(IB1,1)-Q3(1)
           PBQ32=PB(IB1,2)-Q3(2)
           PBQ41=PB(IB1,1)-Q4(1)
           PBQ42=PB(IB1,2)-Q4(2)

           IF ((ABS(PBQ11).LT.0.01).AND.(ABS(PBQ12).LT.0.01)) THEN 
                              PBC(IBC,1)=PB(IB1,1)
                              PBC(IBC,2)=PB(IB1,2)
                       IBC=IBC+1
                         GOTO 201
           ELSEIF ((ABS(PBQ21).LT.0.01).AND.(ABS(PBQ22).LT.0.01)) THEN 
                              PBC(IBC,1)=PB(IB1,1)
                              PBC(IBC,2)=PB(IB1,2)
                         IBC=IBC+1
                         GOTO 201
           ELSEIF ((ABS(PBQ31).LT.0.01).AND.(ABS(PBQ32).LT.0.01)) THEN
                              PBC(IBC,1)=PB(IB1,1)
                              PBC(IBC,2)=PB(IB1,2)
                         IBC=IBC+1
                         GOTO 201
           ELSEIF ((ABS(PBQ41).LT.0.01).AND.(ABS(PBQ42).LT.0.01)) THEN
                              PBC(IBC,1)=PB(IB1,1)
                              PBC(IBC,2)=PB(IB1,2)
                         IBC=IBC+1
                         GOTO 201
           ENDIF

c AGL If none of the four conditions above hold -> then it's an atom
c AGL on the boundary
        
          PBB(IBB,1)=PB(IB1,1)
          PBB(IBB,2)=PB(IB1,2)
          IBB=IBB+1

201        CONTINUE
101        CONTINUE


c AGL Number of atoms at corner
        NBC=IBC-1
c AGL Number of atoms on boundary (not corners) 
        NBB=IBB-1

c AGL So now we have the atoms on the boundary  cell separated into:
c AGL atoms on corners (PBC) and atoms not on corners (PBB)

c AGL Of the four atoms on the corners, there is only 1 inequivalent
c AGL the other 3 are equivalent
c AGL we will take the atom on Q1

        IBNE=1
        PBBNE(IBNE,1)=Q1(1)
        PBBNE(IBNE,2)=Q1(2)
        IBNE=IBNE+1        


c AGL Now we have to check which atoms are equivalent within 
c AGL those atoms not on corners

        DO 102 IB2=1,NBB
        DO 103 IB3=IB2,NBB
         DPB1=PBB(IB2,1)-PBB(IB3,1) 
         DPB2=PBB(IB2,2)-PBB(IB3,2)
c AGL Avoid numerical errors by defining an area of radius RX=RY around 
c AGL atomic positions
         RX=0.01
         RY=0.01

         RXM11=AB1(1)-RX
         RXP11=AB1(1)+RX
         RYM12=AB1(2)-RY
         RYP12=AB1(2)+RY
         RXM21=AB2(1)-RX
         RXP21=AB2(1)+RX
         RYM22=AB2(2)-RY
         RYP22=AB2(2)+RY



c AGL Of each pair of equivalent atoms we will store one in PBBNE
c----------------------------------------------------------------
c----------------------------------------------------------------

c AGL Atoms equivalent through AB1
c --------------------------------

          IF ((DPB1.GE.RXM11).AND.(DPB1.LE.RXP11)) THEN
             IF ((DPB2.GE.RYM12).AND.(DPB2.LE.RYP12)) THEN
                 PBBNE(IBNE,1)=PBB(IB3,1)
                PBBNE(IBNE,2)=PBB(IB3,2)
                 IBNE=IBNE+1
c                 GOTO 401
             ELSE
             ENDIF
          ELSE
          ENDIF
c AGL Atoms equivalent through AB2
c --------------------------------
c          ELSE
            IF ((DPB1.GE.RXM21).AND.(DPB1.LE.RXP21)) THEN
            IF ((DPB2.GE.RYM22).AND.(DPB2.LE.RYP22)) THEN
                 PBBNE(IBNE,1)=PBB(IB3,1)
                 PBBNE(IBNE,2)=PBB(IB3,2)
                 IBNE=IBNE+1
            ELSE
            ENDIF
           ENDIF




401        CONTINUE


103        CONTINUE
102        CONTINUE              

c AGL Total number of non-equivalent atoms on boundary
        NBNE=IBNE-1
        
        end



        SUBROUTINE QUAD_CONTAINS_2D (Q1,Q2,Q3,Q4,P,NYN)
c*******************************************************************************
c
cQUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
c
c  Modified:
c
c    15 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real ( kind = 8 ) Q(2,4), the vertices of the quadrilateral.
c
c    Input, real ( kind = 8 ) P(2), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if the point is in the quadrilateral.
 
c        implicit none


        REAL P, Q1,Q2,Q3,Q4
        DIMENSION P(2), Q1(2),Q2(2),Q3(2),Q4(2)

c
c  This will only handle convex quadrilaterals.
c
        Q12MOD=SQRT((Q1(1)-Q2(1))**2+(Q1(2)-Q2(2))**2)
        Q13MOD=SQRT((Q1(1)-Q3(1))**2+(Q1(2)-Q3(2))**2)
        Q14MOD=SQRT((Q1(1)-Q4(1))**2+(Q1(2)-Q4(2))**2)
        Q23MOD=SQRT((Q2(1)-Q3(1))**2+(Q2(2)-Q3(2))**2)
        Q24MOD=SQRT((Q2(1)-Q4(1))**2+(Q2(2)-Q4(2))**2)
        Q34MOD=SQRT((Q3(1)-Q4(1))**2+(Q3(2)-Q4(2))**2)




        PQ1MOD=SQRT((P(1)-Q1(1))**2+(P(2)-Q1(2))**2)           
        PQ2MOD=SQRT((P(1)-Q2(1))**2+(P(2)-Q2(2))**2)           
        PQ3MOD=SQRT((P(1)-Q3(1))**2+(P(2)-Q3(2))**2)           
        PQ4MOD=SQRT((P(1)-Q4(1))**2+(P(2)-Q4(2))**2)           


c AGL Reference point: Q2
c-------------------------
        CALL ANGLE_2D (Q1,Q2,Q3,ANGLE_1)
        CALL ANGLE_2D (Q1,Q2,P,ANGLE_2)



c AGL If angle_1 is 'almost' equal to angle_2, it may be on a boundary
c AGL To check this we have to compare the projection of P over the
c AGL boundary
        IF(ABS(ANGLE_1-ANGLE_2).LT.0.001) THEN
        IF (PQ2MOD.LE.Q23MOD) THEN
                GOTO 303
           ENDIF 
        ELSE
c AGL If it's not on the boundary, it can be either inside or outside
c AGL If inside, then angle_1 < angle_2
        IF (ANGLE_1.LT.ANGLE_2) GOTO 301
        ENDIF

c AGL Reference point: Q3
c-------------------------
        CALL ANGLE_2D (Q2,Q3,Q4,ANGLE_1)
        CALL ANGLE_2D (Q2,Q3,P,ANGLE_2)


        IF(ABS(ANGLE_1-ANGLE_2).LT.0.001) THEN
           IF (PQ3MOD.LE.Q34MOD) THEN
             GOTO 303
           ENDIF
        ELSE
c AGL If it is not on the boundary, it can be either inside or outside
c AGL If inside, then angle_1 < angle_2
        IF (ANGLE_1.LT.ANGLE_2) GOTO 301
        ENDIF

c AGL Reference point: Q4
c-------------------------
        CALL ANGLE_2D (Q3,Q4,Q1,ANGLE_1)
        CALL ANGLE_2D (Q3,Q4,P,ANGLE_2)

        IF(ABS(ANGLE_1-ANGLE_2).LT.0.001) THEN
           IF (PQ4MOD.LE.Q14MOD) THEN
             GOTO 303
           ENDIF
        ELSE
c AGL If it's not on the boundary, it can be either inside or outside
c AGL If inside, then angle_1 < angle_2
        IF (ANGLE_1.LT.ANGLE_2) GOTO 301
        ENDIF


c AGL Reference point: Q1
c-------------------------
        CALL ANGLE_2D (Q4,Q1,Q2,ANGLE_1)
        CALL ANGLE_2D (Q4,Q1,P,ANGLE_2)


        IF(ABS(ANGLE_1-ANGLE_2).LT.0.001) THEN
           IF (PQ1MOD.LE.Q12MOD) THEN
             GOTO 303
           ENDIF
        ELSE
c AGL If it's not on the boundary, it can be either inside or outside
c AGL If inside, then angle_1 < angle_2
        IF (ANGLE_1.LT.ANGLE_2) GOTO 301
        ENDIF

        NYN=2        
        GOTO 302

301        CONTINUE

        NYN=1
        GOTO 302

303        CONTINUE
        NYN=3
        GOTO 302

        
302        CONTINUE


        end


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
      SUBROUTINE READCT(NLAY,VPOS,CPVPOS,NTAUAW,LPSAW,LMMAX,
     % IPR,LAFLAG,NST1,ASB,VICL,VCL,FRCL,TST,TSTS,ASA,INVECT)
C
      DIMENSION FPOS(60,3),VPOS(NST1,NLAY,3),LPSAW(NST1,NLAY)
      DIMENSION CPVPOS(NST1,NLAY,3),ASB(NST1,3),ASA(10,3)
      DIMENSION LAFLAG(NST1),NTAUAW(NST1)
      DIMENSION VICL(NST1),VCL(NST1),FRCL(NST1)
C
      COMMON /ZMAT/IANZ(60),IZ(60,4),BL(60),ALPHA(60),BETA(60),NZ,IPAR
     & (15,5),NIPAR(5),NPAR,DX(5),NUM,NATOMS,BLS(60),ALPHAS(60),BETAS
     & (60),PHIR,PHIM1,PHIM2
C
160   FORMAT (3F7.4)
161   FORMAT (3F7.2)
180   FORMAT (/4H ASB,1I3,3X,3F7.4)
200   FORMAT (60I3)
215   FORMAT (/26H COMPOSITE LAYER VECTOR   ,3F7.4)
216   FORMAT (8H NTAU = ,60I3)
217   FORMAT (/46H PHASE SHIFT ASSIGNMENT IN COMPOSITE LAYER No.,1I2,
     & 3X,60I3)
444   FORMAT (//21H COMPOSITE LAYER No. ,1I3)
445   FORMAT (21H =================== )
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
c            READ (5,*) (LPSAW(K,I),I=1,LAFLAG(K))
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
cORIG                  READ (5,160) (VPOS(K,I,J),J=1,3)
                  READ (5,*) (VPOS(K,I,J),J=1,3)
cjcm                  WRITE (*,*) (VPOS(K,I,J),J=1,3)
                  IF (IPR.GT.0) WRITE (1,215) (VPOS(K,I,J),J=1,3)
                  DO 480 J=1,3
                     FPOS(I,J)=VPOS(K,I,J)/0.529
                     VPOS(K,I,J)=VPOS(K,I,J)/0.529
480               CONTINUE
1041           CONTINUE
c               READ (5,160) (ASB(K,I),I=1,3)
               READ (5,*) (ASB(K,I),I=1,3)
               IF (IPR.GT.0) WRITE (1,180) K,(ASB(K,I),I=1,3)
               DO 4666 IK=1,3
                  ASB(K,IK)=ASB(K,IK)/0.529
4666           CONTINUE
c               READ (5,161) FRCL(K),VCL(K),VICL(K)
               READ (5,*) FRCL(K),VCL(K),VICL(K)
               VCL(K)=VCL(K)/27.21
               VICL(K)=VICL(K)/27.21
               NLAY2=LAFLAG(K)
               CALL SORT(FPOS,NLAY2)
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
               CALL GEOMV(1,FPOS,NLAY2)
               DO 1043 J=1,NLAY2
                  DO 348 K=1,3
                     VPOS(II,J,K)=FPOS(J,K)
348               CONTINUE
1043           CONTINUE
C CPVPOS corresponds to sorted coordinates needed in subroutine LOOKUP
               CALL SORT(FPOS,NLAY2)
               DO 333 J=1,NLAY2
                  DO 332 K=1,3
                     CPVPOS(II,J,K)=0.529*FPOS(J,K)
332               CONTINUE
333            CONTINUE
1042        CONTINUE
         ENDIF
      ENDIF
      TSTS=TST
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
      SUBROUTINE READE(AE,EE,NBED,NEE,NBEA,BENAME,IPR,IEERG,NBMAX)
C
        DIMENSION AE(NBED,IEERG),EE(NBED,IEERG),NEE(NBED),NBEA(NBED)
        DIMENSION BENAME(5,NBED)
        character(len=20) FMT
        COMMON /REXP/EEINCR
C
30    FORMAT (5(25I3,/))
40    FORMAT (I4,25H EXP. BEAMS TO BE READ IN)
50    FORMAT (30HAVERAGING SCHEME OF EXP. BEAMS,5(25I3,/))
60    FORMAT (20A4)
10    FORMAT (20A4)
70    FORMAT (10HEXP. BEAM ,1I3,2H (,5A4,1H))
35    FORMAT (1I3,1E13.4)
36    FORMAT (5F7.4)
80    FORMAT (31H  EXP. ENERGIES AND INTENSITIES,/,50(5(1F7.2,1E13.4,
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
c         READ (11,FMT) (EE(IB,IE),AE(IB,IE),IE=1,N)
         READ (11,99) (EE(IB,IE),AE(IB,IE),IE=1,N)
 99      FORMAT(2F12.3)

c         READ (11,*) (EE(IB,IE),AE(IB,IE),IE=1,N)
         IF (IPR.GE.0) THEN
            WRITE (1,80) (EE(IB,IE),AE(IB,IE),IE=1,N)
         ENDIF
         DO 86 IE=1,N
            AE(IB,IE)=AE(IB,IE)*FAC
86       CONTINUE
90    CONTINUE
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
      SUBROUTINE READPL(NT0,NSET,PQFEX,NINSET,NDIM,DISP,ICOORD,
     & NSTEP,ANSTEP,NLAY,IPR,ALPHA,BETA,GAMMA,ITMAX,
     & FTOL1,FTOL2,MFLAG,LLFLAG,NGRID)
C
      DIMENSION PQFEX(2,NT0),NINSET(20)
      DIMENSION DISP(NLAY,3),LLFLAG(60)
      COMMON /RPL/DVOPT
      COMMON /WIV2/PERSH,NIV
C
100   FORMAT (20I3)
110   FORMAT (10F7.4)
115   FORMAT (3F7.6,I3)
116   FORMAT (F7.6,I3)
120   FORMAT (/' Calculating for Exit Beams; ')
130   FORMAT (' =========================== ')
200   FORMAT (/' Initial Coordinates in Search ')
210   FORMAT (' ============================= ')
220   FORMAT (/' Line Search Information ')
230   FORMAT (' ======================= ')
240   FORMAT (/' Moving coordinate ',I3,' with ',I3,' moving atoms ')
250   FORMAT (/' Moving Atom Identifiers; ')
260   FORMAT (/' Number of steps in Line Search ',I3,' Step Size',F7.4)
270   FORMAT (/' Search Information; ')
280   FORMAT (/' ALPHA =',F7.4,' BETA =',F7.4,' GAMMA =',F7.4)
290   FORMAT (/' Max. Number of Iterations =',I3)
300   FORMAT (/' Down-hill Simplex Search Requested')
310   FORMAT (/' Powell Direction-set Requested')
320   FORMAT (/' Direction Set with Principal Directions Requested')
350   FORMAT (/' Line Search Requested')
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

C=========================================================================
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
     & TST,TSTS,THETA,FI,LMMAX,NPSI,ES,PHSS,PHSS2,L1,IPR,NEL)
C
      COMPLEX VL(NL,2)
      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),ASA(10,3),ARB1(2)
      DIMENSION ARB2(2),RBR1(2),RBR2(2)
      DIMENSION CPARB1(2),CPARB2(2),PHSS2(NPSI,80)
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
140   FORMAT (/25H PARAMETERS FOR INTERIOR )
160   FORMAT (3F7.4)
166   FORMAT (3F11.8)
161   FORMAT (3F7.2)
170   FORMAT (/10H SURF VECS,2(5X,2F8.4))
171   FORMAT (15X,2F8.4)
172   FORMAT (/' Interlayer Vectors ')
180   FORMAT (/4H ASA,1I3,3X,3F7.4)
200   FORMAT (24I3)
210   FORMAT (/6H FR = ,F7.4,7H ASE = ,F7.4)
280   FORMAT (1H ,1I4,F10.3,F7.3,4X,4I3)
281   FORMAT (7H TST = ,1F7.4)
285   FORMAT (' THETA  FI = ',2F7.2)
290   FORMAT (2HVV,4X,F7.2,4HVPIS,4X,F7.2)
305   FORMAT (7H IT1 = ,3I3)
325   FORMAT (8H THDB = ,1F9.4,6H AM = ,1F9.4,8H FPER = ,1F9.4,
     & 8H FPAR = ,1F9.4,7H DR0 = ,1F9.4)
340   FORMAT (16F7.4)
350   FORMAT (/6X,12HPHASE SHIFTS)
355   FORMAT (8H LMAX = ,1I3)
360   FORMAT (4HE = ,1F7.4,13H  1ST ELEMENT,3X,16F8.4)
361   FORMAT (1H ,11X,11H ELEMENT # ,I2,3X,16F8.4)
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
C  NEL= NO. OF CHEMICAL ELEMENTS IN CALCULATION (.LE.3)
C
      READ (5,200) NEL
      IF (NEL.GT.5) THEN
         WRITE (1,*) ' CURRENT CODE NOT DIMENSIONED FOR NEL.GE.5'
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
           I2=1+I0
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
              WRITE (1,361)(PHSS(I,L),L=1+IO,L1+IO) 
           ENDIF
671   CONTINUE
670   CONTINUE
       ENDIF
C
C  ARA1 AND ARA2 ARE TWO 2-D BASIS VECTORS OF THE SUBSTRATE LAYER
C  LATTICE. THEY SHOULD BE EXPRESSED IN TERMS OF THE PLANAR CARTESIAN
C  Y- AND Z-AXES (X-AXIS IS PERPENDICULAR TO SURFACE)(ANGSTROM)
C
      READ (5,160) (ARA1(I),I=1,2)
      READ (5,160) (ARA2(I),I=1,2)
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
      READ (5,160) (ARB1(I),I=1,2)
      READ (5,160) (ARB2(I),I=1,2)
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
15    FORMAT (16H RENORM. IN BEAM,1I4,7H - AA =,1E15.4,8H, ALPH =,
     & 1E15.4)
60    FORMAT (49H0ENERG. AND INTENS. AFTER RENORMALIZATION IN BEAM,1I3,
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
               CALL SUM1(Y,1,1,1.,1,N1,S1,IEERG)
               CALL SUM1(Y,1,1,1.,N2,NT,S2,IEERG)
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
      SUBROUTINE RFACTF(AT,ETH,INBED,IEERG,AE,EE,NEE,NBEA,BENAME,IPR,
     & XPL,YPL,NNN,AEP,AEPP,YE,TSE,TSE2,TSEP,TSEP2,TSEPP,TSEPP2,TSEY2,
     & WR,WB,IBP,NT0,TSEY,ATP,ATPP,TST,TSTY2,NST1,NST2,RAV,IBK,ROS,R1,
     & R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,EET,NET,AR,BARAV,BV0,YT,NERG)
C
        INTEGER IBP,IBT
c        INTEGER IBK
      DIMENSION IBP(NT0)
      DIMENSION AT(NT0,IEERG),ETH(NT0,IEERG),AE(INBED,IEERG)
      DIMENSION EE(INBED,IEERG),NEE(INBED),NBEA(INBED),YPL(IEERG)
      DIMENSION BENAME(5,INBED),XPL(IEERG),NNN(IEERG),AEP(INBED,IEERG)
      DIMENSION AEPP(INBED,IEERG),YE(INBED,IEERG),TSE(INBED),TSE2(INBED)
      DIMENSION TSEP(INBED),TSEP2(INBED),TSEPP(INBED),TSEPP2(INBED)
      DIMENSION WR(10),WB(NT0),TSEY(INBED),TSEY2(INBED)
      DIMENSION ATP(NT0,IEERG),ATPP(NT0,IEERG),TST(NT0),TSTY2(NT0)
      DIMENSION NST1(NT0),NST2(NT0),RAV(NT0),IBK(NT0),EET(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),NET(NT0)
      DIMENSION AR(11),YT(NT0,IEERG)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
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


        do 33 i=1,NT0
c          IBK(i)=0
33        continue

        CALL RINTAV(AT,NT0,IBP,IEERG,NERG,NBMAX,IBK)

      DO 100 I=1,NBMAX
         NET(I)=NERG
100   CONTINUE
C
C SHIFT THE ENERGY ACCORDING TO THE SHIFT OF THE INNER POTENTIAL
C DETERMINED BY THE SEARCH ALGORITHM
C
      DO 101 I=1,NBMAX
        DO 102 J=1,NERG
         ETH(I,J)=ETH(I,J)-BV0-DVOPT
102     CONTINUE
101   CONTINUE
C
C Now interpolate on the same grid used when smoothing the experimental 
C data (EEINCR)
C and smooth as done with the experimental data
C
         IF (ISMOTH.NE.0) THEN
            CALL INTPOL(AT,ATP,ATPP,NT0,NET,NBMAX,ETH,EEINCR,
     &      IPR,XPL,YPL,IEERG)
            DO 10 I=1,ISMOTH
               CALL SMOOTH(AT,ETH,NT0,NBMAX,NET,IPR,IEERG)
10          CONTINUE
         ENDIF
C
C INTERPOLATE THEORY ONTO WORKING GRID
C
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
c         IF(IFUNC.EQ.0) THEN
C           CALL AVEINT(AT,ATP,NT0,NET,NBMAX,ETH,IEERG,EEAVT)
           CALL AVEINT(AT,ATP,NT0,NST1,NST2,NBMAX,ETH,
     &       IEERG,EEAVT,1)

c         ENDIF
C
C  PRODUCE 1ST AND 2ND DERIVATIVES OF TH. DATA (now done by INTPOL)
C
C      CALL DERL(AT,NET,NT0,NBMAX,ATP,EINCR,IEERG)
C      CALL DERL(ATP,NET,NT0,NBMAX,ATPP,EINCR,IEERG)
C
C  PRODUCE PENDRY Y FUNCTION FOR TH. DATA
C
      CALL YPEND(AT,ATP,ATPP,NT0,NBMAX,NET,ETH,YT,
     %  VPIS,IPR,IEERG,EEAVT)
C      PSH=.0
C      DO 1112 J=1,NBMAX
C      DO 1111 I=1,NET(J)
C         WRITE (*,*) ETH(J,I),zlv
C         WRITE (*,*) ETH(J,I),YT(J,I)
C          IF(J.EQ.5.OR.J.EQ.1) THEN
C             VIEFF=VPIS*27.21
C             TEST= ABS(AT(J,I)/EEAVT(J)*VIEFF)
C             IF(TEST.LT.EINCR/2.) THEN
C               WRITE (*,*) ETH(J,I),0.
C             ELSE
C               WRITE (*,*) ETH(J,I),AT(J,I)+EEAVT(J)*PSH
C             ENDIF
C          ENDIF
C         WRITE (*,*) ETH(J,I),ATPP(J,I)
C1111  CONTINUE
C      WRITE (*,*)
C1112  CONTINUE
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
     &                NV,ROS(IBE))
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
C WT      =  WEIGHTING INFORMATION FOR DIFFERENT TERMINATIONS(SUM WT=1)
C NT     =  NUMBER OF DIFFERENT TERMINATIONS
C
C AUTHOR:WANDER
C
C ==========================================================================
C
      SUBROUTINE RFIN(IBP,NT0,WB,WR,NT,IPR)
C
      DIMENSION IBP(NT0),WR(10),WB(NT0)
c      DIMENSION WT(NT)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
      COMMON /FUN/WT(5)
C
100   FORMAT (40I3)
110   FORMAT (20F7.4)
140   FORMAT (/' Beam Averaging Information; ')
150   FORMAT (' =========================== ')
160   FORMAT (/' R-Factor Information; ')
170   FORMAT (' ===================== ')
180   FORMAT (/' Interpolation Step Size;',1F7.4)
190   FORMAT (' Inner Potential Range;',2F7.4,' Step Size ',1F7.4)
C
C Read in the information needed to link together the experimental and
C theoretical beams. If IBP(i)=IBP(j) then theoretical beams i and j will
C be avveraged and compared with experimental beam IBP(i).
C
      MAXB=0
      READ (12,100) (IBP(J),J=1,NT0)
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
      READ (12,110) EINCR
      IF (IPR.GT.0) THEN
         WRITE (1,160)
         WRITE (1,170)
         WRITE (1,180) EINCR
      ENDIF
C
C Read in weighting information for R-factors.
C WB=Beam weighting in each R-factor
C WR=R-factor weighting in 10 R-factor average
C WT=Termination weighting factor
C
      READ (12,110)(WB(I),I=1,MAXB)
      READ (12,110)(WR(I),I=1,10)
c      READ (12,110)(WT(I),I=1,NT)
C
C Read in remaining parameters;
C ISMOTH=Desired number of experimental smoothing 
C IREN=1(0) Renormalize (or not) the experimental and theoretical curves.
C IRGEXP=0 if experimental data is on a uniform grid,
C        1 if data is on a non-uniform grid. In this case, an initial 
C          interpolation will be done. This is not recommended, as it may 
C          enhance experimental noise.
C
      READ (12,100)ISMOTH,IREN,IRGEXP 
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
C WT      =  WEIGHTING INFORMATION FOR DIFFERENT TERMINATIONS(SUM WT=1)
C NT     =  NUMBER OF DIFFERENT TERMINATIONS
C
C AUTHOR:WANDER
C
C ==========================================================================
C
      SUBROUTINE RFINTF(IBP,NT0,WB,WR,IPR)
C
      DIMENSION IBP(NT0),WR(10),WB(NT0)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
      COMMON /RFACY/MAXB,IREN,ISMOTH,IRGEXP,NBE,NBED
C
100   FORMAT (40I3)
110   FORMAT (20F7.4)
140   FORMAT (/' Beam Averaging Information; ')
150   FORMAT (' =========================== ')
160   FORMAT (/' R-Factor Information; ')
170   FORMAT (' ===================== ')
180   FORMAT (/' Interpolation Step Size;',1F7.4)
190   FORMAT (' Inner Potential Range;',2F7.4,' Step Size ',1F7.4)
C
C Read in the information needed to link together the experimental and
C theoretical beams. If IBP(i)=IBP(j) then theoretical beams i and j will
C be avveraged and compared with experimental beam IBP(i).
C
      MAXB=0
      READ (12,100) (IBP(J),J=1,NT0)
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
      READ (12,110) EINCR
      IF (IPR.GT.0) THEN
         WRITE (1,160)
         WRITE (1,170)
         WRITE (1,180) EINCR
      ENDIF
C
C Read in weighting information for R-factors.
C WB=Beam weighting in each R-factor
C WR=R-factor weighting in 10 R-factor average
C WT=Termination weighting factor
C
      READ (12,110)(WB(I),I=1,MAXB)
      READ (12,110)(WR(I),I=1,10)
C
C Read in remaining parameters;
C ISMOTH=Desired number of experimental smoothing 
C IREN=1(0) Renormalize (or not) the experimental and theoretical curves.
C IRGEXP=0 if experimental data is on a uniform grid,
C        1 if data is on a non-uniform grid. In this case, an initial 
C          interpolation will be done. This is not recommended, as it may 
C          enhance experimental noise.
C
      READ (12,100)ISMOTH,IREN,IRGEXP 
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
c        INTEGER IBP
      DIMENSION ATH(NT0,IEERG),IBP(NT0),TEMP(30,221),ITEMP(30)
      DIMENSION IBK(NT0)
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

C============================================================================
c SLPOS finds the substrate atoms that lie within the superlattice unit
c cell defined by ARB1,ARB2
c
c PARAMETER LIST
c
c Input
c------
c ASB        =  INTERLAYER VECTOR BETWEEN LAST OVERLAYER AND TOP
C                         SUBSTRATE LAYER.
c VPOS        =  Atomic coordinates of the atoms in the composite layers
c LAFLAG = Number of layers in each composite layer
c ASE        = SPACING BETWEEN SURFACE AND OVERLAYER NUCLEI
c ARA1,ARA2 = TWO 2-D BASIS VECTORS OF THE SUBSTRATE LAYER LATTICE
c ARB2,ARB2 = ARE EQUIVALENT TO ARA1,ARA2 BUT FOR AN OVERLAYER
c
c Output
c-------
c SPOSTF = Substrate atoms that lie within the superlattice unit cell
c
c In common blocks
c-----------------
c ASB        =  INTERLAYER VECTOR BETWEEN LAST OVERLAYER AND TOP
C                         SUBSTRATE LAYER.
c ASE        = SPACING BETWEEN SURFACE AND OVERLAYER NUCLEI
c ARA1,ARA2 = TWO 2-D BASIS VECTORS OF THE SUBSTRATE LAYER LATTICE
c ARB2,ARB2 = ARE EQUIVALENT TO ARA1,ARA2 BUT FOR AN OVERLAYER
c NL1,NL2 = SUPERLATTICE CHARACTERIZATION CODES
c RBR1,RBR2,ASL,FR,VPIS,VPIO,VO,VV = not used

        SUBROUTINE SLPOS(ASB,VPOS,LAFLAG,ASE,NST1,NLAYTOT,SPOSTF)
    
        INTEGER NSLBNE,NSLIN,NSLB
 
        DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),ASA(10,3),ARB1(2)
        DIMENSION ARB2(2),RBR1(2),RBR2(2)
        DIMENSION ASB(NST1,3)
        DIMENSION VPOS(NST1,NLAYTOT,3),LAFLAG(NST1)
        DIMENSION SPOSTF(20,1500,3) 
        DIMENSION SPOS2(20,1500,3),ISIN(1),ISB(1)
        DIMENSION PBBNE(200,2),PBX(20)        
        DIMENSION NSLIN(1000),NSLB(1000),NSLBNE(1000)
        DIMENSION SPOS3(20,1500,3),SPOS3B(20,1500,3)
        DIMENSION SPOS3BNE(20,1500,3),PB(200,2)
        DIMENSION Q1(2),Q2(2),Q3(2),Q4(2),P(2)
        DIMENSION Q12(2),Q14(2)

        COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2        

c The 1st substrate layer is referred to the last composite layer, i.e.,
c to NCL=NST1


c 1st atom on 1st SL
            SPOS2(1,1,2)=ASB(NST1,2)
            SPOS2(1,1,3)=ASB(NST1,3)
            SPOS2(1,1,1)=ASB(NST1,1)+VPOS(NST1,LAFLAG(NST1),1)+
     &                                ABS(VPOS(1,1,1))+ASE


c Repeat the substrate configuration of the 1st SL over the yz plane
c using the substrate lattice vectors
c-------------------------------------------------------------------


c MNS1, MNS2 positive integers!!
            MNS1=2
            MNS2=2
230            CONTINUE
            IS2=2
            DO 12 NS1=1,MNS1
            DO 13 NS2=1,MNS2
              MS1=FLOAT(MNS1/2)-NS1
              MS2=FLOAT(MNS2/2)-NS2
              D2=(MS1*ARA1(1))+(MS2*ARA2(1))
              D3=(MS1*ARA1(2))+(MS2*ARA2(2))
              IF ((MS1.EQ.0).AND.(MS2.EQ.0)) GOTO 302 
              SPOS2(1,IS2,2)=SPOS2(1,1,2)+D2
              SPOS2(1,IS2,3)=SPOS2(1,1,3)+D3
              SPOS2(1,IS2,1)=SPOS2(1,1,1)
              IS2=IS2+1
              if (IS2.GT.1000) write(*,*) 'too many substrate atoms'
302        CONTINUE
13        CONTINUE
12        CONTINUE

c AGL Number of atoms in 1st SL 
        NSLXTOT=IS2-1



c We need the substrate atoms of the 1st SL that lie right
c below the superlattice unit cell
c----------------------------------------------------------------------

c Q1,Q2,Q3,Q4 are the corners of the quadrilateral formed by ARB1,ARB2
c the superlattice lattice vectors

        Q1(1)=0.0
        Q1(2)=0.0
        Q2(1)=ARB2(1)
        Q2(2)=ARB2(2)
        Q3(1)=ARB2(1)+ARB1(1)
        Q3(2)=ARB2(2)+ARB1(2)
        Q4(1)=ARB1(1)
        Q4(2)=ARB1(2)

c Moduli of the lattice vectors of superlattice unit cell
        Q1MOD=SQRT((Q1(1)**2)+(Q1(2)**2))
        Q2MOD=SQRT((Q2(1)**2)+(Q2(2)**2))
        Q3MOD=SQRT((Q3(1)**2)+(Q3(2)**2))
        Q4MOD=SQRT((Q4(1)**2)+(Q4(2)**2))


cAGL We must make sure that Q1,Q2,Q4 are taken correctly
cAGL The cross-product between Q1Q2 and Q2Q4 is:

        Q12(1)=Q2(1)-Q1(1)
        Q12(2)=Q2(2)-Q1(2)
        Q14(1)=Q4(1)-Q1(1)
        Q14(2)=Q4(2)-Q1(2)
        
        PROD12=(Q12(1)*Q14(2))-(Q14(1)*Q12(2))
        MOD12=SQRT(PROD12*PROD12)
        SIG12=PROD12/MOD12
        IF (SIG12.LT.0) THEN
        write(*,*) 'ARB1, ARB2 not correctly ordered.'
        write(*,*) 'Exchange them in tleed5.i. STOP'
        STOP
        ELSE
        ENDIF

c AGL Let us separate the substrate lattice (SL) atoms of the 1st SL
c AGL that are inside and on the boundary of the superlattice (SS) unit cell

           ISIN(1)=1
           IS3=1
           ISB(1)=1
           IS4=1
        DO 32 IS=1,NSLXTOT
           P(1)=SPOS2(1,IS,2)
           P(2)=SPOS2(1,IS,3)

 
c Include atoms on boundary
        CALL QUAD_CONTAINS_2D(Q1,Q2,Q3,Q4,P,NYN)




        IF (NYN.EQ.2) THEN
           SPOS3(1,IS3,1)=SPOS2(1,IS,1)
           SPOS3(1,IS3,2)=SPOS2(1,IS,2)
           SPOS3(1,IS3,3)=SPOS2(1,IS,3)
           ISIN(1)=ISIN(1)+1
           IS3=IS3+1
        ELSE
         IF (NYN.EQ.3) THEN
           SPOS3B(1,IS4,1)=SPOS2(1,IS,1)
           SPOS3B(1,IS4,2)=SPOS2(1,IS,2)
           SPOS3B(1,IS4,3)=SPOS2(1,IS,3)
           ISB(1)=ISB(1)+1
           IS4=IS4+1
         ENDIF
        ENDIF        
32      CONTINUE
c AGL Total number of SL atoms inside the SS unit cell
        NSLIN(1)=ISIN(1)-1
c AGL Total number of SL atoms on the boundary of the SS unit cell
        NSLB(1)=ISB(1)-1



c AGL For the 1st SL layer  we need only the
c AGL atoms at the boundaries that are non equivalent

             NBOUND=NSLB(1)
        DO 202 J=1,NSLB(1)
             PBX(1)=SPOS3B(1,J,1)
             PB(J,1)=SPOS3B(1,J,2)
             PB(J,2)=SPOS3B(1,J,3)
202        CONTINUE

c Initialize NSLBNE
        do 980 i=1,1000        
           NSLBNE(i)=0
980        continue


        IF (NBOUND.NE.0) THEN
             CALL QUAD_2D_BOUND (Q1,Q2,Q3,Q4,PB,NBOUND,
     &                  PBBNE,NBNE)
             NSLBNE(1)=NBNE
        ELSE
        ENDIF


c AGL For each layer, the sum of the inequivalent atoms on the boundary
c AGL and the atoms inside the unit cell has to be equal to NL, the
c AGL area of the superlattice unit cell in terms of the area
c AGL of the substrate


        NTOT=NSLIN(1)+NSLBNE(1)


        NL=NL1*NL2

        IF (NTOT.NE.NL) THEN
cAGL Maybe the number of atoms we have considered by translating
cAGL the 1st atom by (ARA1,ARA2) are not enough, so increase
cAGL MNS1,MNS1, i.e., take more SL atoms 'around' the 1st atom
        MNS1=MNS1+5
        MNS2=MNS1+5
cAGL If we have already taken 'many' SL atoms around the 1st atom,
cAGL but we still don't get the correct number of atoms in
cAGL the SS unit cell, maybe there is something else that is wrong..
        IF (MNS1.GT.20) THEN
         write(*,*) 'number of atoms in superlattic
     &                          unit cell
     &                  is not correct.STOP'
        STOP
        ELSE
        ENDIF
        GOTO 230
        ELSE
        ENDIF

        ISL=1

        DO 203 ISLNE=1,NSLBNE(ISL)
             SPOS3BNE(ISL,ISLNE,1)=PBX(ISL)
             SPOS3BNE(ISL,ISLNE,2)=PBBNE(ISLNE,1)
             SPOS3BNE(ISL,ISLNE,3)=PBBNE(ISLNE,2)
203        CONTINUE        

c AGL So, the substrate atoms that need to be included in the 
c AGL toy-function calculation are the following (SPOSTF)

                ISL=1
                ISTF=1
        DO 205 JSL=1,NSLIN(ISL)
                SPOSTF(ISL,ISTF,1)=SPOS3(ISL,JSL,1)
                SPOSTF(ISL,ISTF,2)=SPOS3(ISL,JSL,2)
                SPOSTF(ISL,ISTF,3)=SPOS3(ISL,JSL,3)
                ISTF=ISTF+1
205        CONTINUE
        DO 206 JSL=1,NSLBNE(ISL)
                SPOSTF(ISL,ISTF,1)=SPOS3BNE(ISL,JSL,1)
                SPOSTF(ISL,ISTF,2)=SPOS3BNE(ISL,JSL,2)
                SPOSTF(ISL,ISTF,3)=SPOS3BNE(ISL,JSL,3)
                ISTF=ISTF+1
206        CONTINUE
                

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
60    FORMAT (48H0EXP. ENERG. AND INTENS. AFTER SMOOTHING IN BEAM,1I3,
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
      SUBROUTINE SORT(POS,NLAY)
C
      DIMENSION POS(60,3),POSS(60,3),POSA(3)
C
1000  FORMAT (/,10X,25HCOORDINATES AFTER SORTING,/,13X,1HX,14X,1HY,14X,
     & 1HZ)
1010  FORMAT (7X,F10.7,2(5X,F10.7))
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

C ============================================================================
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
         if (h.eq.0.) then
            write(*,*) 'bad xa input. in SPLINT'
            stop
         endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      yp=(ya(khi)-ya(klo))/h -(3.*a**2 -1.)*h*y2a(klo)/6. +
     *  (3.*b**2 -1.)*h*y2a(khi)/6. 
      return
      end

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
10    TEMP(2, K) = WORY(N)
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
      SUBROUTINE SUM1(Y,NBD,IB,H,I1,I2,S,IEERG)
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
C IT             = NOT USED.
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
700   FORMAT (42H TOO LOW ENERGY FOR AVAILABLE PHASE SHIFTS)
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
C  COMPUTE TEMPERATURE-INDEPENDENT T-MATRIX ELEMENTS
C
         DO 790 L=1,L1
            A=PHS(L)
            AF(L)=A*CI
            AF(L)=CEXP(AF(L))
            A=SIN(A)
            AF(L)=A*AF(L)
            TSF0(IEL,L)=AF(L)
790      CONTINUE
C
C  AVERAGE ANY ANISOTROPY OF RMS VIBRATION AMPLITUDES
C
         DR=SQRT((DRPER(IEL)*DRPER(IEL)+2.0*DRPAR(IEL)*DRPAR(IEL))/3.0)
C
C  COMPUTE TEMPERATURE-DEPENDENT PHASE SHIFTS (DEL)
C
         DR01=DR0(IEL)
         CALL PSTEMP(PPP,NN1,NN2,NN3,DR01,DR,T0,T,E,PHS,DEL)
C
C  PRODUCE TEMPERATURE-DEPENDENT T-MATRIX ELEMENTS
C
         DO 840 L=1,L1
            CA=DEL(L)
            CAF(L)=CA*CI
            CAF(L)=CEXP(CAF(L))
            CAF(L)=-CI*(CAF(L)*CAF(L)-1.0)/2.0
            TSF(IEL,L)=CAF(L)
840      CONTINUE
         RETURN
      ENDIF
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
C IT             = NOT USED.
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
      SUBROUTINE TSCATF_TOY(IEL,L1,ES,PHSS,PHSS2,NPSI,IT,EB,V,PPP,
     & NN1,NN2,NN3,DR0,DRPER,DRPAR,T0,T,TSF0,TSF,AF,CAF,NFLAGINT,
     & PHS,DEL,NERG,IEEV,NEL,PHSSEL)
C
      COMPLEX CI,DEL(10),CA,AF(L1),CAF(L1),TSF0(6,10),TSF(6,10)
      COMPLEX PHSSEL(NERG,NEL,10)
      DIMENSION PHSS(NPSI,50),PHSS2(NPSI,50),ES(NPSI)
      DIMENSION PPP(NN1,NN2,NN3)
      DIMENSION PHS(10),PHSL(90),PHSL2(90)
      DIMENSION DR0(5),DRPER(5),DRPAR(5),IT(5)
C
C
700   FORMAT (42H TOO LOW ENERGY FOR AVAILABLE PHASE SHIFTS)
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
C  COMPUTE TEMPERATURE-INDEPENDENT T-MATRIX ELEMENTS
C
         DO 790 L=1,L1
            A=PHS(L)
            AF(L)=A*CI
            AF(L)=CEXP(AF(L))
            A=SIN(A)
            AF(L)=A*AF(L)
            TSF0(IEL,L)=AF(L)
790      CONTINUE
C
C  AVERAGE ANY ANISOTROPY OF RMS VIBRATION AMPLITUDES
C
         DR=SQRT((DRPER(IEL)*DRPER(IEL)+2.0*DRPAR(IEL)*DRPAR(IEL))/3.0)
C
C  COMPUTE TEMPERATURE-DEPENDENT PHASE SHIFTS (DEL)
C
         DR01=DR0(IEL)
        CALL PSTEMP(PPP,NN1,NN2,NN3,DR01,DR,T0,T,E,PHS,DEL)
C
C  PRODUCE TEMPERATURE-DEPENDENT T-MATRIX ELEMENTS
C
         DO 840 L=1,L1
            CA=DEL(L)
            CAF(L)=CA*CI
            CAF(L)=CEXP(CAF(L))
            CAF(L)=-CI*(CAF(L)*CAF(L)-1.0)/2.0
            TSF(IEL,L)=CAF(L)
            PHSSEL(IEEV,IEL,L)=DEL(L)             
840      CONTINUE
         RETURN
      ENDIF
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
         CALL SUM1(YY,1,1,DE,1,NN,S,IEERG)
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
         CALL SUM1(Y,1,1,EINCR,1,N,S,IEERG)
      ENDIF
      RETURN
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
cNEW            AK2M(IG)=AK2+PSQ(1,IG)
cNEW            AK3M(IG)=AK3+PSQ(2,IG)
100      CONTINUE
      ELSE
C
C NEXIT>0 SO LOCATE THIS BEAM IN MAIN BEAM LIST.
C
C FIRST CHECK FOR EMERGENCE OF THIS BEAM
C
        A=2.0*(E-VV)-AK2M(NEXIT)*AK2M(NEXIT)-AK3M(NEXIT)*AK3M(NEXIT)
cNEW         A=2.0*(E)-AK2M(NEXIT)*AK2M(NEXIT)-AK3M(NEXIT)*AK3M(NEXIT)
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
      SUBROUTINE WRIVTF(AT,ETH,AE,EE,NET,NEE,IEERG,NT0,NBED,VOPT,IBK,
     & WR,WB,ROS,R1,R2,RP1,RP2,RPP1,RPP2,RRZJ,RMZJ,RPE,BENAME,
     & NBE,NST1,NST2)
C
      DIMENSION ETH(NT0,IEERG),AT(NT0,IEERG)
      DIMENSION AE(NBED,IEERG),EE(NBED,IEERG)
      DIMENSION NEE(NBED),NET(NT0),IBK(NT0),NST1(NT0),NST2(NT0)
      DIMENSION ROS(NT0),R1(NT0),R2(NT0),RP1(NT0),RP2(NT0),RPP1(NT0)
      DIMENSION RPP2(NT0),RRZJ(NT0),RMZJ(NT0),RPE(NT0),RFAC(30)
      DIMENSION BENAME(5,NBED),WR(10),WB(NT0)
      CHARACTER NC(9)
      CHARACTER(LEN=5) IV
      CHARACTER(LEN=2) NC2(21)
      CHARACTER(LEN=4) IV2
      CHARACTER(LEN=100) IVNAME
      CHARACTER(LEN=100) IVNAME2
C
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
C110   FORMAT (11H#TitleText: ,5A4,6H Beam ,5A4,5H Rfac,I2,
C     & 1H=,F5.4)
110   FORMAT (11H#TitleText: ,5A4,6H Beam ,5A4,5H Rfac,I2,
     & 1H=,F13.4)
111   FORMAT (5H RFAC,I2)
500   FORMAT (/' THE NUMBER OF THEORETICAL AND EXPERIMENTAL BEAMS 
     & DISAGREE ')
100   FORMAT(20A4)
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
C and record that Rfactor


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
C
C write down the experimental beams(remembering the shift in intensity
C for the Pendry Rfactor)
C
      IV='ivexp'

      NBEAMS=NBE
      NOUT=48
      NBB=NOUT+NBEAMS-1
      DO 1 I=1,NBE
          IF (NOUT.LE.NBB) THEN
c              IF (I.LT.10) IVNAME='/scratch/aran/run1/'//IV//NC(I)
c              IF (I.GE.10) IVNAME='/scratch/aran/run1/'//IV//NC2(I-9)
              IF (I.LT.10) IVNAME='./kleedIV/'//IV//NC(I)
              IF (I.GE.10) IVNAME='./kleedIV/'//IV//NC2(I-9)
              OPEN(UNIT=NOUT,FILE=IVNAME,STATUS='UNKNOWN')
              WRITE (NOUT,110) (TITLE(II),II=1,5),
     &        (BENAME(II,I),II=1,5),NRFAC,RFAC(I) 
              WRITE (NOUT,*) '#"IV exp'
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


      IV2='ivth'
      NOUT2=48
      NBB=NOUT2+NBEAMS-1
      DO 3 I=1,NBE
           IF (NOUT2.LE.NBB) THEN
c              IF (I.LT.10) IVNAME2='/scratch/aran/run1/'//IV2//NC(I)
c              IF (I.GE.10) IVNAME2='/scratch/aran/run1/'//IV2//NC2(I-9)
              IF (I.LT.10) IVNAME2='./kleedIV/' //IV2//NC(I)
              IF (I.GE.10) IVNAME2='./kleedIV/'//IV2//NC2(I-9)
              OPEN(UNIT=NOUT2,FILE=IVNAME2,STATUS='UNKNOWN')
              WRITE (NOUT2,110) (TITLE(II),II=1,5),
     &        (BENAME(II,I),II=1,5),NRFAC,RFAC(I)
              WRITE (NOUT2,*) '#"IV theory'
          ENDIF

        IBT=IBK(I)
C
C SCAL=(I)/(IBT) if we scale each beam separately (NIV=0)
C SCAL=(NBE+1)/(NBE+1)(NIV=1) for only one scale factor
C
c AGL ez dakit zergaitik aspaldi hurrengoa komentatu nuen
c Hurrengoa lerroa originala da
c       IF (NIV.EQ.0.AND.EEAVT(IBT).LE.1.E-8) GOTO 3
c Eta hurrengoa nik idatzi nuen kinematic limit-arentzat aspaldi...
c Baina azkenean limite kinematikoa egiteko ez nuen IF statement hau
c jartzen
c       IF (NIV.EQ.0.AND.EEAVT(IBT).LE.1.E-12) GOTO 3
        SCAL=EEAVE(I)/EEAVT(IBT)
        IF (NIV.EQ.1) SCAL=EEAVE(NBE+1)/EEAVT(NBE+1)
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
            IF (NOUT2.LE.NBB) WRITE (NOUT2,*) E,ZI
c            IF (NOUT2.LE.NBB) WRITE (NOUT2,*) E,AT(IBT,K)
4         CONTINUE
          IF (NOUT2.LE.NBB)WRITE (NOUT2,*)
          NOUT2=NOUT2+1
3     CONTINUE
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

C =========================================================================
C
C Subroutine YPEND calculates the Pendry Y function Y=(A/AP)/((A/AP)**2+VI)
C where AP/A is the logarithmic derivative of the tabulated function A.
C
C =========================================================================
C
      SUBROUTINE YPEND(A,AP,APP,NBD,NB,NE,E,Y,VI,IPR,IEERG,EEAVE)
C
      DIMENSION A(NBD,IEERG),AP(NBD,IEERG),APP(NBD,IEERG),NE(NBD)
      DIMENSION E(NBD,IEERG),Y(NBD,IEERG)
      DIMENSION EEAVE(30)
C
      COMMON /VINY/VMIN,VMAX,DV,EINCR,THETA(1),FI(1)
      COMMON /WIV2/PERSH,NIV,NSE1(30),NSE2(30)
C
40    FORMAT (26H PENDRY Y FUNCTION IN BEAM,1I3,/,50(5(1F7.2,1E13.4,3X)
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

C=======================================================================
C
C  Subroutine ZGT performs Gaussian elimination as the fist step in the 
C  solution of a system of linear equations. This is used to multiply
C  the inverse of a matrix into a vector, the multiplication being done
C  later by subroutine ZST.
C
C  Modified version of routine ZGE from the Van Hove/Tong LEED package.
C  Modifications by WANDER.
C
C=======================================================================
C
      SUBROUTINE ZGE(A,INT,NR,NC,EMACH)
C
      COMPLEX A,YR,DUM
      DIMENSION A(NR,NC),INT(NC)
C
      N=NC
      DO 680 II=2,N
         I=II-1
         YR=A(I,I)
         IN=I
         DO 600 J=II,N
            IF (CABS(YR).LT.CABS(A(J,I))) THEN
               YR=A(J,I)
               IN=J
            ENDIF
600      CONTINUE
         INT(I)=IN
         IF (IN.NE.I) THEN
            DO 620 J=I,N
               DUM=A(I,J)
               A(I,J)=A(IN,J)
               A(IN,J)=DUM
620         CONTINUE
         ENDIF
         IF (CABS(YR).GT.EMACH) THEN
            DO 670 J=II,N
               IF (CABS(A(J,I)).GT.EMACH) THEN
                  A(J,I)=A(J,I)/YR
                  DO 660 K=II,N
                     A(J,K)=A(J,K)-A(I,K)*A(J,I)
660               CONTINUE
               ENDIF
670         CONTINUE
         ENDIF
680   CONTINUE
      RETURN
      END

C=======================================================================
C
C  Subroutine ZST terminates the solution of a system of linear
C  equations initiated by subroutine ZGT, by back-substituting the
C  constant vector.
C
C Modified version of routine ZSU from the Van Hove/Tong LEED package. 
C Modifications by WANDER.
C
C=======================================================================
C
      SUBROUTINE ZSU(A,INT,X,NR,NC,EMACH)
C
      COMPLEX A,X,DUM
      DIMENSION A(NR,NC),X(NC),INT(NC)
C
      N=NC
      DO 730 II=2,N
         I=II-1
         IF (INT(I).NE.I) THEN
            IN=INT(I)
            DUM=X(IN)
            X(IN)=X(I)
            X(I)=DUM
         ENDIF
         DO 720 J=II,N
            IF (CABS(A(J,I)).GT.EMACH) X(J)=X(J)-A(J,I)*X(I)
720      CONTINUE
730   CONTINUE
      DO 780 II=1,N
         I=N-II+1
         IJ=I+1
         IF (I.NE.N) THEN
            DO 750 J=IJ,N
               X(I)=X(I)-A(I,J)*X(J)
750         CONTINUE
         ENDIF
         IF (CABS(A(I,I))-EMACH*1.0E-5.LT.0) A(I,I)=EMACH*1.0E-5*(1.0,
     &    1.0)
         X(I)=X(I)/A(I,I)
780   CONTINUE
      RETURN
      END
