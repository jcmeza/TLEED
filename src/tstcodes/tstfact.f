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
         DFACT=EXP(-X)*(10.0D0**((X-0.5D0)*LOG10(X)-(X-1.0D0)))*
     &    (SQRT(6.283185307179586D0))*(1.0+(1.0/(12.0*X))+(1.0/
     &    (288.0D0*(X**2)))-(139.0D0/(51840.0D0*(X**3)))-(571.0D0/
     &    (2488320.0D0*(X**4))))
!         FACT=SNGL(DFACT)
         FACT=DFACT
      ELSE
         IF (L.EQ.0) FACT=1.0
         IF (L.EQ.1) FACT=0.1
         IF (L.EQ.2) FACT=0.02
         IF (L.EQ.3) FACT=6.0*0.001
         IF (L.EQ.4) FACT=24.0*0.0001
      ENDIF

      RETURN
      END
