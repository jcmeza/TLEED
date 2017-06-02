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
      PARAMETER (MAXC=100 ,tol=3.e-2)
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
99     do 12 j=1,NNDIM
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

C?
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
     & NDIML,DIREC,NLTIN,LPOINT)
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
!300   CONTINUE
13    continue
      FRAC=2.*abs(fp-fret)/(abs(fp)+abs(fret))
C
C are we done?
C
cjcm
      write(*,*) 'powell: ', iter, fret, frac
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
      write(*,*) 'powell: converged ', iter, fret, frac
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
! 152  CONTINUE

      powell = fret

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
