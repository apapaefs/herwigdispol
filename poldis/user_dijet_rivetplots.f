C#######################################################################
C     POLDIS: USER RUTINE
C            . DEFINE MAIN SETTING FOR THE CALCULATION
C            . DEFINE CUTS AND HISTOGRAMS
C#######################################################################
      IMPLICIT NONE
      INTEGER NEV,IPOL,IFRAME,ICH,INNLO,RUNS,IFL
      INTEGER SCHEME,NF,IALG,IBOSON
      DOUBLE PRECISION S,Ee,Ep
      DOUBLE PRECISION R(3),RAD,ALPHAE
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LEPCH,THETA,ZMASS,WMASS,ZGAMMA,WGAMMA,CKM2(3,3)
      EXTERNAL USER,CUTS
      CHARACTER *40 TITLE
      COMMON/TITLEPLOT/TITLE
      COMMON/NEV/NEV
      COMMON/S/S
      COMMON/EeEp/Ee,Ep
      COMMON/FLAGPOL/IPOL
      COMMON/FLAGFRAME/IFRAME
      COMMON/CHANNEL/ICH,IFL
      COMMON/INNLO/INNLO
      DATA R/0,0,0/
      COMMON/JETPARAMS/RAD,IALG
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/ALPHAQED/ALPHAE
      COMMON/EWSETTINGS/LEPCH,THETA,ZMASS,WMASS,ZGAMMA,WGAMMA,CKM2,
     $                  IBOSON

C#######################################################################
C     CHOOSE POLARIZED/UNPOLARIZED XSEC    UNPOL=0, POL=1
C-----------------------------------------------------------------------
      IPOL=1
C-----------------------------------------------------------------------
C     CHOOSE INTERMEDIATE BOSON    ONLY PHOTON = 0 / PHOTON+Z = 1
C                                   ONLY Z = 2     / W = 3
C-----------------------------------------------------------------------
      IBOSON=0
C-----------------------------------------------------------------------
C     CHOOSE CHARGE OF INCIDENT LEPTON
C-----------------------------------------------------------------------
      LEPCH=-1D0
C-----------------------------------------------------------------------
C     CHOOSE LAB FRAME=0/ BREIT FRAME=1 / CMS EP FRAME=2
C-----------------------------------------------------------------------
      IFRAME=1
C-----------------------------------------------------------------------
C     CHOOSE PARTONIC CHANNEL ALL=0/ ONLY QUARKS=1 / ONLY GLUONS=2
C-----------------------------------------------------------------------
      ICH=0
C-----------------------------------------------------------------------
C     CHOOSE INITIAL QUARK FLAVOUR (ONLY FOR FOR ICH=1)
C-----------------------------------------------------------------------
      IFL=0
C-----------------------------------------------------------------------
C     DEFINE PDF SET
      IF (IPOL.EQ.0) THEN
         call InitPDFsetByName("PDF4LHC15_nnlo_100_pdfas")
      ELSEIF (IPOL.EQ.1) THEN
         call InitPDFsetByName("BDSSV24-NNLO")
      ENDIF
      call InitPDF(0)
C-----------------------------------------------------------------------
C     NUMBER OF LIGHT FLAVORS
      NF=5
C-----------------------------------------------------------------------
C     DEFINE LEPTON AND HADRON ENERGIES & CENTER OF MASS ENERGY^2
C     EIC KINEMATICS
      Ee=18
      Ep=275
      S=4*Ee*Ep
C-----------------------------------------------------------------------
C     SET ALPHA_QED VALUE
      ALPHAE=0.00729927d0 ! this is 1/137
C    1/137d0
C-----------------------------------------------------------------------
C     SET EW CONSTANTS
      THETA=0.490194163d0 !0.495674D0 
      ZMASS=91.1876d0
      ZGAMMA=2.4952d0
      WMASS=80.379d0
      WGAMMA=2.085d0

      CKM2(1,1)=(0.9737d0)**2
      CKM2(1,2)=(0.2245d0)**2
      CKM2(1,3)=(0.00382d0)**2
      CKM2(2,1)=(0.2210d0)**2
      CKM2(2,2)=(0.987)**2
      CKM2(2,3)=(0.041d0)**2
      CKM2(3,1)=(0.008d0)**2
      CKM2(3,2)=(0.0388d0)**2
      CKM2(3,3)=(1.013d0)**2
C-----------------------------------------------------------------------
C     DEFINE JET ALGORITHM AND PARAMETERS
      ialg = 1
      RAD  = 1d0
C#######################################################################
C     RUN PARAMETETERS
C-----------------------------------------------------------------------
      INNLO=1
C-----------------------------------------------------------------------
      NEV=  20 000 000
C#######################################################################
C     INITIALISE SEED
      CALL RANDOM_SEED()
      R(1)=14217136
      CALL RANGEN(0,R)
C-----------------------------------------------------------------------
      IF (IPOL .EQ. 1) THEN
        TITLE = "dijets_pol.top"
      ELSE
        TITLE = "dijets_unpol.top"
      ENDIF
C#######################################################################
      CALL DEFINEPLOTS
      CALL POLDIS(NEV,S,NF,USER,CUTS)
      CALL PRINTOUT(NEV)
      WRITE (6,'(/A)') 'Thanks for using POLDIS MF, Please cite:'
      WRITE (6,'(2A/)') '  PhysRevLett.125.082001 (2020) ',
     $ '  (https://arxiv.org/abs/2005.10705)'
      END

C#######################################################################
C     DEFINE CUTS FOR THE EVENT GENERATION
      SUBROUTINE CUTS(S,XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX)
      IMPLICIT NONE
      DOUBLE PRECISION S,XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX
      XMIN=1E-5
      XMAX=1
      YMIN=0.2
      YMAX=0.6
      Q2MIN=49
      Q2MAX=2500
      END

C#######################################################################
C     DEFINE RENORMALIZATION AND FACTORIZATION SCALES
      SUBROUTINE SCALES(SCALEREN,SCALEFACT,Q2)
      COMMON/JETS/POUT,zkt,zeta,zphi,n1,n2,n3
      INTEGER n1,n2,n3
      DOUBLE PRECISION POUT(4,7),zkt(1:3),zeta(1:3),zphi(1:3),Q2
      DOUBLE PRECISION SCALEJET,SCALEREN,SCALEFACT,DOT

      Q2=ABS(DOT(Pout,5,5))
C      SCALEJET=(91.1876d0)**2 ! fixed scale for testing
      SCALEJET=Q2
C      SCALEJET=Q2+((zkt(n1)+zkt(n2))/2)**2)/2 # AP: removed this from line above to try and reproduce Herwig results
      SCALEREN=SCALEJET
      SCALEFACT=SCALEJET
      END

C#######################################################################
      SUBROUTINE DEFINEPLOTS
      IMPLICIT NONE
      INTEGER I,NPLOTS
      COMMON /DEMCOM/CSUM,CSQR
      COMMON /numberofplots/NPLOTS
      DOUBLE PRECISION CSUM(100),CSQR(100)

C-----------------------------------------------------------------------
C  WRITES HISTOGRAM WITH TOPDRAWER PACKAGE (MAX=100)
C  CALL GBOOK1(I,TITLE,BINS,LOW,UP,NEW,HIST,LOGX,LOGY)
C  TITLE is CHARACTER*10 in this codebase.
C-----------------------------------------------------------------------

C number of histograms to show
C (30 original + 6 jet-balance + 3 F2eff(x) + 15 pre-cut observables = 54)
      NPLOTS=54

C INITIALIZE VECTORS TO COMPUTE INTEGRATED X-SECTION PER HISTOGRAM
      DO I=1,100
        CSUM(I)=0
        CSQR(I)=0
      ENDDO

C     pT1 DISTRIBUTION  (per-bin)
      CALL GBOOK1( 1,"pt1 NNLO",15,  5D0, 30D0, 1,1,0,0)
      CALL GBOOK1( 2,"pt1 NLO ",15,  5D0, 30D0, 0,1,0,0)
      CALL GBOOK1( 3,"pt1 LO  ",15,  5D0, 30D0, 0,1,0,0)

C     eta1 (LAB) DISTRIBUTION  (per-bin)
      CALL GBOOK1( 4,"eta1LNNLO",15,-3.5D0, 3.5D0, 1,1,0,0)
      CALL GBOOK1( 5,"eta1L NLO",15,-3.5D0, 3.5D0, 0,1,0,0)
      CALL GBOOK1( 6,"eta1L LO ",15,-3.5D0, 3.5D0, 0,1,0,0)

C     Q2 DEPENDENCE  (per-bin; linear x)
      CALL GBOOK1( 7,"Q2 NNLO ",100, 49D0,2500D0, 1,1,0,0)
      CALL GBOOK1( 8,"Q2 NLO  ",100, 49D0,2500D0, 0,1,0,0)
      CALL GBOOK1( 9,"Q2 LO   ",100, 49D0,2500D0, 0,1,0,0)

C     xBj DEPENDENCE  (per-bin; linear x)
      CALL GBOOK1(10,"X NNLO  ", 20,  0D0,  1D0, 1,1,0,0)
      CALL GBOOK1(11,"X NLO   ", 20,  0D0,  1D0, 0,1,0,0)
      CALL GBOOK1(12,"X LO    ", 20,  0D0,  1D0, 0,1,0,0)

C     Mjj DISTRIBUTION (per-bin; log-x binning)
      CALL GBOOK1(13,"MM NNLO ", 15, 10D0,100D0, 1,1,1,0)
      CALL GBOOK1(14,"MM NLO  ", 15, 10D0,100D0, 0,1,1,0)
      CALL GBOOK1(15,"MM LO   ", 15, 10D0,100D0, 0,1,1,0)

C     <pT> DISTRIBUTION (per-bin)
      CALL GBOOK1(16,"<pt>NNLO", 15,  5D0, 30D0, 1,1,0,0)
      CALL GBOOK1(17,"<pt> NLO", 15,  5D0, 30D0, 0,1,0,0)
      CALL GBOOK1(18,"<pt> LO ", 15,  5D0, 30D0, 0,1,0,0)

C     eta* DISTRIBUTION (per-bin)
      CALL GBOOK1(19,"eta*NNLO", 15,  0D0,2.5D0, 1,1,0,0)
      CALL GBOOK1(20,"eta* NLO", 15,  0D0,2.5D0, 0,1,0,0)
      CALL GBOOK1(21,"eta* LO ", 15,  0D0,2.5D0, 0,1,0,0)

C     Zeta DISTRIBUTION (per-bin)
      CALL GBOOK1(22,"ZetaNNLO", 12,-1.75D0,-0.25D0, 1,1,0,0)
      CALL GBOOK1(23,"Zeta NLO", 12,-1.75D0,-0.25D0, 0,1,0,0)
      CALL GBOOK1(24,"Zeta LO ", 12,-1.75D0,-0.25D0, 0,1,0,0)

C     pT2 DISTRIBUTION  (per-bin)
      CALL GBOOK1(25,"pt2 NNLO", 15,  5D0, 30D0, 1,1,0,0)
      CALL GBOOK1(26,"pt2 NLO ", 15,  5D0, 30D0, 0,1,0,0)
      CALL GBOOK1(27,"pt2 LO  ", 15,  5D0, 30D0, 0,1,0,0)

C     eta2 (LAB) DISTRIBUTION  (per-bin)
      CALL GBOOK1(28,"eta2LNNLO",15,-3.5D0, 3.5D0, 1,1,0,0)
      CALL GBOOK1(29,"eta2L NLO",15,-3.5D0, 3.5D0, 0,1,0,0)
      CALL GBOOK1(30,"eta2L LO ",15,-3.5D0, 3.5D0, 0,1,0,0)

C     pT2/pT1 ratio DISTRIBUTION  (per-bin)
      CALL GBOOK1(31,"r21 NNLO", 15,  0D0,  1D0, 1,1,0,0)
      CALL GBOOK1(32,"r21 NLO ", 15,  0D0,  1D0, 0,1,0,0)
      CALL GBOOK1(33,"r21 LO  ", 15,  0D0,  1D0, 0,1,0,0)

C     pT asymmetry DISTRIBUTION  (per-bin)
      CALL GBOOK1(34,"asymNNLO", 15,  0D0,  1D0, 1,1,0,0)
      CALL GBOOK1(35,"asym NLO", 15,  0D0,  1D0, 0,1,0,0)
      CALL GBOOK1(36,"asym LO ", 15,  0D0,  1D0, 0,1,0,0)

C     F2eff(x) (per-bin)
      CALL GBOOK1(37,"F2 NNLO ", 20,  0D0,  1D0, 1,1,0,0)
      CALL GBOOK1(38,"F2 NLO  ", 20,  0D0,  1D0, 0,1,0,0)
      CALL GBOOK1(39,"F2 LO   ", 20,  0D0,  1D0, 0,1,0,0)

C     Pre-cut Q2 / xBj / y / pT1 / pT2 distributions
      CALL GBOOK1(40,"Q2PreNNLO",100, 49D0,2500D0, 1,1,0,0)
      CALL GBOOK1(41,"Q2PreNLO", 100, 49D0,2500D0, 0,1,0,0)
      CALL GBOOK1(42,"Q2PreLO",  100, 49D0,2500D0, 0,1,0,0)
      CALL GBOOK1(43,"XPreNNLO",  20,  0D0,  1D0, 1,1,0,0)
      CALL GBOOK1(44,"XPreNLO",   20,  0D0,  1D0, 0,1,0,0)
      CALL GBOOK1(45,"XPreLO",    20,  0D0,  1D0, 0,1,0,0)
      CALL GBOOK1(46,"YPreNNLO",  40,0.2D0,0.6D0, 1,1,0,0)
      CALL GBOOK1(47,"YPreNLO",   40,0.2D0,0.6D0, 0,1,0,0)
      CALL GBOOK1(48,"YPreLO",    40,0.2D0,0.6D0, 0,1,0,0)
      CALL GBOOK1(49,"pt1preNNLO",30,  0D0, 30D0, 1,1,0,0)
      CALL GBOOK1(50,"pt1preNLO", 30,  0D0, 30D0, 0,1,0,0)
      CALL GBOOK1(51,"pt1preLO",  30,  0D0, 30D0, 0,1,0,0)
      CALL GBOOK1(52,"pt2preNNLO",30,  0D0, 30D0, 1,1,0,0)
      CALL GBOOK1(53,"pt2preNLO", 30,  0D0, 30D0, 0,1,0,0)
      CALL GBOOK1(54,"pt2preLO",  30,  0D0, 30D0, 0,1,0,0)

      END

C#######################################################################
      SUBROUTINE USER(N,NA,NT,P,S,W)
      IMPLICIT NONE
      INTEGER N,NA,NT,I,J,nplots,IPOL,VEC(4),IFRAME
      DOUBLE PRECISION P(4,7),S,DOT,ETA,Q2,X,Y,TOTWGT,PI,Ee,Ep
      DOUBLE PRECISION CSUM(100),CSQR(100),CS(100),GEV2PB
      COMMON /DEMCOM/CSUM,CSQR
      COMMON /numberofplots/NPLOTS
      INTEGER SCHEME,NF,N1,N2,N3,NEV
      DOUBLE PRECISION zkt(1:3),zeta(1:3),zphi(1:3),W(0:5)
      DOUBLE PRECISION ykt(1:7),yeta(1:7),act(4),Ecurr
      DOUBLE PRECISION F2effLO,F2effNLO,F2effNNLO,Broad,Broadn,Broadd
      DOUBLE PRECISION POUT(4,7),FNORM,POUT2(4,7),QQ,MM,MM2,PTM
      DOUBLE PRECISION POUT3(4,7)
      COMMON/GEV2PB/GEV2PB
      COMMON/JETS/POUT,zkt,zeta,zphi,N1,N2,N3
      COMMON/FLAGPOL/IPOL
      COMMON/NEV/NEV
      COMMON/EeEp/Ee,Ep
      COMMON/FLAGFRAME/IFRAME
      DATA CS/100*0/
      PI=ATAN(1D0)*4

C-----------------------------------------------------------------------
C---FILL THE HISTOGRAMS AT THE END OF EACH EVENT
      IF (N.EQ.0) THEN
        DO I=1,NPLOTS
          CSUM(I)=CSUM(I)+CS(I)
          CSQR(I)=CSQR(I)+CS(I)**2
          CS(I)=0
        ENDDO
        DO J=1,NPLOTS
          CALL GOPERA(200+J,'+',    J,    J,1D0,1D0)
          CALL GOPERA(200+J,'*',200+J,200+J,1D0,1D0)
          CALL GOPERA(200+J,'+',100+J,100+J,1D0,1D0)
          CALL GRESET(200+J)
        ENDDO
        RETURN
      ENDIF
C-----------------------------------------------------------------------
C--- USER ROUTINE STARTS HERE
C-----------------------------------------------------------------------

C     DIS kinematics
      Q2=ABS(DOT(POUT,5,5))
      Y=DOT(POUT,1,5)/DOT(POUT,1,6)
      X=Q2/Y/S

C     Calculate PT's and pseudorapidities in lab frame if Breit chosen
      IF (IFRAME.EQ.1) THEN
        QQ=SQRT(ABS(Q2))
        DO J=2,7
          POUT2(4,J)=POUT(4,J)*(QQ*(2-Y)/4/Ee/Y+Ee*Y/QQ)
     $            -POUT(1,J)*(QQ*SQRT(1-Y)/2/Ee/Y)+POUT(3,J)
     $            *(QQ/4/Ee-Ee*Y/QQ)
          POUT2(1,J)=-POUT(4,J)*(SQRT(1-Y))+POUT(1,J)+
     $               POUT(3,J)*(SQRT(1-Y))
          POUT2(2,J)=POUT(2,J)
          POUT2(3,J)=POUT(4,J)*(QQ*(2-Y)/4/Ee/Y-Ee*Y/QQ)
     $            -POUT(1,J)*(QQ*SQRT(1-Y)/2/Ee/Y)+POUT(3,J)
     $            *(QQ/4/Ee+Ee*Y/QQ)

          YETA(J-1)=0.5D0*DLOG((POUT2(4,J)+POUT2(3,J))
     $              /(POUT2(4,J)-POUT2(3,J)))
          YKT(J-1)=(POUT(1,J)**2+POUT(2,J)**2)**0.5
        ENDDO
      ENDIF

C     Inclusive cross section normalization
      FNORM=2*PI*(1+(1-Y)**2)/Y/Q2**2/137**2*GEV2PB
      IF(IPOL.EQ.1) FNORM=(1-(1-Y)**2)/(1+(1-Y)**2)*FNORM

      F2effLO   = W(0)/FNORM
      F2effNLO  = W(1)/FNORM
      F2effNNLO = W(2)/FNORM

C     Fill F2eff(x) histograms (IDs 37-39)
      CALL GFILL1(37,X,F2effNNLO)
      CALL GFILL1(38,X,F2effNLO)
      CALL GFILL1(39,X,F2effLO)

C     Compute inclusive x-section (kept as in original)
      CS(1)=CS(1)+W(2)
      CS(2)=CS(2)+W(1)
      CS(3)=CS(3)+W(0)

      IF(W(1).LT.W(1))THEN
        WRITE(6,*) 'NaN ALERT'
        STOP
      ENDIF

C---   Di-jet variables
      MM=sqrt((POUT(4,n1+1)+POUT(4,n2+1))**2
     $ -(POUT(1,n1+1)+POUT(1,n2+1))**2
     $ -(POUT(2,n1+1)+POUT(2,n2+1))**2
     $ -(POUT(3,n1+1)+POUT(3,n2+1))**2)

      PTM=(zkt(n1)+zkt(n2))/2

C     Pre-cut inclusive observables: require DIS kinematics and the two
C     leading jets, but not the Breit pT / lab-rapidity dijet acceptance.
      CALL GFILL1(40,Q2,W(2))
      CALL GFILL1(41,Q2,W(1))
      CALL GFILL1(42,Q2,W(0))
      CALL GFILL1(43,X,W(2))
      CALL GFILL1(44,X,W(1))
      CALL GFILL1(45,X,W(0))
      CALL GFILL1(46,Y,W(2))
      CALL GFILL1(47,Y,W(1))
      CALL GFILL1(48,Y,W(0))
      CALL GFILL1(49,zkt(n1),W(2))
      CALL GFILL1(50,zkt(n1),W(1))
      CALL GFILL1(51,zkt(n1),W(0))
      CALL GFILL1(52,zkt(n2),W(2))
      CALL GFILL1(53,zkt(n2),W(1))
      CALL GFILL1(54,zkt(n2),W(0))

      IF (zkt(n1).gt.5d0.and.yeta(n1).gt.-3.5d0
     $ .and.yeta(n1).lt.3.5d0.and.zkt(n2).gt.4d0
     $ .and.yeta(n2).gt.-3.5d0.and.yeta(n2).lt.3.5d0) then

C     pT1/pT2: match Rivet definitions explicitly (n1 hardest, n2 subleading)
        CALL GFILL1(1,zkt(n1),W(2))
        CALL GFILL1(2,zkt(n1),W(1))
        CALL GFILL1(3,zkt(n1),W(0))

        CALL GFILL1(4,yeta(n1),W(2))
        CALL GFILL1(5,yeta(n1),W(1))
        CALL GFILL1(6,yeta(n1),W(0))

        CALL GFILL1(7,Q2,W(2))
        CALL GFILL1(8,Q2,W(1))
        CALL GFILL1(9,Q2,W(0))

        CALL GFILL1(10,X,W(2))
        CALL GFILL1(11,X,W(1))
        CALL GFILL1(12,X,W(0))

        CALL GFILL1(13,MM,W(2))
        CALL GFILL1(14,MM,W(1))
        CALL GFILL1(15,MM,W(0))

        CALL GFILL1(16,PTM,W(2))
        CALL GFILL1(17,PTM,W(1))
        CALL GFILL1(18,PTM,W(0))

        CALL GFILL1(19,ABS(zeta(n1)-zeta(n2))/2,W(2))
        CALL GFILL1(20,ABS(zeta(n1)-zeta(n2))/2,W(1))
        CALL GFILL1(21,ABS(zeta(n1)-zeta(n2))/2,W(0))

        CALL GFILL1(22,DLOG10(x*(1+(MM**2/Q2))),W(2))
        CALL GFILL1(23,DLOG10(x*(1+(MM**2/Q2))),W(1))
        CALL GFILL1(24,DLOG10(x*(1+(MM**2/Q2))),W(0))

        CALL GFILL1(25,zkt(n2),W(2))
        CALL GFILL1(26,zkt(n2),W(1))
        CALL GFILL1(27,zkt(n2),W(0))

        CALL GFILL1(28,yeta(n2),W(2))
        CALL GFILL1(29,yeta(n2),W(1))
        CALL GFILL1(30,yeta(n2),W(0))

        CALL GFILL1(31,zkt(n2)/zkt(n1),W(2))
        CALL GFILL1(32,zkt(n2)/zkt(n1),W(1))
        CALL GFILL1(33,zkt(n2)/zkt(n1),W(0))

        CALL GFILL1(34,(zkt(n1)-zkt(n2))/(zkt(n1)+zkt(n2)),W(2))
        CALL GFILL1(35,(zkt(n1)-zkt(n2))/(zkt(n1)+zkt(n2)),W(1))
        CALL GFILL1(36,(zkt(n1)-zkt(n2))/(zkt(n1)+zkt(n2)),W(0))

C     Compute inclusive dijet x-section (kept as in original)
        CS(7)=CS(7)+W(2)
        CS(8)=CS(8)+W(1)
        CS(9)=CS(9)+W(0)

      ENDIF

 999  END

C#######################################################################
      SUBROUTINE PRINTOUT(NEV)
      IMPLICIT NONE
      INTEGER I,J,NEV
      DOUBLE PRECISION CSUM(100),CSQR(100),SQR,X
      COMMON /DEMCOM/CSUM,CSQR
      COMMON /numberofplots/NPLOTS
      LOGICAL NEW(100),HIST(100),LOGY(100),LOGX(100)
      CHARACTER*10 TITLE(100)
      COMMON /PLOTSCOMT/ NEW,HIST,LOGY,LOGX,TITLE
      INTEGER NPLOTS

C     --- new locals for normalized inclusive totals ---
      DOUBLE PRECISION DN
      DOUBLE PRECISION VARNNLO,VARNLO,VARLO
      DOUBLE PRECISION SIGNNLO,SIGNLO,SIGLO
      DOUBLE PRECISION ERRNNLO,ERRNLO,ERRLO

      SQR(X)=SIGN(SQRT(ABS(X)),X)

C     Build TopDrawer output (as in original)
      DO I=1,NPLOTS
        CALL GOPERA(    I,'*',    I,200+I,1D0,1D0)
        CALL GOPERA(100+I,'-',200+I,100+I,1D0,1D0/NEV)
        CALL GOPERA(100+I,'S',    0,100+I,1D0,0D0)
        CALL GTOPER(I)
      ENDDO

      DO I=1,NPLOTS
        WRITE(*,'(1X,(A,F15.5,A,F15.5))') TITLE(I), CSUM(I),
     $  '    +-', SQR( CSQR(I)-CSUM(I)**2/DBLE(NEV) )
      ENDDO

C-----------------------------------------------------------------------
C  Inclusive total cross sections (normalized by NEV)
C  1: NNLO (W(2)), 2: NLO (W(1)), 3: LO (W(0))
C-----------------------------------------------------------------------
      DN = DBLE(NEV)

      VARNNLO = CSQR(1) - CSUM(1)*CSUM(1)/DN
      VARNLO  = CSQR(2) - CSUM(2)*CSUM(2)/DN
      VARLO   = CSQR(3) - CSUM(3)*CSUM(3)/DN

C     Guard against tiny negative values from rounding
      IF (VARNNLO .LT. 0D0) VARNNLO = 0D0
      IF (VARNLO  .LT. 0D0) VARNLO  = 0D0
      IF (VARLO   .LT. 0D0) VARLO   = 0D0

      SIGNNLO = CSUM(1)
      SIGNLO  = CSUM(2)
      SIGLO   = CSUM(3)

      ERRNNLO = SQRT(VARNNLO)
      ERRNLO  = SQRT(VARNLO)
      ERRLO   = SQRT(VARLO)

      WRITE(*,'(A,F15.6,A,F15.6)') 'NNLO = ', SIGNNLO, '  +- ', ERRNNLO
      WRITE(*,'(A,F15.6,A,F15.6)') 'NLO  = ', SIGNLO,  '  +- ', ERRNLO
      WRITE(*,'(A,F15.6,A,F15.6)') 'LO   = ', SIGLO,   '  +- ', ERRLO

      END
C#######################################################################
