C-----------------------------------------------------------------------
C    CALCULATES STRUCTURES FUNCTIONS (NOT NEEDED IN POLDIS)
C-----------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER NF,IPOL,IORDER,I,IFL,IF3,EWFLAG,IXQ,IGL
       DOUBLE PRECISION X,Q2,MUR,MUF,Q,LEPCHARGE
       COMMON/FLAGPOL/IPOL,IORDER,NF,IFL,IF3,IGL
       COMMON/FLAGEW/LEPCHARGE,EWFLAG
       COMMON/PARAM/X,MUF,MUR,Q2
       DOUBLE PRECISION RES0,RES1,RES2,DINTEG,ER,XI,STR,MIN,MAX,STEP
       EXTERNAL STR
       DOUBLE PRECISION F(-6:6),F2
C-----------------------------------------------------------------------
C    CHOOSE POLARIZED/UNPOLARIZED STRUCTURE FUNCTIONS   UNPOL=0, POL=1
C-----------------------------------------------------------------------
      IPOL=1   ! 2*X*G1
C---- THE LEPCHARGE DETERMINES THE CHARGE OF THE LEPTON
      LEPCHARGE=-1D0
C---- THE ELECTROWEAK FLAG VALUE DETERMINES THE EXCHANGED BOSON:
C         EWFLAG=0 : PHOTON
C         EWFLAG=1 : PHOTON/Z
C         EWFLAG=2 : Z
C         EWFLAG=3 : "W-" (INITIAL ELECTRON - LEPCHARGE = -1)
C                    "W+" (INITIAL POSITRON - LEPCHARGE = +1)
      EWFLAG=0
C-----------------------------------------------------------------------
C    CHOOSE F2, FL OR F3 IN UNPOLARIZED CASE, OR G2 OR G4+G5 IN POLARIZED
C-----------------------------------------------------------------------
      IFL=0  ! F2,G2: IFL=0, FL: IFL=1
      IF3=0  ! F3,G4+GL: IF3=1, F2,FL,G2: IF3=0
      IGL=0  ! ONLY FOR POLARIZED AND IF3=1! G4: IGL=0, GL: IGL=1
c         G2=2*X*G1, G5 = (G4-GL)/(2*X)
C-----------------------------------------------------------------------
C    CHOOSE TO CALCULATE AS A FUNCTION OF X OR Q2 
C-----------------------------------------------------------------------
      IXQ=1  ! X: IXQ=0, Q2: IXQ=1
C-----------------------------------------------------------------------
C---DEFINE PDF SET 
      if(ipol.eq.0) call InitPDFsetByName("PDF4LHC15_nnlo_100_pdfas")
c      if(ipol.eq.0) call InitPDFsetByName("PDF4LHC15_nlo_100_pdfas")
c      if(ipol.eq.0) call InitPDFsetByName("NNPDF30_lo_as_0118")
      if(ipol.eq.1) call InitPDFsetByName("DSSV_REP_LHAPDF6")
      call InitPDF(0) 

C     GAUSSIAN ERROR:      
      ER=1D-4
      NF=3
c      Q=10
c      Q2=Q**2
      IF (IXQ.EQ.0) THEN
C       FIXED Q2 VALUES 
        Q2 = 100
        Q = SQRT(Q2)
C       X LIMITS: 
        MIN = -5 !LOGARITHMIC SCALE!
        MAX = 0
        STEP = 0.5D0
      ELSEIF (IXQ.EQ.1) THEN
C       FIXED X VALUES 
        X = 0.1
C       Q2 LIMITS: 
        MIN = 0!0.301030264 !0.5 !LOGARITHMIC SCALE!
        MAX = 5!4.0001 !3
        STEP = 0.5!0.410996401 !0.3D0
      ELSE
        WRITE(6,*) 'FATAL ERROR'
      ENDIF

C   CALCULATION HEADER DISPLAY
    2 FORMAT (A,$) 
      IF (IPOL.EQ.0) THEN
        IF (IF3.EQ.0) THEN
          IF (IFL.EQ.0) WRITE(6,2) "CALCULATING F2"
          IF (IFL.EQ.1) WRITE(6,2) "CALCULATING FL"
        ELSE
          WRITE(6,2) "CALCULATING F3"
        ENDIF
      ELSE
        IF (IF3.EQ.1) THEN
          IF (IGL.EQ.0) WRITE(6,2) "CALCULATING G4"
          IF (IGL.EQ.1) WRITE(6,2) "CALCULATING GL"
        ELSE
          WRITE(6,2) "CALCULATING 2*X*G1"
        ENDIF
      ENDIF

      IF (IXQ.EQ.0) WRITE(6,*) "AS A FUNCTION OF X AT Q2 =", Q2
      IF (IXQ.EQ.1) WRITE(6,*) "AS A FUNCTION OF Q2 AT X =", X
 

c      do x=0.001,0.5,0.05
      do XI=MIN,MAX,STEP

        IF (IXQ.EQ.0) X=10**XI
        IF (IXQ.EQ.1) THEN
            Q2=10**XI
            Q = SQRT(Q2)
        ENDIF

c        SCALES 
         MUF=Q!**2
         MUR=Q!**2
      
      IORDER=0
      res0=dinteg(str,0d0,1d0,ER)
      IORDER=1
      res1=dinteg(str,0d0,1d0,ER)   
      IORDER=2
      res2=dinteg(str,0d0,1d0,ER)   
        
            WRITE(*,'(7F12.6)')10**XI,RES0,RES1,res2,res1/res0,res2/res1
      ENDDO
      END




      FUNCTION STR(XX)
      IMPLICIT NONE
      DOUBLE PRECISION XX(1)
      DOUBLE PRECISION F(-6:6),F1(-6:6),EQ(-6:6),MUR,MUF,CKM2(3,3),
     $                 CKM(3,3),COMBPDFVEC,COMBPDFAX 
      DOUBLE PRECISION LOGZ,DIMZ,TRIMZ,SMZ,LOG1MZ,DI1MZ,DIZ,TRIZ,
     $       TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,ZETA2,ZETA3   
      DOUBLE PRECISION C2QNSMIN_REG,DLOG,C2QPS_REG,C2G_REG
      DOUBLE PRECISION C2QNSPLUS_SING,C2QNSPLUS_REG,C2QNSPLUS_DELTA
      DOUBLE PRECISION  L3,L2,L1,L0,L3X,L2X,L1X,L0X   
      DOUBLE PRECISION C2QNS_SING2,C2QNS_DELTA2,C2QNS_REG2P,C2QNS_REG2M,
     $        C2QS_SING2,C2QS_DELTA2,C2QS_REG2,C2G_REG2,
     $   G2NNLOq,G2NNLOg,G2NNLO, C2QNS_SING, C2QNS_REG,C2QNS_DELTA  
      DOUBLE PRECISION STR,X,Z,U,PI,CF,CA,LOGMU,TR
      DOUBLE PRECISION A2PI,SIGMA,SIGMA1,DELTA,DELTA1,ZJAC,GL,
     $     F2NLO,FLNLO,F2LO,FLLO,G2NLOg,G2NLOQ,
     $     F2NLOq,FLNLOq,F2NLOg,FLNLOg, SUMC,
     $     G2NLO,G2LO,ALPHAS,ALPHASPDF,ZMAX,CUTOFF,
     $     CLqs1,C2qs_reg1,C2qs_sing1,C2qs_delta1,
     $     GL1, CLqns1,CLg1,C2qns_reg1,
     $     C2qns_sing1,C2qns_delta1,C2g1,KKK,
     $     CLQNSP_delta2,CLQNSM_delta2,C2G_DELTA2,
     $     F2NNLOq,F2NNLOg,F2NNLO,CLQNS_delta2,
     $     C3qns_sing2,C3qns_delta2,C3qns_reg2P,C3qns_reg2M,
     $     F3NNLOq,G4NNLOq,GLNNLOq,
     $     C4QNS_SING2,C4QNS_DELTA2,C4QNS_REG2P,C4QNS_REG2M,C4QPS_REG,
     $     C4QS_SING2,C4QS_DELTA2,C4QS_REG2 
      DOUBLE PRECISION FLNNLOq,FLNNLOG,FLNNLO,CLqns2,CLqs2,CLG2     
      DOUBLE PRECISION SCREN,SCFACT,Q2,LEPCHARGE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,GZ,MZ,THETAW,
     #                 GAMMAZ,PROPZ,INTGZ,SIGMAW,DELTAW,SIGMAW1,DELTAW1
      DOUBLE PRECISION F3LO,CF3,CF3min,CF2min,CF31,C3qns_reg1,
     $     C3qns_sing1,C3qns_delta1,F3NLOq,SUMCEW,DOT,C4qns_reg1,
     $     C4qns_sing1,
     $     C4qns_delta1,C4qs_reg1,C4qs_sing1,C4qs_delta1,GLNLOq,G4NLOq,
     $     G4LO,GLLO,SUMCEWAX,GLNLOg,pisq    
                  
      INTEGER NF,lpol,iorder,I,J,N,IFL,vogt,EWFLAG,IF3,IGL,IND(-6:6)
      COMPLEX*16 WGPLG
      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/PARAM/X,MUF,MUR,Q2
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL
      COMMON/FLAGEW/LEPCHARGE,EWFLAG
      COMMON/EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,CKM2,EQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3  

      vogt=0
      PI=ATAN(1D0)*4
      pisq=pi**2
      ZETA2=PI**2/6.0D0
      ZETA3=1.2020569031D0
      CF=4/3D0
      CA=3
      TR=1/2D0
      CUTOFF=1d-8

C DEFINE CHARGES          
      EQ(0)=0
      EQ(1)=-1D0/3
      EQ(2)=EQ(1)+1
      DO I=1,6
        IF (I.GT.2) EQ(I)=EQ(I-2)
        EQ(-I)=-EQ(I)
      ENDDO
      STR=0D0

C---- CLEAR CHARGES IF THERE IS NO PHOTON
      IF (EWFLAG.GE.2) THEN
        DO I=-6,6
          EQ(I)=0D0
        ENDDO
      ENDIF

C---- DEFINE ELECTROWEAK COUPLINGS
      THETAW=0.501723d0
      IF (EWFLAG.LE.2) THEN
        MZ=91.188d0
        GAMMAZ=2.4952d0
      ELSEIF (EWFLAG.EQ.3) THEN
        MZ=80.385d0
        GAMMAZ=2.085d0
      ENDIF

      CKM(1,1)=(0.9743d0)**2
      CKM(1,2)=(0.2254d0)**2      
      CKM(1,3)=(0.0036d0)**2
      CKM(2,1)=(0.2252d0)**2
      CKM(2,2)=(0.9734)**2
      CKM(2,3)=(0.0414d0)**2
      CKM(3,1)=(0.0089d0)**2
      CKM(3,2)=(0.0405d0)**2
      CKM(3,3)=(0.9991d0)**2
     
      CVZ(0)=0
      CAZ(0)=0
      IF (EWFLAG.LT.3) THEN
        CVZ(1)=-1.D0/2/SIN(2*THETAW)+1.D0/3*TAN(THETAW)
        CVZ(2)=1.D0/2/SIN(2*THETAW)-2.D0/3*TAN(THETAW) 
        CAZ(1)=-1.D0/2/SIN(2*THETAW)
        CAZ(2)=1.D0/2/SIN(2*THETAW) 
        CVZE=(-1.D0/2/SIN(2*THETAW)+TAN(THETAW))
        CAZE=-1.D0/2/SIN(2*THETAW)
      ELSE
C     W BOSON STRUCTURE FUNCTIONS DONT INCLUDE PROPAGATOR OR COUPLING
C     FACTORS, BECAUSE REASONS...        
c        CVZ(1)=1.D0/2/SQRT(2d0)/SIN(THETAW)
c        CVZ(2)=1.D0/2/SQRT(2d0)/SIN(THETAW) 
c        CAZ(1)=1.D0/2/SQRT(2d0)/SIN(THETAW)
c        CAZ(2)=1.D0/2/SQRT(2d0)/SIN(THETAW)
c        CVZE=1.D0/2/SQRT(2d0)/SIN(THETAW)
c        CAZE=1.D0/2/SQRT(2d0)/SIN(THETAW)
        CVZ(1)=1D0/SQRT(2d0)
        CVZ(2)=1D0/SQRT(2d0)
        CAZ(1)=1D0/SQRT(2d0)
        CAZ(2)=1D0/SQRT(2d0)
        CVZE=1D0/SQRT(2d0)
        CAZE=1D0/SQRT(2d0)
      ENDIF

c     write(6,*) CVZ(1), CVZ(2), CAZ(1), CAZ(2)

      DO I=1,6
        IF (I.GT.2) THEN
        CVZ(I) = CVZ(I-2)
        CAZ(I) = CAZ(I-2)
        ENDIF        
        CVZ(-I) = -CVZ(I)
        CAZ(-I) = CAZ(I)
      ENDDO

C---  CREATE EMPTY CKM MATRIX
      DO I=1,3
        DO J=1,3
          CKM2(I,J)=0D0
        ENDDO
      ENDDO

C     ADD NON-EXCLUDED FLAVORS        
      IF (EWFLAG.EQ.3) THEN
        DO I=1,FLOOR(NF/2.d0)
          DO J=1,CEILING(NF/2.d0)
            CKM2(I,J)=CKM(I,J)
          ENDDO
        ENDDO
      ENDIF

C---  Z PROPAGATOR AND WIDTH
      PROPZ=Q2**2/((-Q2-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(-Q2-MZ**2)/(-Q2)
      IF (EWFLAG.EQ.3) THEN
        PROPZ=1D0
        INTGZ=-1d0
      ENDIF

C SCALES 
      SCREN=MUR**2
      SCFACT=MUF**2
      ALPHAS= alphasPDF(MUR)
      A2PI=ALPHAS/2/PI
C CALL LHAPDF AND COMPUTE PDFS
      CALL EVOLVEPDF(X,MUF,F1)

C     ======================= PDF COMBINATIONS =========================

C     CREATE NEW MOMENTUM FRACTION FOR HIGHER ORDERS AND BUILD ITS PDF
       N=1                
       ZMAX=1-CUTOFF/100
       Z=X*(ZMAX/X)**(XX(1)**N)
       ZJAC=LOG(ZMAX/X)*N*(XX(1))**(N-1)*Z
                     
      U=X/Z

      CALL EVOLVEPDF(U,MUF,F)          

C     BUILD THE PDF'S COMBINATIONS:
      SIGMA1=0
      DELTA1=0         
      SIGMA=0
      DELTA=0
      SIGMAW1=0
      DELTAW1=0         
      SIGMAW=0
      DELTAW=0     
      GL1=0
      GL=0
      CF3=0
      CF31=0

      SUMC=0
      SUMCEW=0
      SUMCEWAX=0

C     -- CHARGE SUMS --
      DO I=1,NF
        SUMC = SUMC + EQ(I)**2
        IF ((EWFLAG.GE.1).and.(EWFLAG.LT.3)) THEN
          SUMCEW=SUMCEW+PROPZ*(CVZE**2+CAZE**2)*(CVZ(I)**2+CAZ(I)**2)
     #                 +2*EQ(I)*PROPZ*INTGZ*(CVZ(I)*CVZE)
          SUMCEWAX = SUMCEWAX + PROPZ*(4*CVZ(I)*CAZ(I)*CVZE*CAZE)
     #               +2*EQ(I)*PROPZ*INTGZ*(CAZ(I)*CAZE)
        ENDIF
      ENDDO

C     W BOSON WORKS DIFFERENTLY (WITH CKM)
      IF (EWFLAG.EQ.3) THEN
        IF (LEPCHARGE.LT.0) THEN
          DO I=1,3
            DO J=1,3
              SUMCEW = SUMCEW + PROPZ*(CVZE**2+CAZE**2)
     $                          *(CVZ(2*I)**2+CAZ(2*I)**2)*CKM2(I,J)
              SUMCEWAX = SUMCEWAX + PROPZ*(4*CVZ(2*I)*CAZ(2*I)
     $                              *CVZE*CAZE)*CKM2(I,J)
            ENDDO
          ENDDO
        ELSE
          DO I=1,3
            DO J=1,3
              SUMCEW = SUMCEW + PROPZ*(CVZE**2+CAZE**2)
     $                          *(CVZ(2*I-1)**2+CAZ(2*I-1)**2)*CKM2(J,I)
              SUMCEWAX = SUMCEWAX + PROPZ*(4*CVZ(2*I-1)*CAZ(2*I-1)
     $                              *CVZE*CAZE)*CKM2(J,I)
            ENDDO
          ENDDO
        ENDIF
      ENDIF 

C     -- SIGMAS --
C     INDEXES OF FLAVOR TO USE
c      IF ((EWFLAG.EQ.3).AND.(LEPCHARGE.LT.0)) THEN
c        IND = (/-5,-3,-1,2,4,6,0,0,0,0,0,0,0/)
c      ELSE IF ((EWFLAG.EQ.3).AND.(LEPCHARGE.GT.0)) THEN
c        IND = (/-6,-4,-2,1,3,5,0,0,0,0,0,0,0/)
c      ELSE
c        IND = (/-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6/)
c      ENDIF

      DO I=-NF,NF
c        IF ((I.NE.0).AND.(ANY(IND==I))) THEN
        IF (I.NE.0) THEN
          SIGMA = SIGMA+ F(I)/U * (SUMC+SUMCEW)/NF
          SIGMA1 = SIGMA1+ F1(I)/X * (SUMC+SUMCEW)/NF
          
          IF (EWFLAG.GE.1) THEN
            SIGMAW = SIGMAW + (SUMCEWAX/NF)*F(I)/U
            SIGMAW1 = SIGMAW1 + (SUMCEWAX/NF)*F1(I)/X
          ENDIF
        ENDIF
      ENDDO

C     -- DELTAS --
      DELTA = DELTA     + COMBPDFVEC(NF,Q2,F,0)/U - SIGMA
      DELTA1 = DELTA1   + COMBPDFVEC(NF,Q2,F1,0)/X - SIGMA1
      DELTAW = DELTAW   + COMBPDFAX(NF,Q2,F,0)/U - SIGMAW
      DELTAW1 = DELTAW1 + COMBPDFAX(NF,Q2,F1,0)/X - SIGMAW1

C     COMBINATIONS FOR MINUS COEFFICIENTS (ONLY RELEVANT FOR W)      
      CF2min = COMBPDFVEC(NF,Q2,F,1)/U
      CF3min = COMBPDFAX(NF,Q2,F,1)/U

C     POR FAVOR CAMBIAR A FUTURO      
      CF3 = DELTAW + SIGMAW
      CF31 = DELTAW1 + SIGMAW1

C     -- GLUON -- (SUMCEWAX DOESNT CONTRIBUTE HERE!)
      GL1 = GL1+ F1(0)/X * (SUMC+SUMCEW)/NF
      GL  = GL +  F(0)/U * (SUMC+SUMCEW)/NF 

C     ====================== STRUCTURE FNCTIONS ========================
C LO - ZERO ORDER IN ALPHAS

C COMPUTE STRUCTURE FUNCTIONS
      IF (LPOL.EQ.0) THEN
        F2LO=(DELTA1+SIGMA1)*X
        FLLO=0D0
        F3LO=(DELTAW1+SIGMAW1)*X

       STR=F2LO                
       IF (IFL.EQ.1) STR=FLLO
       IF (IF3.EQ.1) STR=F3LO

C -- POLARIZED CASE
      ELSEIF (LPOL.Eq.1) THEN
        G2LO=(DELTA1+SIGMA1)*X ! G2 MEANS 2*X*G1
        GLLO=0D0
        G4LO=(DELTAW1+SIGMAW1)*X

        STR=G2LO
        IF (IF3.EQ.1) THEN
            STR=G4LO
            IF (IGL.EQ.1) STR=0 ! GLLO=FLLO=0
        ENDIF

      ENDIF
 

C FIRST ORDER IN ALPHAS     
       IF (IORDER.gt.0) THEN 
         
c     (in units of alphas/4/Pi)
       LOGMU=LOG(Q2/SCFACT)   
       
       IF (LPOL.EQ.0) THEN

C LONGITUDINAL       
       CLqns1=CF*4*z
       CLqs1=CLqns1
       CLg1=TR*16*z*(1-z)*NF
C NON-SINGLET  F2     
       C2qns_reg1=CF*(-2*(1+z)*log(1-z)-2*(1+z**2)/(1-z)*log(z)+
     $             6+4*z-2*(1+z)*LOGMU  )        
       C2qns_sing1=CF*(4*log(1-z)/(1-z)-3/(1-z)+4/(1-z)*LOGMU )
       C2qns_delta1=CF*( -4*PI**2/6-9    
     $         +3*LOGMU -3*LOG(1-X)+2*LOG(1-X)**2+4*LOG(1-X)*LOGMU)
C SINGLET F2
       C2qs_reg1=C2qns_reg1
       C2qs_sing1=C2qns_sing1      
       C2qs_delta1=C2qns_delta1
C NON-SINGLET  F3
       C3qns_reg1=CF*(-2*(1+z)*(log(1-z)-log(z))-4*log(z)/(1-z)+4+2*z
     #                -2*(1+z)*LOGMU) !debug(chequear término LOGMU)        
       C3qns_sing1=CF*(4*log(1-z)/(1-z)-3/(1-z)
     #                 +4/(1-z)*LOGMU)
       C3qns_delta1=CF*(-9-4*ZETA2
     #                  -3*log(1-X)+2*log(1-X)**2
     #                  +4*LOG(1-X)*LOGMU+3*LOGMU)    
C GLUON   F2         
       C2g1=NF*TR*(4*(1-2*z+2*z**2)*(Log(1-z)-Log(z))-4+32*z*(1-z)
     $       + 4*(1-2*z+2*z**2)*LOGMU )
C
C QUARK       
       FLNLOq=X*DELTA*CLqns1/Z*A2PI/2 !ONLY REGULAR PART NS 
     $        + X*SIGMA*CLqs1/Z*A2PI/2 !ONLY REGULAR PART S
c     
       F2NLOq=X*A2PI/2*   (C2qns_reg1*DELTA/Z  !REGULAR PART NS
     $   +    C2qns_sing1*(DELTA/Z-DELTA1)! PLUS DISTRIBUTION NS
     $   + DELTA1* C2qns_delta1/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*A2PI/2*   (C2qs_reg1*SIGMA/Z  !REGULAR PART  S
     $   +    C2qs_sing1*(SIGMA/Z-SIGMA1)! PLUS DISTRIBUTION  S
     $   + SIGMA1* C2qs_delta1/ZJAC )   ! DELTA TERM   S

       F3NLOq=X*A2PI/2*   (C3qns_reg1*CF3/Z  !REGULAR PART NS
     $   +    C3qns_sing1*(CF3/Z-CF31)! PLUS DISTRIBUTION NS
     $   + CF31* C3qns_delta1/ZJAC )   ! DELTA TERM  NS 

C GLUON
       F2NLOg=X*(GL)*C2g1/Z*A2PI/2 !ONLY REGULAR PART  G
       FLNLOg=X*(GL)*CLg1/Z*A2PI/2 !ONLY REGULAR PART  G
C             
       FLNLO=(FLNLOq+FLNLOg)*zjac
       F2NLO=(F2NLOq+F2NLOg)*zjac


       STR=F2LO+F2NLO
       IF (IFL.EQ.1) STR=FLLO+FLNLO 
       IF (IF3.EQ.1) STR=F3LO+F3NLOq*zjac 
       
       ELSEIF (LPOL.EQ.1) THEN
C NON-SINGLET  G2     
       C2qns_reg1=CF*(-2*(1+z)*log(1-z)-2*(1+z**2)/(1-z)*log(z)+
     $             4+2*z-2*(1+z)*LOGMU  )        
       C2qns_sing1=CF*(4*log(1-z)/(1-z)-3/(1-z)+4/(1-z)*LOGMU )
       C2qns_delta1=CF*( -4*PI**2/6-9    
     $         +3*LOGMU -3*LOG(1-X)+2*LOG(1-X)**2+4*LOG(1-X)*LOGMU)
C SINGLET G2
       C2qs_reg1=C2qns_reg1
       C2qs_sing1=C2qns_sing1      
       C2qs_delta1=C2qns_delta1     
C GLUON  G2         
       C2g1=NF*TR*(4*(2*z-1)*(Log(1-z)-Log(z))+4*(3-4*z)
     $       + 4*(2*z-1)*LOGMU )

C EW - PART
C G - LONGITUDINAL (=FL)      
       CLqns1=CF*4*z
       CLqs1=CLqns1
       CLg1=0d0
C NON-SINGLET  G4 (=F2)    
       C4qns_reg1=CF*(-2*(1+z)*log(1-z)-2*(1+z**2)/(1-z)*log(z)+
     $             6+4*z-2*(1+z)*LOGMU  )        
       C4qns_sing1=CF*(4*log(1-z)/(1-z)-3/(1-z)+4/(1-z)*LOGMU )
       C4qns_delta1=CF*( -4*PI**2/6-9    
     $         +3*LOGMU -3*LOG(1-X)+2*LOG(1-X)**2+4*LOG(1-X)*LOGMU)
C SINGLET G4 (=F2)
       C4qs_reg1=C4qns_reg1
       C4qs_sing1=C4qns_sing1      
       C4qs_delta1=C4qns_delta1 

C QUARK            
       G2NLOq=X*A2PI/2*   (C2qns_reg1*DELTA/Z  !REGULAR PART NS
     $   +    C2qns_sing1*(DELTA/Z-DELTA1)! PLUS DISTRIBUTION NS
     $   + DELTA1* C2qns_delta1/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*A2PI/2*   (C2qs_reg1*SIGMA/Z  !REGULAR PART  S
     $   +    C2qs_sing1*(SIGMA/Z-SIGMA1)! PLUS DISTRIBUTION  S
     $   + SIGMA1* C2qs_delta1/ZJAC )   ! DELTA TERM   S

C EW PART 
       GLNLOq=X*DELTAW*CLqns1/Z*A2PI/2  !ONLY REGULAR PART NS 
     $        + X*SIGMAW*CLqs1/Z*A2PI/2  !ONLY REGULAR PART S
c
       G4NLOq=X*A2PI/2*    (C4qns_reg1*deltaw/Z  !REGULAR PART NS
     $   +    C4qns_sing1*(deltaw/Z-deltaw1)! PLUS DISTRIBUTION NS
     $   + deltaw1* C4qns_delta1/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*A2PI/2*    (C4qs_reg1*sigmaw/Z  !REGULAR PART  S
     $   +    C4qs_sing1*(sigmaw/Z-sigmaw1)! PLUS DISTRIBUTION  S
     $   + sigmaw1* C4qs_delta1/ZJAC )   ! DELTA TERM   S

C GLUON
       G2NLOg=X*(GL)*C2g1/Z*A2PI/2  !ONLY REGULAR PART  G
C             
       G2NLO=(G2NLOq+G2NLOg)*zjac

           
       STR=G2LO+G2NLO
       IF (IF3.EQ.1) THEN
         STR=G4LO+G4NLOq*zjac
         IF (IGL.EQ.1) STR=GLNLOq*zjac
       ENDIF

       ENDIF
       ENDIF


C SECOND ORDER IN ALPHAS

       IF (IORDER.eq.2) THEN

      LOGZ=DLOG(Z)  
C      
      DIMZ=WGPLG(1,1,-Z)
      TRIMZ=WGPLG(2,1,-Z)
      SMZ=WGPLG(1,2,-Z)
      DIZ=WGPLG(1,1,Z)
      TRIZ=WGPLG(2,1,Z)   
C           
      LOG1MZ=DLOG(1.0D0-Z)
      DI1MZ=WGPLG(1,1,1.0D0-Z)
      TRI1MZ=WGPLG(2,1,1.0D0-Z)      
      S1MZ=WGPLG(1,2,1.0D0-Z)      
C           
      LOG1PZ=DLOG(1.0D0+Z)
C           
      TRIPCO=WGPLG(2,1,(1.0D0-Z)/(1.0D0+Z))
      TRIMCO=WGPLG(2,1,-(1.0D0-Z)/(1.0D0+Z))

      L3=LOG1MZ**3/(1-Z)
      L2=LOG1MZ**2/(1-Z)
      L1=LOG1MZ**1/(1-Z)
      L0=1D0/(1-Z)

      L0X=DLOG(1-X)    
      L1X=DLOG(1-X)**2/2    
      L2X=DLOG(1-X)**3/3    
      L3X=DLOG(1-X)**4/4     

       IF (LPOL.EQ.0) THEN


C NONSINGLET PART    
      CALL C2NS(Z,VOGT,C2QNS_SING2,C2QNS_DELTA2,C2QNS_REG2P,C2QNS_REG2M)    

C SINGLET PART
      CALL C2PS(Z,VOGT,C2QPS_REG)

      C2QS_SING2 = C2QNS_SING2
      C2QS_DELTA2 = C2QNS_DELTA2
      C2QS_REG2 = C2QNS_REG2P + C2QPS_REG

C GLUON PART
      CALL C2G(Z,VOGT,C2G_REG2,C2G_DELTA2)


C F2 NNLO STRUCTURE FUNCTION
       F2NNLOq=X*(A2PI/2)**2 * (C2qns_reg2p*DELTA/Z  !REGULAR PART NS
     $   + C2qns_reg2m*CF2min/Z          !REGULAR PART MINUS
     $   + C2qns_sing2*(DELTA/Z-DELTA1)  ! PLUS DISTRIBUTION NS
     $   + DELTA1* C2qns_delta2/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*(A2PI/2)**2 *  (C2qs_reg2*SIGMA/Z  !REGULAR PART  S
     $   +    C2qs_sing2*(SIGMA/Z-SIGMA1)! PLUS DISTRIBUTION  S
     $   + SIGMA1* C2qs_delta2/ZJAC )   ! DELTA TERM   S

       F2NNLOg=X*((GL)*C2G_REG2/Z  !ONLY REGULAR PART  G
     $   + GL1* C2G_delta2/ZJAC ) *(A2PI/2)**2  ! DELTA TERM  GLU VOGT

       F2NNLO=(F2NNLOq+F2NNLOg)*zjac


C LONGITUDINAL STRUCTURE FUNCTION
      CALL CFL(Z,VOGT,CLQS2,CLQNS2,CLQNSP_DELTA2,CLG2)

         FLNNLOq=X*(DELTA*CLqns2/Z+ !REGULAR PART NS
     $     DELTA1*(CLQNSP_delta2)/ZJAC )*(A2PI/2)**2 !DELTA PART 
     $        + X*(SIGMA*(CLqs2+CLqns2)/Z+ !REGULAR PART S
     $     SIGMA1*(CLQNSP_delta2)/ZJAC )*(A2PI/2)**2 !DELTA PART S

C GLUON
       FLNNLOg=X*((GL)*CLg2/Z) *(A2PI/2)**2   !ONLY REGULAR PART  G
             
       FLNNLO=(FLNNLOq +FLNNLOg )*zjac

C3qns_sing2
C F3 STRUCTURE FUNCTION
      CALL C2POLNS(Z,C3qns_sing2,C3qns_delta2,C3qns_reg2p,C3qns_reg2m) 

       F3NNLOq=X*(A2PI/2)**2 * (C3qns_reg2p*CF3/Z  !REGULAR PART NS
     $   + C3qns_reg2m*CF3min/Z        !REGULAR PART MINUS  
     $   + C3qns_sing2*(CF3/Z-CF31)    ! PLUS DISTRIBUTION NS
     $   + CF31* C3qns_delta2/ZJAC )   ! DELTA TERM  NS 


C STR=F2 BY DEFAULT          
      STR=F2LO+F2NLO*(1+LOG(MUR**2/MUF**2)*(11*CA-2*NF)/(6D0)*(A2PI))
     $          +F2NNLO
C adding term related to change of renormalization scale from van neerven          

C STR=FL IF SELECTED        
       IF (IFL.EQ.1) STR=FLLO+FLNLO*(1+LOG(MUR**2/MUF**2)
     $       *(11*CA-2*NF)/(6D0)*(A2PI))   +FLNNLO
C STR=F3 IF SELECTED        
       IF (IF3.EQ.1) STR=F3LO+F3NLOq*zjac*(1+LOG(MUR**2/MUF**2)
     $       *(11*CA-2*NF)/(6D0)*(A2PI)) + F3NNLOq*zjac



C -- POLARIZED NNLO -- 
       ELSEIF (LPOL.EQ.1) THEN
      
c NON SINGLET G1 (WE KEEP CALLING THE COEFFICIENT C2)
       CALL C2POLNS(Z,C2QNS_SING2,C2QNS_DELTA2,C2QNS_REG2P,C2QNS_REG2M)

c PURE SINGLET G1
       CALL C2POLPS(Z,C2QPS_REG)
        
       C2QS_SING2 = C2QNS_SING2
       C2QS_DELTA2 = C2QNS_DELTA2
       C2QS_REG2 = C2QNS_REG2P + C2QPS_REG
        
C GLUON G1
       CALL C2POLG(Z,C2G_REG2)
        
       G2NNLOq=X*(A2PI/2)**2 * (C2qns_reg2p*DELTA/Z  !REGULAR PART NS
     $   + C2qns_reg2m*CF2min/Z   
     $   + C2qns_sing2*(DELTA/Z-DELTA1)! PLUS DISTRIBUTION NS
     $   + DELTA1* C2qns_delta2/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*(A2PI/2)**2 *  (C2qs_reg2*SIGMA/Z  !REGULAR PART  S
     $   +    C2qs_sing2*(SIGMA/Z-SIGMA1)! PLUS DISTRIBUTION  S
     $   + SIGMA1* C2qs_delta2/ZJAC )   ! DELTA TERM   S

       G2NNLOg=X*(GL)*C2G_REG2/Z*(A2PI/2)**2  !ONLY REGULAR PART  G

       G2NNLO=(G2NNLOq+G2NNLOg)*zjac


C EW PART (NO GLUONS NOR TRIANGLES)
C     NON SINGLET G4 (=F2)    
      CALL C2NS(Z,VOGT,C4QNS_SING2,C4QNS_DELTA2,C4QNS_REG2P,C4QNS_REG2M)

C     PURE SINGLET G4 (=F2)
c       CALL C2PS(Z,VOGT,C4QPS_REG)

       C4QS_SING2 = C4QNS_SING2
       C4QS_DELTA2 = C4QNS_DELTA2
       C4QS_REG2 = C4QNS_REG2P !+ C4QPS_REG

C     GL (=FL)
       CALL CFL(Z,VOGT,CLQS2,CLQNS2,CLQNSP_DELTA2,CLG2)  

       GLNNLOq=X*DELTAW*CLqns2/Z*(A2PI/2)**2 !ONLY REGULAR PART NS 
     $        + X*SIGMAW*CLqns2/Z*(A2PI/2)**2 !ONLY REGULAR PART S
c
       G4NNLOq=X*(A2PI/2)**2* (C4qns_reg2p*deltaw/Z  !REGULAR PART NS
     $   + C4qns_reg2m*CF3min/Z   
     $   + C4qns_sing2*(deltaw/Z-deltaw1)! PLUS DISTRIBUTION NS
     $   + deltaw1* C4qns_delta2/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*(A2PI/2)**2* (C4qs_reg2*sigmaw/Z  !REGULAR PART  S
     $   +    C4qs_sing2*(sigmaw/Z-sigmaw1)! PLUS DISTRIBUTION  S
     $   + sigmaw1* C4qs_delta2/ZJAC )   ! DELTA TERM   S


      STR=G2LO+G2NLO*(1+LOG(MUR**2/MUF**2)*(11*CA-2*NF)/(6D0)*(A2PI))
     $          +G2NNLO
C adding term related to change of renormalization scale from van neerven 
C STR=G4/G5 IF SELECTED        
       IF (IF3.EQ.1) THEN 
        STR=G4LO+(G4NLOq+GLNLOq)*zjac*(1+LOG(MUR**2/MUF**2)
     $       *(11*CA-2*NF)/(6D0)*(A2PI)) + (G4NNLOq+GLNNLOq)*zjac
        IF (IGL.EQ.1) STR=(GLNLOq)*zjac*(1+LOG(MUR**2/MUF**2)
     $       *(11*CA-2*NF)/(6D0)*(A2PI)) + (GLNNLOq)*zjac
       ENDIF

      endif
C 
       endif

      
C     W STRUCTURE FUNCTIONS INCLUDE A *2 FACTOR        
       IF (EWFLAG.EQ.3) STR=2*STR

 100       RETURN
       END





C...DOUBLE PRECISE GAUSS INTEGRATION :
       FUNCTION DINTEG (F, ALFA, BETA, EPS)
       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       DIMENSION W(12), X(12)
       DATA CONST / 1.0 D-12 /
       DATA W
     1  /0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2   0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38647,
     3   0.09515 85116 82492, 0.12462 89712 55533, 0.14959 59888 16576,
     4   0.16915 65193 95002, 0.18260 34150 44923, 0.18945 06104 55068/
       DATA X
     1  /0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2   0.18343 46424 95650, 0.98940 09349 91649, 0.94457 50230 73232,
     3   0.86563 12023 87831, 0.75540 44083 55003, 0.61787 62444 02643,
     4   0.45801 67776 57227, 0.28160 35507 79258, 0.09501 25098 37637/
       DINTEG = 0.0 D0
       IF ( ALFA . EQ. BETA ) RETURN
       A = ALFA
       B = BETA
       DELTA = CONST * (DABS(A-B))
       AA = A
    1  Y = B - AA
       IF( DABS(Y) .LE. DELTA ) RETURN
    2  BB = AA + Y
       C1 = 0.5 D0 * (AA + BB)
       C2 = C1 - AA
       S8 = 0.0 D0
       S16 = 0.0 D0
       DO 15 I = 1, 4
          C3 = X(I) * C2
          S8 = S8 + W(I) * (F(C1+C3) + F(C1-C3))
   15  CONTINUE
       DO 16 I = 5, 12
          C3 = X(I) * C2
          S16 = S16 + W(I) * (F(C1+C3) + F(C1-C3))
  16   CONTINUE
       S8 = S8 * C2
       S16= S16 * C2
       IF( DABS(S16-S8) .GT. EPS * DABS(S8)) THEN
          Y = 0.5 * Y
          IF ( DABS(Y) .LE. DELTA ) THEN
             DINTEG = 0.0
             WRITE (*,10)
  10         FORMAT (1X,' DINTEG : TOO HIGH ACCURACY ')
          ELSE
             GOTO 2
          END IF
       ELSE
          DINTEG = DINTEG + S16
          AA = BB
          GOTO 1
       END IF
       RETURN
       END
c
C    

C   24/08/89 101231638  MEMBER NAME  WGPLG    (ZWPROD.S)    F77
      COMPLEX*16 FUNCTION WGPLG(N,P,X)
 
      INTEGER P,P1,NC(10),INDEX(31)
      DOUBLE PRECISION FCT(0:4),SGN(0:4),U(0:4),S1(4,4),C(4,4)
      DOUBLE PRECISION A(0:30,10)
      DOUBLE PRECISION X,X1,H,ALFA,R,Q,C1,C2,B0,B1,B2,ZERO,HALF
 
      COMPLEX*16 V(0:5),SK,SM
 
      DATA FCT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0/
      DATA SGN /1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0/
      DATA ZERO /0.0D0/, HALF /0.5D0/
      DATA C1 /1.33333 33333 333D0/, C2 /0.33333 33333 3333D0/
 
      DATA S1(1,1) /1.64493 40668 482D0/
      DATA S1(1,2) /1.20205 69031 596D0/
      DATA S1(1,3) /1.08232 32337 111D0/
      DATA S1(1,4) /1.03692 77551 434D0/
      DATA S1(2,1) /1.20205 69031 596D0/
      DATA S1(2,2) /2.70580 80842 778D-1/
      DATA S1(2,3) /9.65511 59989 444D-2/
      DATA S1(3,1) /1.08232 32337 111D0/
      DATA S1(3,2) /9.65511 59989 444D-2/
      DATA S1(4,1) /1.03692 77551 434D0/
 
      DATA C(1,1) / 1.64493 40668 482D0/
      DATA C(1,2) / 1.20205 69031 596D0/
      DATA C(1,3) / 1.08232 32337 111D0/
      DATA C(1,4) / 1.03692 77551 434D0/
      DATA C(2,1) / 0.00000 00000 000D0/
      DATA C(2,2) /-1.89406 56589 945D0/
      DATA C(2,3) /-3.01423 21054 407D0/
      DATA C(3,1) / 1.89406 56589 945D0/
      DATA C(3,2) / 3.01423 21054 407D0/
      DATA C(4,1) / 0.00000 00000 000D0/
 
      DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/
 
      DATA NC /24,26,28,30,22,24,26,19,22,17/
 
      DATA A( 0,1) / .96753 21504 3498D0/
      DATA A( 1,1) / .16607 30329 2785D0/
      DATA A( 2,1) / .02487 93229 2423D0/
      DATA A( 3,1) / .00468 63619 5945D0/
      DATA A( 4,1) / .00100 16274 9616D0/
      DATA A( 5,1) / .00023 20021 9609D0/
      DATA A( 6,1) / .00005 68178 2272D0/
      DATA A( 7,1) / .00001 44963 0056D0/
      DATA A( 8,1) / .00000 38163 2946D0/
      DATA A( 9,1) / .00000 10299 0426D0/
      DATA A(10,1) / .00000 02835 7538D0/
      DATA A(11,1) / .00000 00793 8705D0/
      DATA A(12,1) / .00000 00225 3670D0/
      DATA A(13,1) / .00000 00064 7434D0/
      DATA A(14,1) / .00000 00018 7912D0/
      DATA A(15,1) / .00000 00005 5029D0/
      DATA A(16,1) / .00000 00001 6242D0/
      DATA A(17,1) / .00000 00000 4827D0/
      DATA A(18,1) / .00000 00000 1444D0/
      DATA A(19,1) / .00000 00000 0434D0/
      DATA A(20,1) / .00000 00000 0131D0/
      DATA A(21,1) / .00000 00000 0040D0/
      DATA A(22,1) / .00000 00000 0012D0/
      DATA A(23,1) / .00000 00000 0004D0/
      DATA A(24,1) / .00000 00000 0001D0/
 
      DATA A( 0,2) / .95180 88912 7832D0/
      DATA A( 1,2) / .43131 13184 6532D0/
      DATA A( 2,2) / .10002 25071 4905D0/
      DATA A( 3,2) / .02442 41559 5220D0/
      DATA A( 4,2) / .00622 51246 3724D0/
      DATA A( 5,2) / .00164 07883 1235D0/
      DATA A( 6,2) / .00044 40792 0265D0/
      DATA A( 7,2) / .00012 27749 4168D0/
      DATA A( 8,2) / .00003 45398 1284D0/
      DATA A( 9,2) / .00000 98586 9565D0/
      DATA A(10,2) / .00000 28485 6995D0/
      DATA A(11,2) / .00000 08317 0847D0/
      DATA A(12,2) / .00000 02450 3950D0/
      DATA A(13,2) / .00000 00727 6496D0/
      DATA A(14,2) / .00000 00217 5802D0/
      DATA A(15,2) / .00000 00065 4616D0/
      DATA A(16,2) / .00000 00019 8033D0/
      DATA A(17,2) / .00000 00006 0204D0/
      DATA A(18,2) / .00000 00001 8385D0/
      DATA A(19,2) / .00000 00000 5637D0/
      DATA A(20,2) / .00000 00000 1735D0/
      DATA A(21,2) / .00000 00000 0536D0/
      DATA A(22,2) / .00000 00000 0166D0/
      DATA A(23,2) / .00000 00000 0052D0/
      DATA A(24,2) / .00000 00000 0016D0/
      DATA A(25,2) / .00000 00000 0005D0/
      DATA A(26,2) / .00000 00000 0002D0/
 
      DATA A( 0,3) / .98161 02799 1365D0/
      DATA A( 1,3) / .72926 80632 0726D0/
      DATA A( 2,3) / .22774 71490 9321D0/
      DATA A( 3,3) / .06809 08329 6197D0/
      DATA A( 4,3) / .02013 70118 3064D0/
      DATA A( 5,3) / .00595 47848 0197D0/
      DATA A( 6,3) / .00176 76901 3959D0/
      DATA A( 7,3) / .00052 74821 8502D0/
      DATA A( 8,3) / .00015 82746 1460D0/
      DATA A( 9,3) / .00004 77492 2076D0/
      DATA A(10,3) / .00001 44792 0408D0/
      DATA A(11,3) / .00000 44115 4886D0/
      DATA A(12,3) / .00000 13500 3870D0/
      DATA A(13,3) / .00000 04148 1779D0/
      DATA A(14,3) / .00000 01279 3307D0/
      DATA A(15,3) / .00000 00395 9070D0/
      DATA A(16,3) / .00000 00122 9055D0/
      DATA A(17,3) / .00000 00038 2658D0/
      DATA A(18,3) / .00000 00011 9459D0/
      DATA A(19,3) / .00000 00003 7386D0/
      DATA A(20,3) / .00000 00001 1727D0/
      DATA A(21,3) / .00000 00000 3687D0/
      DATA A(22,3) / .00000 00000 1161D0/
      DATA A(23,3) / .00000 00000 0366D0/
      DATA A(24,3) / .00000 00000 0116D0/
      DATA A(25,3) / .00000 00000 0037D0/
      DATA A(26,3) / .00000 00000 0012D0/
      DATA A(27,3) / .00000 00000 0004D0/
      DATA A(28,3) / .00000 00000 0001D0/
 
      DATA A( 0,4) /1.06405 21184 614 D0/
      DATA A( 1,4) /1.06917 20744 981 D0/
      DATA A( 2,4) / .41527 19325 1768D0/
      DATA A( 3,4) / .14610 33293 6222D0/
      DATA A( 4,4) / .04904 73264 8784D0/
      DATA A( 5,4) / .01606 34086 0396D0/
      DATA A( 6,4) / .00518 88935 0790D0/
      DATA A( 7,4) / .00166 29871 7324D0/
      DATA A( 8,4) / .00053 05827 9969D0/
      DATA A( 9,4) / .00016 88702 9251D0/
      DATA A(10,4) / .00005 36832 8059D0/
      DATA A(11,4) / .00001 70592 3313D0/
      DATA A(12,4) / .00000 54217 4374D0/
      DATA A(13,4) / .00000 17239 4082D0/
      DATA A(14,4) / .00000 05485 3275D0/
      DATA A(15,4) / .00000 01746 7795D0/
      DATA A(16,4) / .00000 00556 7550D0/
      DATA A(17,4) / .00000 00177 6234D0/
      DATA A(18,4) / .00000 00056 7224D0/
      DATA A(19,4) / .00000 00018 1313D0/
      DATA A(20,4) / .00000 00005 8012D0/
      DATA A(21,4) / .00000 00001 8579D0/
      DATA A(22,4) / .00000 00000 5955D0/
      DATA A(23,4) / .00000 00000 1911D0/
      DATA A(24,4) / .00000 00000 0614D0/
      DATA A(25,4) / .00000 00000 0197D0/
      DATA A(26,4) / .00000 00000 0063D0/
      DATA A(27,4) / .00000 00000 0020D0/
      DATA A(28,4) / .00000 00000 0007D0/
      DATA A(29,4) / .00000 00000 0002D0/
      DATA A(30,4) / .00000 00000 0001D0/
 
      DATA A( 0,5) / .97920 86066 9175D0/
      DATA A( 1,5) / .08518 81314 8683D0/
      DATA A( 2,5) / .00855 98522 2013D0/
      DATA A( 3,5) / .00121 17721 4413D0/
      DATA A( 4,5) / .00020 72276 8531D0/
      DATA A( 5,5) / .00003 99695 8691D0/
      DATA A( 6,5) / .00000 83806 4065D0/
      DATA A( 7,5) / .00000 18684 8945D0/
      DATA A( 8,5) / .00000 04366 6087D0/
      DATA A( 9,5) / .00000 01059 1733D0/
      DATA A(10,5) / .00000 00264 7892D0/
      DATA A(11,5) / .00000 00067 8700D0/
      DATA A(12,5) / .00000 00017 7654D0/
      DATA A(13,5) / .00000 00004 7342D0/
      DATA A(14,5) / .00000 00001 2812D0/
      DATA A(15,5) / .00000 00000 3514D0/
      DATA A(16,5) / .00000 00000 0975D0/
      DATA A(17,5) / .00000 00000 0274D0/
      DATA A(18,5) / .00000 00000 0077D0/
      DATA A(19,5) / .00000 00000 0022D0/
      DATA A(20,5) / .00000 00000 0006D0/
      DATA A(21,5) / .00000 00000 0002D0/
      DATA A(22,5) / .00000 00000 0001D0/
 
      DATA A( 0,6) / .95021 85196 3952D0/
      DATA A( 1,6) / .29052 52916 1433D0/
      DATA A( 2,6) / .05081 77406 1716D0/
      DATA A( 3,6) / .00995 54376 7280D0/
      DATA A( 4,6) / .00211 73389 5031D0/
      DATA A( 5,6) / .00047 85947 0550D0/
      DATA A( 6,6) / .00011 33432 1308D0/
      DATA A( 7,6) / .00002 78473 3104D0/
      DATA A( 8,6) / .00000 70478 8108D0/
      DATA A( 9,6) / .00000 18278 8740D0/
      DATA A(10,6) / .00000 04838 7492D0/
      DATA A(11,6) / .00000 01303 3842D0/
      DATA A(12,6) / .00000 00356 3769D0/
      DATA A(13,6) / .00000 00098 7174D0/
      DATA A(14,6) / .00000 00027 6586D0/
      DATA A(15,6) / .00000 00007 8279D0/
      DATA A(16,6) / .00000 00002 2354D0/
      DATA A(17,6) / .00000 00000 6435D0/
      DATA A(18,6) / .00000 00000 1866D0/
      DATA A(19,6) / .00000 00000 0545D0/
      DATA A(20,6) / .00000 00000 0160D0/
      DATA A(21,6) / .00000 00000 0047D0/
      DATA A(22,6) / .00000 00000 0014D0/
      DATA A(23,6) / .00000 00000 0004D0/
      DATA A(24,6) / .00000 00000 0001D0/
 
      DATA A( 0,7) / .95064 03218 6777D0/
      DATA A( 1,7) / .54138 28546 5171D0/
      DATA A( 2,7) / .13649 97959 0321D0/
      DATA A( 3,7) / .03417 94232 8207D0/
      DATA A( 4,7) / .00869 02788 3583D0/
      DATA A( 5,7) / .00225 28408 4155D0/
      DATA A( 6,7) / .00059 51608 9806D0/
      DATA A( 7,7) / .00015 99561 7766D0/
      DATA A( 8,7) / .00004 36521 3096D0/
      DATA A( 9,7) / .00001 20747 4688D0/
      DATA A(10,7) / .00000 33801 8176D0/
      DATA A(11,7) / .00000 09563 2476D0/
      DATA A(12,7) / .00000 02731 3129D0/
      DATA A(13,7) / .00000 00786 6968D0/
      DATA A(14,7) / .00000 00228 3195D0/
      DATA A(15,7) / .00000 00066 7205D0/
      DATA A(16,7) / .00000 00019 6191D0/
      DATA A(17,7) / .00000 00005 8018D0/
      DATA A(18,7) / .00000 00001 7246D0/
      DATA A(19,7) / .00000 00000 5151D0/
      DATA A(20,7) / .00000 00000 1545D0/
      DATA A(21,7) / .00000 00000 0465D0/
      DATA A(22,7) / .00000 00000 0141D0/
      DATA A(23,7) / .00000 00000 0043D0/
      DATA A(24,7) / .00000 00000 0013D0/
      DATA A(25,7) / .00000 00000 0004D0/
      DATA A(26,7) / .00000 00000 0001D0/
 
      DATA A( 0,8) / .98800 01167 2229D0/
      DATA A( 1,8) / .04364 06760 9601D0/
      DATA A( 2,8) / .00295 09117 8278D0/
      DATA A( 3,8) / .00031 47780 9720D0/
      DATA A( 4,8) / .00004 31484 6029D0/
      DATA A( 5,8) / .00000 69381 8230D0/
      DATA A( 6,8) / .00000 12464 0350D0/
      DATA A( 7,8) / .00000 02429 3628D0/
      DATA A( 8,8) / .00000 00504 0827D0/
      DATA A( 9,8) / .00000 00109 9075D0/
      DATA A(10,8) / .00000 00024 9467D0/
      DATA A(11,8) / .00000 00005 8540D0/
      DATA A(12,8) / .00000 00001 4127D0/
      DATA A(13,8) / .00000 00000 3492D0/
      DATA A(14,8) / .00000 00000 0881D0/
      DATA A(15,8) / .00000 00000 0226D0/
      DATA A(16,8) / .00000 00000 0059D0/
      DATA A(17,8) / .00000 00000 0016D0/
      DATA A(18,8) / .00000 00000 0004D0/
      DATA A(19,8) / .00000 00000 0001D0/
 
      DATA A( 0,9) / .95768 50654 6350D0/
      DATA A( 1,9) / .19725 24967 9534D0/
      DATA A( 2,9) / .02603 37031 3918D0/
      DATA A( 3,9) / .00409 38216 8261D0/
      DATA A( 4,9) / .00072 68170 7110D0/
      DATA A( 5,9) / .00014 09187 9261D0/
      DATA A( 6,9) / .00002 92045 8914D0/
      DATA A( 7,9) / .00000 63763 1144D0/
      DATA A( 8,9) / .00000 14516 7850D0/
      DATA A( 9,9) / .00000 03420 5281D0/
      DATA A(10,9) / .00000 00829 4302D0/
      DATA A(11,9) / .00000 00206 0784D0/
      DATA A(12,9) / .00000 00052 2823D0/
      DATA A(13,9) / .00000 00013 5066D0/
      DATA A(14,9) / .00000 00003 5451D0/
      DATA A(15,9) / .00000 00000 9436D0/
      DATA A(16,9) / .00000 00000 2543D0/
      DATA A(17,9) / .00000 00000 0693D0/
      DATA A(18,9) / .00000 00000 0191D0/
      DATA A(19,9) / .00000 00000 0053D0/
      DATA A(20,9) / .00000 00000 0015D0/
      DATA A(21,9) / .00000 00000 0004D0/
      DATA A(22,9) / .00000 00000 0001D0/
 
      DATA A( 0,10) / .99343 65167 1347D0/
      DATA A( 1,10) / .02225 77012 6826D0/
      DATA A( 2,10) / .00101 47557 4703D0/
      DATA A( 3,10) / .00008 17515 6250D0/
      DATA A( 4,10) / .00000 89997 3547D0/
      DATA A( 5,10) / .00000 12082 3987D0/
      DATA A( 6,10) / .00000 01861 6913D0/
      DATA A( 7,10) / .00000 00317 4723D0/
      DATA A( 8,10) / .00000 00058 5215D0/
      DATA A( 9,10) / .00000 00011 4739D0/
      DATA A(10,10) / .00000 00002 3652D0/
      DATA A(11,10) / .00000 00000 5082D0/
      DATA A(12,10) / .00000 00000 1131D0/
      DATA A(13,10) / .00000 00000 0259D0/
      DATA A(14,10) / .00000 00000 0061D0/
      DATA A(15,10) / .00000 00000 0015D0/
      DATA A(16,10) / .00000 00000 0004D0/
      DATA A(17,10) / .00000 00000 0001D0/
 
      IF(N .LT. 1 .OR. N .GT. 4 .OR. P .LT. 1 .OR. P .GT. 4 .OR.
     1   N+P .GT. 5) THEN
       WGPLG=ZERO
       PRINT 1000, N,P
       RETURN
      END IF
      IF(X .EQ. SGN(0)) THEN
       WGPLG=S1(N,P)
       RETURN
      END IF
 
      IF(X .GT. FCT(2) .OR. X .LT. SGN(1)) THEN
       X1=SGN(0)/X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       V(1)=LOG(DCMPLX(-X,ZERO))
       DO 33 L = 2,N+P
   33  V(L)=V(1)*V(L-1)/L
       SK=ZERO
       DO 34 K = 0,P-1
       P1=P-K
       R=X1**P1/(FCT(P1)*FCT(N-1))
       SM=ZERO
       DO 35 M = 0,K
       N1=N+K-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 31 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   31  B1=B0
       Q=(FCT(N1-1)/FCT(K-M))*(B0-H*B2)*R/P1**N1
   35  SM=SM+V(M)*Q
   34  SK=SK+SGN(K)*SM
       SM=ZERO
       DO 36 M = 0,N-1
   36  SM=SM+V(M)*C(N-M,P)
       WGPLG=SGN(N)*SK+SGN(P)*(SM+V(N+P))
       RETURN
      END IF
 
      IF(X .GT. HALF) THEN
       X1=SGN(0)-X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       U(0)=SGN(0)
       V(1)=LOG(DCMPLX(X1,ZERO))
       U(1)=LOG(X)
       DO 23 L = 2,P
   23  V(L)=V(1)*V(L-1)/L
       DO 26 L = 2,N
   26  U(L)=U(1)*U(L-1)/L
       SK=ZERO
       DO 24 K = 0,N-1
       P1=N-K
       R=X1**P1/FCT(P1)
       SM=ZERO
       DO 25 M = 0,P-1
       N1=P-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 12 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   12  B1=B0
       Q=SGN(M)*(B0-H*B2)*R/P1**N1
   25  SM=SM+V(M)*Q
   24  SK=SK+U(K)*(S1(P1,P)-SM)
       WGPLG=SK+SGN(P)*U(N)*V(P)
       RETURN
      END IF
 
      L=INDEX(10*N+P-10)
      H=C1*X+C2
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 11 I = NC(L),0,-1
      B0=A(I,L)+ALFA*B1-B2
      B2=B1
   11 B1=B0
      WGPLG=(B0-H*B2)*X**P/(FCT(P)*P**N)
      RETURN
 1000 FORMAT(/' ***** CERN SUBROUTINE WGPLG ... ILLEGAL VALUES',
     1        '   N = ',I3,'   P = ',I3)
      END
C-----------------------------------------------------------------------
      DOUBLEPRECISION FUNCTION COMBPDFVEC(NF,Q2,F,IMINUS)
C---COMBINE THE QUARK PDFS FOR THE F2 LIKE STRUCTURE FUNCTIONS.  
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF,IMINUS
      DOUBLE PRECISION P(4,7),F(-6:6),Q2
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CKM2(3,3),CVZE,CAZE,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CHARGE
      COMMON/FLAGEW/LEPCHARGE,EWFLAG
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,CKM2,EQ

      IF (IMINUS.EQ.1) THEN
      CHARGE = - LEPCHARGE
      ELSE
      CHARGE = LEPCHARGE
      ENDIF          

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=(-Q2)**2/(((-Q2)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*((-Q2)-MZ**2)/(-Q2)
      IF (EWFLAG.EQ.3) THEN
        PROPZ=1D0
        INTGZ=-1d0
      ENDIF

      COMBPDFVEC = 0D0
      IF (EWFLAG.LT.3) THEN
        DO I=-NF,NF
C         Photon Part: 
          COMBPDFVEC=COMBPDFVEC+(EQ(I)**2)*F(I)
          IF (EWFLAG.GE.1) THEN
C           Z Part: VVVV-like coupligs:        
            COMBPDFVEC=COMBPDFVEC+F(I)*PROPZ*(CVZ(I)**2+CAZ(I)**2)
     $                *(CVZE**2+CAZE**2)

C           Photon/Z interference terms (null for W)
            COMBPDFVEC=COMBPDFVEC+2*EQ(I)*F(I)*PROPZ*INTGZ*(CVZ(I)*CVZE)
          ENDIF
        ENDDO

      ELSEIF (CHARGE.LT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: VVVV-like coupligs: 
            COMBPDFVEC=COMBPDFVEC+F(2*I)*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(2*I)**2+CAZ(2*I)**2)*CKM2(I,J)
            COMBPDFVEC=COMBPDFVEC+F(-2*I+1)*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(-2*I+1)**2+CAZ(-2*I+1)**2)*CKM2(J,I)
          ENDDO
        ENDDO

      ELSEIF (CHARGE.GT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: VVVV-like coupligs: 
            COMBPDFVEC=COMBPDFVEC+F(2*I-1)*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(2*I-1)**2+CAZ(2*I-1)**2)*CKM2(J,I)
            COMBPDFVEC=COMBPDFVEC+F(-2*I)*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(-2*I)**2+CAZ(-2*I)**2)*CKM2(I,J)
          ENDDO
        ENDDO
      ENDIF
      END
C-----------------------------------------------------------------------
      DOUBLEPRECISION FUNCTION COMBPDFAX(NF,Q2,F,IMINUS)
C---COMBINE THE QUARK PDFS FOR THE F3 LIKE STRUCTURE FUNCTIONS. 
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF,IMINUS
      DOUBLE PRECISION P(4,7),F(-6:6),Q2
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CKM2(3,3),CVZE,CAZE,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CHARGE
      COMMON/FLAGEW/LEPCHARGE,EWFLAG
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,CKM2,EQ

      IF (IMINUS.EQ.1) THEN
      CHARGE = -LEPCHARGE
      ELSE
      CHARGE = LEPCHARGE
      ENDIF

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=(-Q2)**2/(((-Q2)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*((-Q2)-MZ**2)/(-Q2)
      IF (EWFLAG.EQ.3) THEN
        PROPZ=1D0
        INTGZ=-1d0
      ENDIF

      COMBPDFAX = 0D0
      IF (EWFLAG.LT.3) THEN
        DO I=-NF,NF
          IF (EWFLAG.GE.1) THEN
C           Z Part: AVAV-like coupligs:
            COMBPDFAX=COMBPDFAX+F(I)*PROPZ*(4*CVZ(I)*CAZ(I)*CVZE*CAZE)

C           Photon/Z interference terms (null for W)
            COMBPDFAX=COMBPDFAX+2*EQ(I)*F(I)*PROPZ*INTGZ*(CAZ(I)*CAZE)
          ENDIF
        ENDDO
      ELSEIF (CHARGE.LT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: AVAV-like coupligs:
            COMBPDFAX=COMBPDFAX+F(2*I)*PROPZ*(4*CVZE*CAZE)
     $                    *(CVZ(2*I)*CAZ(2*I)*CKM2(I,J))
            COMBPDFAX=COMBPDFAX+F(-2*I+1)*PROPZ*(4*CVZE*CAZE)
     $                     *(CVZ(-2*I+1)*CAZ(-2*I+1)*CKM2(J,I))
          ENDDO
        ENDDO
      ELSEIF (CHARGE.GT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: AVAV-like coupligs:
            COMBPDFAX=COMBPDFAX+F(2*I-1)*PROPZ*(4*CVZE*CAZE)
     $                    *(CVZ(2*I-1)*CAZ(2*I-1)*CKM2(J,I))
            COMBPDFAX=COMBPDFAX+F(-2*I)*PROPZ*(4*CVZE*CAZE)
     $                     *(CVZ(-2*I)*CAZ(-2*I)*CKM2(I,J))
          ENDDO
        ENDDO
      ENDIF
      END
C-----------------------------------------------------------------------
      SUBROUTINE C2NS(Z,VOGT,C2QNS_SING2,C2QNS_DELTA2,
     $                C2QNS_REG2P,C2QNS_REG2M)
      IMPLICIT NONE
C---CALCULATE THE C2NS COEFFICIENTS
      DOUBLE PRECISION Z,C2QNS_SING2,C2QNS_DELTA2,
     $                 C2QNS_REG2P,C2QNS_REG2M
      DOUBLE PRECISION C2QNS_SING,C2QNS_DELTA,C2QNS_REG,C2QNSPLUS_SING,
     $    C2QNSPLUS_DELTA,C2QNSPLUS_REG,C2QNSMIN_REG
      INTEGER VOGT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3
      
      INTEGER LPOL,IORDER,NF,IFL,IF3,IGL
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL

      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

C SCALE DEPENDENT PART      
      C2QNS_SING=( CF**2*(16*L1+12*L0)+CA*CF*(-22D0/3*L0)
     $        +NF*CF*4D0/3*L0)*LOGMU**2  
     $    + (  CF**2*(24*L2-12*L1-(32*ZETA2+45)*L0   )       
     $    + CA*CF*(    -44D0/3*L1+(367D0/9-8*ZETA2)*L0 )
     $    +NF*CF*(8D0/3*L1-58D0/9*L0)             
     $      ) *LOGMU    

      C2QNS_REG=( CF**2*(-8*(1+Z)*LOG1MZ-4*(1+Z**2)/(1-Z)*LOGZ 
     $         +2*(1+Z)*LOGZ-2*(5+Z))
     $        +CA*CF*(11D0/3*(1+Z) )
     $        +NF*CF*(-2D0/3*(1+Z) ))*LOGMU**2  
     $ +( CF**2*((1+Z**2)/(1-Z)*(4*LOGZ**2-24*LOGZ*LOG1MZ-6*LOGZ)
     $ +(1+Z**2)/(1+Z)*(4*LOGZ**2-16*DIMZ-16*LOGZ*LOG1PZ-8*ZETA2) 
     $ +(1+Z)*(4*DI1MZ+4*LOGZ*LOG1MZ-12*LOG1MZ**2-4*LOGZ**2+16*ZETA2)
     $ +8*(2+3*Z)*LOG1MZ-2*(3+11*Z)*LOGZ+2*(19+14*Z))     
     $ +CA*CF*((1+Z**2)/(1-Z)*(44D0/3*LOGZ+2*LOGZ**2)    
     $      +(1+Z**2)/(1+Z)*(8*DIMZ-2*LOGZ**2+8*LOGZ*LOG1PZ+4*ZETA2)
     $    +(1+Z)*(22D0/3*LOG1MZ+4*ZETA2)-1D0/9*(164+434*Z))
     $  +NF*CF*((-8D0/3)*(1+Z**2)/(1-Z)*LOGZ-4D0/3*(1+Z)*LOG1MZ
     $  +1D0/9*(32+68*Z) ) )*LOGMU

      C2QNS_DELTA=( CF**2*(-8*ZETA2+9D0/2)+CA*CF*(-11D0/2)
     $        +NF*CF*1)*LOGMU**2  +
     $     (CF**2* (40*ZETA3-12*ZETA2-51D0/2)
     $       +CA*CF*(-12*ZETA3+88D0/3*ZETA2+215D0/6)
     $        +NF*CF*(-16D0/3*ZETA2-19D0/3 )          
     $    )*LOGMU
     $   + ( CF**2*(16*L1X+12*L0X)+CA*CF*(-22D0/3*L0X)
     $        +NF*CF*4D0/3*L0X)*LOGMU**2  
     $    + (  CF**2*(24*L2X-12*L1X-(32*ZETA2+45)*L0X   )       
     $    + CA*CF*(    -44D0/3*L1X+(367D0/9-8*ZETA2)*L0X )
     $    +NF*CF*(8D0/3*L1X-58D0/9*L0X)             
     $      ) *LOGMU     


       IF (VOGT.eq.1) then

C NON-SCALE PART FROM VOGT
       C2QNSPLUS_SING=
     1          + 14.2222 * L3 - 61.3333 * L2 - 31.105 * L1 
     2          + 188.64 *L0
     3        + NF * ( 1.77778 * L2 - 8.5926 * L1 + 6.3489 *L0) 
       C2QNSPLUS_DELTA= - 338.531 + 0.485 + NF * (46.844 - 0.0035)
     1          + 14.2222 * L3X - 61.3333 * L2X - 31.105 * L1X 
     2          + 188.64 *L0X
     3        + NF * ( 1.77778 * L2X - 8.5926 * L1X + 6.3489 *L0X)           
       C2QNSPLUS_REG =   ! ACTUALLY PLUS+MINUS (VN TERMINOLOGY)
     1          - 69.59 - 1008.* Z
     2          - 2.835 * LOGZ**3 - 17.08 * LOGZ**2 + 5.986 * LOGZ 
     3      - 17.19 * LOG1MZ**3 + 71.08 * LOG1MZ**2 - 660.7 * LOG1MZ
     4         - 174.8 * LOGZ * LOG1MZ**2 + 95.09 * LOGZ**2 * LOG1MZ
     5        + NF * ( - 5.691 - 37.91 * Z 
     6          + 2.244 * LOGZ**2 + 5.770 * LOGZ 
     7          - 1.707 * LOG1MZ**2  + 22.95 * LOG1MZ
     8          + 3.036 * LOGZ**2 * LOG1MZ + 17.97 * LOGZ * LOG1MZ )        
c            
       C2QNSMIN_REG=0 ! MINUS ALREADY ADDED TO C2QNSPLUS_REG

       C2QNS_SING2 = C2QNSPLUS_SING + C2QNS_SING
       C2QNS_DELTA2 = C2QNSPLUS_DELTA + C2QNS_DELTA
       C2QNS_REG2P = (C2QNSPLUS_REG + C2QNS_REG) 
       C2QNS_REG2M = C2QNSMIN_REG

       ELSEIF (vogt.eq.0) then
C NON-SCALE PART
       C2QNSPLUS_SING=(CF*NF*(-247 + 174*LOG1MZ - 36*LOG1MZ**2 
     &        + 72*ZETA2))/(27d0*(-1 + Z)) + (CA*CF*
     &     (3155 + 396*LOG1MZ**2 - 792*ZETA2 + 
     &       6*LOG1MZ*(-367 + 72*ZETA2) - 2160*ZETA3))/(54.*(-1 + Z)) + 
     &  (CF**2*(51 - 54*LOG1MZ - 36*LOG1MZ**2 + 16*LOG1MZ**3 + 
     &       72*ZETA2 - 64*LOG1MZ*ZETA2 - 16*ZETA3))/(2 - 2*Z)
       C2QNSPLUS_DELTA= CF**2*(331d0/8+69*ZETA2+6*ZETA2**2-78*ZETA3+
     &        (51*L0X)/2. - 27*L1X - 18*L2X + 8*L3X + 36*L0X*ZETA2 - 
     &                           32*L1X*ZETA2 - 8*L0X*ZETA3) + 
     &  CA*CF*(-5465/72. - 251*ZETA2/3. + 71*ZETA2**2/5. + 140*ZETA3/3.
     &    +(L0X*(-3155 + 792*ZETA2 + 2160*ZETA3 + 
     &      (1101 - 216*ZETA2)*L0X - 132*L0X**2))/54.)
     &      + CF*NF*((457 + 456*ZETA2 + 48*ZETA3)/36d0 + 
     &           (L0X*(247 - 87*L0X + 24*L1X - 72*ZETA2))/27d0)

       C2QNSPLUS_REG = CF*NF*(-23d0/18 + 
     -     (2*((-4*DI1MZ)/3. + (19*LOGZ)/3. - (8*LOG1MZ*LOGZ)/3. + 
     -          (5*LOGZ**2)/3.))/(1 - Z) - (27*Z)/2. + 
     -     (LOG1MZ*(1 + 13*Z))/3. - (LOGZ*(7 + 19*Z))/3. - 
     -     (1 + Z)*(247d0/54 - (4*DI1MZ)/3. + 
     -        (2*LOG1MZ**2)/3. + (19*LOGZ)/3. + (5*LOGZ**2)/3. - 
     -        LOG1MZ*(29d0/9 + (8*LOGZ)/3.) - (4*ZETA2)/3.))
     -    + CA*CF*(-3229d0/180 + 
     -     LOG1MZ*(133/6d0 - (371*Z)/6.) - 4/(5.*Z) + 
     -     (2191*Z)/20. - (36*Z**2)/5. + 
     -     4*(DI1MZ + LOG1MZ*LOGZ)*(1 + Z) + 
     -     (LOGZ*(13 + 24/Z + 1753*Z - 216*Z**2))/30. + 
     -     LOGZ**2*(-2 + 2*Z - 12*Z**2 - (18*Z**3)/5.) + 
     -     (DIMZ + LOG1PZ*LOGZ)*
     -      (-20 - 4/(5.*Z**2) - 4*Z + 24*Z**2 + (36*Z**3)/5.) + 
     -     (-2 - 10*Z + 24*Z**2 + (36*Z**3)/5.)*ZETA2 + 
     -     (4 - 20*Z)*(DIMZ*LOGZ + S1MZ - 2*TRIMZ - LOG1MZ*ZETA2) + 
     -     (2*((22*DI1MZ)/3. + 
     -          (-239/6d0 + 4*DI1MZ + 12*DIMZ)*LOGZ - 
     -          (55*LOGZ**2)/6. - LOGZ**3 + 
     -          LOG1MZ*(4*DI1MZ + (44*LOGZ)/3. + 2*LOGZ**2) + 12*S1MZ - 
     -          12*TRI1MZ - 24*TRIMZ - 18*ZETA3))/(1 - Z) - 
     -     (1 + Z)*(-3155/108d0 + (22*DI1MZ)/3. - 
     -        (11*LOG1MZ**2)/3. + 
     -        (-239/6d0 + 4*DI1MZ + 12*DIMZ)*LOGZ - 
     -        (55*LOGZ**2)/6. - LOGZ**3 + 12*S1MZ - 12*TRI1MZ - 
     -        24*TRIMZ + LOG1MZ*
     -         (367/18d0 + 4*DI1MZ + (44*LOGZ)/3. + 
     -           2*LOGZ**2 - 4*ZETA2) + (22*ZETA2)/3. + 2*ZETA3))
      C2QNSPLUS_REG = C2QNSPLUS_REG + 
     -  CF**2*(407/20d0 + 8/(5.*Z) - (1917*Z)/20. + (72*Z**2)/5. + 
     -     LOG1MZ**2*(5 + 9*Z) - DI1MZ*(14 + 30*Z) - 
     -     LOG1MZ*LOGZ*(28 + 44*Z) + (LOG1MZ*(-91 + 141*Z))/2. + 
     -     (LOGZ*(13 - 16/Z - 407*Z + 144*Z**2))/10. + 
     -     (DIMZ + LOG1PZ*LOGZ)*
     -      (40 + 8/(5.*Z**2) + 8*Z - 48*Z**2 - (72*Z**3)/5.) + 
     -     LOGZ**2*(29/2d0 + (25*Z)/2. + 24*Z**2 + (36*Z**3)/5.) + 
     -     (-10 + 6*Z - 48*Z**2 - (72*Z**3)/5.)*ZETA2 + 
     -     (-8 + 40*Z)*(DIMZ*LOGZ + S1MZ - 2*TRIMZ - LOG1MZ*ZETA2) + 
     -     (1 + Z)*(2*LOG1MZ**2*LOGZ + (5*LOGZ**3)/3. + 
     -        4*LOG1MZ*(DI1MZ - LOGZ**2) - 4*TRI1MZ - 
     -        4*LOGZ*(DI1MZ + ZETA2)) - 
     -     (1 + Z)*(51/4d0 - 6*DI1MZ + 4*LOG1MZ**3 - (3*LOGZ**2)/2. - 
     -        (4*LOGZ**3)/3. - LOG1MZ**2*(9 + 14*LOGZ) - 12*S1MZ + 
     -        12*TRI1MZ + 48*TRIMZ + 18*ZETA2 - 
     -        LOG1MZ*(27/2d0 + 4*DI1MZ - 12*LOGZ - 12*LOGZ**2 + 
     -        16*ZETA2) + LOGZ*(61/2d0 - 24*DIMZ + 24*ZETA2) + 32*ZETA3)
     -        + (2*(-6*DI1MZ - 14*LOG1MZ**2*LOGZ - (3*LOGZ**2)/2. - 
     -          (4*LOGZ**3)/3. - 
     -          LOG1MZ*(4*DI1MZ - 12*LOGZ - 12*LOGZ**2) - 12*S1MZ + 
     -          12*TRI1MZ + 48*TRIMZ + 
     -          LOGZ*(61/2d0 - 24*DIMZ + 24*ZETA2) + 36*ZETA3))/(1 - Z))
          
       C2QNSMIN_REG =(CF*(CF-CA/2d0))*((DIMZ + LOG1PZ*LOGZ)*
     -   (32 + 8/(5.*Z**2) + 32*Z + 48*Z**2 - (72*Z**3)/5.) + 
     -  8*(DI1MZ + LOG1MZ*LOGZ)*(1 + Z)+
     -   (16*LOG1MZ*(1 - Z) + (LOGZ*(-26 - 8/Z - 106*Z + 72*Z**2))/5. + 
     -     (-162 + 8/Z + 82*Z + 72*Z**2)/5. + 
     -     LOGZ**2*(-4 - 16*Z - 24*Z**2 + (36*Z**3)/5.) + 
     -     (-4 + 20*Z + 48*Z**2 - (72*Z**3)/5.)*ZETA2) + 
     -  (4 + 20*Z)*(-4*DIMZ*LOG1PZ - 2*LOG1PZ**2*LOGZ + 
     -     LOG1PZ*LOGZ**2 - 4*SMZ + 2*TRIMZ - 2*LOG1PZ*ZETA2 + 2*ZETA3)
     -   + ((1 + Z**2)*(-16*DIMZ*LOG1PZ + 
     -       LOGZ*(-8 + 8*DI1MZ + 16*DIMZ - 8*LOG1PZ**2 + 
     -          20*LOG1PZ*LOGZ - 2*LOGZ**2) + 8*S1MZ - 16*SMZ - 
     -       16*TRI1MZ - 16*TRIMCO + 8*TRIMZ + 16*TRIPCO + 
     -       LOG1MZ*(-16*DIMZ - 16*LOG1PZ*LOGZ + 4*LOGZ**2 - 8*ZETA2) - 
     -       8*LOG1PZ*ZETA2 + 8*ZETA3))/(1 + Z))

       C2QNS_SING2 = C2QNSPLUS_SING + C2QNS_SING
       C2QNS_DELTA2 = C2QNSPLUS_DELTA + C2QNS_DELTA
       C2QNS_REG2P = (C2QNSPLUS_REG + C2QNS_REG)
       C2QNS_REG2M = C2QNSMIN_REG       

       ENDIF

      END
C-----------------------------------------------------------------------
      SUBROUTINE C2PS(Z,VOGT,C2QPS_REG)
      IMPLICIT NONE
C---CALCULATE THE C2PS COEFFICIENTS
      DOUBLE PRECISION Z,C2QPS_REG
      INTEGER VOGT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      INTEGER LPOL,IORDER,NF,IFL,IF3,IGL
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL

      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3


       IF (VOGT.eq.1) then

C SAME STUFF BY VOGT + SCALE FROM VN
        C2QPS_REG= NF * ( 5.290 * (1./Z-1.) + 4.310 * LOGZ**3   
     1   - 2.086 * LOGZ**2 + 39.78 * LOGZ - 0.101 * (1.-Z) *LOG1MZ**3 
     2   - (24.75 - 13.80 * Z) * LOGZ**2 * LOG1MZ+30.23*LOGZ*LOG1MZ )
        C2QPS_REG = C2QPS_REG + (CF*NF)*(
     $   (( 8*(1+Z)*LOGZ+4D0/3*(3-4*Z**2-3*Z+4D0/Z))*LOGMU**2
     $   +( 16*(1+Z)*(DI1MZ+LOGZ*LOG1MZ-LOGZ**2)+32*Z**2*LOGZ
     $   +8D0/3*(3-4*Z**2-3*Z+4D0/Z)*LOG1MZ
     $    -16D0/9*(39+4*Z**2-30*Z-13D0/Z))*LOGMU)/2)

       ELSEIF (vogt.eq.0) then
C SAME STUFF + SCALE FROM VN
       C2QPS_REG = CF*NF*(158/9d0 + 344/(27.*Z) - (422*Z)/9. + 
     -    (448*Z**2)/27. + 16*LOG1MZ*LOGZ*Z**2 + 
     -    LOGZ*(56 - (88*Z)/3. - (128*Z**2)/9.) + 
     -    LOGZ**2*(-1 + 15*Z - (32*Z**2)/3.) + 
     -    (DIMZ + LOG1PZ*LOGZ)*
     -     (-16 - 16/(3.*Z) - 16*Z - (16*Z**2)/3.) + 
     -    LOG1MZ*(-104/3d0 + 104/(9.*Z) + (80*Z)/3. - 
     -       (32*Z**2)/9.) + 
     -    LOG1MZ**2*(2 + 8/(3.*Z) - 2*Z - (8*Z**2)/3.) + 
     -    DI1MZ*(4 + 16/(3.*Z) - 4*Z + (32*Z**2)/3.) + 
     -    (-4 - 32/(3.*Z) - 12*Z + (16*Z**2)/3.)*ZETA2 + 
     -    (1 + Z)*(8*DI1MZ*LOG1MZ - 8*DI1MZ*LOGZ + 4*LOG1MZ**2*LOGZ - 
     -       8*LOG1MZ*LOGZ**2 + (10*LOGZ**3)/3. - 8*TRI1MZ - 
     -       8*LOGZ*ZETA2))
        C2QPS_REG = C2QPS_REG + (CF*NF)*(
     $   (( 8*(1+Z)*LOGZ+4D0/3*(3-4*Z**2-3*Z+4D0/Z))*LOGMU**2
     $   +( 16*(1+Z)*(DI1MZ+LOGZ*LOG1MZ-LOGZ**2)+32*Z**2*LOGZ
     $   +8D0/3*(3-4*Z**2-3*Z+4D0/Z)*LOG1MZ
     $    -16D0/9*(39+4*Z**2-30*Z-13D0/Z))*LOGMU)/2)
       ENDIF

      END
C-----------------------------------------------------------------------
      SUBROUTINE C2G(Z,VOGT,C2G_REG2,C2G_DELTA2)
      IMPLICIT NONE
C---CALCULATE THE C2G COEFFICIENTS
      DOUBLE PRECISION Z,C2G_REG2,C2G_DELTA2
      INTEGER VOGT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      INTEGER LPOL,IORDER,NF,IFL,IF3,IGL
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL

      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3


       IF (VOGT.eq.1) THEN

C GLUON FROM VOGT + SCALE FROM VN
      C2G_REG2 = NF *(1./Z * (11.90 + 1494.* LOG1MZ) + 5.319 * LOGZ**3  
     1        - 59.48 * LOGZ**2 - 284.8 * LOGZ + 392.4 - 1483.* LOG1MZ
     2      + (6.445 + 209.4 * (1.-Z)) * LOG1MZ**3 - 24.00 * LOG1MZ**2
     3         - 724.1 * LOGZ**2 * LOG1MZ - 871.8 * LOGZ * LOG1MZ**2 )
C SCALE TERMS      
      C2G_REG2=C2G_REG2
     $   +NF*( CF*TR*( -4*(1+4*Z**2-2*Z)*LOGZ +8*(1+2*Z**2-2*Z)*LOG1MZ
     $   -2+8*Z) +CA*TR*(8*(1+4*Z)*LOGZ+8*(1+2*Z**2-2*Z)*LOG1MZ
     $         +4D0/3*(3-31*Z**2+24*Z+4D0/Z))) *LOGMU**2        
     $   +NF*( CF*TR*(16*(1+2*Z**2-2*Z)*(LOG1MZ**2-2*ZETA2)+
     $  8*(1+4*Z**2-2*Z)*LOGZ**2-8*(3+8*Z**2-6*Z)*LOGZ*LOG1MZ
     $+8*(1+10*Z**2-6*Z)*LOGZ+8*(1-2*Z)*DI1MZ-4*(7+20*Z**2-24*Z)*LOG1MZ
     $                                            +4*(9+4*z**2-17*z))
     $ +CA*TR*(32*Z*(3-Z)*LOGZ*LOG1MZ-16*(1+2*Z**2+2*Z)*(LOGZ*LOG1PZ
     $ +DIMZ)-16*(1+3*Z)*LOGZ**2+16*(1+4*Z)*DI1MZ-16*(1+2*Z**2)*ZETA2
     $ +8*Z*(25*Z-24)*LOGZ+8*(1+2*Z**2-2*Z)*LOG1MZ**2-8D0/3*(3+67*Z**2
     $  -60*Z-4D0/Z)*LOG1MZ-4D0/9*(165-407*Z**2+276*Z-52D0/Z)))*LOGMU
     
C "DELTA TERM" FOR GLUONS. NOT REAL, JUST ADDED BY VOGT TO IMPROVE ACCURACY 
       C2G_DELTA2 = - NF * 0.28 
         
         ELSEIF (vogt.eq.0) then

C GLUON COMPLETE
      C2G_REG2 =CA*NF*(239d0/9 + 344d0/(27*Z) + (1072d0*Z)/9 - 
     -     (4493d0*Z**2)/27 + 
     -     (8*DI1MZ*LOGZ + 8*LOG1PZ*LOGZ**2 - 4*S1MZ + 16*TRIMZ)*Z**2 + 
     -     (LOGZ**3*(10 + 28*Z))/3 + 
     -     LOGZ*(58 + (584d0*Z)/3. - (2090d0*Z**2)/9.) + 
     -     LOGZ**2*(-1 + 88*Z - (194*Z**2)/3.) + 
     -     LOG1MZ**2*(-2 + 8/(3.*Z) + 36*Z - (122*Z**2)/3.) + 
     -     LOG1MZ**2*LOGZ*(24*Z - 8*Z**2) + 
     -     DI1MZ*LOG1MZ*(4 + 40*Z - 8*Z**2) + 
     -     (2*LOG1MZ**3*(1 - 2*Z + 2*Z**2))/3. + 
     -     8*(-(DIMZ*LOG1MZ) - LOG1MZ*LOG1PZ*LOGZ - TRIMCO + TRIPCO)*
     -      (1 + 2*Z + 2*Z**2) + TRI1MZ*(-4 - 72*Z + 8*Z**2) + 
     -     LOG1MZ*LOGZ**2*(-4 - 32*Z + 8*Z**2) + 
     -     (DIMZ + LOG1PZ*LOGZ)*(-24 - 16/(3.*Z) + (80*Z**2)/3.) + 
     -     LOG1MZ*LOGZ*(8 - 144*Z + 148*Z**2) + 
     -     (DI1MZ*(12 + 16/Z - 192*Z + 176*Z**2))/3. + 
     -     (LOG1MZ*(-186 + 104/Z - 1362*Z + 1570*Z**2))/9. + 
     -     LOG1MZ*(-20 + 24*Z - 32*Z**2)*ZETA2 + 
     -     LOGZ*(-48*Z + 16*Z**2)*ZETA2 + 
     -     ((12 - 32/Z - 240*Z + 268*Z**2)*ZETA2)/3. + 
     -     4*(1 + Z)**2*(4*DIMZ*LOG1PZ - 2*DI1MZ*LOGZ + 2*DIMZ*LOGZ + 
     -        2*LOG1PZ**2*LOGZ + LOG1PZ*LOGZ**2 + S1MZ + 4*SMZ - 
     -        2*TRIMZ + 2*LOG1PZ*ZETA2) - (10 + 12*Z + 12*Z**2)*ZETA3)
      C2G_REG2 = C2G_REG2 +
     -    CF*NF*(-647/15d0 + 8/(15.*Z) + (239*Z)/5. + 
     -     64*TRIMZ*Z - (36*Z**2)/5. + DI1MZ*(-10 + 24*Z) + 
     -     (LOGZ*(-236 - 8/Z + 339*Z - 648*Z**2))/15. + 
     -     LOG1MZ**2*(14*Z - 23*Z**2) + LOG1MZ*(-12*Z + 10*Z**2) + 
     -     LOG1MZ*LOGZ*(-24*Z + 56*Z**2) + 
     -     LOGZ**2*(-3/2d0 + (22*Z)/3. - 36*Z**2 - (48*Z**3)/5.) + 
     -     (DIMZ + LOG1PZ*LOGZ)*
     -      (48 + 8/(15.*Z**2) + (64*Z)/3. + (96*Z**3)/5.) + 
     -     ((-20*Z)/3. + 46*Z**2 + (96*Z**3)/5.)*ZETA2 + 
     -     8*(1 + Z)**2*(-4*DIMZ*LOG1PZ - 2*LOG1PZ**2*LOGZ + 
     -        LOG1PZ*LOGZ**2 - 4*SMZ - 2*LOG1PZ*ZETA2) + 
     -     Z**2*((10*LOG1MZ**3)/3. - 12*LOG1MZ**2*LOGZ - 5*LOGZ**3 + 
     -        12*S1MZ - 8*TRI1MZ + LOG1MZ*(16*LOGZ**2 - 16*ZETA2) + 
     -        LOGZ*(12*DI1MZ + 20*ZETA2)) + (64*Z + 36*Z**2)*ZETA3 + 
     -     4*(1 - Z)**2*((5*LOG1MZ**3)/6. - (5*LOGZ**3)/12. - 
     -        LOG1MZ**2*(13/4d0 + 2*LOGZ) + 
     -        LOG1MZ*(7/2d0 + 2*DI1MZ + 4*LOGZ + 2*LOGZ**2) - S1MZ - 
     -        4*TRI1MZ + 12*TRIMZ + (13*ZETA2)/2. + 
     -        LOGZ*(DI1MZ - 4*DIMZ + 3*ZETA2) + 13*ZETA3))
C SCALE TERMS      
      C2G_REG2=C2G_REG2
     $   +NF*( CF*TR*( -4*(1+4*Z**2-2*Z)*LOGZ +8*(1+2*Z**2-2*Z)*LOG1MZ
     $   -2+8*Z) +CA*TR*(8*(1+4*Z)*LOGZ+8*(1+2*Z**2-2*Z)*LOG1MZ
     $         +4D0/3*(3-31*Z**2+24*Z+4D0/Z))) *LOGMU**2        
     $   +NF*( CF*TR*(16*(1+2*Z**2-2*Z)*(LOG1MZ**2-2*ZETA2)+
     $  8*(1+4*Z**2-2*Z)*LOGZ**2-8*(3+8*Z**2-6*Z)*LOGZ*LOG1MZ
     $+8*(1+10*Z**2-6*Z)*LOGZ+8*(1-2*Z)*DI1MZ-4*(7+20*Z**2-24*Z)*LOG1MZ
     $                                            +4*(9+4*z**2-17*z))
     $ +CA*TR*(32*Z*(3-Z)*LOGZ*LOG1MZ-16*(1+2*Z**2+2*Z)*(LOGZ*LOG1PZ
     $ +DIMZ)-16*(1+3*Z)*LOGZ**2+16*(1+4*Z)*DI1MZ-16*(1+2*Z**2)*ZETA2
     $ +8*Z*(25*Z-24)*LOGZ+8*(1+2*Z**2-2*Z)*LOG1MZ**2-8D0/3*(3+67*Z**2
     $  -60*Z-4D0/Z)*LOG1MZ-4D0/9*(165-407*Z**2+276*Z-52D0/Z)))*LOGMU
     
C "DELTA TERM" FOR GLUONS. NOT REAL, JUST ADDED BY VOGT TO IMPROVE ACCURACY 
       C2G_DELTA2 = 0

       ENDIF

      END
C-----------------------------------------------------------------------
      SUBROUTINE CFL(Z,VOGT,CLQS2,CLQNS2,CLQNSP_DELTA2,CLG2)
      IMPLICIT NONE
C---CALCULATE THE FL COEFFICIENTS
      DOUBLE PRECISION Z,CLQS2,CLQNS2,CLQNSP_DELTA2,CLG2
      INTEGER VOGT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      INTEGER LPOL,IORDER,NF,IFL,IF3,IGL
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL

      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3


      IF (vogt.eq.1) then
C LONGITUDINAL STRUCTURE FUNCTION FROM VOGT + SCALE FROM VN

       CLqs2 = NF*CF*TR*(-32*Z*LOGZ-32D0/3*(3-2*Z**2-1D0/Z))*LOGMU
     &      +  NF * ( (15.94 - 5.212 * Z) * (1.-Z)**2 * LOG1MZ
     1         + (0.421 + 1.520 * Z) * LOGZ**2 + 28.09 * (1.-Z) * LOGZ
     2         - (2.370/Z - 19.27) * (1.-Z)**3 )
C
       CLqns2 = (CF**2*(8*Z*(2*LOG1MZ-LOGZ)+4*(2+Z)) +CA*CF*(-44D0/3*Z)
     $      +NF*CF*8D0/3*Z) *LOGMU
C       plus PART
     1    - 40.41 + 97.48 * Z
     2    + (26.56 * Z - 0.031) * LOGZ**2 - 14.85 * LOGZ 
     3    + 13.62 * LOG1MZ**2 - 55.79 * LOG1MZ - 150.5 * LOGZ * LOG1MZ 
     4    + NF * 16./27.D0 * ( 6.* Z*LOG1MZ - 12.* Z*LOGZ - 25.* Z + 6.)

C "DELTA TERM" FOR NS QUARK. NOT REAL, JUST ADDED BY VOGT TO IMPROVE ACCURACY 
       CLQNSP_delta2 = -0.164

               
       CLg2 =NF*( ( CF*TR*(32*Z*LOGZ+16*(1-2*Z**2+Z)) 
     $   +CA*TR*(64*Z*(1-Z)*LOG1MZ-128*Z*LOGZ 
     $  -32D0/3*(3-17*Z**2+15*Z-1D0/Z)))*LOGMU )
C VOGT      
     $ +  NF * ( (94.74 - 49.20 * Z) * (1.-Z) * LOG1MZ**2 
     1         + 864.8 * (1.-Z) * LOG1MZ + 1161.* Z * LOGZ * LOG1MZ 
     2         + 60.06 * Z * LOGZ**2 + 39.66 * (1.-Z) * LOGZ 
     3         - 5.333 * (1./Z - 1.) )    

       ELSEIF (vogt.eq.0) then

C LONGITUDINAL STRUCTURE FUNCTION + SCALE FROM VN
       CLqs2 = NF*CF*TR*(-32*Z*LOGZ-32D0/3*(3-2*Z**2-1D0/Z))*LOGMU 
     &   + (16/9d0)/Z*CF*NF*(-(1 - Z)**3 - 9*(1 - Z)*Z**2 + 
     &      3*LOG1MZ*(1 - Z)*(1 - 2*Z - 2*Z**2) + 
     &      9*LOGZ*Z*(1 - Z - 2*Z**2) + 9*Z**2*(DIZ + LOGZ**2 - ZETA2))
C
       CLqns2 = (CF**2*(8*Z*(2*LOG1MZ-LOGZ)+4*(2+Z)) +CA*CF*(-44D0/3*Z)
     $      +NF*CF*8D0/3*Z) *LOGMU
     &   + 4*(CA - 2*CF)*CF*Z*
     &     ((-23*LOG1MZ)/3. - 4*DIZ*LOGZ + 4*LOG1PZ**2*LOGZ - 
     &    4*DIMZ*(-2*LOG1PZ + LOGZ) + 8*SMZ + 4*TRIMZ + 4*TRIZ + 
     &    (2*LOGZ**2*(5 - 3*Z**2))/5. + 
     &    (4*LOGZ*(6 - 3*Z + 47*Z**2 - 9*Z**3))/(15.*Z**2) - 
     &    (144 + 294*Z - 1729*Z**2 + 216*Z**3)/(90.*Z**2) - 
     &    (4*(DIMZ + LOG1PZ*LOGZ)*(2 + 10*Z**2 + 5*Z**3 - 3*Z**5))/
     &     (5.*Z**3) - 8*ZETA3 - 2*LOGZ**2*DLOG(1 - Z**2) + 
     &    4*ZETA2*((-5 + 3*Z**2)/5. + DLOG(1 - Z**2)))
     &   + 8*CF**2*Z*(DIZ - (78 - 355*Z)/(36.*Z) + 
     &    (LOG1MZ*(6 - 25*Z))/(6.*Z) - (LOGZ*(3 - 22*Z))/(3.*Z) - 
     &    3*ZETA2 + DLOG(Z/(1 - Z))**2)
     &   - (8*CF*NF*Z*(-(6 - 25*Z)/(6.*Z) + DLOG(Z**2/(1 - Z))))/3.

C "DELTA TERM" FOR NS QUARK. NOT REAL, JUST ADDED BY VOGT TO IMPROVE ACCURACY 
       CLQNSP_delta2 = 0
               
       CLg2 =NF*( ( CF*TR*(32*Z*LOGZ+16*(1-2*Z**2+Z)) 
     $   +CA*TR*(64*Z*(1-Z)*LOG1MZ-128*Z*LOGZ 
     $  -32D0/3*(3-17*Z**2+15*Z-1D0/Z)))*LOGMU )
C NON-SCALE      
     &   +CA*NF*(16/3d0 - 16/(9.*Z) + (272*Z)/3. - 
     &     64*DI1MZ*Z + 48*LOGZ**2*Z - (848*Z**2)/9. + 
     &     LOGZ*(16 + 128*Z - 208*Z**2) + LOG1MZ**2*(16*Z - 16*Z**2) + 
     &     LOG1MZ*LOGZ*(-96*Z + 32*Z**2) + 
     &     (DIMZ + LOG1PZ*LOGZ)*(32*Z + 32*Z**2) + 
     &     LOG1MZ*(-16 + 16/(3.*Z) - 144*Z + (464*Z**2)/3.) + 
     &     32*Z**2*ZETA2) + CF*NF*
     &   (-128/15d0 + 32/(15.*Z) - (304*Z)/5. + 
     &     16*(DI1MZ + LOG1MZ*LOGZ)*Z + (336*Z**2)/5. + 
     &     LOG1MZ*(8 + 24*Z - 32*Z**2) + 
     &     (LOGZ*(-104 - 32/Z - 624*Z + 288*Z**2))/15. - 
     &     LOGZ**2*((32*Z)/3. + (32*Z**3)/5.) + 
     &     (DIMZ + LOG1PZ*LOGZ)*
     &      (32/(15.*Z**2) - (32*Z)/3. + (64*Z**3)/5.) + 
     &     ((-32*Z)/3. + (64*Z**3)/5.)*ZETA2)

       ENDIF

      END
C-----------------------------------------------------------------------
      SUBROUTINE C2POLNS(Z,C2QNS_SING2,C2QNS_DELTA2,
     $                   C2QNS_REG2P,C2QNS_REG2M)
      IMPLICIT NONE
C---CALCULATE THE C2NS COEFFICIENTS
      DOUBLE PRECISION Z,C2QNS_SING2,C2QNS_DELTA2,
     $                 C2QNS_REG2P,C2QNS_REG2M
      DOUBLE PRECISION C2QNS_SING,C2QNS_DELTA,C2QNS_REG,C2QNSPLUS_SING,
     $    C2QNSPLUS_DELTA,C2QNSPLUS_REG,C2QNSMIN_REG
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3
      DOUBLE PRECISION C3NSP2A,C3NSM2A,C3NS2B,C3NSP2C,C3NSM2C,DL,DL1

      INTEGER LPOL,IORDER,NF,IFL,IF3,IGL,vogt
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL

      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      vogt = 0

      if (vogt.eq.0) then
      C2QNSMIN_REG=(1+Z**2)/(1+Z)*(16*TRIPCO-16*TRIMCO
     $         +8*S1MZ-16*TRI1MZ-16*SMZ +8*TRIMZ
     $         -16*LOG1MZ*DIMZ-16*LOG1PZ*DIMZ+8*LOGZ*DI1MZ
     $   +16*LOGZ*DIMZ-16*LOGZ*LOG1PZ*LOG1MZ+20*LOGZ**2*LOG1PZ
     $   +4*LOGZ**2*LOG1MZ-8*LOGZ*LOG1PZ**2-8*ZETA2*LOG1MZ
     $   -8*ZETA2*LOG1PZ-2*LOGZ**3-8*LOGZ+8*ZETA3)
     $   +(1+Z)*(16*SMZ-8*TRIMZ+16*LOG1PZ*DIMZ
     $   +8*ZETA2*LOG1PZ+8*LOGZ*LOG1PZ**2-4*LOGZ**2*LOG1PZ+8*DI1MZ
     $   +8*LOGZ*LOG1MZ-8*ZETA3)+(1-Z)*(16*LOG1MZ+30)
     $   +8*(Z**2+1D0/Z)*(DIMZ+LOGZ*LOG1PZ)-4*(2+Z**2+Z)*LOGZ**2
     $   +4*(1+2*Z**2-Z)*ZETA2+(6+38*Z)*LOGZ 
   
      C2QNSMIN_REG=(CF**2-CA*CF/2)*C2QNSMIN_REG
      
      C2QNS_SING=( CF**2*(16*L1+12*L0)+CA*CF*(-22D0/3*L0)
     $        +NF*CF*4D0/3*L0)*LOGMU**2  
     $    + (  CF**2*(24*L2-12*L1-(32*ZETA2+45)*L0   )       
     $    + CA*CF*(    -44D0/3*L1+(367D0/9-8*ZETA2)*L0 )
     $    +NF*CF*(8D0/3*L1-58D0/9*L0)             
     $      ) *LOGMU    
      
      C2QNS_REG=( CF**2*(-8*(1+Z)*LOG1MZ-4*(1+Z**2)/(1-Z)*LOGZ 
     $         +2*(1+Z)*LOGZ-2*(5+Z))
     $        +CA*CF*(11D0/3*(1+Z) )
     $        +NF*CF*(-2D0/3*(1+Z) ))*LOGMU**2  
     $ +( CF**2*((1+Z**2)/(1-Z)*(4*LOGZ**2-24*LOGZ*LOG1MZ-6*LOGZ)
     $ +(1+Z**2)/(1+Z)*(4*LOGZ**2-16*DIMZ-16*LOGZ*LOG1PZ-8*ZETA2) 
     $ +(1+Z)*(4*DI1MZ+4*LOGZ*LOG1MZ-12*LOG1MZ**2-4*LOGZ**2+16*ZETA2)
     $ +8*(1+2*Z)*LOG1MZ-2*(1+9*Z)*LOGZ+2*(16+11*Z))
     $ +CA*CF*((1+Z**2)/(1-Z)*(44D0/3*LOGZ+2*LOGZ**2)    
     $      +(1+Z**2)/(1+Z)*(8*DIMZ-2*LOGZ**2+8*LOGZ*LOG1PZ+4*ZETA2)
     $    +(1+Z)*(22D0/3*LOG1MZ+4*ZETA2)-1D0/9*(98+368*Z))
     $  +NF*CF*((-8D0/3)*(1+Z**2)/(1-Z)*LOGZ-4D0/3*(1+Z)*LOG1MZ
     $  +1D0/9*(20+56*Z) ) )*LOGMU


      C2QNS_DELTA=( CF**2*(-8*ZETA2+9D0/2)+CA*CF*(-11D0/2)
     $        +NF*CF*1)*LOGMU**2  +
     $     (CF**2* (40*ZETA3-12*ZETA2-51D0/2)
     $       +CA*CF*(-12*ZETA3+88D0/3*ZETA2+215D0/6)
     $        +NF*CF*(-16D0/3*ZETA2-19D0/3 )          
     $    )*LOGMU
     $   + ( CF**2*(16*L1X+12*L0X)+CA*CF*(-22D0/3*L0X)
     $        +NF*CF*4D0/3*L0X)*LOGMU**2  
     $    + (  CF**2*(24*L2X-12*L1X-(32*ZETA2+45)*L0X   )       
     $    + CA*CF*(    -44D0/3*L1X+(367D0/9-8*ZETA2)*L0X )
     $    +NF*CF*(8D0/3*L1X-58D0/9*L0X)             
     $      ) *LOGMU     
    
      C2QNSPLUS_SING= CF**2*(8*L3-18*L2
     $  -(32*ZETA2+27)*L1+(-8*ZETA3+36*ZETA2+51D0/2)*L0)
     $ +CA*CF*(-22D0/3*L2+(-8*ZETA2+367D0/9)*L1 
     $  +(40*ZETA3+44D0/3*ZETA2-3155D0/54)*L0   )
     $   +NF*CF*(4D0/3*L2-58D0/9*L1
     $           +(-8D0/3*ZETA2+247D0/27)*L0)
     
       C2QNSPLUS_DELTA=CF**2*( 6*ZETA2**2-78*ZETA3+69*ZETA2+331D0/8)   
     $ +CA*CF*(71d0/5*ZETA2**2+140D0/3*ZETA3-251D0/3*ZETA2-5465D0/72)     
     $   +NF*CF*(4D0/3*ZETA3+38D0/3*ZETA2+457D0/36)
     $  +    CF**2*(8*L3X-18*L2X
     $  -(32*ZETA2+27)*L1X+(-8*ZETA3+36*ZETA2+51D0/2)*L0X)
     $ +CA*CF*(-22D0/3*L2X+(-8*ZETA2+367D0/9)*L1X 
     $  +(40*ZETA3+44D0/3*ZETA2-3155D0/54)*L0X   )
     $   +NF*CF*(4D0/3*L2X-58D0/9*L1X
     $           +(-8D0/3*ZETA2+247D0/27)*L0X)
    
      C2QNSPLUS_REG=( CF**2*(
     $  (1+Z**2)/(1-Z)*(-12*S1MZ+12*TRI1MZ+48*TRIMZ+36*ZETA3
     $  -6*DI1MZ-24*LOGZ*DIMZ+24*ZETA2*LOGZ-4*LOG1MZ*DI1MZ
     $ +12*LOGZ**2*LOG1MZ-14*LOGZ*LOG1MZ**2-4D0/3*LOGZ**3-3D0/2*LOGZ**2
     $ +12*LOGZ*LOG1MZ
     $ +61D0/2*LOGZ)+(1+Z)*(-4*TRI1MZ+4*LOG1MZ*DI1MZ
     $ -4*LOG1MZ**3-4*LOGZ*DI1MZ-4*ZETA2*LOGZ+2*LOGZ*LOG1MZ**2
     $ -4*LOGZ**2*LOG1MZ+5D0/3*LOGZ**3+4*ZETA3)
     $ +(1-Z)*(8*S1MZ-16*TRIMZ+8*LOGZ*DIMZ)
     $ -8*(1+Z+Z**2+1D0/Z)*(DIMZ+LOGZ*LOG1PZ)-2*(5+13*Z)*DI1MZ
     $-4*(7+2*Z**2+7*Z)*ZETA2+2*(5+7*Z)*LOG1MZ**2+8*(1+3*Z)*ZETA2*LOG1MZ
     $ +2*(5+3*Z)*LOG1MZ+(29D0/2+4*Z**2+41D0/2*Z)*LOGZ**2
     $ -16*(1+2*Z)*LOGZ*LOG1MZ
     $ +3D0/2*(3-Z)*LOGZ-41-10*Z)
     $     + CA*CF*(     (1+Z**2)/(1-Z)*(12*S1MZ
     $   -12*TRI1MZ-24*TRIMZ-18*ZETA3+22D0/3*DI1MZ+12*LOGZ*DIMZ  
     $   +4*LOG1MZ*DI1MZ+4*LOGZ*DI1MZ+2*LOGZ**2*LOG1MZ
     $   +44D0/3*LOGZ*LOG1MZ-LOGZ**3-55D0/6*LOGZ**2-239D0/6*LOGZ)
     $   +(1+Z)*(4*DI1MZ
     $   +11D0/3*LOG1MZ**2+4*LOGZ*LOG1MZ-20*ZETA3)+(1-Z)*(-4*S1MZ
     $   +8*TRIMZ-4*LOGZ*DIMZ)+(8*ZETA2-134D0/9-314D0/9*Z)*LOG1MZ
     $ +4*(1+Z+Z**2+1D0/Z)*(DIMZ+LOGZ*LOG1PZ)+4D0/3*ZETA2*(3*Z**2-4*Z-4)
     $ -2*(2+Z**2+2*Z)*LOGZ**2+1D0/6*(157*Z-47)*LOGZ+5D0/27*(145+364*Z) 
     $   )
     $   +NF*CF*(  (1+Z**2)/(1-Z)*
     $    (-4D0/3*DI1MZ-8D0/3*LOGZ*LOG1MZ+5D0/3*LOGZ**2+19D0/3*LOGZ)
     $   +(1+Z)*(4D0/3*ZETA2-2D0/3*LOG1MZ**2)+4D0/9*(5+14*Z)*LOG1MZ
     $   +1D0/3*(1-11*Z)*LOGZ-2D0/27*(58+151*Z) )
     $   )

      C2QNS_SING2 = C2QNSPLUS_SING + C2QNS_SING
      C2QNS_DELTA2 = (C2QNSPLUS_DELTA + C2QNS_DELTA)
      C2QNS_REG2P = C2QNSPLUS_REG + C2QNS_REG
      C2QNS_REG2M = - C2QNSMIN_REG     

      elseif (vogt.eq.1) then

      DL  = LOGZ
      DL1 = LOG1MZ

c     F3: quark non-singlet plus - regular term (A)
       C3NSP2A = 
     1          - 242.9 - 467.2 * Z
     2          - 3.049 * DL**3 - 30.14 * DL**2 - 79.14 * DL 
     3          - 15.20 * DL1**3 + 94.61 * DL1**2 - 396.1 * DL1
     4          - 92.43 * DL * DL1**2 
     5        + NF * ( - 6.337 - 14.97 * Z 
     6          + 2.207 * DL**2 + 8.683 * DL 
     7          + 0.042 * DL1**3 - 0.808 * DL1**2  + 25.00 * DL1
     8          + 9.684 * DL * DL1 )

c     F3: quark non-singlet minus - regular term (A)
      C3NSM2A = 
     1          - 206.1 - 576.8 * Z
     2          - 3.922 * DL**3 - 33.31 * DL**2 - 67.60 * DL 
     3          - 15.20 * DL1**3 + 94.61 * DL1**2 - 409.6 * DL1
     4          - 147.9 * DL * DL1**2 
     5        + NF * ( - 6.337 - 14.97 * Z 
     6          + 2.207 * DL**2 + 8.683 * DL 
     7          + 0.042 * DL1**3 - 0.808 * DL1**2 + 25.00 * DL1
     8          + 9.684 * DL * DL1 )

c     F3: quark non-singlet - singular term (B) - equal for +-
      C3NS2B = 
     1          + 14.2222 * L3 - 61.3333 * L2 - 31.105 * L1 
     2          + 188.64 * L0
     3        + NF * ( 1.77778 * L2 - 8.5926 * L1 + 6.3489 * L0 ) 
c      C3NS2B = DM * C3NS2B

c     F3: quark non-singlet plus - local term (C)
      C3NSP2C = 
c     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
c     2          + 188.64 * DL1 - 338.531  - 0.152 
c     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
c     4          + 6.3489 * DL1 + 46.844 + 0.013)
     1           - 338.046 + NF * 46.8405
c     additional terms from plus distribution
     5          + 14.2222 * L3X - 61.3333 * L2X - 31.105 * L1X 
     6          + 188.64 * L0X
     7        + NF * ( 1.77778 * L2X - 8.5926 * L1X + 6.3489 * L0X ) 

c     F3: quark non-singlet minus - local term (C)
      C3NSM2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 - 0.104 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 + 0.013)

      C2QNS_SING2 = C3NS2B
      C2QNS_DELTA2 = C3NSP2C !C3NSM2C
      C2QNS_REG2P = (C3NSP2A + C3NSM2A)/2.0D0
      C2QNS_REG2M = -(C3NSP2A - C3NSM2A)/2.0D0

      endif

      END
C-----------------------------------------------------------------------
      SUBROUTINE C2POLPS(Z,C2QPS_REG)
      IMPLICIT NONE
C---CALCULATE THE C2PS COEFFICIENTS
      DOUBLE PRECISION Z,C2QPS_REG
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      INTEGER LPOL,IORDER,NF,IFL,IF3,IGL
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL

      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3


      C2QPS_REG= (20*(1-Z)+8*(1+Z)*LOGZ)*LOGMU**2
     $  +(16*(1+Z)*(DI1MZ+LOGZ*LOG1MZ-LOGZ**2)+8*(1-Z)*(5*LOG1MZ
     $  -11)-32*(2-Z)*LOGZ)*LOGMU
     $  +(1+Z)*(-16*TRI1MZ
     $  +16*LOG1MZ*DI1MZ-16*LOGZ*DI1MZ-16*ZETA2*LOGZ
     $  +8*LOGZ*LOG1MZ**2-16*LOGZ**2*LOG1MZ+20D0/3*LOGZ**3)
     $  +(1-Z)*(20*LOG1MZ**2-88*LOG1MZ+760D0/3)-32*(1+Z**2/3+Z+1D0/3/Z)
     $    *(DIMZ+LOGZ*LOG1PZ)+(50+16D0/3*Z**2-10*Z)*LOGZ**2
     $  -32*(2-Z)*LOGZ*LOG1MZ
     $  +4D0/3*(119-13*Z)*LOGZ-(72+32D0/3*Z**2-40*Z)*ZETA2
     $   -8*(3+Z)*DI1MZ 

      C2QPS_REG=(CF*TR*NF)*C2QPS_REG

      END
C-----------------------------------------------------------------------
      SUBROUTINE C2POLG(Z,C2G_REG2)
      IMPLICIT NONE
C---CALCULATE THE C2G COEFFICIENTS
      DOUBLE PRECISION Z,C2G_REG2
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      INTEGER LPOL,IORDER,NF,IFL,IF3,IGL
      COMMON/FLAGPOL/LPOL,IORDER,NF,IFL,IF3,IGL

      COMMON/COLFAC/CF,CA,TR,PI,PISQ
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3


        C2G_REG2 = CF*TR*NF*(
     $   (4*(1-2*Z)*(-2*LOG1MZ+LOGZ)+6)*LOGMU**2
     $ +(8*(1-2*Z)*(-DI1MZ-2*LOG1MZ**2+3*LOGZ*LOG1MZ-LOGZ**2+4*ZETA2) 
     $  +4*(17-20*Z)*LOG1MZ-4*(12-8*Z)*LOGZ-68+52*Z)*LOGMU 
     $  +(1-2*Z)*(32*TRI1MZ-16*LOG1MZ*DI1MZ-8*LOGZ*DI1MZ
     $   - 24*ZETA2*LOGZ-20D0/3*LOG1MZ**3+16*LOGZ*LOG1MZ**2
     $   -16*LOGZ**2*LOG1MZ
     $ +10D0/3*LOGZ**3)-16*(1+Z**2+2*Z)*(4*SMZ+4*LOG1PZ*DIMZ 
     $  +2*LOGZ*dLOG(1+Z)**2-LOGZ**2*LOG1PZ+2*ZETA2*LOG1PZ)-32*(1+Z**2
     $  -6*Z)*TRIMZ+8*(1+4*Z**2-2*Z)*S1MZ
     $  +16D0/3*(13*Z**2+12*Z+4D0/Z)*(DIMZ+LOGZ*LOG1PZ)
C from here Errata 2007     
     $  +4*(5-12*Z)*DI1MZ+32*(1+Z**2-2*Z)*LOGZ*DIMZ
     $  +1D0/3*(123-104*Z**2-48*Z)*LOGZ**2
     $ -(88-96*Z)*LOGZ*LOG1MZ+6*(9-12*Z)*LOG1MZ**2-32*ZETA2*Z**2*LOG1MZ
     $ -4*(31-4*Z**2-26*Z)*LOG1MZ+1D0/3*(416-48*Z**2-274*Z)*LOGZ
     $  -8*(5-4*Z**2-26*Z)*ZETA3
     $  -4D0/3*(81-52*Z**2-108*Z)*ZETA2+2D0/3*(233-239*Z) 
C from here Errata 1997 add term     
     $  +16*(1+2*Z)*(2*DI1MZ+2*LOGZ*LOG1MZ-LOGZ**2)
     $  +96*(1-Z)*LOG1MZ-(144+64*Z)*LOGZ-304*(1-Z)
     $  )

     $ +CA*TR*NF*((-8*(1-2*Z)*LOG1MZ+16*(1+Z)*LOGZ+48*(1-Z))*LOGMU**2  
     $ +(-16*(1+2*Z)*(DIMZ+LOGZ*LOG1PZ)+32*(1+Z)*DI1MZ 
     $ +48*LOGZ*LOG1MZ-32*Z*ZETA2-8*(1-2*Z)*LOG1MZ**2+16*(7-8*Z)*LOG1MZ
     $ -8*(3+4*Z)*LOGZ**2-24*(5-4*Z)*LOGZ-8*(20-21*Z))*LOGMU
     $ +16*(1+2*Z)*(TRIPCO-TRIMCO -LOG1MZ*DIMZ-LOGZ*DI1MZ
     $ -LOGZ*LOG1MZ*LOG1PZ) 
     $  +16*(1+Z**2+2*Z)*(2*SMZ+LOGZ*LOG1PZ**2+2*LOG1PZ*DIMZ
     $  +ZETA2*LOG1PZ)+8*(1-2*Z**2+2*Z)*S1MZ-8*(9+2*Z)*TRI1MZ
     $  -8*(1-Z**2+2*Z)*(2*TRIMZ-LOGZ**2*LOG1PZ-2*LOGZ*DIMZ)
     $  -16*(2+Z)*LOGZ**2*LOG1MZ+24*(LOGZ*LOG1MZ**2-2*ZETA2*LOGZ
     $  -DI1MZ)
     $  -16D0/3*(6+11*Z**2+12*Z+2D0/Z)*(DIMZ+LOGZ*LOG1PZ)
     $  -4D0/3*(1-2*Z)*LOG1MZ**3+8*(6-7*Z)*LOG1MZ**2
     $  +8*(3+2*Z**2-10*Z)*ZETA2*LOG1MZ+4D0/3*(7+10*Z)*LOGZ**3
     $  +2D0/3*(135+44*Z**2-48*Z)*LOGZ**2-8*(17-16*Z)*LOGZ*LOG1MZ
     $  +8D0/3*(118+3*Z**2-26*Z)*LOGZ-4*(44+2*Z**2-53*Z)*LOG1MZ
     $  -4*(3+4*Z**2+10*Z)*ZETA3-16D0/3*(27+11*Z**2-24*Z)*ZETA2
     $ +8*(5+2*Z)*LOG1MZ*DI1MZ+4D0/3*(355-367*Z) )

      END
C-----------------------------------------------------------------------
