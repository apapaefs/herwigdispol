c     ##################################################################
c     # #                 _____      _     _ _                       # #                    
c     # #                |  __ \    | |   | (_)                      # #
c     # #                | |__) |__ | | __| |_ ___                   # #
c     # #                |  ___/ _ \| |/ _` | / __|                  # #
c     # #                | |  | (_) | | (_| | \__ \                  # #
c     # #                |_|   \___/|_|\__,_|_|___/                  # #
c     # #                                                            # #
c     ##################################################################
c     # #       by  I. Borsa, D.de Florian & I. Pedron (2020)        # #
c     ##################################################################                 
c     # Significant extension of DISENT v0.1 by Mike Seymour, with     #
c     # minor modififications by Gavin Salam. It includes polarization #
c     # of initial estate particles and the Projection-to-Born method  #
c     # to obtain exclusive NNLO results.                              #
c     #                                                                #
c     # NOTE: a bug was found in one of the dipoles of the gluon       #
c     #       channel in DISENT v0.1!!! This accounts for the          # 
c     #       discrepancies observed with respect to DISASTER++ in     #
c     #       G. McCance, DESY-PROC-1999-02, also hep-ph/9912481, and  #
c     #       V. Antonelli, M. Dasgupta and G.P. Salam,                #
c     #       JHEP 0002 (2000) 001, among others.                      #
c     ##################################################################
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
      IPOL=0
C-----------------------------------------------------------------------
C     CHOOSE INTERMEDIATE BOSON    ONLY PHOTON = 0 / PHOTON+Z = 1 
C                                   ONLY Z = 2     / W = 3
C-----------------------------------------------------------------------
      IBOSON=1
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
      call InitPDF(0) !debug replicas
C-----------------------------------------------------------------------        
C     NUMBER OF LIGHT FLAVORS   
      NF=4
C-----------------------------------------------------------------------      
C     DEFINE LEPTON AND HADRON ENERGIES & CENTER OF MASS ENERGY^2 
C     EIC KINEMATICS      
      Ee=18
      Ep=275     
      S=4*Ee*Ep
C-----------------------------------------------------------------------
C     SET ALPHA_QED VALUE
      ALPHAE=1/137d0
C-----------------------------------------------------------------------
C     SET EW CONSTANTS
      THETA=0.492334
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
c         . ialg=1  !antikt with E-recombination
c         . ialg=2  !antikt with ET weighted recombination   
c         . ialg=3  !kt with ET weighted recombination 
      ialg = 2  ! JET ALGORITHM DEFINITION
      RAD = 0.8d0  ! JET RADIUS
C####################################################################### 
C     RUN PARAMETETERS
C-----------------------------------------------------------------------      
C     CHOOSE RUN  NNLO=1 / JUST UP TO NLO=0  
      INNLO=1
C-----------------------------------------------------------------------
C     NUMBER OF EVENTS      
      NEV=  1 00 000
C####################################################################### 
c     INITIALISE SEED
      CALL RANDOM_SEED()
      R(1)=14217136  ! standard seed
      CALL RANGEN(0,R) 
C-----------------------------------------------------------------------
      TITLE = 'inclusivejet.top'
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
C     SET KINEMATICAL CUTS      
      XMIN=0
      XMAX=1
c     LIMITS FOR EIC       
      YMIN=0.2
      YMAX=0.6
      Q2MIN=25
      Q2MAX=2500
      END      
C#######################################################################
C     DEFINE RENORMALIZATION AND FACTORIZATION SCALES
      SUBROUTINE SCALES(SCALEREN,SCALEFACT,Q2)
      COMMON/JETS/POUT,zkt,zeta,zphi,n1,n2,n3
      INTEGER n1,n2,n3
      DOUBLE PRECISION POUT(4,7),zkt(1:3),zeta(1:3),zphi(1:3),Q2
      DOUBLE PRECISION SCALEJET,SCALEREN,SCALEFACT,DOT  

C     POUT(:,1): Incoming Parton four-momentum
C     POUT(:,2-4): Outgoing Jets four-momentum 
C     POUT(:,5): Photon four-momentum
C     POUT(:,6): Incoming lepton four-momentum 
C     POUT(:,7): Outgoing lepton four-momentum

C     ZKT(n=1-3): Jet transverse momenta (n1 being the hardest)
C     ZETA(n=1-3): Jet pseudorapidities (n1 being the hardest)
C     ZPHI(n=1-3): Jet phi (n1 being the hardest)

C     SET FACTORIZATION AND RENORMALIZATION SCALES USING KINEMATICS OF
C     LEPTONS, PHOTON AND JETS
C
C     NOTICE THAT IF YOU ARE GOING TO COMPUTE AN INCLUSIVE QUANTITY,
C     YOU SHOULD CHOOSE AN "INCLUSIVE" SCALE DEFINITION 
C     (NOT THE KT OF EACH JET)
C
C     FOR THE TOTAL CROSS SECTION IT ONLY MAKES SENSE TO USE 
C     SOMETHING PROP. TO Q2 IF YOU HAVE SOMETHING LIKE THE CROSS 
C     SECTION AFTER A JET VETO (WITH KTMAX) OR THE INCLUSIVE  
C     CROSS SECTION FOR THE PRODUCTION OF JETS WITH KT>KTMIN
C
C     YOU CAN CHOSE A COMBINATION OF Q2 AND KTMIN OR KTMAX
C     IF YOU LOOK FOR A DIFFERENTIAL JET DISTRIBUTION THE KT OF EACH 
C     JET (COMBINED WITH Q2) MAKES SENSE AS A SCALE

      Q2=ABS(DOT(Pout,5,5))
      SCALEJET=Q2
c      SCALEJET=(zkt(n1)+zkt(n2))**2
c      SCALEJET=(Q2+((zkt(n1)+zkt(n2))/2)**2)/2 
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
C     CALL GBOOK1(I,TITLE,BINS,LOW,UP,NEW,HIST,LOGX,LOGY)
C           INT CHAR  INT DP  DP LOG LOG  LOG  LOG   variable type
C
C     I   IS HISTOGRAM NUMBER (INTEGER)
C
C     TITLE IS HISTOGRAM TITLE (CHARACTER*10)
C
C     BINS  NUMBER OF BINS  (INTEGER)
C
C     LOW  LOWER VALUE FOR X AXIS (DOUBLE PRECISION)
C
C     UP   UPPER VALUE FOR X AXIS (DOUBLE PRECISION)
C
C     NEW  IS .TRUE.  TO PUT HISTOGRAM ON A NEW PAGE
C          .FALSE. TO PUT IT ON TOP OF THE LAST ONE
C
C     HIST IS .TRUE.  FOR HISTOGRAM (ACCUMULATED X-SECTION PER BIN)
C          .FALSE. FOR GRAPH (DIFFERENTIAL X-SECTION, DIVIDED BY X BIN-SIZE)
C
C     LOGX IS .TRUE.  FOR LOGARITHMIC X-SCALE
C          .FALSE. FOR LINEAR X-SCALE
C
C     LOGY IS .TRUE.  FOR LOGARITHMIC Y-SCALE
C          .FALSE. FOR LINEAR Y-SCALE
C
C     HERE CAN USE 1(0) INSTEAD OF .TRUE.(.FALSE.)
C----------------------------------------------------------------------- 
C number of histograms to show
      NPLOTS=30
C INITIALIZE VECTORS TO COMPUTE INTEGRATED X-SECTION PER HISTOGRAM  
      DO I=1,100
        CSUM(I)=0
        CSQR(I)=0
      ENDDO

C     Q^2 DEPENDENCE
       CALL GBOOK1(1,"Q2 NNLO",10,25d0,1000D0,1,0,1,1)
       CALL GBOOK1(2,"Q2 NLO ",10,25d0,1000D0,0,0,1,1)
       CALL GBOOK1(3,"Q2 LO  ",10,25d0,1000D0,0,0,1,1)
C     X DEPENDENCE      
       CALL GBOOK1(4,"X NNLO",10,1d-3,1d0,1,0,1,1)
       CALL GBOOK1(5,"X NLO ",10,1d-3,1d0,0,0,1,1)
       CALL GBOOK1(6,"X LO  ",10,1d-3,1d0,0,0,1,1)
C INCLUSIVE JET PRODUCTION X-SECTION
       CALL GBOOK1(7,"pt 30 NNLO",30,5d0,35d0,1,0,1,1)
       CALL GBOOK1(8,"pt 30 NLO ",30,5d0,35d0,0,0,1,1)
       CALL GBOOK1(9,"pt 30 LO  ",30,5d0,35d0,0,0,1,1)
       
       CALL GBOOK1(10,"pt 15 NNLO",15,5d0,35d0,1,0,1,1)
       CALL GBOOK1(11,"pt 15 NLO ",15,5d0,35d0,0,0,1,1)
       CALL GBOOK1(12,"pt 15 LO  ",15,5d0,35d0,0,0,1,1)

       CALL GBOOK1(13,"pt 10 NNLO",10,5d0,35d0,1,0,1,1)
       CALL GBOOK1(14,"pt 10 NLO ",10,5d0,35d0,0,0,1,1)
       CALL GBOOK1(15,"pt 10 LO  ",10,5d0,35d0,0,0,1,1)

       CALL GBOOK1(16,"eta 15 NNLO",15,-1.8d0,3.2d0,1,0,0,0)
       CALL GBOOK1(17,"eta 15 NLO ",15,-1.8d0,3.2d0,0,0,0,0)
       CALL GBOOK1(18,"eta 15 LO  ",15,-1.8d0,3.2d0,0,0,0,0)

       CALL GBOOK1(19,"eta 10 NNLO",10,-1.8d0,3.2d0,1,0,0,0)
       CALL GBOOK1(20,"eta 10 NLO ",10,-1.8d0,3.2d0,0,0,0,0)
       CALL GBOOK1(21,"eta 10 LO  ",10,-1.8d0,3.2d0,0,0,0,0)

       CALL GBOOK1(22,"eta 5 NNLO",5,-1.8d0,3.2d0,1,0,0,0)
       CALL GBOOK1(23,"eta 5 NLO ",5,-1.8d0,3.2d0,0,0,0,0)
       CALL GBOOK1(24,"eta 5 LO  ",5,-1.8d0,3.2d0,0,0,0,0)

C X DEPENDENCE of full effective structure function (no cuts)   
       CALL GBOOK1(25,"FX NNLO",10,5d-5,9d-1,1,0,1,1)
       CALL GBOOK1(26,"FX NLO ",10,5d-5,9d-1,0,0,1,1)
       CALL GBOOK1(27,"FX LO  ",10,5d-5,9d-1,0,0,1,1)
       
       CALL GBOOK1(28,"pt y NNLO",10,5d0,35d0,1,0,1,1)
       CALL GBOOK1(29,"pt y NLO ",10,5d0,35d0,0,0,1,1)
       CALL GBOOK1(30,"pt y LO  ",10,5d0,35d0,0,0,1,1)
                   
       CALL GBOOK1(31,"eta y NNLO",10,-1.8d0,3.2d0,1,0,0,0)
       CALL GBOOK1(32,"eta y NLO ",10,-1.8d0,3.2d0,0,0,0,0)
       CALL GBOOK1(33,"eta y LO  ",10,-1.8d0,3.2d0,0,0,0,0)
       
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
C  (THIS IS SO THAT PAIRS OF LARGE POSITIVE AND NEGATIVE WEIGHTS
C   GET CANCELLED BEFORE CALCULATING THE VARIANCE)
C  N=0 means all the events/counterevents have been computed (end of event)
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
C DO NOT CHANGE ANYTHING ABOVE
C-----------------------------------------------------------------------
C--- USER ROUTINE STARTS HERE
C     
C      YOU CAN USE THE MOMENTA OF THE PARTICLES 
C     (AFTER PASSING THE JET ALGORITHM)
C     IN POUT(I,J)  WITH I=1,2,3,4 AS  PX,PY,PZ,E
C     AND I=1 INCOMING PARTON, I=2,3,4 OUTGOING JETS (MIGHT BE MASSIVE!)
C     I=5 BOSON, I=6 INCOMING LEPTON, I=7 OUTGOING LEPTON
C     SOME INFORMATION IS AVAILABLE IN VECTOR P(I,J), SAME CONVENTION 
C     BUT PARTONS INSTEAD OF JETS IN THE FINAL STATE (I=2,3,4)
C     ALSO JET INFORMATION IS AVAILABLE IN TERMS OF THE VARIABLES
C     ZKT(ni),ZETA(ni),ZPHI(ni). with ni=1,2,3 
C     HERE et(n3)<=et(n2)<=et(n1)  HARDEST JET IS N1, SECOND HARDEST JET IS N2
C     THEREFORE  
C            ZKT(N1) IS THE TRANSVERSE MOMENTUM OF THE HARDEST JET
C            ZETA(N1) IS THE RAPIDITY OF THE HARDEST JET 
C            ZPHI(N1) IS THE AZIMUTHAL ANGLE OF THE HARDEST JET
C     AND SO ON..
C-----------------------------------------------------------------------
C
C---DEFINE WEIGHTS ACCORDING TO POWER OF ALPHAS^(0,1,2) (NA)
C
C---W(0) WEIGHT FOR ORDER ALPHAS^0 CONTRIBUTION 
C---W(1) WEIGHT FOR ORDER ALPHAS^0+ALPHAS^1 CONTRIBUTION 
C---W(2) WEIGHT FOR ORDER ALPHAS^0+ALPHAS^1+ALPHAS^2 CONTRIBUTION 
C  
C---FILL HISTOGRAMS
C
C---SIMPLY CALL GFILL1(plotnumber, X variable, WEIGHT at a given order)
C
C---LO AND NLO EFFECTIVE STRUCTURE FUNCTIONS
C-----------------------------------------------------------------------
C     CALCULATE DIS KINEMATICS
      Q2=ABS(DOT(POUT,5,5))
      Y=DOT(POUT,1,5)/DOT(POUT,1,6)
      X=Q2/Y/S
C-----------------------------------------------------------------------
C     CALCULATE PT'S AND PSEUDORAPIDITIES IN THE LAB FRAME IN CASE 
C     BREIT FRAME IS CHOSEN

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
c
          YETA(J-1)=0.5D0*DLOG((POUT2(4,J)+POUT2(3,J))
     $              /(POUT2(4,J)-POUT2(3,J)))
          YKT(J-1)=(POUT(1,J)**2+POUT(2,J)**2)**0.5
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Inclusive Cross Section
        FNORM=2*PI*(1+(1-Y)**2)/Y/Q2**2/137**2*GEV2PB
        IF(IPOL.EQ.1) FNORM=(1-(1-Y)**2)/(1+(1-Y)**2)*FNORM
   
        F2effLO=W(0)/FNORM
        F2effNLO=W(1)/FNORM
        F2effNNLO=W(2)/FNORM
          
C     FILL HISTOGRAMS        
        CALL GFILL1(25,X,F2effNNLO)
        CALL GFILL1(26,X,F2effNLO)
        CALL GFILL1(27,X,F2effLO)

C     COMPUTE INCLUSIVE X-SECTION             
             CS(1)=CS(1)+W(2)           
             CS(2)=CS(2)+W(1)  
             CS(3)=CS(3)+W(0)

        IF(W(1).LT.W(1))THEN
        WRITE(6,*) 'NaN ALERT'
        STOP
        ENDIF  
                 
C     Inclusive jet
        DO I=1,3

           IF (zkt(i).gt.5.and.zkt(i).lt.36 
     $     .and. abs(zeta(i)).lt.3) then 


        CALL GFILL1(1,Q2,W(2))
        CALL GFILL1(2,Q2,W(1))
        CALL GFILL1(3,Q2,W(0))
        
        CALL GFILL1(4,X,W(2))
        CALL GFILL1(5,X,W(1))
        CALL GFILL1(6,X,W(0))

           CALL GFILL1(7,zkt(i),W(2))               
           CALL GFILL1(8,zkt(i),W(1))
           CALL GFILL1(9,zkt(i),W(0))

           CALL GFILL1(10,zkt(i),W(2))               
           CALL GFILL1(11,zkt(i),W(1))
           CALL GFILL1(12,zkt(i),W(0))
           
           CALL GFILL1(13,zkt(i),W(2))               
           CALL GFILL1(14,zkt(i),W(1))
           CALL GFILL1(15,zkt(i),W(0))
           
           CALL GFILL1(16,zeta(i),W(2))               
           CALL GFILL1(17,zeta(i),W(1))
           CALL GFILL1(18,zeta(i),W(0))

           CALL GFILL1(19,zeta(i),W(2))               
           CALL GFILL1(20,zeta(i),W(1))
           CALL GFILL1(21,zeta(i),W(0))

           CALL GFILL1(22,zeta(i),W(2))               
           CALL GFILL1(23,zeta(i),W(1))
           CALL GFILL1(24,zeta(i),W(0))
                   
C COMPUTE INCLUSIVE JET X-SECTION             
             CS(7)=CS(7)+W(2)
             CS(8)=CS(8)+W(1)
             CS(9)=CS(9)+W(0)           
c select y>0.04 (HERA like)              
                IF (Y.gt.0.04d0) THEN                
                    CALL GFILL1(28,zkt(i),W(2))               
                  CALL GFILL1(29,zkt(i),W(1))
                  CALL GFILL1(30,zkt(i),W(0))           
                  CALL GFILL1(31,zeta(i),W(2))               
                  CALL GFILL1(32,zeta(i),W(1))
                  CALL GFILL1(33,zeta(i),W(0))
C COMPUTE INCLUSIVE JET X-SECTION WITH Y>0.04             
                 CS(31)=CS(31)+W(2)
                 CS(32)=CS(32)+W(1)
                 CS(33)=CS(33)+W(0)                   
                ENDIF  
            ENDIF        
       ENDDO

   50 format(8f10.5)      
   51 format(4f11.5)            
 999  END
C#######################################################################
      SUBROUTINE PRINTOUT(NEV)
      IMPLICIT NONE
      INTEGER I,J,NEV,kk
      DOUBLE PRECISION CSUM(100),CSQR(100),SQR,X,EQSQ
      COMMON /DEMCOM/CSUM,CSQR
      COMMON /numberofplots/nplots
      LOGICAL NEW(100),HIST(100),LOGY(100),LOGX(100)
      CHARACTER*10 TITLE(100)         
      COMMON / PLOTSCOMT/ NEW,HIST,LOGY,LOGX,TITLE
      INTEGER NPLOTS
      SQR(X)=SIGN(SQRT(ABS(X)),X)
        
      DO I=1,NPLOTS
c        IF (I.GT.1.AND.CSUM(1).NE.0.AND.CSUM(2).NE.0) THEN
c          CALL GSCALE(    I,1/CSUM(2))
c          CALL GSCALE(100+I,1/CSUM(2)**2)
c        ENDIF
        CALL GOPERA(    I,'*',    I,200+I,1D0,1D0)
        CALL GOPERA(100+I,'-',200+I,100+I,1D0,1D0/NEV)
        CALL GOPERA(100+I,'S',    0,100+I,1D0,0D0)         
        CALL GTOPER(I)
      ENDDO

      DO  I=1,NPLOTS
           WRITE(*,'(1X,(A,F15.5,A,F15.5))')TITLE(I),CSUM(I),
     $     '    +-',SQR( CSQR(I)-CSUM(I)**2/NEV)
      ENDDO

      END
C#######################################################################
