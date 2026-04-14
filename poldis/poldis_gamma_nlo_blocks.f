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
C-----------------------------------------------------------------------
      SUBROUTINE POLDIS(NEV,S,NFL,USER,CUTS)
      IMPLICIT NONE
C---CALCULATE DIS EVENT FEATURES TO NEXT-TO-LEADING ORDER
C   ACCORDING TO THE METHOD OF CATANI AND SEYMOUR NPB485 (1997) 291
C
C   VERSION 0.0  5TH May 2020
C
C - NEV IS THE NUMBER OF EVENTS TO GENERATE
C - S IS THE TOTAL LEPTON-HADRON CENTRE-OF-MASS ENERGY SQUARED
C - NFL IS THE NUMBER OF FLAVOURS
C - USER IS THE NAME OF ROUTINE TO ANALYSE THE EVENTS
C   IT TAKES SIX ARGUMENTS:
C     SUBROUTINE USER(N,NA,ITYPE,P,S,WEIGHT)
C     N = THE NUMBER OF PARTONS (FIRST IS INCOMING, N-1 OUTGOING)
C     NA = THE ORDER IN ALPHA
C     ITYPE = THE TYPE: 0=TREE LEVEL, 1=SUBTRACTION, 2=FINITE VIRTUAL,
C                       3=FINITE COLLINEAR
C     P(4,7) = THEIR MOMENTA (P(1-4,I) IS 4-MOMENTUM OF PARTICLE I,
C              P(1-4,5) IS TOTAL EXCHANGED MOMENTUM (IE FIRST-REST)
C              P(1-4,6/7) IS THE 4-MOM OF THE INCOMING/OUTGOING LEPTON)
C     S = THE TOTAL LEPTON-HADRON CENTRE-OF-MASS ENERGY SQUARED
C     WEIGHT(I) = THE WEIGHT IN THE INTEGRAL INITIATED BY PARTONS OF
C                 TYPE I (0=G,1=D,2=U,3=S,...,6=T,-1=DBAR,...,-6=TBAR)
C   IT IS ALSO CALLED AT THE END OF EACH FULL EVENT WITH N=0
C   A SIMPLE EXAMPLE (CALLED USER) IS GIVEN BELOW
C - CUTS IS THE NAME OF ROUTINE THAT PROVIDES CUTS ON X-Q**2-Y
C   IT TAKES ONE ARGUMENT AND PROVIDES SIX OUTPUTS:
C     SUBROUTINE CUTS(S,XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX)
C   IF XMIN=XMAX, THE RESULTING CROSS-SECTION IS DIFFERENTIAL IN X
C   IF XMAX>XMIN, IT IS INTEGRATED OVER THE GIVEN X RANGE
C   IF XMAX=0, THE FULL ALLOWED X RANGE IS INTEGRATED OVER
C   AND LIKEWISE FOR Q2 AND Y.
C   AT MOST TWO OF THE VARIABLES CAN BE FIXED, BUT ALL CAN HAVE LIMITS
C   A LOWER LIMIT ON Q2 MUST EITHER BE EXPLICITLY GIVEN,
C   OR IMPLIED BY THE OTHER TWO LIMITS
C
C   IT IS IMPORTANT TO THINK ABOUT THE EFFICIENCY OF THE USER ROUTINE,
C   AS IT IS CALLED 14 TIMES PER EVENT
C  (THIS IS 8 3-PARTON CALLS, 4 2-PARTON CALLS AND 2 1-PARTON CALLS)
C
C
      INTEGER NEV,NFL,NPERM3,NPERM4,I,J,INNLO,IBOSON,KNLOBLK,
     $     NNLOBLK
      PARAMETER (NPERM3=1,NPERM4=6)
      PARAMETER (NNLOBLK=8)
      DOUBLE PRECISION S,CA_in, CF_in, TR_in,LALA,
     $     WTWO,WTHR,WFOR,JTHR,JFOR,JTMP,MTWO(-6:6),MTHR(-6:6),
     $     MFOR(-6:6),VTWO(-6:6),VTHR(-6:6),CTHR(-6:6),CFOR(-6:6),
     $     STHR(-6:6,NPERM3),SFOR(-6:6,NPERM4),WEIGHT(-6:6),ZERO(-6:6),
     $     P(4,7),Q(4,7,NPERM4),NRM,NPOW1,NPOW2,CUTOFF_IN,SCALE_in,
     $     NLOBLOCKSUM(NNLOBLK),NLOBLOCKQR(NNLOBLK)
      INTEGER SCHEME,NF,ORDER
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION THETAW,CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,MZ,GAMMAZ
      DOUBLE PRECISION DOT,LEPCHARGE,CKM(3,3)
      DOUBLE PRECISION LEPCH,THETA,ZMASS,WMASS,ZGAMMA,WGAMMA,CKM2(3,3)       
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG
      INTEGER NPOW(2)
      INTEGER ISUB,IPOL,IMERGE,ICH,IFL,EWFLAG
      DOUBLE PRECISION XPOW(2)
      DOUBLE PRECISION SCALEREN,SCALEFACT,Q2,WWW(0:5)
      DOUBLE PRECISION PBORNP(4,7),PBORN(4,7),WBORN(0:5)
      COMMON/PBORNV/PBORNP 
      COMMON  /SAMPLE/ XPOW,NPOW
      COMMON/JETSCALES/SCALEREN,SCALEFACT,Q2
      COMMON/GEV2PB/NRM
      COMMON/FLAGPOL/IPOL
      COMMON/CHANNEL/ICH,IFL
      common/imerge/imerge
      COMMON /NLOBLOCKACC/ NLOBLOCKSUM,NLOBLOCKQR
      EXTERNAL USER,CUTS
      COMMON/INNLO/INNLO
      COMMON/EWSETTINGS/LEPCH,THETA,ZMASS,WMASS,ZGAMMA,WGAMMA,CKM,
     $                  IBOSON          
      DATA ZERO/13*0/
C---PRINT OPENING MESSAGE
      WRITE (6,'(A)')'#################################################'
      WRITE (6,'(A)')'# #         _____      _     _ _ ______       # #'
      WRITE (6,'(A)')'# #        |  __ \    | |   | (_)___  /       # #'
      WRITE (6,'(A)')'# #        | |__) |__ | | __| |_   / /        # #'
      WRITE (6,'(A)')'# #        |  ___/ _ \| |/ _` | | / /         # #'
      WRITE (6,'(A)')'# #        | |  | (_) | | (_| | |/ /__        # #'
      WRITE (6,'(A)')'# #        |_|   \___/|_|\__,_|_/_____|       # #'
      WRITE (6,'(A)')'# #                                           # #'
      WRITE (6,'(A)')'#################################################'
      WRITE (6,'(A)')'## by I. Borsa, D.de Florian & I. Pedron (2020)##'       
      WRITE (6,'(A)')'#################################################'
      WRITE (6,'(/2A)')' This is POLDIS, a program for calculating',
     $     ' jet quantities in'
      WRITE (6,'(1A)')   ' polarized deep inelastic scattering to'
      WRITE (6,'(1A)') 'next-to-next-to-leading order in alpha_s'
      WRITE (6,'(A)')    ' It is based on the original DISENT',
     $  ' written by Mike Seymour'
      WRITE (6,'(2A)')   ' S.Catani & M.H.Seymour,',
     $     ' Nucl. Phys. B485 (1997) 291'
      WRITE (6,'(/A)') ' Extensively modified to include polarization'
      WRITE (6,'(A)')' by I.Borsa, D.de Florian & I.Pedron, March 2020'
      WRITE (6,'(A/)')' PhysRevLett.125.082001 (2020)'
      WRITE (6,'(A,I10)') '  Number of Events=',NEV
      IF (IPOL.eq.0)  WRITE (6,'(A)')  '  Polarization: no'
      IF (IPOL.eq.1)  WRITE (6,'(A)')  '  Polarization: yes'
      IF (ICH.eq.1)  WRITE (6,'(A)')  '  Channel: quarks only'
      IF (ICH.eq.2)  WRITE (6,'(A)')  '  Channel: gluons only'

C---- THE ELECTROWEAK FLAG VALUE DETERMINES THE EXCHANGED BOSON:
C         EWFLAG=0 : PHOTON
C         EWFLAG=1 : PHOTON/Z
C         EWFLAG=2 : Z
C         EWFLAG=3 : "W-" (INITIAL ELECTRON - LEPCHARGE = -1)
C                    "W+" (INITIAL POSITRON - LEPCHARGE = +1)
      EWFLAG = IBOSON
      IF (EWFLAG.eq.0) WRITE (6,'(A)')  '  Exchanged boson: Photon'
      IF (EWFLAG.eq.1) WRITE (6,'(A)')  '  Exchanged boson: Photon/Z'
      IF (EWFLAG.eq.2) WRITE (6,'(A)')  '  Exchanged boson: Z'
      IF (EWFLAG.eq.3) WRITE (6,'(A)')  '  Exchanged boson: W'


C---INITIALIZE COLOUR FACTORS AND OTHER CONSTANTS
      CF=4/3D0
      CA=3 
      TR=1/2D0
      PI=ATAN(1D0)*4
      PISQ=PI**2
      HF=0.5D0
      EQ(0)=0
      EQ(1)=-1D0/3
      EQ(2)=EQ(1)+1
      DO I=1,6
        IF (I.GT.2) EQ(I)=EQ(I-2)
        EQ(-I)=-EQ(I)
      ENDDO
      LEPCHARGE = LEPCH      

C---- CLEAR CHARGES IF THERE IS NO PHOTON
      IF (EWFLAG.GE.2) THEN
        DO I=-6,6
          EQ(I)=0D0
        ENDDO
      ENDIF

C---- DEFINE ELECTROWEAK COUPLINGS
      THETAW = THETA
      IF (EWFLAG.LE.2) THEN
        MZ = ZMASS  
        GAMMAZ = ZGAMMA
      ELSEIF (EWFLAG.EQ.3) THEN
        MZ = WMASS
        GAMMAZ = WGAMMA 
      ENDIF

c      CVZ(0)=0
c      CVZ(1)=(-1.D0/2+2.D0/3d0*SIN(THETAW)**2)/SIN(2*THETAW)
c      CVZ(2)=(1.D0/2-4.D0/3d0*SIN(THETAW)**2)/SIN(2*THETAW) 
c      CAZ(0)=0
c      CAZ(1)=-1.D0/2/SIN(2*THETAW)
c      CAZ(2)=1.D0/2 /SIN(2*THETAW)
c      CVZE=(-1.D0/2+2.D0*SIN(THETAW)**2)/SIN(2*THETAW)
c      CAZE=-1.D0/2/SIN(2*THETAW)
     

      CVZ(0)=0
      CAZ(0)=0
      IF (EWFLAG.LT.3) THEN
        CVZ(1)=-1.D0/2/SIN(2*THETAW)+1.D0/3*TAN(THETAW)
        CVZ(2)=1.D0/2/SIN(2*THETAW)-2.D0/3*TAN(THETAW) 
        CAZ(1)=-1.D0/2/SIN(2*THETAW)
        CAZ(2)=1.D0/2/SIN(2*THETAW) 
        CVZE=(-1.D0/2/SIN(2*THETAW)+TAN(THETAW))*(-LEPCHARGE)
        CAZE=-1.D0/2/SIN(2*THETAW)
      ELSE
        CVZ(1)=1.D0/2/SQRT(2d0)/SIN(THETAW)
        CVZ(2)=1.D0/2/SQRT(2d0)/SIN(THETAW) 
        CAZ(1)=1.D0/2/SQRT(2d0)/SIN(THETAW)
        CAZ(2)=1.D0/2/SQRT(2d0)/SIN(THETAW)
        CVZE=(-LEPCHARGE)*1.D0/2/SQRT(2d0)/SIN(THETAW)
        CAZE=1.D0/2/SQRT(2d0)/SIN(THETAW)
      ENDIF

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
C-----------------------------------------------------------------------
C     SUBTRACTION FLAGS ! 1=BORN-PROJECTED 0=DIPOLES
      ISUB=1 !wip
      IF (ISUB.eq.1) WRITE (6,'(A)')  '  Subtraction: P2B'
      IF (ISUB.eq.0) WRITE (6,'(A)')  '  Subtraction: Dipoles only'
C-----------------------------------------------------------------------                    
C     CHANGE AT YOUR OWN RISK IF DARING TO FACE THE CONSEQUENCES
C     CHOOSE CUTOFF       
      CUTOFF=1D-8
      WRITE (6,*) ' CUTOFF=',CUTOFF 
C     PARAMETERS RELATED TO IMPORTANCE SAMPLING      
      NPOW(1)=2
      NPOW(2)=4
        
      XPOW(1)=1-1D0/NPOW(1)
      XPOW(2)=1-1D0/NPOW(2)
C-----------------------------------------------------------------------
C---SCHEME IS 0 FOR MSbar AND OTHERWISE FOR DIS (OBSOLETE!!!)
      SCHEME=0
C      write(6,*) 'SCHEME =', SCHEME
C---SCALE IS FACTORIZATION SCALE**2/Q**2
C X-SECTION IN PB Gev to nb to pb
      NRM=389385*1000
      DO KNLOBLK=1,NNLOBLK
        NLOBLOCKSUM(KNLOBLK)=0D0
        NLOBLOCKQR(KNLOBLK)=0D0
      ENDDO
      
C FROM HERE
C JETALG IS THE ROUTINE THAT APPLIES THE JET ALGORITHM TO FINAL STATE PARTONS
C IT HAS TO BE CALLED EVERYTIME THE KINEMATICS IS MODIFIED
C TO ACCOUNT FOR PROPER JET COMBINATION
C AND FACTORIZATION AND RENORMALIZATION SCALES       
      
      
C---START MAIN LOOP
      DO I=1,NEV
        IF (MOD(I,100 000).EQ.1) CALL RANGEN(-I,SFOR)

C--- STARTS ORDER ALPHAS^0  --------------------

C---GENERATE A TWO-PARTON STATE
        CALL GENTWO(S,P,WTWO,CUTS,*1000)
        CALL JETALG(P)
        SCALE=SCALEFACT/Q2        
C---EVALUATE THE TWO-PARTON TREE-LEVEL MATRIX ELEMENT
        CALL MATTWO(P,MTWO)
C---GIVE IT TO THE USER
        CALL VECMUL(13,NRM*WTWO/NEV,MTWO,WEIGHT)
        CALL GETWEIGHT(2,0,0,P,S,WEIGHT,SCALEREN,SCALEFACT,Q2,NF,WWW,1)             
        IF (ISUB.EQ.0) CALL USER(2,0,0,P,S,WWW)   
C--- COMPUTE BORNPROJ NORMALIZATION (STRUCTURE FUNCTION)         
        IF (ISUB.EQ.1) THEN
          CALL BORNPROJ(S,PBORNP,WBORN,*1000)
          CALL USER(5,0,0,PBORNP,S,WBORN)
        ENDIF 
         
C--- STARTS ORDER ALPHAS^1  --------------------

C---KEEP USING TWO-PARTON KINEMATICS (ALSO TO EVALUATE SCALE FOR COLTHR)
        CALL JETALG(P)
        SCALE=SCALEFACT/Q2            
C---EVALUATE THE TWO-PARTON ONE-LOOP MATRIX ELEMENT
        CALL VIRTWO(S,P,VTWO,*1000)
C---GIVE IT TO THE USER
        CALL VECMUL(13,NRM*WTWO/NEV,VTWO,WEIGHT)
        CALL GETWEIGHT(2,1,2,P,S,WEIGHT,SCALEREN,SCALEFACT,Q2,NF,WWW,0)                
        IF (ISUB.EQ.0) CALL USER(2,1,2,P,S,WWW)

C---EVALUATE THE THREE-PARTON COLLINEAR SUBTRACTION (CHANGES P)
        CALL COLTHR(S,P,CTHR,*1000)
        CALL JETALG(P)
        SCALE=SCALEFACT/Q2
C---GIVE IT TO THE USER
        CALL VECMUL(13,NRM*WTWO/NEV,CTHR,WEIGHT)
        CALL GETWEIGHT(3,1,3,P,S,WEIGHT,SCALEREN,SCALEFACT,Q2,NF,WWW,1)                
        IF (ISUB.EQ.0)     CALL USER(3,1,3,P,S,WWW)

C---GENERATE A THREE-PARTON STATE
        CALL GENTHR(P,WTHR,*1000)
C---CALCULATE THE JACOBIAN FACTOR AND SUBTRACTION CONFIGURATIONS
        JTHR=0
        DO J=1,NPERM3
          CALL SUBTHR(J,S,P,Q(1,1,J),STHR(-6,J),JTMP,*1000)
          JTHR=JTHR+JTMP
        ENDDO
c        IF (DOT(P,1,2).LT.1D-6) WRITE(6,*) 'SUB3:', -STHR(0,1), Q2       
C---GIVE THE SUBTRACTION CONFIGURATIONS TO THE USER
        DO J=1,NPERM3
          CALL VECMUL(13,-NRM*WTWO*WTHR/JTHR/NEV,STHR(-6,J),WEIGHT)
          CALL JETALG(Q(1,1,J))
          SCALE=SCALEFACT/Q2
          CALL GETWEIGHT(3,1,1,Q(1,1,J),S,WEIGHT,
     $                                 SCALEREN,SCALEFACT,Q2,NF,WWW,0)                
          IF (ISUB.EQ.0) CALL USER(3,1,1,Q(1,1,J),S,WWW)
        ENDDO

C---EVALUATE THE THREE-PARTON TREE-LEVEL MATRIX ELEMENT
        CALL MATTHR(P,MTHR)
C---GIVE IT TO THE USER
c        IF (DOT(P,1,2).LT.1D-6) WRITE(6,*) 'MAT3:', MTHR(0)  
        CALL VECMUL(13,NRM*WTWO*WTHR/JTHR/NEV,MTHR,WEIGHT)
        CALL JETALG(P)
        SCALE=SCALEFACT/Q2
        CALL GETWEIGHT(3,1,0,P,S,WEIGHT,SCALEREN,SCALEFACT,Q2,NF,WWW,1)                
        CALL USER(3,1,0,P,S,WWW)

C---GIVE BORNPROJECTED COUNTER-EVENT FOR THREE-PARTON TREE-LEVEL MATRIX ELEMENT
C        CALL PSBORN(P,PBORNP)
        CALL JETALG(PBORNP)
        SCALE=SCALEFACT/Q2
        IF (ISUB.EQ.1)   CALL USER(3,1,5,PBORNP,S,-WWW)

C--- (HIGHER ORDER SKIPPER)      
        IF (INNLO.NE.1) GOTO 1000               
 
C--- STARTS ORDER ALPHAS^2  --------------------
C---  EVALUATE THE THREE-PARTON ONE-LOOP MATRIX ELEMENT
          CALL JETALG(P)
          SCALE=SCALEFACT/Q2
          CALL VIRTHR(S,P,VTHR,*1000)
C---  GIVE IT TO THE USER
          CALL VECMUL(13,NRM*WTWO*WTHR/JTHR/NEV,VTHR,WEIGHT)
          CALL GETWEIGHT(3,2,2,P,S,WEIGHT,SCALEREN,SCALEFACT,
     $                                                    Q2,NF,WWW,0)                
          CALL USER(3,2,2,P,S,WWW)

C---GIVE BORNPROJECTED COUNTER-EVENT FOR THREE-PARTON ONE-LOOP MATRIX ELEMENT
c          CALL PSBORN(P,PBORNP)
          CALL JETALG(PBORNP)
          SCALE=SCALEFACT/Q2
          IF (ISUB.EQ.1)   CALL USER(3,2,5,PBORNP,S,-WWW)                        
C---KEEP USING THREE-PARTON KINEMATICS FROM P (ALSO EVALUATE SCALE FOR COLFOR)
C--- CALL JETALG TO MAKE SURE SCALE IS RIGHT (not modified by BORNPROJ) USED By COLFOR
          CALL JETALG(P)
          SCALE=SCALEFACT/Q2

C---  EVALUATE THE FOUR-PARTON COLLINEAR SUBTRACTION (and change phase space P)
          CALL COLFOR(S,P,CFOR,*1000)
          CALL JETALG(P)
          SCALE=SCALEFACT/Q2          
C---  GIVE IT TO THE USER
          CALL VECMUL(13,NRM*WTWO*WTHR/JTHR/NEV,CFOR,WEIGHT)
          CALL GETWEIGHT(4,2,3,P,S,WEIGHT,SCALEREN,SCALEFACT,
     $                                                    Q2,NF,WWW,2)                
          CALL USER(4,2,3,P,S,WWW)

C---GIVE BORNPROJECTED COUNTER-EVENT FOR FOUR-PARTON COLLINEAR SUBTRACTION
c          CALL PSBORN(P,PBORNP)
          CALL JETALG(PBORNP)
          SCALE=SCALEFACT/Q2
          IF (ISUB.EQ.1)    CALL USER(4,2,5,PBORNP,S,-WWW)

C---  GENERATE A FOUR-PARTON STATE
          CALL GENFOR(P,WFOR,*1000)
          CALL JETALG(P)
          SCALE=SCALEFACT/Q2
C---  CALCULATE THE JACOBIAN FACTOR AND SUBTRACTION CONFIGURATIONS
          JFOR=0
          DO J=1,NPERM4
            CALL SUBFOR(J,S,P,Q(1,1,J),SFOR(-6,J),JTMP,*1000)
            JFOR=JFOR+JTMP
          ENDDO
c        IF (DOT(P,3,4).LT.1D-6) WRITE(6,*) 'SUBFOR_5:', SFOR(1,5)
c        IF (DOT(P,3,4).LT.1D-6) WRITE(6,*) 'SUBFOR_6:', SFOR(1,6)          

C---  EVALUATE THE FOUR-PARTON TREE-LEVEL MATRIX ELEMENT
c          CALL JETALG(P)  !debug
          SCALE=SCALEFACT/Q2         
          CALL MATFOR(P,MFOR)
c        IF (DOT(P,1,2).LT.1D-6) WRITE(6,*) 'MATFOR:', MFOR(1)
C---  GIVE IT TO THE USER
          CALL VECMUL(13,NRM*WTWO*WTHR*WFOR/JFOR/NEV,MFOR,WEIGHT)
          CALL JETALG(P)

          CALL GETWEIGHT(4,2,0,P,S,WEIGHT,SCALEREN,SCALEFACT,
     $                                                    Q2,NF,WWW,1)

          CALL USER(4,2,0,P,S,WWW)

C---GIVE BORNPROJECTED COUNTER-EVENT FOR FOUR-PARTON TREE-LEVEL MATRIX ELEMENT
c          CALL PSBORN(P,PBORNP)
          CALL JETALG(PBORNP)
          SCALE=SCALEFACT/Q2
          IF (ISUB.EQ.1)    CALL USER(4,2,5,PBORNP,S,-WWW)

C---  GIVE THE SUBTRACTION CONFIGURATIONS TO THE USER
          DO J=1,NPERM4
            CALL VECMUL(13,
     $                 -NRM*WTWO*WTHR*WFOR/JFOR/NEV,SFOR(-6,J),WEIGHT)
            CALL JETALG(Q(1,1,J))
            SCALE=SCALEFACT/Q2
            IF (SCALE.NE.SCALE) WRITE(6,*) J,Q(:,:,J)
            CALL GETWEIGHT(4,2,1,Q(1,1,J),S,WEIGHT,SCALEREN,
     $                                          SCALEFACT,Q2,NF,WWW,2)

c     (for even faster results use the following if the scale does not
c      depend on the event, eg scale=Q2):
c            CALL GETWEIGHT(4,2,1,Q(1,1,J),S,WEIGHT,SCALEREN,
c     $                                          SCALEFACT,Q2,NF,WWW,0)                
c            IF (ABS(WWW(2)-WWW(1)).GT.5.D0) WRITE(6,*) 'PERM',J
            CALL USER(4,2,1,Q(1,1,J),S,WWW)

C---GIVE BORNPROJECTED COUNTER-EVENT FOR SUBTRACTION CONFIGURATIONS 
c            CALL PSBORN(Q(1,1,J),PBORNP)
            CALL JETALG(PBORNP)
            SCALE=SCALEFACT/Q2
            IF (ISUB.EQ.1)   CALL USER(4,2,5,PBORNP,S,-WWW)                                
          ENDDO
C---TELL THE USER THAT THE EVENT IS COMPLETE
 1000   CALL USER(0,0,0,P,S,WWW)
      ENDDO
      IF (IPOL.EQ.1) CALL REPORT_NLO_BLOCKS(NEV)
      WRITE (6,'(/A)') 'Thanks for using POLDIS MF, Please cite:'
      WRITE (6,'(2A/)') '  PhysRevLett.125.082001 (2020) ',
     $ '  (https://arxiv.org/abs/2005.10705)'
      END
C-----------------------------------------------------------------------
      SUBROUTINE ACCUMULATE_NLO_BLOCKS(PREF,G2NLOQ,G2NLOG,G4NLOQ,
     $ GLNLOQ,Y)
      IMPLICIT NONE
      INTEGER I,NNLOBLK
      PARAMETER (NNLOBLK=8)
      DOUBLE PRECISION PREF,G2NLOQ,G2NLOG,G4NLOQ,GLNLOQ,Y
      DOUBLE PRECISION NLOBLOCKSUM(NNLOBLK),NLOBLOCKQR(NNLOBLK)
      DOUBLE PRECISION VALUES(NNLOBLK),YMINUS,YPLUS,MINUSY2
      COMMON /NLOBLOCKACC/ NLOBLOCKSUM,NLOBLOCKQR

      YMINUS=1D0-(1D0-Y)**2
      YPLUS=1D0+(1D0-Y)**2
      MINUSY2=-Y**2

      VALUES(1)=PREF*G2NLOQ
      VALUES(2)=PREF*G2NLOG
      VALUES(3)=PREF*G4NLOQ
      VALUES(4)=PREF*GLNLOQ
      VALUES(5)=PREF*YMINUS*(G2NLOQ+G2NLOG)
      VALUES(6)=PREF*YPLUS*G4NLOQ
      VALUES(7)=PREF*MINUSY2*GLNLOQ
      VALUES(8)=VALUES(5)+VALUES(6)+VALUES(7)

      DO I=1,NNLOBLK
        NLOBLOCKSUM(I)=NLOBLOCKSUM(I)+VALUES(I)
        NLOBLOCKQR(I)=NLOBLOCKQR(I)+VALUES(I)**2
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE REPORT_NLO_BLOCKS(NEV)
      IMPLICIT NONE
      INTEGER I,NEV,NNLOBLK
      PARAMETER (NNLOBLK=8)
      DOUBLE PRECISION DN,ERR,VAR
      DOUBLE PRECISION NLOBLOCKSUM(NNLOBLK),NLOBLOCKQR(NNLOBLK)
      CHARACTER*24 NLOBLOCKLABEL(NNLOBLK)
      COMMON /NLOBLOCKACC/ NLOBLOCKSUM,NLOBLOCKQR
      DATA NLOBLOCKLABEL /
     $ 'G2NLOq_pb',
     $ 'G2NLOg_pb',
     $ 'G4NLOq_pb',
     $ 'GLNLOq_pb',
     $ 'Yminus_times_G2_pb',
     $ 'Yplus_times_G4_pb',
     $ 'minusY2_times_GL_pb',
     $ 'NLO_pol_pb' /

      DN=DBLE(NEV)
      WRITE(*,'(/A)') 'POLDIS polarized NLO block summary'
      WRITE(*,'(A)') '================================='
      DO I=1,NNLOBLK
        VAR=NLOBLOCKQR(I)-NLOBLOCKSUM(I)*NLOBLOCKSUM(I)/DN
        IF (VAR.LT.0D0) VAR=0D0
        ERR=SQRT(VAR)
        WRITE(*,'(A,1X,A,1X,A,ES24.16,1X,A,ES24.16)')
     $   'POLDIS_NLO_BLOCK',
     $   'label='//NLOBLOCKLABEL(I),
     $   'value=', NLOBLOCKSUM(I),
     $   'error=', ERR
      ENDDO

      END
C-----------------------------------------------------------------------
      FUNCTION DOT(P,I,J)
      IMPLICIT NONE
C---RETURN THE DOT PRODUCT OF P(*,I) AND P(*,J)
      INTEGER I,J
      DOUBLE PRECISION DOT,P(4,7)
      DOT=P(4,I)*P(4,J)-P(3,I)*P(3,J)-P(2,I)*P(2,J)-P(1,I)*P(1,J)
      END
C-----------------------------------------------------------------------
      SUBROUTINE VECMUL(N,S,A,B)
      IMPLICIT NONE
C---MULTIPLY THE VECTOR A BY THE SCALAR S TO GIVE THE VECTOR B
      INTEGER N,I
      DOUBLE PRECISION S,A(N),B(N)
      DO I=1,N
        B(I)=S*A(I)
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENTWO(S,P,W,CUTS,*)
      IMPLICIT NONE
C---GENERATE A TWO-PARTON CONFIGURATION
C   WITHIN THE PHASE-SPACE LIMITS GIVEN IN CUTS()
C   THE WEIGHT GIVES THE TOTAL VOLUME OF PHASE-SPACE
      DOUBLE PRECISION S,P(4,7),W,R(2),XMIN,XMAX,QMIN,QMAX,YMIN,YMAX,
     $     Q2,QJAC,Y,YJAC,E,Ee,Ep,X,eta,dot,PBORNV(4,7)
      INTEGER SCHEME,NF,IFRAME,I,J,IPOL
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGFRAME/IFRAME
      COMMON/EeEp/Ee,Ep
      COMMON/PBORNV/PBORNV
      COMMON/PBORNPAR/X,Y,Q2,YJAC,QJAC
      COMMON/FLAGPOL/IPOL
      EXTERNAL CUTS
      CALL RANGEN(2,R)
C---FIRST GET THE USER'S CUTS
      CALL CUTS(S,XMIN,XMAX,QMIN,QMAX,YMIN,YMAX)
C---CHECK THAT SYSTEM IS NOT 0VER-CONSTRAINED
      IF (XMIN.EQ.XMAX.AND.XMIN.NE.0.AND.
     $    YMIN.EQ.YMAX.AND.YMIN.NE.0.AND.
     $    QMIN.EQ.QMAX.AND.QMIN.NE.0)
     $     STOP 'AT MOST TWO KINEMATIC VARIABLES CAN BE CONSTRAINED!'
C---INSERT KINEMATIC LIMITS WHERE NONE ARE GIVEN
      IF (XMAX.LE.0) XMAX=1
      IF (YMAX.LE.0) YMAX=1
      IF (QMAX.LE.0) QMAX=S
C---CHECK IF Q2 IS CONSTRAINED BY X AND Y
      IF (XMAX*YMAX*S.LT.QMAX) QMAX=XMAX*YMAX*S
      IF (XMIN*YMIN*S.GT.QMIN) QMIN=XMIN*YMIN*S
C---CHECK THERE IS SOME RANGE AVAILABLE IN Q2
      IF (QMAX.LT.QMIN) STOP 'NO PHASE-SPACE AVAILABLE IN GENTWO!'
      IF (QMIN.EQ.0) STOP 'Q2 MUST BE CONSTRAINED LARGER THAN ZERO!'
C---GENERATE A Q2 VALUE
      IF (QMAX.EQ.QMIN) THEN
        Q2=QMAX
        QJAC=1
      ELSE
C---CONVENTIONAL GENERATION      
        Q2=(QMAX/QMIN)**R(1)*QMIN  !debug
        QJAC=LOG(QMAX/QMIN)*Q2
c        IF (IPOL.EQ.0) THEN
C---MORE EVENTS AT SMALL Q2 GENERATION        
        Q2=(QMAX/QMIN)**(R(1)**2)*QMIN
        QJAC=LOG(QMAX/QMIN)*Q2*2*R(1)
c        ENDIF
c        Q2=(QMAX-QMIN)*R(1)+QMIN
c        QJAC=QMAX-QMIN 
      ENDIF
C---CHECK IF Y IS CONSTRAINED BY X AND Q2
      IF (XMIN*YMAX*S.GT.Q2) YMAX=Q2/(XMIN*S)
      IF (XMAX*YMIN*S.LT.Q2) YMIN=Q2/(XMAX*S)
C---CHECK THERE IS SOME RANGE AVAILABLE IN Y
      IF (YMAX.LT.YMIN) STOP 'NO PHASE-SPACE AVAILABLE IN GENTWO!'
C---GENERATE A Y VALUE
      IF (YMAX.EQ.YMIN) THEN
        Y=YMAX
        YJAC=1
      ELSE
C---CONVENTIONAL GENERATION      
        Y=(YMAX/YMIN)**R(2)*YMIN
        YJAC=LOG(YMAX/YMIN)*Y
C---MORE EVENTS AT SMALL Y GENERATION   R(2)^1.5     
        Y=(YMAX/YMIN)**(R(2)**1.5)*YMIN
        YJAC=LOG(YMAX/YMIN)*Y*1.5*R(2)**0.5
c        IF (IPOL.EQ.0) THEN
C---EVEN MORE EVENTS AT SMALL Y GENERATION   R(2)^2      
        Y=(YMAX/YMIN)**(R(2)**2)*YMIN
        YJAC=LOG(YMAX/YMIN)*Y*2*R(2)             
c        ENDIF 
      ENDIF
C.......................................................................
      IF (IFRAME.EQ.1) THEN
C     CONSTRUCT MOMENTA IN THE BREIT FRAME
C.......................................................................
      E=SQRT(Q2)/2
      X=Q2/Y/S
      P(1,1)=   0
      P(2,1)=   0
      P(3,1)=   E
      P(4,1)=   E
      
      P(1,2)=   0
      P(2,2)=   0
      P(3,2)=  -E
      P(4,2)=   E
      
      P(1,5)=   0
      P(2,5)=   0
      P(3,5)=-2*E
      P(4,5)=   0
      
      P(1,6)=   E/Y*2*SQRT(1-Y)
      P(2,6)=   0
      P(3,6)=  -E
      P(4,6)=   E/Y*(2-Y)
      
      P(1,7)=   E/Y*2*SQRT(1-Y)
      P(2,7)=   0
      P(3,7)=   E
      P(4,7)=   E/Y*(2-Y)
C Correction by DdeF March 2020     
      P(1,3)=   0
      P(2,3)=   0
      P(3,3)=   0
      P(4,3)=   0  
         
      P(1,4)=   0
      P(2,4)=   0
      P(3,4)=   0
      P(4,4)=   0
C.......................................................................
      ELSEIF (IFRAME.EQ.0) THEN
C     LAB FRAME      
C.......................................................................
      X=Q2/Y/S
      P(1,1)=   0
      P(2,1)=   0
      P(3,1)=   X*Ep
      P(4,1)=   X*Ep
             
      P(1,5)= -(Q2*(1-Y))**0.5
      P(2,5)=   0
      P(3,5)=-Y*(Ee+X*Ep)
      P(4,5)= Y*(Ee-X*Ep)
      
      P(1,6)=   0
      P(2,6)=   0
      P(3,6)=  -Ee
      P(4,6)=   Ee
      
      P(1,7)=  P(1,6)-P(1,5)
      P(2,7)=  P(2,6)-P(2,5)
      P(3,7)=  P(3,6)-P(3,5)
      P(4,7)=  P(4,6)-P(4,5)    
        
      P(1,2)=  P(1,1)+P(1,5)
      P(2,2)=  P(2,1)+P(2,5)
      P(3,2)=  P(3,1)+P(3,5)
      P(4,2)=  P(4,1)+P(4,5)  

      P(1,3)=   0
      P(2,3)=   0
      P(3,3)=   0
      P(4,3)=   0   
        
      P(1,4)=   0
      P(2,4)=   0
      P(3,4)=   0
      P(4,4)=   0    
C.......................................................................  
      ELSEIF (IFRAME.EQ.2) THEN
C     eP CM FRAME      
C.......................................................................
      X=Q2/Y/S
      E=S**0.5/2
      P(1,1)=   0
      P(2,1)=   0
      P(3,1)=   X*E
      P(4,1)=   X*E
             
      P(1,5)= -(Q2*(1-Y))**0.5
      P(2,5)=   0
      P(3,5)=-Y*(E+X*E)
      P(4,5)= Y*(E-X*E)
      
      P(1,6)=   0
      P(2,6)=   0
      P(3,6)=  -E
      P(4,6)=   E
      
      P(1,7)=  P(1,6)-P(1,5)
      P(2,7)=  P(2,6)-P(2,5)
      P(3,7)=  P(3,6)-P(3,5)
      P(4,7)=  P(4,6)-P(4,5)    
        
      P(1,2)=  P(1,1)+P(1,5)
      P(2,2)=  P(2,1)+P(2,5)
      P(3,2)=  P(3,1)+P(3,5)
      P(4,2)=  P(4,1)+P(4,5)  

      P(1,3)=   0
      P(2,3)=   0
      P(3,3)=   0
      P(4,3)=   0   
        
      P(1,4)=   0
      P(2,4)=   0
      P(3,4)=   0
      P(4,4)=   0    
C.......................................................................  
      ENDIF
      DO J=1,7
        DO I=1,4
          PBORNV(I,J)=P(I,J)
        ENDDO
      ENDDO
C---CALCULATE WEIGHT
      W=Y/(16*PI*Q2**2)*QJAC*YJAC
      END    
C-----------------------------------------------------------------------
      SUBROUTINE GENTHR(P,W,*)
      IMPLICIT NONE
C---GENERATE A THREE-PARTON CONFIGURATION FROM A GIVEN TWO-PARTON ONE
      DOUBLE PRECISION P(4,7),W,R(2),X,XJAC,Z,EMSQ,ABS,DOT
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      CALL RANGEN(2,R)
C---FIND OUT WHAT X VALUE WAS GENERATED AND THE JACOBIAN FACTOR
      CALL GETCOL(X,XJAC)
C---GENERATE A Z VALUE
      Z=1-R(1)**NPOW(1)
      IF (R(2).GT.0.5) Z=1-Z
C---ENFORCE INVARIANT MASS CUTOFF
      IF (Z.LT.CUTOFF.OR.1-Z.LT.CUTOFF) RETURN 1
C---GENERATE THEIR MOMENTA
      CALL GENDEC(P,2,3,Z,*999)
C---CALCULATE WEIGHT
      EMSQ=P(4,5)**2-P(3,5)**2-P(2,5)**2-P(1,5)**2
      W=-EMSQ/(16*PISQ)
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENFOR(P,W,*)
      IMPLICIT NONE
C---GENERATE A FOUR-PARTON CONFIGURATION FROM A GIVEN THREE-PARTON ONE
      INTEGER EMIT,OMIT,I
      DOUBLE PRECISION P(4,7),W,R(6),X,XJAC,Z,EMSQ,XX(2),T
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      CALL RANGEN(6,R)
C---CHOOSE WHICH PARTON TO SPLIT UNIFORMLY
      EMIT=INT(R(1)*2)+2
C---AND WHICH OF THE OTHERS TO BALANCE ITS MOMENTUM WITH
      OMIT=INT(R(2)*2)+1
      IF (OMIT.EQ.EMIT) OMIT=3
C---IF INCOMING PARTON IS THE SPECTATOR GENERATE AN INITIAL-STATE DIPOLE
      IF (OMIT.EQ.1) THEN
C---FIND OUT WHAT X VALUE WAS GENERATED AND THE JACOBIAN FACTOR
        CALL GETCOL(X,XJAC)
C---GENERATE A Z VALUE
        Z=1-R(3)**NPOW(2)
        IF (R(4).GT.0.5) Z=1-Z
C---ENFORCE INVARIANT MASS CUTOFF
        IF (Z.LT.CUTOFF.OR.1-Z.LT.CUTOFF) RETURN 1
C---GENERATE THEIR MOMENTA
        CALL GENDEC(P,EMIT,4,Z,*999)
C---OTHERWISE GENERATE A FINAL-STATE DIPOLE
      ELSE
C---AND FORGET THE PROPOSED COLLINEAR SPLITTING
        DO I=1,4
          P(I,1)=P(I,1)-P(I,4)
        ENDDO
C---GENERATE X1,X2 VALUES
        XX(1)=1-MIN(R(3),R(4))**NPOW(2)
        XX(2)=2-XX(1)-(R(5)+(1-R(5))*(1-XX(1))**(1D0/NPOW(2)))**NPOW(2)
C---ENFORCE INVARIANT MASS CUTOFFS
        IF (1-XX(1).LT.CUTOFF.OR.1-XX(2).LT.CUTOFF
     $       .OR.XX(1)+XX(2)-1.LT.CUTOFF) RETURN 1
C---GENERATE THEIR MOMENTA
        CALL GENDIP(P,EMIT,OMIT,4,XX,*999)
      ENDIF
C---HALF THE TIME, SWAP PARTICLES 3 AND 4 SO THEY ARE SYMMETRIC
      IF (R(6).GT.0.5) THEN
        DO I=1,4
          T=P(I,3)
          P(I,3)=P(I,4)
          P(I,4)=T
        ENDDO
      ENDIF
C---CALCULATE WEIGHT
      EMSQ=P(4,5)**2-P(3,5)**2-P(2,5)**2-P(1,5)**2
      W=-EMSQ/(16*PISQ)
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENDEC(P,I1,I2,Z,*)
      IMPLICIT NONE
C---GENERATE A 2->2 DECAY
      DOUBLE PRECISION P(4,7),Z,R(1),Q(4,2),PTDQ,C(4),D(4)
      INTEGER I1,I2,I
      LOGICAL FIRST
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DATA FIRST/.TRUE./
      CALL RANGEN(1,R)
      PTDQ=Z*(1-Z)
      IF (PTDQ.LT.0) THEN
        IF (FIRST) THEN
          FIRST=.FALSE.
          WRITE (*,*) 'Numerical error in GENDEC!'
          WRITE (*,*) 'Please report this to seymour@surya11.cern.ch,'
          WRITE (*,*) 'giving the date of this version,'
          WRITE (*,*) 'and the latest value of ISEED.'
          WRITE (*,*) 'Event generation has not been affected.'
        ELSE
          WRITE (*,*) 'Another numerical error in GENDEC!'
        ENDIF
        RETURN 1
      ENDIF
      PTDQ=SQRT(PTDQ)
      CALL GTPERP(PTDQ,P,I2,I1,6,C,D)
      DO I=1,4
        Q(I,1)=Z*P(I,I1)+(1-Z)*P(I,I2)
     $       +COS(R(1)*2*PI)*C(I)+SIN(R(1)*2*PI)*D(I)
        Q(I,2)=(1-Z)*P(I,I1)+Z*P(I,I2)
     $       -COS(R(1)*2*PI)*C(I)-SIN(R(1)*2*PI)*D(I)
      ENDDO
      DO I=1,4
        P(I,I1)=Q(I,1)
        P(I,I2)=Q(I,2)
      ENDDO
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENDIP(P,EMIT,OMIT,THRD,X,*)
      IMPLICIT NONE
C---GENERATE A 2->3 DIPOLE EMISSION
      DOUBLE PRECISION P(4,7),X(2),R(1),Q(4,2),PTDQ,C(4),D(4)
      INTEGER EMIT,OMIT,THRD,I
      LOGICAL FIRST
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DATA FIRST/.TRUE./
      CALL RANGEN(1,R)
      PTDQ=(X(1)+X(2)-1)*(1-X(1))*(1-X(2))/X(1)**2
      IF (PTDQ.LT.0) THEN
        IF (FIRST) THEN
          FIRST=.FALSE.
          WRITE (*,*) 'Numerical error in GENDIP!'
          WRITE (*,*) 'Please report this to seymour@surya11.cern.ch,'
          WRITE (*,*) 'giving the date of this version,'
          WRITE (*,*) 'and the latest value of ISEED.'
          WRITE (*,*) 'Event generation has not been affected.'
        ELSE
          WRITE (*,*) 'Another numerical error in GENDIP!'
        ENDIF
        RETURN 1
      ENDIF
      PTDQ=SQRT(PTDQ)
      CALL GTPERP(PTDQ,P,EMIT,OMIT,6,C,D)
      DO I=1,4
        Q(I,1)=X(1)*P(I,OMIT)
        Q(I,2)=(1-(1-X(2))/X(1))*P(I,EMIT)
     $       +(1-X(1))*(1-X(2))/X(1)*P(I,OMIT)
     $       +COS(R(1)*2*PI)*C(I)+SIN(R(1)*2*PI)*D(I)
      ENDDO
      DO I=1,4
        P(I,THRD)=P(I,EMIT)+P(I,OMIT)-Q(I,1)-Q(I,2)
        P(I,OMIT)=Q(I,1)
        P(I,EMIT)=Q(I,2)
      ENDDO
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENCOL(I,X,XJAC,XMIN)
      IMPLICIT NONE
C---GENERATE AN X VALUE AND STORE IT FOR LATER RETRIEVAL
      INTEGER I
      DOUBLE PRECISION X,XJAC,XMIN,XL,XLJAC,R(2)
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      SAVE XL,XLJAC
      CALL RANGEN(2,R)
      IF (R(1).LT.0.5) THEN
        X=1-(1-XMIN)*R(2)**NPOW(I)
      ELSE
        X=XMIN**R(2)
      ENDIF
      IF (X.EQ.1) THEN
        XJAC=0
      ELSE
        XJAC=1/(0.5/(-X*LOG(XMIN))
     $       +0.5*((1-XMIN)/(1-X))**XPOW(I)/(NPOW(I)*(1-XMIN)))
      ENDIF
      XL=X
      XLJAC=XJAC
      RETURN
      ENTRY GETCOL(X,XJAC)
C---RETURN THE SAME VALUES AS LAST TIME
      X=XL
      XJAC=XLJAC
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTPERP(PTDQ,P,I,J,K,C,D)
      IMPLICIT NONE
C---FIND THE VECTORS PERPENDICULAR TO P(I) AND P(J)
C   C AND D ARE PURELY SPACE-LIKE VECTORS IN THE P(I)+P(J) CMF,
C   WITH C IN THE SAME PLANE AS P(K) AND D PERPENDICULAR TO IT,
C   BOTH HAVING LENGTH PTDQ*SQRT(2*DOT(P,I,J))
      DOUBLE PRECISION PTDQ,P(4,7),C(4),D(4),PTF,DIJ,DIK,DJK,DOT,EPS4
      INTEGER I,J,K,L
      DIJ=DOT(P,I,J)
      DIK=DOT(P,I,K)
      DJK=DOT(P,J,K)
      PTF=PTDQ/SQRT(DIK*DJK)
      DO L=1,4
        C(L)=PTF*(DIJ*P(L,K)-DJK*P(L,I)-DIK*P(L,J))
      ENDDO
      DO L=1,4
        D(L)=EPS4(L,P(1,I),P(1,J),C)/DIJ
      ENDDO
      END
C-----------------------------------------------------------------------
      FUNCTION EPS4(I,A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION EPS4,EPS3,A(4),B(4),C(4),AA(3),BB(3),CC(3)
      INTEGER I,J,K,S(4)
      DATA S/+1,-1,+1,+1/
      J=1
      DO K=1,3
        IF (I.EQ.J) J=J+1
        AA(K)=A(J)
        BB(K)=B(J)
        CC(K)=C(J)
        J=J+1
      ENDDO
      EPS4=0
      DO J=1,3
        EPS4=EPS4+CC(J)*EPS3(J,AA,BB)
      ENDDO
      EPS4=S(I)*EPS4
      END
C-----------------------------------------------------------------------
      FUNCTION EPS3(I,A,B)
      IMPLICIT NONE
      DOUBLE PRECISION EPS3,A(3),B(3),AA(2),BB(2)
      INTEGER I,J,K,S(3)
      DATA S/+1,-1,+1/
      J=1
      DO K=1,2
        IF (I.EQ.J) J=J+1
        AA(K)=A(J)
        BB(K)=B(J)
        J=J+1
      ENDDO
      EPS3=S(I)*(AA(1)*BB(2)-AA(2)*BB(1))
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATTWO(P,M)
      IMPLICIT NONE
C---EVALUATE THE TWO-PARTON MATRIX ELEMENT SQUARED FOR THE GIVEN
C   CONFIGURATION.
      INTEGER I,IPOL,EWFLAG
      DOUBLE PRECISION P(4,7),M(-6:6),Q,DOT
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,GZ,MZ,
     #                 GAMMAZ,LEPCHARGE,CKM2(3,3) 
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG       
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL     
C-----SELECT POLARIZED/UNPOLARIZED MATRIX ELEMENT
C#######################################################################
      IF (IPOL.EQ.0) THEN
      Q=4*(4*PI)**2/DOT(P,5,5)**2*
     $     (DOT(P,1,6)**2+DOT(P,1,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)
C.......................................................................
      IF (EWFLAG.EQ.0) GOTO 79 
      QZ=4*(4*PI)**2/DOT(P,5,5)**2*
     $     (DOT(P,1,6)**2-DOT(P,1,7)**2+DOT(P,2,7)**2-DOT(P,2,6)**2)                
C#######################################################################
      ELSEIF(IPOL.EQ.1) THEN
      Q=4*(4*PI)**2/DOT(P,5,5)**2*
     $     (DOT(P,1,6)**2-DOT(P,1,7)**2+DOT(P,2,7)**2-DOT(P,2,6)**2)
C.......................................................................
      IF (EWFLAG.EQ.0) GOTO 79 
      QZ=4*(4*PI)**2/DOT(P,5,5)**2*
     $     (DOT(P,1,6)**2+DOT(P,1,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)     
C.......................................................................
      ENDIF
 79   CONTINUE

      DO I=-6,6
        M(I)=0
      ENDDO

C     COMBINE THE CROSS SECTION WITH COUPLINGS
      CALL COMBINEQI(P,M,Q,QZ)

      END
C-----------------------------------------------------------------------
      SUBROUTINE MATTHR(P,M)
      IMPLICIT NONE
C---EVALUATE THE THREE-PARTON MATRIX ELEMENT SQUARED FOR THE GIVEN
C   CONFIGURATION.
      INTEGER I,IPOL,EWFLAG
      DOUBLE PRECISION P(4,7),M(-6:6),QQ,GQ,DOT
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QQZ,GQZ,MZ,
     #                 GAMMAZ,LEPCHARGE,CKM2(3,3)
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL
      
C-----SELECT POLARIZED/UNPOLARIZED MATRIX ELEMENT
C#######################################################################
      IF (IPOL.EQ.0) THEN
      QQ=8*(4*PI)**2*
     $     (DOT(P,1,6)**2+DOT(P,1,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
      GQ=8*(4*PI)**2*
     $     (DOT(P,3,6)**2+DOT(P,3,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*TR/(-4*DOT(P,2,1)*DOT(P,3,1)*DOT(P,5,5))
C.......................................................................
C     Agrego las estructuras mixtas
      IF (EWFLAG.EQ.0) GOTO 89
C     Para quarks uso que AVAV(Q)=VVVV_pol(Q)
      QQZ=-8*(4*PI)**2*
     $     (-DOT(P,1,6)**2+DOT(P,1,7)**2-DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
C     Para gluones uso que AVAV(G)=-VVVV_pol(Q)
      GQZ=-8*(4*PI)**2*
     $     (-DOT(P,3,6)**2+DOT(P,3,7)**2-DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*TR/(-4*DOT(P,2,1)*DOT(P,3,1)*DOT(P,5,5))
C#######################################################################
      ELSEIF (IPOL.EQ.1) THEN
      QQ=-8*(4*PI)**2*
     $     (-DOT(P,1,6)**2+DOT(P,1,7)**2-DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
      GQ=8*(4*PI)**2*
     $     (-DOT(P,3,6)**2+DOT(P,3,7)**2+DOT(P,2,7)**2-DOT(P,2,6)**2)
     $     *16*PISQ*TR/(-4*DOT(P,2,1)*DOT(P,3,1)*DOT(P,5,5))
C.......................................................................
C     Agrego las estructuras mixtas
      IF (EWFLAG.EQ.0) GOTO 89
C     Para quarks uso que AVAV_pol(Q)=VVVV(Q)
      QQZ=8*(4*PI)**2*
     $     (DOT(P,1,6)**2+DOT(P,1,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
C     AVAV_pol(G)
      GQZ=-8*(4*PI)**2*
     $     (DOT(P,3,6)**2+DOT(P,3,7)**2-DOT(P,2,7)**2-DOT(P,2,6)**2)
     $     *16*PISQ*TR/(-4*DOT(P,2,1)*DOT(P,3,1)*DOT(P,5,5))
C#######################################################################  
      ENDIF
 89   CONTINUE
 
c      QQ=0.d0 !debug
c      GQ=0.d0 !debug
      DO I=-6,6
        M(I)=0
      ENDDO

C     COMBINE THE CROSS SECTION WITH COUPLINGS
      CALL COMBINEQI(P,M,QQ,QQZ)
      CALL COMBINEGLUON(P,M,GQ,GQZ)

      END
C-----------------------------------------------------------------------
      SUBROUTINE CONTHR(P,V,VV,Q,QBAR,G,M)
      IMPLICIT NONE
C---EVALUATE THE CONTRACTION OF THE THREE-PARTON MATRIX ELEMENT SQUARED
C   FOR THE GIVEN CONFIGURATION WITH THE VECTOR V.  IT IS ASSUMED THAT
C   V IS PERPENDICULAR TO THE GLUON, V.P3.EQ.0, AND THAT V.V.NE.0 .
C   THE NORMALIZATION IS SUCH THAT THE AVERAGE OVER V'S AZIMUTH AROUND
C   P3 IS EQUAL TO HALF OF MATTHR
C
C GPS: added VV as one of the arguments, to allow a better calculation 
C      of it by the calling routines
      INTEGER I,J,Q,QBAR,G
      DOUBLE PRECISION P(4,7),V(4),M,DOT,VDOT,EMSQ,OMSQ,VV,Y1,Y2,L1,L2,
     $     DOTS,T1,T2,T3,T4,T,V1,V2,V6
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      VDOT(I)=P(4,I)*V(4)-P(3,I)*V(3)-P(2,I)*V(2)-P(1,I)*V(1)
      EMSQ=DOT(P,5,5)
      OMSQ=1/EMSQ
C---  let VV be calculated by the calling routine
C      VV=V(4)**2-V(3)**2-V(2)**2-V(1)**2
      Y1=2*DOT(P,ABS(QBAR),ABS(G))*SIGN(1,QBAR*G)
      Y2=2*DOT(P,ABS(Q),ABS(G))*SIGN(1,Q*G)
      L1=2*DOT(P,ABS(Q),6)*SIGN(1,Q)
      L2=2*DOT(P,ABS(QBAR),6)*SIGN(1,QBAR)
      V1=VDOT(ABS(Q))*SIGN(1,Q)
      V2=VDOT(ABS(QBAR))*SIGN(1,QBAR)
      V6=VDOT(6)
      DOTS=1/(VV*(Y1+Y2)**2)
      T=CF*16*PISQ/(Y1*Y2)**2
      T1=T*4*(Y1**2*(L2**2+(EMSQ-L2)**2)+Y2**2*(L1**2+(EMSQ-L1)**2)
     $     +2*Y1*Y2*(L1*EMSQ+L2*EMSQ-2*L1*L2))*OMSQ**2*DOTS
      T2=T*4*(Y1*(EMSQ-2*L2)-Y2*(EMSQ-2*L1))*Y1*Y2*OMSQ**2*DOTS
      T3=T*8*(Y1*Y2*OMSQ)**2*DOTS
      T4=T*((EMSQ-Y1-L1-L2)**2+(EMSQ-Y2-L1-L2)**2)
     $     *Y1*Y2*OMSQ
      M=-T1*(V1*Y1-V2*Y2)**2
     $  -T2*2*(V1*Y1-V2*Y2)*(-V6*(Y2+Y1)+(EMSQ-L1-L2)*(V1+V2))
     $  -T3*(-V6*(Y2+Y1)+(EMSQ-L1-L2)*(V1+V2))**2
     $  +T4
      m=m*(4*pi)**2
      END
C-----------------------------------------------------------------------
      SUBROUTINE CONTHRPOLQ(P,V,VV,Q,QBAR,G,M)
      IMPLICIT NONE
C---EVALUATE THE CONTRACTION OF THE THREE-PARTON MATRIX ELEMENT SQUARED
C   FOR THE GIVEN CONFIGURATION WITH THE VECTOR V.  IT IS ASSUMED THAT
C   V IS PERPENDICULAR TO THE GLUON, V.P3.EQ.0, AND THAT V.V.NE.0 .
C   THE NORMALIZATION IS SUCH THAT THE AVERAGE OVER V'S AZIMUTH AROUND
C   P3 IS EQUAL TO HALF OF MATTHR
      INTEGER I,J,Q,QBAR,G
      DOUBLE PRECISION P(4,7),V(4),M,DOT,VDOT,EMSQ,OMSQ,VV,Y1,Y2,L1,L2,
     $     DOTS,V1,V2,V3,V6,V7,p13,p12,p16,p17,p23,p36,p37,T
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      VDOT(I)=P(4,I)*V(4)-P(3,I)*V(3)-P(2,I)*V(2)-P(1,I)*V(1)
      EMSQ=DOT(P,5,5)
      OMSQ=1/EMSQ
C---  let VV be calculated by the calling routine
C      VV=V(4)**2-V(3)**2-V(2)**2-V(1)**2
      Y1=2*DOT(P,ABS(QBAR),ABS(G))*SIGN(1,QBAR*G)
      Y2=2*DOT(P,ABS(Q),ABS(G))*SIGN(1,Q*G)

      p12=DOT(P,ABS(QBAR),ABS(G))*SIGN(1,QBAR*G)
      p13=DOT(P,ABS(QBAR),ABS(Q))*SIGN(1,QBAR*Q)
      p16=DOT(P,ABS(QBAR),6)*SIGN(1,QBAR)
      p17=DOT(P,ABS(QBAR),7)*SIGN(1,QBAR)
      p23=DOT(P,ABS(Q),ABS(G))*SIGN(1,Q*G)
      p36=DOT(P,ABS(Q),6)*SIGN(1,Q)
      p37=DOT(P,ABS(Q),7)*SIGN(1,Q)
      V1=VDOT(ABS(QBAR))*SIGN(1,QBAR)
      V3=VDOT(ABS(Q))*SIGN(1,Q)
      V6=VDOT(6)
      V7=VDOT(7)

c      DOTS=1/(VV*(Y1+Y2)**2)
      DOTS=1/VV
      T=CF*8*PISQ!/(Y1*Y2)**2

      M= (2*(-2*(p13 - p16 + p17)**2*(p16 + p17)*V3**2 - 
     -      2*(p16 - p17 + p36 - p37)*
     -       ((p13 + p36 - p37)*(p36 + p37)*V1 + 
     -         (p13 - p16 + p17)*(p16 + p17)*V3)*(-V1 + V3 - V6 + V7) + 
     -      2*(p13 - p16 + p17)*(p13 + p36 - p37)*V3*
     -       (-((p16 + p17 + p36 + p37)*V1) + (p13 - p16 + p17)*V6 + 
     -         (p13 - p16 + p17)*V7) - 
     -      (p13 + p36 - p37)*
     -       (2*p13*p36*V1**2 + 2*p36**2*V1**2 + 2*p13*p37*V1**2 - 
     -         2*p37**2*V1**2 - 
     -         2*(p13 - p16 + p17)*(p13 + p36 - p37)*V1*V6 - 
     -         2*(p13 - p16 + p17)*(p13 + p36 - p37)*V1*V7 + 
     -         p13*p16**2*VV - p16**3*VV + p16**2*p17*VV - 
     -         p13*p17**2*VV + p16*p17**2*VV - p17**3*VV + 
     -         2*p13*p17*p36*VV - 2*p16*p17*p36*VV + 
     -         2*p17**2*p36*VV - p13*p36**2*VV + p16*p36**2*VV - 
     -         p17*p36**2*VV - 2*p13*p16*p37*VV + 2*p16**2*p37*VV - 
     -         2*p16*p17*p37*VV + p13*p37**2*VV - p16*p37**2*VV + 
     -         p17*p37**2*VV)))/(p12**2*p23**2*EMSQ)

      m=m*(-DOTS*T)*(4*pi)**2 !!!!OJO CON EL SIGNO!!!
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATFOR(P,M)
      IMPLICIT NONE
C---EVALUATE THE FOUR-PARTON MATRIX ELEMENT SQUARED FOR THE GIVEN
C   CONFIGURATION.
      INTEGER I,J,EWFLAG
      DOUBLE PRECISION P(4,7),M(-6:6),ERTA,ERTB,ERTC,ERTD,ERTE,
     $     BDPQA,BDPQB,BDPQC,BDPGAXP,
     $     BDPQD1,BDPQD2,BDPQE,BDPG,
     $     LEIA,LEIB,LEIC,LEID,LEIE,LEIQAXD,BDPQD2AX,
     $     BDPQFVV,BDPQFAA,BDPQFAV,BDPQFVA,
     $     A,B,C,DS,D1,D2,E,Q,G,QQ,EMSQ,DOT,F2,
     $     FVV,FAA,FAV,FVA,QQFVV,QQFAA,QQFAV,QQFVA,QW,QWZ,E2,EXTFACT

      DOUBLE PRECISION ERTGAXA,ERTGAXB,ERTGAXC,LEIGAXA,LEIGAXB,
     $     LEIGAXC,ERTEW,BDPQEW 
      INTEGER SCHEME,NF,IPOL
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     #                 GAMMAZ,LEPCHARGE,CKM2(3,3)
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL
      EMSQ=-DOT(P,5,5)

C#######################################################################      
      IF (IPOL.EQ.0) THEN        
C#######################################################################
C.......................................................................
C     SOLO FOTON
C.......................................................................
      A=2*(LEIA(P,P(1,6),-1,2,3,4)+LEIA(P,P(1,6),2,-1,3,4)
     $    +LEIA(P,P(1,6),-1,2,4,3)+LEIA(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTA(P,-1,2,3,4)+ERTA(P,2,-1,3,4)
     $    +ERTA(P,-1,2,4,3)+ERTA(P,2,-1,4,3))
      B=2*(LEIB(P,P(1,6),-1,2,3,4)+LEIB(P,P(1,6),2,-1,3,4)
     $    +LEIB(P,P(1,6),-1,2,4,3)+LEIB(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTB(P,-1,2,3,4)+ERTB(P,2,-1,3,4)
     $    +ERTB(P,-1,2,4,3)+ERTB(P,2,-1,4,3))     
      C=2*(LEIC(P,P(1,6),-1,2,3,4)+LEIC(P,P(1,6),2,-1,3,4)
     $    +LEIC(P,P(1,6),-1,2,4,3)+LEIC(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTC(P,-1,2,3,4)+ERTC(P,2,-1,3,4)
     $    +ERTC(P,-1,2,4,3)+ERTC(P,2,-1,4,3))
      Q=HF*(CF*A+(CF-CA/2)*B+CA*C)


      A=2*(LEIA(P,P(1,6),2,3,-1,4)+LEIA(P,P(1,6),3,2,-1,4)
     $    +LEIA(P,P(1,6),2,3,4,-1)+LEIA(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $     ERTA(P,2,3,-1,4)+ERTA(P,3,2,-1,4)
     $    +ERTA(P,2,3,4,-1)+ERTA(P,3,2,4,-1))
      B=2*(LEIB(P,P(1,6),2,3,-1,4)+LEIB(P,P(1,6),3,2,-1,4)
     $    +LEIB(P,P(1,6),2,3,4,-1)+LEIB(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $     ERTB(P,2,3,-1,4)+ERTB(P,3,2,-1,4)
     $    +ERTB(P,2,3,4,-1)+ERTB(P,3,2,4,-1))
      C=2*(LEIC(P,P(1,6),2,3,-1,4)+LEIC(P,P(1,6),3,2,-1,4)
     $    +LEIC(P,P(1,6),2,3,4,-1)+LEIC(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $     ERTC(P,2,3,-1,4)+ERTC(P,3,2,-1,4)
     $    +ERTC(P,2,3,4,-1)+ERTC(P,3,2,4,-1))
      G=-(CF*A+(CF-CA/2)*B+CA*C)


      D1=2*(LEID(P,P(1,6),4,-1,3,2)+LEID(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $      ERTD(P,4,-1,3,2)+ERTD(P,3,2,4,-1))


      Q=Q+NF*TR*D1
      D2=2*(LEID(P,P(1,6),-1,4,2,3)+LEID(P,P(1,6),2,3,-1,4))-EMSQ/2*(
     $      ERTD(P,-1,4,2,3)+ERTD(P,2,3,-1,4))
      QQ=TR*D2

      IF (EWFLAG.EQ.3) THEN
C     W-BOSON CASE:        
        E=-ERTEW(P,4,-1,2,3)-ERTEW(P,2,3,4,-1)-ERTEW(P,3,2,4,-1)
        E2=-ERTEW(P,-1,4,2,3)
        Q=Q+HF*(CF-CA/2)*E
        QW=HF*(CF-CA/2)*E2
      ELSE
        E=2*(LEIE(P,P(1,6),4,-1,3,2)+LEIE(P,P(1,6),3,2,4,-1)
     $      +LEIE(P,P(1,6),4,-1,2,3)+LEIE(P,P(1,6),2,3,4,-1)
     $      +LEIE(P,P(1,6),-1,4,2,3)+LEIE(P,P(1,6),2,3,-1,4)
     $      +LEIE(P,P(1,6),-1,4,3,2)+LEIE(P,P(1,6),3,2,-1,4))-EMSQ/2*(
     $       ERTE(P,4,-1,3,2)+ERTE(P,3,2,4,-1)
     $      +ERTE(P,4,-1,2,3)+ERTE(P,2,3,4,-1)
     $      +ERTE(P,-1,4,2,3)+ERTE(P,2,3,-1,4)
     $      +ERTE(P,-1,4,3,2)+ERTE(P,3,2,-1,4))
        Q=Q+HF*(CF-CA/2)*E
        QW=0D0
      ENDIF

      FVV=-BDPQFVV(P,-1,4,2,3)
      QQFVV=HF*FVV
C.......................................................................
C     Z/W: Agrego las estructuras con acoplamientos mixtos faltantes
C          Para el canal de quarks uso que AVAV(Q)=VVVVpol(Q)
C          Para el canal de gluones, uso que AVAV(G)=-VVVVpol(Q)
C.......................................................................
      IF (EWFLAG.EQ.0) GOTO 99

      A=-(BDPQA(P,-1,2,3,4)-BDPQA(P,2,-1,3,4)
     $    +BDPQA(P,-1,2,4,3)-BDPQA(P,2,-1,4,3))
      B=-(BDPQB(P,-1,2,3,4)-BDPQB(P,2,-1,3,4)
     $    +BDPQB(P,-1,2,4,3)-BDPQB(P,2,-1,4,3))
      C=-(BDPQC(P,-1,2,3,4)-BDPQC(P,2,-1,3,4)
     $    +BDPQC(P,-1,2,4,3)-BDPQC(P,2,-1,4,3))
      
      QZ=HF*(CF*A+(CF-CA/2)*B+CA*C)

      A=-(BDPQA(P,2,3,-1,4)-BDPQA(P,3,2,-1,4)
     $    +BDPQA(P,2,3,4,-1)-BDPQA(P,3,2,4,-1))
      B=-(BDPQB(P,2,3,-1,4)-BDPQB(P,3,2,-1,4)
     $    +BDPQB(P,2,3,4,-1)-BDPQB(P,3,2,4,-1))
      C=-(BDPQC(P,2,3,-1,4)-BDPQC(P,3,2,-1,4)
     $    +BDPQC(P,2,3,4,-1)-BDPQC(P,3,2,4,-1))
      GZ=(CF*A+(CF-CA/2)*B+CA*C)

      D1=BDPQD1(P,-1,4,3,2)
      QZ=QZ+NF*TR*D1

      D2=BDPQD1(P,4,-1,2,3)     
      QQZ=TR*D2

      IF (EWFLAG.EQ.3) THEN
C     W-BOSON CASE:        
        E=-BDPQEW(P,4,-1,2,3)-BDPQEW(P,2,3,4,-1)-BDPQEW(P,3,2,4,-1)
        E2=-BDPQEW(P,-1,4,2,3)
        QZ=QZ+HF*(CF-CA/2)*E
        QWZ=HF*(CF-CA/2)*E2
      ELSE
        E=-(BDPQE(P,4,-1,3,2) + BDPQE(P,4,-1,2,3)
     $      + BDPQE(P,-1,4,2,3) + BDPQE(P,-1,4,3,2))
        QZ=QZ+HF*(CF-CA/2)*E
        QWZ=0D0
      ENDIF

      FAA=-BDPQFAA(P,-1,4,2,3)
      QQFAA=HF*FAA
      FAV=-BDPQFAV(P,-1,4,2,3)
      QQFAV=HF*FAV
      FVA=-BDPQFVA(P,-1,4,2,3)
      QQFVA=HF*FVA
C#######################################################################
      ELSEIF (IPOL.EQ.1) THEN
C#######################################################################                
      A=-(BDPQA(P,-1,2,3,4)-BDPQA(P,2,-1,3,4)
     $    +BDPQA(P,-1,2,4,3)-BDPQA(P,2,-1,4,3))
      B=-(BDPQB(P,-1,2,3,4)-BDPQB(P,2,-1,3,4)
     $    +BDPQB(P,-1,2,4,3)-BDPQB(P,2,-1,4,3))
      C=-(BDPQC(P,-1,2,3,4)-BDPQC(P,2,-1,3,4)
     $    +BDPQC(P,-1,2,4,3)-BDPQC(P,2,-1,4,3))

      Q=HF*(CF*A+(CF-CA/2)*B+CA*C)

      G=(BDPG(P,1,4,2,3)+BDPG(P,1,4,3,2))/16D0

      D1=BDPQD1(P,-1,4,3,2)      
      Q=Q+NF*TR*D1

      D2=BDPQD2(P,-1,4,2,3)   
      QQ=TR*D2

      IF (EWFLAG.EQ.3) THEN
C     W-BOSON CASE:        
        E=-BDPQEW(P,4,-1,2,3)-BDPQEW(P,2,3,4,-1)-BDPQEW(P,3,2,4,-1)
        E2=-BDPQEW(P,-1,4,2,3)
        Q=Q+HF*(CF-CA/2)*E
        QW=HF*(CF-CA/2)*E2
      ELSE
        E=-(BDPQE(P,4,-1,3,2) + BDPQE(P,4,-1,2,3)
     $      + BDPQE(P,-1,4,2,3) + BDPQE(P,-1,4,3,2))

        Q=Q+HF*(CF-CA/2)*E
        QW=0D0
      ENDIF

      FVV=-(BDPQFAV(P,-1,4,2,3))
      QQFVV=HF*FVV
C.......................................................................
C     Z: Agrego las estructuras con acoplamientos mixtos faltantes
C        Para el canal de quarks uso que AVAVpol(Q)=VVVV(Q)
C        Para el canal de gluones,  AVAV(G)pol--> funciones ERTGAX
C.......................................................................
      IF (EWFLAG.EQ.0) GOTO 99

      A=2*(LEIA(P,P(1,6),-1,2,3,4)+LEIA(P,P(1,6),2,-1,3,4)
     $    +LEIA(P,P(1,6),-1,2,4,3)+LEIA(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTA(P,-1,2,3,4)+ERTA(P,2,-1,3,4)
     $    +ERTA(P,-1,2,4,3)+ERTA(P,2,-1,4,3))
      B=2*(LEIB(P,P(1,6),-1,2,3,4)+LEIB(P,P(1,6),2,-1,3,4)
     $    +LEIB(P,P(1,6),-1,2,4,3)+LEIB(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTB(P,-1,2,3,4)+ERTB(P,2,-1,3,4)
     $    +ERTB(P,-1,2,4,3)+ERTB(P,2,-1,4,3))    
      C=2*(LEIC(P,P(1,6),-1,2,3,4)+LEIC(P,P(1,6),2,-1,3,4)
     $    +LEIC(P,P(1,6),-1,2,4,3)+LEIC(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTC(P,-1,2,3,4)+ERTC(P,2,-1,3,4)
     $    +ERTC(P,-1,2,4,3)+ERTC(P,2,-1,4,3))
    
      QZ=HF*(CF*A+(CF-CA/2)*B+CA*C)

c       A=(-ERTGAXA(P,2,3,-1,4)+ERTGAXA(P,3,2,-1,4)+
c     $   -LEIGAXA(P,2,3,-1,4)+LEIGAXA(P,3,2,-1,4))
c       B=(-ERTGAXB(P,2,3,-1,4)+ERTGAXB(P,3,2,-1,4)+
c     $   -LEIGAXB(P,2,3,-1,4)+LEIGAXB(P,3,2,-1,4))
c       C=(-ERTGAXC(P,2,3,-1,4)+ERTGAXC(P,3,2,-1,4)+
c     $   -LEIGAXC(P,2,3,-1,4)+LEIGAXC(P,3,2,-1,4))
c      GZ=(CF*A+(CF-CA/2)*B+CA*C)

      GZ=(BDPGAXP(P,1,4,2,3)-BDPGAXP(P,1,4,3,2))/16D0

      D1=2*(LEID(P,P(1,6),4,-1,3,2)+LEID(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $      ERTD(P,4,-1,3,2)+ERTD(P,3,2,4,-1))
      QZ=QZ+NF*TR*D1


      D2=-BDPQD2AX(P,-1,4,2,3)      
      QQZ=TR*D2

      IF (EWFLAG.EQ.3) THEN
C     W-BOSON CASE:        
        E=-ERTEW(P,4,-1,2,3)-ERTEW(P,2,3,4,-1)-ERTEW(P,3,2,4,-1)
        E2=-ERTEW(P,-1,4,2,3)
        QZ=QZ+HF*(CF-CA/2)*E
        QWZ=HF*(CF-CA/2)*E2
      ELSE
        E=2*(LEIE(P,P(1,6),4,-1,3,2)+LEIE(P,P(1,6),3,2,4,-1)
     $      +LEIE(P,P(1,6),4,-1,2,3)+LEIE(P,P(1,6),2,3,4,-1)
     $      +LEIE(P,P(1,6),-1,4,2,3)+LEIE(P,P(1,6),2,3,-1,4)
     $      +LEIE(P,P(1,6),-1,4,3,2)+LEIE(P,P(1,6),3,2,-1,4))-EMSQ/2*(
     $       ERTE(P,4,-1,3,2)+ERTE(P,3,2,4,-1)
     $      +ERTE(P,4,-1,2,3)+ERTE(P,2,3,4,-1)
     $      +ERTE(P,-1,4,2,3)+ERTE(P,2,3,-1,4)
     $      +ERTE(P,-1,4,3,2)+ERTE(P,3,2,-1,4))
        QZ=QZ+HF*(CF-CA/2)*E
        QWZ=0D0
      ENDIF

      FAA=-BDPQFVA(P,-1,4,2,3)
      QQFAA=HF*FAA
      FAV=-BDPQFVV(P,-1,4,2,3)
      QQFAV=HF*FAV
      FVA=-BDPQFAA(P,-1,4,2,3)
      QQFVA=HF*FVA      


C#######################################################################      
      ENDIF
 99   CONTINUE
C.......................................................................      
C---INCLUDE EXTERNAL FACTORS
      EXTFACT = 256*PI**4/EMSQ *(4*PI)**2*4/EMSQ

      Q=Q*CF*EXTFACT
      G=G*TR*EXTFACT
      QQ=QQ*CF*EXTFACT
      QW=QW*CF*EXTFACT

      QZ=QZ*CF*EXTFACT
      GZ=GZ*TR*EXTFACT
      QQZ=QQZ*CF*EXTFACT
      QWZ=QWZ*CF*EXTFACT

      QQFVV=QQFVV*CF*EXTFACT
      QQFAA=QQFAA*CF*EXTFACT
      QQFAV=QQFAV*CF*EXTFACT
      QQFVA=QQFVA*CF*EXTFACT                        

      DO I=-6,6
        M(I)=0
      ENDDO

C     COMBINE THE CROSS SECTION WITH COUPLINGS
      CALL COMBINEQI(P,M,Q,QZ)
      CALL COMBINEGLUON(P,M,G,GZ)
      CALL COMBINEQF(P,M,QQ,QQZ)
      
      IF (EWFLAG.LE.2) THEN
C     ADD ERTF DIAGRAMS ONLY FOR Z (AND TECHNICALLY PHOTON)
        CALL COMBINEQIF(P,M,QQFVV,0,0)
        CALL COMBINEQIF(P,M,QQFAA,1,1)
        CALL COMBINEQIF(P,M,QQFVA,0,1)
        CALL COMBINEQIF(P,M,QQFAV,1,0)
      ELSE
C     ADD THE SPECIAL COMBINATION OF IDENTICAL PARTICLE TERMS FOR W
        CALL COMBINEQIW(P,M,QW,QWZ)
      ENDIF

      END
C-----------------------------------------------------------------------
      SUBROUTINE VIRTWO(S,P,V,*)
      IMPLICIT NONE
C---CALCULATE THE TWO-PARTON MATRIX-ELEMENT AT NEXT-TO-LEADING ORDER
      INTEGER I,IPOL
      DOUBLE PRECISION S,P(4,7),V(-6:6),M(-6:6),O,X,XJAC,XMIN,
     $     QQ,GQ,QG,GG,KQF,DOT
      PARAMETER (O=0)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTWO(P,M)
C---SUM OF FACTORIZING VIRTUAL CROSS-SECTION AND SUBTRACTION COUNTERTERM
      QQ=CF*(2-PISQ)
      GG=0
C---THE NON-FACTORIZING VIRTUAL CROSS-SECTION
      DO I=-6,6
        V(I)=0
      ENDDO
C---GENERATE A COLLINEAR EMISSION
      XMIN=2*DOT(P,1,6)/S
      CALL GENCOL(1,X,XJAC,XMIN)
C---ENFORCE INVARIANT MASS CUTOFF
      IF (1-X.LT.CUTOFF) RETURN 1
C---CALCULATE THE COLLINEAR COUNTERTERM
      GQ=0
      QG=0
      KQF=1.5
      CALL KPFUNS(-X,XJAC,XMIN,KQF,O,O,O,QQ,GQ,QG,GG)
C---THE TOTAL
      V(0)=V(0)+GG*M(0)
      DO I=-6,6
        IF (I.NE.0) THEN
          V(I)=V(I)+QG*M(0)+QQ*M(I)
          IF (ABS(I).LE.NF) V(0)=V(0)+GQ*M(I)
        ENDIF
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE VIRTHR(S,P,V,*)
      IMPLICIT NONE
C---CALCULATE THE THREE-PARTON MATRIX-ELEMENT AT NEXT-TO-LEADING ORDER
      INTEGER I,EWFLAG
      DOUBLE PRECISION S,P(4,7),V(-6:6),M(-6:6),X,XJAC,XMIN,
     $     QQ,GQ,QG,GG,KQF,KGF,PQF,PGF,L12,L13,L23,DOT,ERTV,LEIV,EMSQ,
     $     BDPV,QQQ,GGG,BDPVGAX,BDPVTRI,QQQZ,GGGZ
      INTEGER SCHEME,NF,IPOL
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QQZ,GGZ,MZ,
     #                 GAMMAZ,LEPCHARGE,CKM2(3,3)
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL

      EMSQ=-DOT(P,5,5)
      L12=LOG(2*DOT(P,1,2)/EMSQ)
      L13=LOG(2*DOT(P,1,3)/EMSQ)
      L23=LOG(2*DOT(P,2,3)/EMSQ)
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTHR(P,M)
C---THE NON-FACTORIZING VIRTUAL CROSS-SECTION
C--------------------------------------------
      IF (IPOL.EQ.0) THEN
      QQ=-((4*PI)**2*4/EMSQ)*
     $     (2*LEIV(P,P(1,6),2,-1,3)-EMSQ/2*ERTV(P,2,-1,3))
        IF (EWFLAG.GE.1) THEN
          QQZ = -((4*PI)**2*4/EMSQ)*(BDPV(P,P(1,6),1,3,2,0))
        ENDIF

      GG=TR/CF*((4*PI)**2*4/EMSQ)*
     $     (2*LEIV(P,P(1,6),2,3,-1)-EMSQ/2*ERTV(P,2,3,-1))
        IF (EWFLAG.GE.1) THEN
          GGZ = -TR/CF*((4*PI)**2*4/EMSQ)*(BDPV(P,P(1,6),-2,3,-1,0))
        ENDIF
C--------------------------------------------     
      ELSEIF (IPOL.EQ.1) THEN
c     (minus compensates the one in EMSQ)
      QQ=-((4*PI)**2*4/EMSQ)*(BDPV(P,P(1,6),1,3,2,0))
        IF (EWFLAG.GE.1) THEN
          QQZ = -((4*PI)**2*4/EMSQ)*
     $         (2*LEIV(P,P(1,6),2,-1,3)-EMSQ/2*ERTV(P,2,-1,3))
        ENDIF
      
      GG=-((4*PI)**2*4/EMSQ)*(BDPV(P,P(1,6),1,3,2,1))
        IF (EWFLAG.GE.1) THEN
          GGZ = ((4*PI)**2*4/EMSQ)*(BDPVGAX(P,P(1,6),1,3,2)
     $                             -BDPVGAX(P,P(1,6),1,2,3))
        ENDIF
      ENDIF
C--------------------------------------------
c     TRIANGLE CONTRIBUTIONS (ONLY FOR Z) 
      IF (EWFLAG.GE.1) THEN
        QQQ=-((4*PI)**2*4/EMSQ**2)*BDPVTRI(P,P(1,6),1,3,2,0,0)
        QQQZ=-((4*PI)**2*4/EMSQ**2)*BDPVTRI(P,P(1,6),1,3,2,1,0)
        GGG=-((4*PI)**2*4/EMSQ**2)*BDPVTRI(P,P(1,6),1,3,2,1,1)
        GGGZ=-((4*PI)**2*4/EMSQ**2)*BDPVTRI(P,P(1,6),1,3,2,1,1)
      ENDIF
C--------------------------------------------
C      QQQ=0.d0
C      QQQZ=0.d0
C      GGG=0.d0
C      GGGZ=0.d0

      DO I=-6,6
        V(I)=0
      ENDDO

C     COMBINE THE CROSS SECTION WITH COUPLINGS
      CALL COMBINEQI(P,V,QQ,QQZ)
      CALL COMBINEGLUON(P,V,GG,GGZ)

C     ADD TRAINGLE DIAGRAMS ONLY FOR Z     
      IF ((EWFLAG.EQ.1).OR.(EWFLAG.EQ.2)) THEN
        CALL COMBINEQITRIANG(P,V,QQQ,QQQZ)
        CALL COMBINEGTRIANG(P,V,GGG,GGGZ)
      ENDIF

C---SUM OF FACTORIZING VIRTUAL CROSS-SECTION AND SUBTRACTION COUNTERTERM
      QQ=CF*2+CA*50D0/9-TR*NF*16D0/9-CF*PISQ
     $     -3*(CF-CA/2)*L12-(5*CA-TR*NF)/3*(L13+L23)
      GG=CF*2+CA*50D0/9-TR*NF*16D0/9-CA*PISQ
     $     -3*(CF-CA/2)*L23-(5*CA-TR*NF)/3*(L12+L13)
C---GENERATE A COLLINEAR EMISSION
      XMIN=2*DOT(P,1,6)/S
      CALL GENCOL(2,X,XJAC,XMIN)
C---ENFORCE INVARIANT MASS CUTOFF
      IF (1-X.LT.CUTOFF) RETURN 1
C---CALCULATE THE COLLINEAR COUNTERTERM
      GQ=0
      QG=0
      KQF=(1.5*(CF-CA/2)+0.5*(11D0/6*CA-2D0/3*NF*TR))/CF
      KGF=1.5
      PQF=-((CF-CA/2)*L12+CA/2*L13)/CF
      PGF=-(L12+L13)/2
      CALL KPFUNS(-X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
C---THE TOTAL
      V(0)=V(0)+GG*M(0)
      DO I=-6,6
        IF (I.NE.0) THEN
          V(I)=V(I)+QG*M(0)+QQ*M(I)
          IF (ABS(I).LE.NF) V(0)=V(0)+GQ*M(I)
        ENDIF
      ENDDO

      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ERTV(P,I,J,K)
      IMPLICIT NONE
C---RETURN THE ERT F FUNCTION FOR VIRTUAL TERMS THAT ARE NOT
C   TRIVIALLY PROPORTIONAL TO TREE-LEVEL
      INTEGER I,J,K
      DOUBLE PRECISION P(4,7),DOT,R,RR,X,Y,DILOG,EMSQ,
     $     Y12,Y13,Y23,L12,L13,L23,R1312,R2312,R2313
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON  /SUBCOM/ CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      RR(X,Y)=DILOG(X)+DILOG(Y)-DILOG(X*Y)-PISQ/6
      R(X,Y)=RR((1-X)/Y,(1-Y)/X)
      EMSQ=DOT(P,5,5)
      Y12=2*DOT(P,ABS(I),ABS(J))/EMSQ*SIGN(1,I*J)
      Y13=2*DOT(P,ABS(I),ABS(K))/EMSQ*SIGN(1,I*K)
      Y23=2*DOT(P,ABS(J),ABS(K))/EMSQ*SIGN(1,J*K)
      L12=LOG(ABS(Y12))
      L13=LOG(ABS(Y13))
      L23=LOG(ABS(Y23))
      R1312=R(Y13,Y12)
      R2312=R(Y23,Y12)
      R2313=R(Y23,Y13)
      CFSUB=Y12/(Y12+Y13)+Y12/(Y12+Y23)+(Y12+Y23)/Y13+(Y12+Y13)/Y23
     $     +L13*(4*Y12**2+2*Y12*Y13+4*Y12*Y23+Y13*Y23)/(Y12+Y23)**2
     $     +L23*(4*Y12**2+2*Y12*Y23+4*Y12*Y13+Y13*Y23)/(Y12+Y13)**2
     $     -2*((Y12**2+(Y12+Y13)**2)/(Y13*Y23)*R2312
     $        +(Y12**2+(Y12+Y23)**2)/(Y13*Y23)*R1312
     $        +(Y13**2+Y23**2)/(Y13*Y23*(Y13+Y23))
     $     -2*L12*(Y12**2/(Y13+Y23)**2+2*Y12/(Y13+Y23)))
      CASUB=L13*Y13/(Y12+Y23)+L23*Y23/(Y12+Y13)
     $     +((Y12**2+(Y12+Y13)**2)/(Y13*Y23)*R2312
     $     +(Y12**2+(Y12+Y23)**2)/(Y13*Y23)*R1312
     $     +(Y13**2+Y23**2)/(Y13*Y23*(Y13+Y23))
     $     -2*L12*(Y12**2/(Y13+Y23)**2+2*Y12/(Y13+Y23)))
     $     -R2313*(Y13/Y23+Y23/Y13+2*Y12/(Y13*Y23))
      TFSUB=0
      ERTV=CF*(CF*CFSUB+CA*CASUB)
      CFSUB=CFSUB/ERTV
      CASUB=CASUB/ERTV
      ERTV=16*PISQ/EMSQ*ERTV
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION LEIV(P,Q,I,J,K)
      IMPLICIT NONE
C---RETURN THE PART OF THE ONE-LOOP HADRONIC TENSOR NOT TRIVIALLY
C   PROPORTIONAL TO TREE-LEVEL ONE CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
      INTEGER I,J,K
      DOUBLE PRECISION P(4,7),Q(4),DOT,R,RR,X,Y,Z,DILOG,EMSQ,
     $     Y12,Y13,Y23,L12,L13,L23,A,B,C,D,P1Q,P2Q,P3Q,QQ,
     $     CF1,CF2,CF3,CF4,CA1,CA2,CA3,CA4,R1312,R2312,R2313
      PARAMETER (Z=0)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON  /SUBCOM/ CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      RR(X,Y)=DILOG(X)+DILOG(Y)-DILOG(X*Y)-PISQ/6
      R(X,Y)=RR((1-X)/Y,(1-Y)/X)
      EMSQ=DOT(P,5,5)
      Y12=2*DOT(P,ABS(I),ABS(J))/EMSQ*SIGN(1,I*J)
      Y13=2*DOT(P,ABS(I),ABS(K))/EMSQ*SIGN(1,I*K)
      Y23=2*DOT(P,ABS(J),ABS(K))/EMSQ*SIGN(1,J*K)
      L12=LOG(ABS(Y12))
      L13=LOG(ABS(Y13))
      L23=LOG(ABS(Y23))
C---FIRST DECOMPOSE Q AS A*P1+B*P2+C*P3+D*E_(MU,P1,P2,P3)/EMSQ
      P1Q=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)-P(2,ABS(I))*Q(2)
     $     -P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2Q=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)-P(2,ABS(J))*Q(2)
     $     -P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3Q=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)-P(2,ABS(K))*Q(2)
     $     -P(1,ABS(K))*Q(1))*SIGN(1,K)
      QQ=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      A=(P3Q*Y12+P2Q*Y13-P1Q*Y23)/(Y12*Y13*EMSQ)
      B=(P1Q*Y23+P3Q*Y12-P2Q*Y13)/(Y12*Y23*EMSQ)
      C=(P2Q*Y13+P1Q*Y23-P3Q*Y12)/(Y23*Y13*EMSQ)
      D=SQRT(MAX(Z,4*(QQ/EMSQ-A*B*Y12-A*C*Y13-B*C*Y23)/(-Y12*Y13*Y23)))
C---THEN CALCULATE THE CONTRACTIONS WITH -P1(MU)P1(NU) ETC
      R1312=R(Y13,Y12)
      R2312=R(Y23,Y12)
      R2313=R(Y23,Y13)
      CF1=R1312*Y12-R2312*Y12**2/Y13
     $     -L12*Y12*Y13/(Y13+Y23)**2-L23*Y12+Y12/2*(Y23-Y13)/(Y23+Y13)
      CA1=-R1312*Y12/2+R2312*Y12**2/Y13/2+R2313*Y12/2
     $     +L12*Y12*Y13/(Y13+Y23)**2/2-Y12*Y23/(Y23+Y13)/2
      CF2=R2312*Y12-R1312*Y12**2/Y23
     $     -L12*Y12*Y23/(Y13+Y23)**2-L13*Y12+Y12/2*(Y13-Y23)/(Y13+Y23)
      CA2=-R2312*Y12/2+R1312*Y12**2/Y23/2+R2313*Y12/2
     $     +L12*Y12*Y23/(Y13+Y23)**2/2-Y12*Y13/(Y13+Y23)/2
      CF3=R1312*Y12+R2312*Y12-2*L12*Y12
     $     -L13*Y12*Y23/2/(Y12+Y23)**2*(1+Y12+Y23)
     $     -L23*Y12*Y13/2/(Y12+Y13)**2*(1+Y12+Y13)
     $     -Y12/2*(Y13/(Y12+Y13)+Y23/(Y12+Y23))
      CA3=-R1312*Y12/2-R2312*Y12/2+R2313*Y12
     $     +L12*Y12+L13*Y12*Y13/(Y12+Y23)/2+L23*Y12*Y23/(Y12+Y13)/2
      CF4=R1312*Y12/4/Y23*(Y12-Y12**2-Y12*Y23+2*Y12**2*Y23
     $                         +2*Y12*Y23**2+Y23**3)
     $     +R2312*Y12/4/Y13*(Y12-Y12**2-Y12*Y13+2*Y12**2*Y13
     $                           +2*Y12*Y13**2+Y13**3)
     $     +L12*Y12**2/4/(Y13+Y23)*(Y13+Y23-2*Y13*Y23)
     $     -L13*Y12/8/(Y12+Y23)*Y13*(Y12+Y13)*(Y23-2*Y12)
     $     -L23*Y12/8/(Y12+Y13)*Y23*(Y12+Y23)*(Y13-2*Y12)
     $     +Y12/8*(Y13+Y23-2*Y13*Y23)
      CA4=-R1312*Y12/8/Y23*(Y12-Y12**2-Y12*Y23+2*Y12**2*Y23
     $                          +2*Y12*Y23**2+Y23**3)
     $     -R2312*Y12/8/Y13*(Y12-Y12**2-Y12*Y13+2*Y12**2*Y13
     $                           +2*Y12*Y13**2+Y13**3)
     $     +R2313*Y12/8*((1-Y13)**2+(1-Y23)**2)
     $     -L12*Y12**2/8/(Y13+Y23)*(Y13+Y23-2*Y13*Y23)
     $     -L13*Y12/8*Y13*(Y12+Y13)-L23*Y12/8*Y23*(Y12+Y23)
     $     -Y12/8*(Y13+Y23-2*Y13*Y23)
C---AND COMBINE THEM WITH THE APPROPRIATE WEIGHTS
      CFSUB=(A-B)*(A-C)*CF1+(B-C)*(B-A)*CF2+(C-A)*(C-B)*CF3+D*D*CF4
      CASUB=(A-B)*(A-C)*CA1+(B-C)*(B-A)*CA2+(C-A)*(C-B)*CA3+D*D*CA4
      TFSUB=0
      LEIV=CF*(CF*CFSUB+CA*CASUB)
      CFSUB=CFSUB/LEIV
      CASUB=CASUB/LEIV
      LEIV=16*PISQ*LEIV
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BDPV(P,Q,I,J,K,CHAN)
      IMPLICIT NONE
C---RETURN THE PART OF THE ONE-LOOP HADRONIC TENSOR NOT TRIVIALLY
C   PROPORTIONAL TO TREE-LEVEL ONE CONTRACTED WITH THE POLARIZED PART
C   OF THE LEPTONIC TENSOR.
      INTEGER I,J,K,CHAN
      DOUBLE PRECISION P(4,7),Q(4),DOT,R,RR,X,Y,Z,DILOG,EMSQ,
     $     Y12,Y13,Y23,L12,L13,L23,A,B,C,D,P1Q,P2Q,P3Q,QQ,
     $     PC1,PB1,R1312,R2312,R2313,R1213,R1223,R1323
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON  /SUBCOM/ CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB

      RR(X,Y)=DILOG(X)+DILOG(Y)-DILOG(X*Y)-PISQ/6
      R(X,Y)=RR((1-X)/Y,(1-Y)/X)
      EMSQ=DOT(P,5,5)
      Y12=2*DOT(P,ABS(I),ABS(J))/EMSQ*SIGN(1,I*J)
      Y13=2*DOT(P,ABS(I),ABS(K))/EMSQ*SIGN(1,I*K)
      Y23=2*DOT(P,ABS(J),ABS(K))/EMSQ*SIGN(1,J*K)
      L12=LOG(ABS(Y12))
      L13=LOG(ABS(Y13))
      L23=LOG(ABS(Y23))

C---WE DECOMPOSE LEPTON Q AS A*P1+B*P2+C*P3+D*E_(MU,P1,P2,P3)/EMSQ
      P1Q=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)-P(2,ABS(I))*Q(2)
     $     -P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2Q=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)-P(2,ABS(J))*Q(2)
     $     -P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3Q=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)-P(2,ABS(K))*Q(2)
     $     -P(1,ABS(K))*Q(1))*SIGN(1,K)
      QQ=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)

      A=(P3Q*Y12+P2Q*Y13-P1Q*Y23)/(Y12*Y13*EMSQ)
      B=(P1Q*Y23+P3Q*Y12-P2Q*Y13)/(Y12*Y23*EMSQ)
      C=(P2Q*Y13+P1Q*Y23-P3Q*Y12)/(Y23*Y13*EMSQ)
C     D=SQRT(MAX(Z,4*(QQ/EMSQ-A*B*Y12-A*C*Y13-B*C*Y23)/(-Y12*Y13*Y23)))

C---THEN CALCULATE THE CONTRACTION WITH THE LEVI-CIVITA TENSOR
      R1213=R(-Y13,-Y12)
      R1223=R(Y23,-Y12)
      R1323=R(Y23,-Y13)
      
C---FIRST WE DO THE QUARK CHANNEL
      IF (CHAN.EQ.0) THEN

      PB1=(CF*L12*Y12*Y13*(Y12+Y13))/(Y13-Y23)**2
     $    + (L12*(Y12+Y13)*(CA*Y12+CF*(Y12+4*Y13)))/(Y13-Y23)
     $    - (2*(CA-2*CF)*L13*Y13*(Y12+Y13-Y23)**2)/(Y12-Y23)**2
     $    + (L23*(CF*(Y12-2*Y13)+CA*(Y12+Y13))*Y23)/(Y12+Y13)
     $    - (2*(CA-2*CF)*Y13**2*Y23)/Y12/(Y12-Y23) 
     $    + (CF*Y23*(Y12+Y23))/(Y13-Y23)
     $    + (CF*(Y12**3+Y12**2*(2*Y13+Y23)+Y13*Y23*(3*Y13+Y23)
     $      +Y12*Y13*(Y13+3*Y23)))/Y12/Y23
     $    - (CA*(Y12**3+2*Y12**2*Y13+Y13*Y23*(Y13+Y23) 
     $      +Y12*(Y13**2+Y13*Y23-Y23**2)))/Y12/Y23
     $    - ((CA-2*CF)*R1213*(2*Y13**2*(Y12+Y13)-2*Y13**2*Y23 
     $      -(Y12-Y13)*Y23**2))/Y12/Y23
     $    + (CA*R1223*(Y12**3+3*Y12**2*Y13+Y12*(4*Y13**2-Y23**2)
     $      +Y13*(2*Y13**2-2*Y13*Y23+Y23**2)))/Y12/Y23
     $    - ((CA-2*CF)*R1323*(Y12**3+3*Y12**2*Y13+4*Y12*Y13**2
     $      +2*Y13**2*(Y13-Y23)))/Y12/Y23
      PC1=(L12*(CA*Y12*(Y12+Y13) + CF*(Y12**2+2*Y12*Y13+4*Y13**2)))
     $      /(Y13-Y23) + (CF*L12*Y12*Y13*(Y12+Y13))/(Y13-Y23)**2
     $    + (CA+CF)*(L12*Y12+L23*Y23)
     $    + (CF*L23*Y13*(Y13-Y23)*Y23)/(Y12+Y13)**2
     $    - (L23*(CA*Y23*(Y23-Y13)+CF*(4*Y13**2-2*Y13*Y23+Y23**2)))
     $      /(Y12+Y13) + (CF*Y13**2*(Y13-Y23))/Y12/(Y12+Y13)
     $    + (CF*Y13**2*(Y12+Y13))/((Y13-Y23)*Y23) - ((Y12+Y13-Y23)
     $      *(Y12+Y23)*(CA*(Y12-Y23)+CF*(Y13-Y12+Y23)))/(Y12*Y23)
     $    - ((CA-2*CF)*R1213*(2*Y13**2-2*Y13*Y23-(Y12-Y23)*Y23))/Y12
     $    + (CA*R1223*(Y12+Y23)*(Y12**2+2*Y12*Y13+2*Y13**2
     $      -2*(Y12+Y13)*Y23+Y23**2))/(Y12*Y23)
     $    - ((CA-2*CF)*R1323*(Y12**2+2*Y12*Y13+2*Y13**2-Y12*Y23))/Y23
      
C---AND WE COMBINE IT WITH THE APPROPRIATE WEIGHTS
      BDPV=CF*((-A-C)*PB1+(C-B)*PC1)

C---NOW WE DO THE GLUON CHANNEL
      ELSEIF (CHAN.EQ.1) THEN
      PB1=(2*(CA-2*CF)*L23*Y23**2)/(Y12+Y13) 
     $    + CF*L12*Y12*(Y12-Y23)*Y23/(Y13-Y23)**2
     $    - L12*(CA*Y12*(Y12-Y23)+CF*((Y12-Y23)**2+3*Y23**2))/(Y13-Y23)
     $    - (CA+CF)*L12*Y12 + (CF*L13*Y13*(Y13-Y23)*Y23)/(Y12-Y23)**2
     $    - L13*(CA*Y13*(Y13-Y23)+CF*((Y13-Y23)**2+3*Y23**2))/(Y12-Y23)
     $    - (CA+CF)*L13*Y13 + (CF*(Y12-Y23)*Y23**2)/Y13/(Y13-Y23)
     $    + (CF*(Y13-Y23)*Y23**2)/Y12/(Y12-Y23)
     $    - (CA*(Y12**2+Y13**2)*(Y12+Y13-Y23)
     $      -CF*(Y12+Y13)*(Y12**2+Y13**2-Y23**2))/Y12/Y13
     $    + CA*R1213*(Y12**2+Y13**2)*(Y12+Y13-2*Y23)/Y12/Y13
     $    - (CA-2*CF)*R1223*(Y13**2*(Y12+Y13-2*Y23)-2*Y12*Y23**2)
     $      /Y12/Y13
     $    - (CA-2*CF)*R1323*(Y12**2*(Y12+Y13-2*Y23)-2*Y13*Y23**2)
     $      /Y12/Y13
      PC1=2*(CA-2*CF)*L23*(Y13**2-Y12*(Y12-Y23))*Y23/(Y12+Y13)**2
     $    - L13*Y13*(CA*(Y12-Y23)+CF*(Y12+2*Y23))/(Y12-Y23)
     $    - L12*(Y12-Y23)*((CA+CF)*Y12*Y13
     $      -(CA*Y12+2*CF*(Y12+2*Y13))*Y23+4*CF*Y23**2)/(Y13-Y23)**2
     $    + CF*(Y12-Y13)*Y13/(Y13-Y23)
     $    - 2*(CA-2*CF)*Y13**2*Y23/Y12/(Y12+Y13)
     $    - (CA*(Y12-Y23)*(Y12**2-Y12*Y23+Y13*(Y13+Y23)) 
     $      -CF*(Y12**3-Y13*Y23*(3*Y13+Y23)-Y12**2*(Y13+2*Y23)
     $       +Y12*(2*Y13**2+3*Y13*Y23+Y23**2)))/Y12/Y13
     $    + CA*R1213*(Y12**3-3*Y12**2*Y23+Y13*(Y13-2*Y23)*Y23
     $      +Y12*(Y13**2+2*Y23**2))/Y12/Y13
     $    - (CA-2*CF)*R1223*(Y12*Y13+(Y13-2*Y23)*Y23)/Y12
     $    - (CA-2*CF)*R1323*(Y12**3-3*Y12**2*Y23+2*(Y12-Y13)*Y23**2)
     $      /(Y12*Y13)
      
C---AND WE COMBINE IT WITH THE APPROPRIATE WEIGHTS
      BDPV=TR*((-A-C)*PB1+(C-B)*PC1)
      ENDIF

C---  Hay un fator  1/2 adicional porque cómo comvertimos p12->Y12
      BDPV=16*PISQ*BDPV/2d0
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BDPVGAX(P,Q,I,J,K)
      IMPLICIT NONE      
C---RETURN THE PART OF THE ONE-LOOP HADRONIC TENSOR NOT TRIVIALLY
C   PROPORTIONAL TO TREE-LEVEL ONE CONTRACTED WITH THE POLARIZED PART
C   OF THE LEPTONIC TENSOR.

C   DOS LLAMADAS SON NECESARIAS PARA EL EM COMPLETO (HACER CROSSING 2<->3
C   CON UN SIGNO - EN LA SEGUNDA LLAMADA)    
      INTEGER I,J,K
      DOUBLE PRECISION P(4,7),Q(4),DOT,R,RR,X,Y,Z,DILOG,EMSQ,
     $     Y12,Y13,Y23,L12,L13,L23,A,B,C,D,P1Q,P2Q,P3Q,QQ,
     $     R1312,R2312,R2313,R1213,R1223,R1323,
     $     PA2,PBC2,PABC,PBC,PERT   
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON  /SUBCOM/ CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB

      RR(X,Y)=DILOG(X)+DILOG(Y)-DILOG(X*Y)-PISQ/6
      R(X,Y)=RR((1-X)/Y,(1-Y)/X)
      EMSQ=DOT(P,5,5)
      Y12=2*DOT(P,ABS(I),ABS(J))/EMSQ*SIGN(1,I*J)
      Y13=2*DOT(P,ABS(I),ABS(K))/EMSQ*SIGN(1,I*K)
      Y23=2*DOT(P,ABS(J),ABS(K))/EMSQ*SIGN(1,J*K)
      L12=LOG(ABS(Y12))
      L13=LOG(ABS(Y13))
      L23=LOG(ABS(Y23))

C---WE DECOMPOSE LEPTON Q AS A*P1+B*P2+C*P3+D*E_(MU,P1,P2,P3)/EMSQ
      P1Q=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)-P(2,ABS(I))*Q(2)
     $     -P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2Q=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)-P(2,ABS(J))*Q(2)
     $     -P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3Q=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)-P(2,ABS(K))*Q(2)
     $     -P(1,ABS(K))*Q(1))*SIGN(1,K)
      QQ=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)

      A=(P3Q*Y12+P2Q*Y13-P1Q*Y23)/(Y12*Y13*EMSQ)
      B=(P1Q*Y23+P3Q*Y12-P2Q*Y13)/(Y12*Y23*EMSQ)
      C=(P2Q*Y13+P1Q*Y23-P3Q*Y12)/(Y23*Y13*EMSQ)
C     D=SQRT(MAX(Z,4*(QQ/EMSQ-A*B*Y12-A*C*Y13-B*C*Y23)/(-Y12*Y13*Y23)))

C---THEN CALCULATE THE CONTRACTION WITH THE LEVI-CIVITA TENSOR
      R1213=R(-Y13,-Y12)
      R1223=R(Y23,-Y12)
      R1323=R(Y23,-Y13)

C---  PARTE A^2
      PA2=  (CA - 2*CF)*R1223*Y23 + 
     -  (CF*Y12*Y23**2)/((Y12 - Y23)*(-Y13 + Y23)) - 
     -  (L12*Y23*(CF*Y13*(Y12 + 2*Y13 - 2*Y23) + CA*Y12*(-Y13 + Y23)))/
     -   (Y13 - Y23)**2


C---  PARTE (B^2+C^2)
      PBC2= 2*CF*L13*Y23 - CA*R1213*Y23 + 
     -  ((CF*(Y12 - Y13) + CA*Y13)*Y23)/(Y12 + Y13) - 
     -  ((CA - 2*CF)*L23*Y12*(Y12 + Y13 - Y23)*Y23)/(Y12 + Y13)**2 + 
     -  ((CA - 2*CF)*Y23*(R1223*Y12 + R1323*Y23))/Y12

C---  PARTE A*(B+C)
      PABC= ((CA - 2*CF)*L23*(Y12 - 2*Y13)*(Y12 + Y13 - Y23)*Y23)/
     -   (Y13*(Y12 + Y13)) + 
     -  (CA*R1213*(Y12**2 - Y13**2 - 2*Y12*Y23))/Y13 - 
     -  ((CA - 2*CF)*R1323*(Y12**3 - 2*Y12**2*Y23 + Y12*Y13*Y23 + 
     -       Y13*(Y13 - 2*Y23)*Y23))/(Y12*Y13) + 
     -  (CA*(Y13**2 + Y12*(-Y12 + Y23)))/Y13 + 
     -  ((CA - 2*CF)*R1223*(Y13**3 + Y12*(Y12 - Y23)*Y23 + 
     -       Y13*Y23*(Y12 + Y23)))/Y13**2 + 
     -  CF*(Y12**2/(Y12 - Y23) - Y13**2/(Y13 - Y23) + Y23 - 
     -     (-Y12**2 + Y12*(Y13 + Y23))/Y13) + 
     -  (L13*(CA*(Y12 - Y23)*(Y12*(Y13 - Y23) + Y23*(-2*Y13 + Y23)) + 
     -       CF*(Y12*(2*Y13 - 7*Y23)*Y23 - 2*(Y13 - 2*Y23)*Y23**2 + 
     -          Y12**2*(Y13 + 3*Y23))))/(Y12 - Y23)**2 - 
     -  (L12*(CA*Y12*(Y13 - Y23)*
     -        (Y12*(Y13 - Y23) + Y23*(-2*Y13 + Y23)) + 
     -       CF*(2*Y13*(Y13 - Y23)*Y23**2 + Y12*Y23**2*(-Y13 + 2*Y23) + 
     -          Y12**2*(Y13**2 + Y13*Y23 - 2*Y23**2))))/
     -   (Y13*(Y13 - Y23)**2)


C---  PARTE B*C
      PBC = -(CF*Y12) + (CF*Y12**2)/(Y12 - Y23) + 
     -  ((-(CA*Y12**2) + CF*(3*Y12**2 + 4*Y12*Y13 + 3*Y13**2))*Y23)/
     -   (Y13*(Y12 + Y13)) + (CA*R1213*(Y12 + Y13 - 2*Y23)*Y23)/Y13 + 
     -  ((CA - 2*CF)*R1223*Y23*
     -     (Y13**3 + 2*Y12*Y13*Y23 - 2*Y13**2*Y23 + Y12*(Y12 - Y23)*Y23)
     -     )/(Y12*Y13**2) + ((CA - 2*CF)*L23*Y12*Y23*
     -     (-Y23**2 + Y12*(Y13 + Y23)))/(Y13*(Y12 + Y13)**2) + 
     -  (Y23*(CA*Y23 - CF*(2*Y12 + 2*Y13 + Y23)))/Y13 + 
     -  (L12*Y23*(CA*(Y13 - Y23)*(Y12 + Y13 - Y23)*Y23 + 
     -       CF*(-2*Y12*Y13**2 - Y13*(Y12 + Y13)*Y23 + 
     -          (2*Y12 + 3*Y13)*Y23**2 - 2*Y23**3)))/
     -   (Y13*(Y13 - Y23)**2)


C---  PARTE ERT (Contracción con g_{mu nu} en el tensor leptónico)

      PERT = (-(CA*(Y12 + Y13 - Y23)) + CF*(Y12 + Y13 - Y23))/Y13 + 
     -  ((CA - 2*CF)*R1223*(Y13 - 2*Y23))/Y12 + 
     -  (CA*R1213*(Y12 + Y13 - 2*Y23))/Y13 - (CF*Y12)/(Y12 - Y23) - 
     -  (L12*(CA*Y12*(Y13 - Y23) + 
     -       CF*(Y12*(Y13 - 2*Y23) + 4*Y23*(-Y13 + Y23))))/
     -   (Y13 - Y23)**2
     

C---AND WE COMBINE IT WITH THE APPROPRIATE WEIGHTS
      BDPVGAX = A**2*PA2+B**2*PBC2+ A*B*PABC+B*C*PBC+PERT*(-1/2.D0)

C---  Hay un fator  1/2 adicional porque cómo comvertimos p12->Y12
      BDPVGAX=16*PISQ*BDPVGAX/2.D0
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BDPVTRI(P,Q,I,J,K,AX,CH)
      IMPLICIT NONE      
C---RETURN THE PART OF THE ONE-LOOP CONTRIBUTION COMING FROM TRIANGLE
C     DIAGRAMS (ONLY FOR Z).
C     CH: CHANNEL, 1 FOR GLUON - 0 FOR QUARKS
C     AX: 0 FOR AAAA PART, 1 FOR AVAVA PART
  
      INTEGER I,J,K,AX,CH,IPOL
      DOUBLE PRECISION P(4,7),Q(4),DOT,R,RR,X,Y,Z,DILOG,EMSQ,
     $     Y12,Y13,Y23,L12,L13,L23,A,B,C,D,P1Q,P2Q,P3Q,QQ,
     $     R1312,R2312,R2313,R1213,R1223,R1323,
     $     PA2,PBC2,PABC,PBC,PERT  
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON  /SUBCOM/ CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON/FLAGPOL/IPOL

      RR(X,Y)=DILOG(X)+DILOG(Y)-DILOG(X*Y)-PISQ/6
      R(X,Y)=RR((1-X)/Y,(1-Y)/X)
      EMSQ=DOT(P,5,5)
      Y12=2*DOT(P,ABS(I),ABS(J))/EMSQ*SIGN(1,I*J)
      Y13=2*DOT(P,ABS(I),ABS(K))/EMSQ*SIGN(1,I*K)
      Y23=2*DOT(P,ABS(J),ABS(K))/EMSQ*SIGN(1,J*K)
      L12=LOG(ABS(Y12))
      L13=LOG(ABS(Y13))
      L23=LOG(ABS(Y23))
C---WE DECOMPOSE LEPTON Q AS A*P1+B*P2+C*P3+D*E_(MU,P1,P2,P3)/EMSQ
      P1Q=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)-P(2,ABS(I))*Q(2)
     $     -P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2Q=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)-P(2,ABS(J))*Q(2)
     $     -P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3Q=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)-P(2,ABS(K))*Q(2)
     $     -P(1,ABS(K))*Q(1))*SIGN(1,K)
      QQ=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)

      A=(P3Q*Y12+P2Q*Y13-P1Q*Y23)/(Y12*Y13*EMSQ)
      B=(P1Q*Y23+P3Q*Y12-P2Q*Y13)/(Y12*Y23*EMSQ)
      C=(P2Q*Y13+P1Q*Y23-P3Q*Y12)/(Y23*Y13*EMSQ)
C     D=SQRT(MAX(Z,4*(QQ/EMSQ-A*B*Y12-A*C*Y13-B*C*Y23)/(-Y12*Y13*Y23)))

C---THEN CALCULATE THE CONTRACTION WITH THE LEVI-CIVITA TENSOR
      R1213=R(-Y13,-Y12)
      R1223=R(Y23,-Y12)
      R1323=R(Y23,-Y13)
C-----------------------------------------------------------------------
C---  GLUON CHANNEL
      IF (CH.EQ.1) THEN
        IF (AX.EQ.0) THEN
          IF (IPOL.EQ.0) THEN
C---  PARTE AAAA
      BDPVTRI =  TR*((-L23 + Y12 + Y13)*(-8*P3Q**2*Y12*(-1 + Y13) - 
     -      8*P2Q**2*(-1 + Y12)*Y13 + 
     -      2*EMSQ*P3Q*(Y12*(-3 + Y13) + Y12**2*(-3 + Y13) - 
     -         Y13*(1 + Y13)**2) + 
     -      EMSQ**2*(Y12**3 - Y12**2*(-2 + Y13) + Y13*(1 + Y13)**2 - 
     -         Y12*(-1 + 2*Y13 + Y13**2)) + 
     -      2*P2Q*(4*P3Q*(Y12 + Y12**2 + Y13 + Y13**2) - 
     -         EMSQ*(2*Y12**2 + Y12**3 + 3*Y13*(1 + Y13) - 
     -            Y12*(-1 + Y13 + Y13**2)))))/
     -  (2.*EMSQ*Y12*Y13*(Y12 + Y13)**2)
          ELSE
C---  PARTE AAAApol
      BDPVTRI = -TR*(((1 + L23)*(Y12 + Y13) - L23*Y23)*
     -     (2*(EMSQ - 4*(P2Q + P3Q))*Y12*Y13*(Y12 + Y13) + 
     -       2*(P3Q*Y12*(Y12 + 7*Y13) + 
     -          Y13*(-3*EMSQ*Y12 + P2Q*(7*Y12 + Y13)))*Y23 + 
     -       (EMSQ - 2*(P2Q + P3Q))*(Y12 + Y13)*Y23**2))/
     -  (Y12*Y13*(Y12 + Y13)**2)
          ENDIF

        ELSEIF (AX.EQ.1) THEN
          IF (IPOL.EQ.0) THEN
C---  PARTE AVAV
      BDPVTRI =  TR*((Y12 + L23*Y12 + Y13 + L23*Y13 - L23*Y23)*
     -    (2*P2Q*(4*Y12**2*Y13 + Y13*Y23*(-Y13 + Y23) + 
     -         Y12*(2*Y13**2 - 3*Y13*Y23 - Y23**2)) + 
     -      EMSQ*(-2*Y12**2*Y13 - Y13*Y23**2 + 
     -         Y12*(2*Y13**2 + Y23**2)) + 
     -      P3Q*(2*Y13*Y23**2 + Y12**2*(-4*Y13 + 2*Y23) - 
     -         2*Y12*(4*Y13**2 - 3*Y13*Y23 + Y23**2))))/
     -  (Y12*Y13*(Y12 + Y13)**2)
          ELSE
C---  PARTE AVAVpol
      BDPVTRI = -TR*((-L23 + Y12 + Y13)*
     -     (-8*P3Q**2*Y12*(-1 + Y13) + 8*P2Q**2*(-1 + Y12)*Y13 + 
     -       8*P2Q*P3Q*(Y12 - Y13)*(1 + Y12 + Y13) + 
     -       2*EMSQ*(-(Y12*(1 + Y12)*(P2Q + 3*P3Q + P2Q*Y12)) + 
     -          (P3Q - P2Q*(-3 + Y12) + P3Q*Y12*(1 + Y12))*Y13 + 
     -          (2*P3Q - P2Q*(-3 + Y12))*Y13**2 + P3Q*Y13**3) + 
     -       EMSQ**2*(Y12 - Y13)*(Y12*(2 + Y12) + (1 + Y13)**2)))/
     -  (EMSQ*Y12*Y13*(Y12 + Y13)**2)
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
C---  QUARK CHANNEL
      ELSE!IF (CH.EQ.0) THEN
C---  PARTE AAAA (o AVAVpol)

       IF ((IPOL.EQ.0.AND.AX.EQ.0) .OR. (IPOL.EQ.1.AND.AX.EQ.1)) THEN 
      BDPVTRI = CF*((1 + L13 + Y13)*(8*P3Q**2*(1 + Y13)**2 - 
     -      8*P2Q**2*(-1 + Y12)*(1 + Y12 + Y13) - 
     -      2*EMSQ*P3Q*(1 + Y13)*
     -       (4 + 5*Y13 + Y13**2 + 2*Y12*(2 + Y13)) + 
     -      EMSQ**2*(2*Y12**2*Y13 + (1 + Y13)**2*(2 + Y13) + 
     -         Y12*(2 + 5*Y13 + 3*Y13**2)) + 
     -      2*P2Q*(4*P3Q*(2 + 2*Y12 + 3*Y13 + Y13**2) + 
     -         EMSQ*(-4 - 7*Y13 + 2*Y12**2*Y13 - 3*Y13**2 + 
     -            Y12*(-4 - Y13 + Y13**2)))))/
     -  (4.*EMSQ*Y12*(1 + Y13)**2*(1 + Y12 + Y13))

        ELSEIF ((IPOL.EQ.1.AND.AX.EQ.0).OR.(IPOL.EQ.0.AND.AX.EQ.1)) THEN
C---  PARTE AVAV (o AAAApol)               
      BDPVTRI = CF*((Y12 + L13*Y12 + L13*Y13 - (1 + L13)*Y23)*
     -    (2*Y12*Y13*(-(P3Q*Y12) + P2Q*Y13) + 
     -      (12*P3Q*Y12*(Y12 + Y13) - EMSQ*Y12*(2*Y12 + 3*Y13) + 
     -         2*P2Q*(4*Y12**2 + 3*Y12*Y13 + Y13**2))*Y23 + 
     -      (-4*(P2Q + 3*P3Q)*Y12 - 2*(P2Q + P3Q)*Y13 + 
     -         EMSQ*(4*Y12 + Y13))*Y23**2))/(4.*Y12*(Y12 - Y23)**2*Y23)
        ENDIF
C-----------------------------------------------------------------------
      ENDIF

      BDPVTRI=16*PISQ*BDPVTRI
      END
C-----------------------------------------------------------------------
      SUBROUTINE KPFUNS(X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
      IMPLICIT NONE
C---EVALUATE THE SUM OF THE K AND P FUNCTIONS.
C   IF X<0, RETURN THE DELTA-FUNCTION AND `PLUS' SUBTRACTIONS FOR -X
C   KQF=-SUM_I T_I.T_Q/T_I.T_I GAMMA_I/C_Q
C   PQF=-SUM_I T_I.T_Q/T_Q.T_Q*LOG(Q^2/2P_I.P_Q)
C   WHERE Q IS AN EXTERNALLY AGREED RENORMALIZATION POINT
C   AND LIKEWISE KGF AND PGF
      DOUBLE PRECISION X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG,Z,L,S,
     $     DIS,DILOG,D,LM
      INTEGER SCHEME,NF,IPOL
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL
      Z=ABS(X)
      S=1
      IF (X.LE.0) S=-1
      L=LOG((1-Z)/Z)
C---THE PLUS DISTRIBUTIONS
      QQ=QQ+S*XJAC*CF*2/(1-Z)*(L-KQF/2)
     $     -S*XJAC*CF*(1+Z**2)/(1-Z)*(LOG(SCALE)+PQF)
      GG=GG+S*XJAC*CA*2/(1-Z)*(L-KGF/2)
     $     -S*XJAC*CA*2/(1-Z)*(LOG(SCALE)+PGF)
      IF (SCHEME.NE.0) THEN
        DIS=S*XJAC*CF*((1+Z**2)/(1-Z)*(L-0.75)+0.25*(9+5*Z))
        QQ=QQ-DIS
        QG=QG+DIS
      ENDIF
      IF (X.LE.0) THEN
C---THE DELTA FUNCTIONS
c     (the XMIN terms appear because Z is integrated between Xmin and 1)
        D=DILOG(1-XMIN)
        LM=LOG(1-XMIN)
        QQ=QQ-CF*(5-PISQ+KQF+PISQ/3-LM**2-2*D+KQF*LM
     $       +(2*LM+XMIN+XMIN**2/2)*(LOG(SCALE)+PQF))
        GG=GG-CA*(50D0/9-PISQ+KGF+PISQ/3-LM**2-2*D+KGF*LM
     $       +2*LM*(LOG(SCALE)+PGF))+TR*NF*16D0/9
     $       -(11D0/6*CA-2D0/3*NF*TR)*(LOG(SCALE)+PGF)
      ELSE
C---THE SMOOTH FUNCTIONS
C.......................................................................
      IF (IPOL.EQ.0) THEN
        QQ=QQ+XJAC*CF*(-(1+Z)*L+(1-Z))
        GQ=GQ+XJAC*TR*((Z**2+(1-Z)**2)*(L-LOG(SCALE)-PQF)+2*Z*(1-Z))
        QG=QG+XJAC*CF*((1+(1-Z)**2)/Z*(L-LOG(SCALE)-PGF)+Z)
        GG=GG+XJAC*CA*((1-Z)/Z-1+Z*(1-Z))*2*(L-LOG(SCALE)-PGF)
C.......................................................................
      ELSEIF (IPOL.EQ.1) THEN
        QQ=QQ+XJAC*CF*(-(1+Z)*L+(1-Z))
        GQ=GQ+XJAC*TR*((Z**2-(1-Z)**2)*(L-LOG(SCALE)-PQF)+2*(1-Z))
        QG=QG+XJAC*CF*(((1-Z)**2-1)/Z*(L-LOG(SCALE)-PGF)+2*(Z-1))
        GG=GG+XJAC*CA*(((1-Z)/Z-1+Z*(1-Z)-(1-Z)**3/Z)*2*(L-LOG(SCALE)-
     $     PGF)-4*(1-Z))
c     (we work in the scheme were the polarized Pqq is equal to Pqq)
C.......................................................................
      ENDIF

        IF (SCHEME.NE.0) THEN
          DIS=XJAC*TR*((Z**2+(1-Z)**2)*L+8*Z*(1-Z)-1)
          GQ=GQ-DIS
          GG=GG+2*NF*DIS
        ENDIF
      ENDIF
      END
C-----------------------------------------------------------------------
      SUBROUTINE COLTHR(S,P,W,*)
      IMPLICIT NONE
C---GENERATE A COLLINEAR SPLITTING TO GIVE THREE PARTONS
C   AND EVALUATE THE WEIGHT FOR IT
      INTEGER I,IPOL
      DOUBLE PRECISION S,P(4,7),W(-6:6),M(-6:6),O,X,XJAC,XMIN,
     $     QQ,GQ,QG,GG,KQF,DOT
      PARAMETER (O=0)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTWO(P,M)
C---IN FACT THE GENERATION WAS ALREADY DONE EARLIER
      CALL GETCOL(X,XJAC)
C---SO WE JUST HAVE TO CALCULATE THE WEIGHT
      XMIN=2*DOT(P,1,6)/S
      QQ=0
      GG=0
      GQ=0
      QG=0
      KQF=1.5
      CALL KPFUNS(X,XJAC,XMIN,KQF,O,O,O,QQ,GQ,QG,GG)
C---THE TOTAL
      W(0)=GG*M(0)
      DO I=-6,6
        IF (I.NE.0) THEN
          W(I)=QG*M(0)+QQ*M(I)
          IF (ABS(I).LE.NF) W(0)=W(0)+GQ*M(I)
        ENDIF
      ENDDO
C---AND THE KINEMATICS
      DO I=1,4
        P(I,1)=P(I,1)/X
        P(I,3)=P(I,1)*(1-X)
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE COLFOR(S,P,W,*)
      IMPLICIT NONE
C---GENERATE A COLLINEAR SPLITTING TO GIVE FOUR PARTONS
C   AND EVALUATE THE WEIGHT FOR IT
      INTEGER I
      DOUBLE PRECISION S,P(4,7),W(-6:6),M(-6:6),X,XJAC,XMIN,
     $     QQ,GQ,QG,GG,KQF,KGF,PQF,PGF,L12,L13,DOT,EMSQ
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTHR(P,M)
C---IN FACT THE GENERATION WAS ALREADY DONE EARLIER
      CALL GETCOL(X,XJAC)
C---SO WE JUST HAVE TO CALCULATE THE WEIGHT
      XMIN=2*DOT(P,1,6)/S
      QQ=0
      QG=0
      GQ=0
      GG=0
      KQF=(1.5*(CF-CA/2)+0.5*(11D0/6*CA-2D0/3*NF*TR))/CF
      KGF=1.5
      EMSQ=-DOT(P,5,5)
      L12=LOG(2*DOT(P,1,2)/EMSQ)
      L13=LOG(2*DOT(P,1,3)/EMSQ)
      PQF=-((CF-CA/2)/CF*L12+CA/2/CF*L13)
      PGF=-(L12+L13)/2
      CALL KPFUNS(X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
C---THE TOTAL
      W(0)=GG*M(0)
      DO I=-6,6
        IF (I.NE.0) THEN
          W(I)=QG*M(0)+QQ*M(I)
          IF (ABS(I).LE.NF) W(0)=W(0)+GQ*M(I)
        ENDIF
      ENDDO
C---AND THE KINEMATICS
      DO I=1,4
        P(I,1)=P(I,1)/X
        P(I,4)=P(I,1)*(1-X)
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE SUBTHR(PERM,SS,P,Q,S,JAC,*)
      IMPLICIT NONE
C---GENERATE A TWO-PARTON STATE FROM A THREE-PARTON STATE,
C   CALCULATE THE JACOBIAN FACTOR FOR THE CORRESPONDING CHANNEL,
C   AND (IF PERM.GT.0) THE APPROXIMATE MATRIX-ELEMENT.
      INTEGER NPERM,PERM,IPERM,IJF,KF,I,J,K,M
      DOUBLE PRECISION SS,P(4,7),Q(4,7),S(-6:6),JAC,
     $     X,Z,DEN,F,QQ,GQ,EMSQ,DOT,XMIN,GQ1,GQ2
      PARAMETER (NPERM=1)
      DIMENSION IPERM(3,NPERM),IJF(NPERM),KF(NPERM)
      INTEGER SCHEME,NF,IPOL
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      COMMON/FLAGPOL/IPOL
      DATA IPERM/2,3,1/
      DATA IJF/2/
      DATA  KF/1/
      IF (ABS(PERM).GT.NPERM.OR.PERM.EQ.0)
     $     STOP 'PERM TOO BIG IN SUBTHR!'
C---FIND WHICH PARTONS TO TREAT
      I=IPERM(1,ABS(PERM))
      J=IPERM(2,ABS(PERM))
      K=IPERM(3,ABS(PERM))
C---FIND THE KINEMATIC VARIABLES
      EMSQ=2*(-DOT(P,I,J)+DOT(P,J,K)+DOT(P,K,I))
      DEN=2/EMSQ
      X=1/(1+DOT(P,I,J)*DEN)
      Z=DOT(P,I,K)*X*DEN
C--- GPS -- avoid crashes later on
      if (z .lt. cutoff .or. (1d0-z) .lt. cutoff 
     $     .or. (1d0-x) .lt. cutoff .or. x .lt. cutoff) return 1
C---COPY INTO Q, REPLACING K BY KTILDE AND I BY IJTILDE
      DO M=1,4
        Q(M,IJF(ABS(PERM)))=P(M,I)+P(M,J)-(1-X)*P(M,K)
        Q(M,KF(ABS(PERM)))=X*P(M,K)
        Q(M,5)=P(M,5)
        Q(M,6)=P(M,6)
        Q(M,7)=P(M,7)
      ENDDO
C---CALCULATE THE CORRESPONDING JACOBIAN FACTOR
      XMIN=2*DOT(Q,1,6)/SS
C--- GPS to avoid crashes -------
      if (xmin .lt. cutoff) return 1
      JAC=1/(2*NPOW(1)*(Z*(1-Z))**XPOW(1))*(Z**XPOW(1)+(1-Z)**XPOW(1))
     $     *(0.5/(-X*LOG(XMIN))
     $     +0.5*((1-XMIN)/(1-X))**XPOW(1)/(NPOW(1)*(1-XMIN)))
      IF (PERM.LT.0) RETURN
C---OVERALL FACTOR
      F=16*PISQ/EMSQ

      IF (IPOL.EQ.0) THEN
C-----------------------------------------------------------------------
C---CALCULATE WEIGHT FOR FINAL-STATE SPLITTING FUNCTION
      QQ=F*CF*(2/(2-Z-X)-(1+Z))/(1-X)
C---CALCULATE WEIGHTS FOR INITIAL-STATE SPLITTING FUNCTIONS
      QQ=QQ+F*CF*(2/(2-Z-X)-(1+X))/(1-Z)
c      GQ=F*TR*(X**2+(1-X)**2)/(Z*(1-Z))
      GQ2=F*TR*(X**2+(1-X)**2)/(Z)
      GQ1=F*TR*(X**2+(1-X)**2)/((1-Z))
c     (total of 4 dipoles)  
C-----------------------------------------------------------------------
      ELSEIF (IPOL.EQ.1) THEN
C-----------------------------------------------------------------------
C---CALCULATE WEIGHT FOR FINAL-STATE SPLITTING FUNCTION
      QQ=F*CF*(2/(2-Z-X)-(1+Z))/(1-X)
C---CALCULATE WEIGHTS FOR INITIAL-STATE SPLITTING FUNCTIONS
      QQ=QQ+F*CF*(2/(2-Z-X)-(1+X))/(1-Z)
      GQ2=F*TR*(2*X-1)/(Z)
      GQ1=F*TR*(2*X-1)/(1-Z) 
C-----------------------------------------------------------------------
      ENDIF

C---MULTIPLY WITH THE LOWEST-ORDER MATRIX ELEMENT
      CALL MATTWO(Q,S)
      S(0)=0
C---COMBINE THE GLUON DIPOLES   
      DO I=1,NF
        S(0)=S(0)+(GQ1*S(I)+GQ2*S(-I))
      ENDDO

      DO I=-6,6
        IF (I.NE.0) S(I)=QQ*S(I)
      ENDDO

C---READJUST THE MOMENTA A BIT
      DO M=1,4
        Q(M,3)=P(M,1)-Q(M,1)
        Q(M,1)=P(M,1)
C just to make sure DdeF where P(M,4)=0d0  4 does not exist here
        Q(M,4)=P(M,4)            
      ENDDO

 999  END
C-----------------------------------------------------------------------
      SUBROUTINE SUBFOR(PERM,SS,P,Q,S,JAC,*)
      IMPLICIT NONE
C---GENERATE A THREE-PARTON STATE FROM A FOUR-PARTON STATE,
C   CALCULATE THE JACOBIAN FACTOR FOR THE CORRESPONDING CHANNEL,
C   AND (IF PERM.GT.0) THE APPROXIMATE MATRIX-ELEMENT.
      INTEGER NPERM,PERM,IPERM,IJF,KF,LF,NPERM3,I,J,K,L,M,n
      INTEGER EWFLAG
      DOUBLE PRECISION SS,P(4,7),Q(4,7),S(-6:6),JAC,
     $     X,Z,DEN,QQ,GQ,EMSQ,DOT,XMIN,JTMP,STMP(-6:6),QTMP(4,7),QUSQ,
     $     Y,ZI,ZJ,OY,ZTI,ZTJ,V(4),VV,S1(-6:6),S2(-6:6),SC,S3(-6:6),CUT
     $     ,s4(-6:6),s5(-6:6),s6(-6:6),s7(-6:6),s8(-6:6),gg,qg,gq2
     $     ,s9(-6:6),s10(-6:6),s11(-6:6),temp,VDOTI,LEPCHARGE,CKM2(3,3)
      PARAMETER (NPERM=6,NPERM3=1)
      DIMENSION IPERM(4,NPERM),IJF(NPERM),KF(NPERM),LF(NPERM)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE,VDOT(4)
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/FLAGPOL/IPOL
      INTEGER NPOW(2),IPOL
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,SCZ,MZ,GAMMAZ
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG  

      DATA IPERM/2,3,4,1, 2,4,3,1, 2,3,1,4, 2,4,1,3, 3,4,2,1, 3,4,1,2/
      DATA IJF/2,2,2,2,3,3/
      DATA  KF/3,3,1,1,2,1/
      DATA  LF/1,1,3,3,1,2/
      IF (ABS(PERM).GT.NPERM.OR.PERM.EQ.0)
     $     STOP 'PERM TOO BIG IN SUBFOR!'

      do i=-6,6
        s(i)=0
        s1(i)=0
        s2(i)=0
        s3(i)=0
        s4(i)=0
        s5(i)=0
        s6(i)=0
        s7(i)=0
        s8(i)=0
        s9(i)=0
        s10(i)=0
        s11(i)=0
      enddo

C---FIND WHICH PARTONS TO TREAT
      I=IPERM(1,ABS(PERM))
      J=IPERM(2,ABS(PERM))
      K=IPERM(3,ABS(PERM))
      L=IPERM(4,ABS(PERM))
C---FIND THE KINEMATIC VARIABLES
      IF (K.EQ.1) THEN
        EMSQ=2*(-DOT(P,I,J)+DOT(P,J,K)+DOT(P,K,I))
        DEN=2/EMSQ
        !write(0,*) DOT(P,I,J)*DEN
C---GPS avoid division by zero: but should one perhaps be in any case
C   worried if DOT(P,I,J)*DEN<0?
        X = (1+DOT(P,I,J)*DEN)
        if (X .LE. CUTOFF) return 1
        X = 1/X
c$$$        X=1/(1+DOT(P,I,J)*DEN)
        Z=DOT(P,I,K)*X*DEN
C---GPS --------------
        if (Z .LE. CUTOFF .OR. 1-Z .LE. CUTOFF .or.(1-x).le.cutoff) then
           RETURN 1
        end if
      ELSE
        EMSQ=2*(DOT(P,I,J)+DOT(P,J,K)+DOT(P,K,I))
        DEN=2/EMSQ
        Y=DOT(P,I,J)*DEN
        ZI=DOT(P,I,K)*DEN
        ZJ=DOT(P,J,K)*DEN
C---GPS-------------------
        if (y .le. cutoff .or. 1-y .le. cutoff .or.
     $       1-zj .le. cutoff .or. 1-zi .le. cutoff) then
           return 1
        end if
        OY=1/(1-Y)
        ZTI=ZI*OY
        ZTJ=ZJ*OY
      ENDIF
C---COPY INTO Q, REPLACING K BY KTILDE AND I BY IJTILDE
      DO M=1,4
        IF (K.EQ.1) THEN
          Q(M,IJF(ABS(PERM)))=P(M,I)+P(M,J)-(1-X)*P(M,K)
          Q(M,KF(ABS(PERM)))=X*P(M,K)
        ELSE
          Q(M,IJF(ABS(PERM)))=P(M,I)+P(M,J)-Y*OY*P(M,K)
          Q(M,KF(ABS(PERM)))=OY*P(M,K)
        ENDIF
        Q(M,LF(ABS(PERM)))=P(M,L)
        Q(M,5)=P(M,5)
        Q(M,6)=P(M,6)
        Q(M,7)=P(M,7)
c        Q(M,4)=0.d0 !debug
      ENDDO
C---REENFORCE THE MINIMUM CUTOFF ON ALL PAIR MASSES
      CUT=CUTOFF*(DOT(Q,1,2)+DOT(Q,1,3))
      IF (DOT(Q,1,2).LT.CUT.OR.DOT(Q,1,3).LT.CUT.OR.DOT(Q,2,3).LT.CUT)
     $     RETURN 1
C---CALCULATE THE CORRESPONDING JACOBIAN FACTOR
      JAC=0
      DO M=1,NPERM3
        CALL SUBTHR(-M,SS,Q,QTMP,STMP,JTMP,*999)
        JAC=JAC+JTMP
      ENDDO
      IF (K.EQ.1) THEN
        XMIN=2*DOT(Q,1,6)/SS
        JAC=JAC/(2*NPOW(2)*(Z*(1-Z))**XPOW(2))
     $       *(Z**XPOW(2)+(1-Z)**XPOW(2))
     $       *(0.5/(-X*LOG(XMIN))
     $       +0.5*((1-XMIN)/(1-X))**XPOW(2)/(NPOW(2)*(1-XMIN)))
      ELSE
        IF (ABS(PERM).LE.4) THEN
          JAC=JAC/(NPOW(2)**2*(Y*(1-ZI))**XPOW(2))
        ELSE
          JAC=JAC/(2*NPOW(2)**2*(Y*(1-ZI)*(1-ZJ))**XPOW(2)
     $         /((1-ZI)**XPOW(2)+(1-ZJ)**XPOW(2)))
        ENDIF
        JAC=JAC*2
      ENDIF
      QUSQ=-DOT(P,5,5)
      JAC=JAC*QUSQ/EMSQ
C---INCLUDE A PRIORI CHANNEL WEIGHTS
      JAC=JAC/8
      IF (ABS(PERM).GT.4) JAC=JAC*2
      IF (PERM.LT.0) RETURN

C#######################################################################
C     QUARK CHANNEL Q->QGG
C#######################################################################
C.......................................................................
C---CALCULATE WEIGHT FOR QUARK-GLUON SPLITTING FUNCTION
C.......................................................................
      IF (PERM.LE.4) THEN
        CALL MATTHR(Q,S1)
        IF (K.EQ.1) THEN
C     DIPOLO DE ESTADO FINAL QG (1 ESPECTADOR) + DIPOLO DE 
C     ESTADO INICIAL QG (+ U<->1-Z) (Vqg_pol=Vqg_nopol en 4D)          
          QQ=16*PISQ/EMSQ*(X**2+Z**2)/((1-X)*(1-Z))         
        ELSE
C     DIPOLO DE ESTADO FINAL QG SIN PARTON INICIAL      
          QQ=16*PISQ/EMSQ*(2/(1-ZTI*(1-Y))-(1+ZTI))/Y

        ENDIF
        IF (KF(PERM).EQ.3) THEN
          QQ=QQ*HF*CA
        ELSE
          QQ=QQ*(CF-HF*CA)
        ENDIF
C---SYMMETRY FACTORS
        QQ=QQ/2
        S1(0)=0
        DO M=-6,6
          IF (M.NE.0) S1(M)=QQ*S1(M)
        ENDDO
        DO M=-6,6
          S2(M)=0
          S3(M)=0
        ENDDO
      ELSE
C.......................................................................        
C     CALCULATE WEIGHT FOR GLUON-GLUON SPLITTING FUNCTION
C.......................................................................
        DO M=-6,6
          S1(M)=0
        ENDDO
        CALL MATTHR(Q,S2)
        IF (K.EQ.1) THEN
C     DIPOLOS QG DE ESTADO INCIAL QG (1-3 Y 1-4 CON 4/3 DE ESPECATDOR)+
C     DIPOLO DE ESTADO FINAL GG (CON 1 PARTON ININCIAL).
C     (LOS DIPOLOS QG POLARIZADOS SON IGUALES A LOS NO POLARIZADOS EN 4D)
          DO M=1,4
            V(M)=Z*P(M,I)-(1-Z)*P(M,J)
          ENDDO
          VV = -2*Z*(1-Z)*DOT(P,I,J)
          IF (IPOL.EQ.0) THEN
            CALL CONTHR(Q,V,VV,2,-1,3,SC)
            IF (EWFLAG.GE.1) CALL CONTHRPOLQ(Q,V,VV,2,1,3,SCZ)
          
          ELSEIF (IPOL.EQ.1) THEN
            CALL CONTHRPOLQ(Q,V,VV,2,1,3,SC)
            IF (EWFLAG.GE.1) CALL CONTHR(Q,V,VV,2,-1,3,SCZ)
          ENDIF

          do m=-6,6
            s3(m)=0
          enddo

          SC=16*pisq/emsq*SC*4*z*(1-z)*hf*ca/2/(1-x)
          SCZ=16*pisq/emsq*SCZ*4*z*(1-z)*hf*ca/2/(1-x)
          CALL COMBINEQI(P,s3,SC,SCZ)

          QQ=16*PISQ/EMSQ*
     $         ((2/(2-Z-X)+2/(Z+1-X)-4)/(1-X)
     $         +(2/(2-Z-X)-(1+X))/(1-Z)+(2/(Z+1-X)-(1+X))/Z)
          QQ=QQ*HF*CA
          QQ=QQ/2
        ELSE
C---  DIPOLO DE ESTADO FINAL GG (IGUAL AL NO POLARIZADO)          
          DO M=1,4
            V(M)=ZTI*P(M,I)-ZTJ*P(M,J)
          ENDDO
          VV = -2*ZTI*ZTJ*DOT(P,I,J)
          IF(IPOL.EQ.0) THEN
            CALL CONTHR(Q,V,VV,2,-1,3,SC)
            IF (EWFLAG.GE.1) CALL CONTHRPOLQ(Q,V,VV,2,1,3,SCZ)
          ELSEIF(IPOL.EQ.1) THEN
            CALL CONTHRPOLQ(Q,V,VV,2,1,3,SC)
            IF (EWFLAG.GE.1) CALL CONTHR(Q,V,VV,2,-1,3,SCZ)
          ENDIF

          do m=-6,6
            s3(m)=0
          enddo

          SC=16*pisq/emsq*SC*4*zti*ztj*hf*ca/2/y
          SCZ=16*pisq/emsq*SCZ*4*zti*ztj*hf*ca/2/y
          CALL COMBINEQI(P,s3,SC,SCZ)

          QQ=16*PISQ/EMSQ*(2/(1-ZTI*(1-Y))+2/(1-ZTJ*(1-Y))-4)/Y
          QQ=QQ*HF*CA
          QQ=QQ/2
        ENDIF
        S2(0)=0
        DO M=-6,6
          IF (M.NE.0) S2(M)=QQ*S2(M)
        ENDDO
      ENDIF

C#######################################################################
C     GLUON CHANNEL G->QQbG
C#######################################################################
C.......................................................................
c     calculate weight for quark-gluon splitting function (all fsr)
C.......................................................................
C     DIP DE ESTADO FINAL QG CON PARTON ININCIAL (1) (= NO POL)
      if (perm.ne.1.and.perm.ne.3) then
        call matthr(q,s4)
        if (k.eq.1) then
          gg=16*pisq/emsq*(2/(2-z-x)-(1+z))/(1-x)
          gg=gg*hf*ca
        else
C     DIP DE ESTADO FINAL QG SIN PARTON ININCIAL (1) (= NO POL)          
          gg=16*pisq/emsq*(2/(1-zti*(1-y))-(1+zti))/y
          gg=gg*(cf-hf*ca)
        endif
        s4(0)=s4(0)*gg
        do m=-6,6
          if (m.ne.0) s4(m)=0
        enddo
      else
        do m=-6,6
          s4(m)=0
        enddo
      endif
C.......................................................................      
c     calculate weight for gluon-quark splitting function (isr)
C.......................................................................
      if (perm.eq.3.or.perm.eq.4.or.perm.eq.6) then
        if (perm.eq.3) then
C     DIV DE ESTADO FINAL GQ CON PARTON ININCIAL (1) 
          call matthr(q,s5)
          IF (IPOL.EQ.0) gq=16*pisq/emsq*(x**2+(1-x)**2)
          IF (IPOL.EQ.1) gq=16*pisq/emsq*(x**2-(1-x)**2)
          gq2=gq*(cf-hf*ca)/cf*tr/z !initial antiquark in born          
          gq=gq*(cf-hf*ca)/cf*tr/(1-z) !initial quark in born
        else
CC-   EN EL BORN TIENE LAS ETIQUETEAS AL REVES (2=Q, 3=QB), 
CC    ASI QUE CAMBIA LOS MOMENTOS ANTES DE LLAMAR AL BORN
        IF(perm.eq.4) THEN !debug_new
          do m=1,4
            temp=q(m,2)
            q(m,2)=q(m,3)
            q(m,3)=temp
          enddo
          call matthr(q,s5)
          do m=1,4
            temp=q(m,2)
            q(m,2)=q(m,3)
            q(m,3)=temp
          enddo
          IF (IPOL.EQ.0) gq=16*pisq/emsq*(x**2+(1-x)**2)/z
          IF (IPOL.EQ.1) gq=16*pisq/emsq*(x**2-(1-x)**2)/z
c          gq=gq*hf*ca/cf*tr
c          gq2=0.d0          
          gq2=gq*hf*ca/cf*tr
          gq=0.d0          
        ELSE
            call matthr(q,s5)
          IF (IPOL.EQ.0) gq=16*pisq/emsq*(x**2+(1-x)**2)/z
          IF (IPOL.EQ.1) gq=16*pisq/emsq*(x**2-(1-x)**2)/z
c          gq2=gq*hf*ca/cf*tr
c          gq=0.d0            
          gq=gq*hf*ca/cf*tr
          gq2=0.d0            
        ENDIF
        endif
        s5(0)=0
c       Q channel born is not summed nf flavours, unlike the G one.
        do m=1,nf
c          s5(0)=s5(0)+gq*s5(m)
          s5(0)=s5(0)+gq*s5(m)+gq2*s5(-m)
        enddo
        do m=-6,6
          if (m.ne.0) s5(m)=0
        enddo
      else
        do m=-6,6
          s5(m)=0
        enddo
      endif
C.......................................................................      
c---calculate weight for gluon-gluon splitting function (isr)
C.......................................................................
      if (perm.eq.4.or.perm.eq.6) then
        call matthr(q,s6)
        do m=1,4
          v(m)=p(m,i)/z-p(m,j)/(1-z)
        enddo
        VV = -2*DOT(P,I,J)/(Z*(1-Z))
        call conthr(q,v,VV,2,3,-1,sc)
        if (ewflag.ge.1) call conthrpolq(q,v,VV,2,-3,-1,scz)

        s7(0)=0
        IF (IPOL.EQ.0) THEN
C         COMBINE THE (UNPOLARIZED) CONTRACTED DIPOLE            
          sc = -16*pisq/emsq*sc*4*(1-x)/x*hf*ca/(1-z) *tr/cf
          scz = -16*pisq/emsq*scz*4*(1-x)/x*hf*ca/(1-z) *tr/cf
          CALL COMBINEGLUON(P,s7,SC,SCZ)
        ENDIF

        do m=-6,6
          if (m.ne.0) s7(m)=0
        enddo

C     CONTRACCIÓN CON G_{MU NU} PARA EL NO POLARIZADO  
        IF (IPOL.EQ.0) gg=16*pisq/emsq*
     $       (2/(2-x-z)-2+2*x*(1-x))/(1-z)
C     CONTRACCIÓN CON VGG_POL PARA EL POLARIZADO     
        IF (IPOL.EQ.1) gg=16*pisq/emsq*
     $       (2/(2-z-x)-2+4*(1-x))/(1-z)

        gg=gg*hf*ca
        s6(0)=s6(0)*gg
        do m=-6,6
          if (m.ne.0) s6(m)=0
        enddo
      else
        do m=-6,6
          s6(m)=0
          s7(m)=0
        enddo
      endif

C#######################################################################
C     QUARK CHANNEL Q->Q QQb
C#######################################################################
C.......................................................................
c     calculate weight for quark-antiquark splitting function (fsr)
C.......................................................................
      if (perm.ge.5) then
        call matthr(q,s8)
        if (k.eq.1) then
          do m=1,4
            v(m)=z*p(m,i)-(1-z)*p(m,j)
          enddo
          VV = -2*Z*(1-Z)*DOT(P,I,J)
          IF (IPOL.EQ.0) THEN
            call CONTHR(q,v,VV,2,-1,3,sc)
            IF (EWFLAG.GE.1) CALL CONTHRPOLQ(Q,V,VV,2,1,3,SCZ)
          ELSEIF (IPOL.EQ.1) THEN
            call CONTHRPOLQ(q,v,VV,2,1,3,sc)
            IF (EWFLAG.GE.1) CALL CONTHR(Q,V,VV,2,-1,3,SCZ)
          ENDIF

          do m=-6,6
            s9(m)=0D0
          enddo

C         COMBINE THE CONTRACTED DIPOLE            
          sc = -16*pisq/emsq*sc*4*z*(1-z)*hf*tr*nf/(1-x)
          scz = -16*pisq/emsq*scz*4*z*(1-z)*hf*tr*nf/(1-x)
          CALL COMBINEQI(P,s9,SC,SCZ)

          qq=16*pisq/emsq/(1-x)
          qq=qq*hf*tr*nf

        else
          do m=1,4
            v(m)=zti*p(m,i)-ztj*p(m,j)
          enddo
          VV = -2*ZTI*ZTJ*DOT(P,I,J)
          IF (IPOL.EQ.0) THEN
            call CONTHR(q,v,VV,2,-1,3,sc)
            IF (EWFLAG.GE.1) CALL CONTHRPOLQ(Q,V,VV,2,1,3,SCZ)
          ELSEIF (IPOL.EQ.1) THEN
            call CONTHRPOLQ(q,v,VV,2,1,3,sc)
            IF (EWFLAG.GE.1) CALL CONTHR(Q,V,VV,2,-1,3,SCZ)
          ENDIF
            
          do m=-6,6
            s9(m)=0D0
          enddo

C         COMBINE THE CONTRACTED DIPOLE            
          sc = -16*pisq/emsq*sc*4*zti*ztj*hf*tr*nf/y
          scz = -16*pisq/emsq*scz*4*zti*ztj*hf*tr*nf/y
          CALL COMBINEQI(P,s9,SC,SCZ)

          qq=16*pisq/emsq/y
          qq=qq*hf*tr*nf
        endif
        s8(0)=0
        do m=-6,6
          if (m.ne.0) s8(m)=qq*s8(m)
        enddo
      endif
C.......................................................................      
c     calculate weight for quark-antiquark splitting function (isr)
C.......................................................................
      if (perm.eq.3.or.perm.eq.4) then

        if (perm.eq.4) then
          do m=1,4
            temp=q(m,2)
            q(m,2)=q(m,3)
            q(m,3)=temp
          enddo
        endif

        call matthr(q,s10)
        do m=1,4
          v(m)=p(m,i)/z-p(m,j)/(1-z)
        enddo
        VV = -2*DOT(P,I,J)/(Z*(1-Z))
        call conthr(q,v,VV,2,3,-1,sc)
        IF (EWFLAG.GE.1) CALL conthrpolq(q,v,VV,2,-3,-1,scz)

        do m=-6,6
          s11(m)=0D0
        enddo          

        IF (IPOL.EQ.0) THEN
C         COMBINE THE (UNPOLARIZED) CONTRACTED DIPOLE            
          sc = -16*pisq/emsq*sc*4*(1-x)/x*hf*cf/z *tr/cf
          scz = -16*pisq/emsq*scz*4*(1-x)/x*hf*cf/z *tr/cf
          CALL COMBINEQF(P,s11,SC,SCZ)
        ENDIF
        s11(0)=0

C     CONTRACCIÓN CON G_{mu nu} PARA EL NO POLARIZADO         
        IF (IPOL.EQ.0) qg=16*pisq/emsq*x/z
C     CONTRACCIÓN CON VGG_POL PARA EL POLARIZADO        
        IF (IPOL.EQ.1) qg=16*pisq/emsq*(2-x)/z
        qg=qg*hf*cf
        do m=-6,6
          if (m.ne.0) s10(m)=s10(0)*qg     
        enddo
        s10(0)=0
      if (perm.eq.4) then
          do m=1,4
            temp=q(m,2)
            q(m,2)=q(m,3)
            q(m,3)=temp
          enddo
      endif        
      endif    
C.......................................................................
c      DO M=-6,6 !debug dipolos
c     canal q->qgg 
c      	S1(M)=0d0
c      	S2(M)=0d0
c      	S3(M)=0d0
c     canal g->qqbg       	
c      	S4(M)=0d0
c      	S5(M)=0d0
c      	S6(M)=0d0
c      	S7(M)=0d0
c     canal q->qqbqb 	
c      	S8(M)=0d0
c      	S9(M)=0d0
c      	S10(M)=0d0
c      	S11(M)=0d0
c      ENDDO
C.......................................................................
C---ADD THEM TOGETHER
      DO M=-6,6
        S(M)=S1(M)+S2(M)+S3(M)+s4(m)+s5(m)+s6(m)+s7(m)
     $       +s8(m)+s9(m)+s10(m)+s11(m)
      ENDDO
C---READJUST THE MOMENTA A BIT
      DO M=1,4
        Q(M,4)=P(M,1)-Q(M,1)
        Q(M,1)=P(M,1)
      ENDDO   !debug

c      DO M=1,4
c             write(6,*)Q(M,4)/Q(4,4)

c        write(6,*)"r",(Q(M,4)+Q(M,3)+Q(M,2))/(P(M,4)+P(M,3)+P(M,2))
c        write(6,*)(Q(M,4)+Q(M,3)+Q(M,2)),(P(M,4)+P(M,3)+P(M,2))
c        Q(M,1)=P(M,1)
c      ENDDO      

 999  END
C-----------------------------------------------------------------------
      FUNCTION ERTA(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT A FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S=S12+S13+S14+S23+S24+S34
      ERTA=(S12*S34**2-S13*S24*S34+S14*S23*S34+3*S12*S23*S34+
     $     3*S12*S14*S34+4*S12**2*S34-S13*S23*S24+2*S12*S23*S24-
     $     S13*S14*S24-2*S12*S13*S24+2*S12**2*S24+S14*S23**2+
     $     2*S12*S23**2+S14**2*S23+4*S12*S14*S23+4*S12**2*S23+
     $     2*S12*S14**2+2*S12*S13*S14+4*S12**2*S14+2*S12**2*S13+
     $     2*S12**3)/(2*S13*S134*S234*S24)+
     $     (S24*S34+S12*S34+S13*S24-S14*S23+S12*S13)/(S13*S134**2)+
     $     2*S23*(S-S13)/(S13*S134*S24)+
     $     S34/(2*S13*S24)
      END
C-----------------------------------------------------------------------
      FUNCTION BDPQA(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=2*DOT(P,6,7)

      BDPQA=(-(S37*S46) + S36*S47)/(4.*S13*S24)

      BDPQA = BDPQA +
     -   (-(S17*S26*S34) + S16*S27*S34-S14*S27*S36+S14*S26*S37 + 
     -    S27*S34*S46 - S26*S34*S47 + 
     -    S13*(-(S17*S26) + S16*S27 + S27*S46 - S26*S47))/
     -  (2.*S13*(S13 + S14 + S34)**2)

      BDPQA = BDPQA + 
     -   (S12*S27*S36+S14*S27*S36+ 2*S23*S27*S36 + S24*S27*S36 + 
     -    S27*S34*S36 - S12*S26*S37 - S14*S26*S37 - 2*S23*S26*S37 - 
     -    S24*S26*S37 - S26*S34*S37 + S23*S27*S46 - S23*S37*S46 - 
     -    S17*S23*(S26 + S46) - S23*S26*S47 + S23*S36*S47 + 
     -    S16*S23*(S27 + S47))/(2.*S13*S24*(S13 + S14 + S34))

      BDPQA = BDPQA + 
     -    (-(S14**2*S17*S26)-2*S14*S17*S23*S26 - S17*S23**2*S26 - 
     -    S17*S23*S24*S26 + S14**2*S16*S27 + 2*S14*S16*S23*S27 + 
     -    S16*S23**2*S27 + S16*S23*S24*S27 - 2*S14*S17*S26*S34 - 
     -    2*S17*S23*S26*S34 + 2*S14*S16*S27*S34 + 2*S16*S23*S27*S34 - 
     -    S17*S26*S34**2 + S16*S27*S34**2 + S14*S17*S24*S36 + 
     -    S14*S23*S27*S36 + S17*S24*S34*S36 - S14*S16*S24*S37 - 
     -    S14*S23*S26*S37 - S16*S24*S34*S37 - S14*S17*S23*S46 - 
     -    S14*S23*S37*S46 + S14*S16*S23*S47 + S14*S23*S36*S47 + 
     -    S13*(-(S14*S17*S26) + S14*S16*S27 - S23*S27*S46 - 
     -       S27*S34*S46 + S24*S37*S46 + S23*S26*S47 + S26*S34*S47 - 
     -       S24*S36*S47) + S12**2*
     -     (S27*S36 - S26*S37 + S27*S46 - S17*(2*S26 + S36 + S46) - 
     -       S26*S47 + S16*(2*S27 + S37 + S47)) + 
     -    S12*(-3*S17*S23*S26 - S17*S24*S26 + 3*S16*S23*S27 + 
     -       S16*S24*S27 - 2*S17*S26*S34 + 2*S16*S27*S34 - 
     -       S17*S23*S36 + S23*S27*S36 - S17*S34*S36 + S16*S23*S37 - 
     -       S23*S26*S37 + S16*S34*S37 - 2*S17*S23*S46 - S17*S24*S46 + 
     -       S23*S27*S46 + S27*S34*S46 - S23*S37*S46 - S24*S37*S46 + 
     -       2*S16*S23*S47 + S16*S24*S47 - S23*S26*S47 - S26*S34*S47 + 
     -       S23*S36*S47 + S24*S36*S47 + 
     -       S13*(-(S17*S26) + S16*S27 + S27*S36 - S26*S37 - S37*S46 + 
     -          S36*S47) + S14*
     -        (2*S27*S36 - 2*S26*S37 + S27*S46 - S37*S46 - 
     -          S17*(3*S26 + S36 + S46) - S26*S47 + S36*S47 + 
     -          S16*(3*S27 + S37 + S47))))/
     -  (4.*S13*S24*(S13 + S14 + S34)*(S23 + S24 + S34)) 
      END
C-----------------------------------------------------------------------
      FUNCTION ERTGAXA(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTGAXA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=2*DOT(P,6,7)

      ERTGAXA=((-(S14**2*S23*S36)-S14*S23*(S16*S34+S34*S36 - S13*S46) + 
     -      S13*S24*(S26*S34 + S24*S36 + S34*S36 - S23*S46))*S67)/
     -  (4.*S13*S14*S23*S24*S36)

      ERTGAXA = ERTGAXA -
     -         ((2*S12**2*S14*S23*S36 + 
     -       2*S14**2*S23*(S16*S23 - S13*S26 + S23*S36) + 
     -       2*S14*S23*(S23 + S24 + S34)*
     -        (S16*S23 - S13*S26 + S23*S36) + 
     -       S12*(2*S14**2*S23*S36 + 
     -          2*S14*S23*(S16*S23 - S13*S26 + 
     -             (2*S23 + S24 + S34)*S36) + 
     -          S13*S24*(S26*S34 - S23*S46)) - 
     -       S13*S24*(-(S16*S23*S24) - S23*S26*S34 - S26*S34**2 - 
     -          S23*S24*S36 + S24*S34*S36 + S23**2*S46 + S23*S34*S46 + 
     -          S13*(-(S26*S34) + S24*(S26 + S36) + S23*S46)))*S67)/
     -  (2.*S13*S14*S23*S24*(S13 + S14 + S34)*S36)

      ERTGAXA = ERTGAXA -
     -         ((S14*(S12*S16*S34 + S24*S34*(S16 + S36) + 
     -          S14*(-(S16*S23) + S24*S36)) + 
     -       S13**2*(S14*S26 - S26*S34 - (S12 + S24)*S46) + 
     -       S13*(S14**2*S26 + 
     -          S34*(S16*(S23 + S24) - S26*S34 + S23*S36 + S24*S36 + 
     -             S12*(S16 + S36) + S23*S46) - 
     -          S14*(S16*S23 + S12*S46 + S24*(-S36 + S46))))*S67)/
     -  (2.*S13*S14*(S13 + S14 + S34)**2*S36)

      ERTGAXA = ERTGAXA +
     -         ((-(S14**3*S23**2*S36) + 
     -      S12**2*(S13*S14*S23 + S13**2*S24 + S14*S23*(-S23 + S34) - 
     -         S13*S24*(S23 + S34))*S36 + 
     -      S14**2*S23*(S16*S23*(S23 + S34) - S23*(S23 + S34)*S36 - 
     -         S13*(S23*S26 + S26*S34 - S24*S36)) - 
     -      S13*S24*(S13 + S24 + S34)*
     -       (S13**2*S26 - S16*S23*S34 - 
     -         S13*(S16*S23 - S26*S34 + S24*S36)) + 
     -      S14*(S16*S23*(S23**3 + S13**2*S24 + S23**2*(S24 + 2*S34) + 
     -            S23*(-2*S13*S24 + S34**2)) - 
     -         S13*(S23**3*S26 + S13**2*S24*S26 + 
     -            S23**2*(2*S26*S34 + S24*(S26 - S36)) + 
     -            S23*(S26*S34**2 + S24**2*S36 + S13*S24*(-2*S26 + S36))
     -            )) + S12*(S13**3*S24*(-S26 + S36) - 
     -         S14*S23**2*(-(S16*(S23 + S34)) + 
     -            (3*S14 + S23 + S24 + S34)*S36) + 
     -         S13*S23*(S14**2*S36 + 
     -            S14*(-(S23*S26) - 2*S26*S34 + S16*(S23 + S34) + 
     -               S23*S46) + 
     -            S24*(-(S26*S34) + S16*(S23 + 2*S34) - S24*S36 + 
     -               S23*S46)) + 
     -         S13**2*(S16*S23*S24 + S14*S24*S36 - 
     -            S14*S23*(S26 + S46) + 
     -            S24*(-(S26*S34) + 3*S24*S36 + S34*S36 - 
     -               S23*(S26 + S46)))))*S67)/
     -  (4.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34)*S36)

      ERTGAXA = -ERTGAXA          

      END
C-----------------------------------------------------------------------
      FUNCTION ERTB(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT B FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTB,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S123,S124,S134,S234,S
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S123=S12+S13+S23
      S124=S12+S14+S24
      S134=S13+S14+S34
      S234=S23+S24+S34
      S=S12+S13+S14+S23+S24+S34
      ERTB=(S12*S24*S34+S12*S14*S34-S13*S24**2+S13*S14*S24+
     $     2*S12*S14*S24)/(S13*S134*S23*S14)+
     $     S12*(S+S34)*S124/(S134*S234*S14*S24)-
     $     (2*S13*S24+S14**2+S13*S23+2*S12*S13)/(S13*S134*S14)+
     $     S12*S123*S124/(2*S13*S14*S23*S24)
      END
C-----------------------------------------------------------------------
      FUNCTION BDPQB(P,I,J,K,L)
C     EVALUATE THE POLARIZED (QUARK)  ERT B FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQB,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=2*DOT(P,6,7)

      BDPQB=(-(S16*S23*S27*S34) + S14*S23*S27*S36+S14*S16*S23*S37 - 
     -    S14*S23*S26*S37 - S12*S16*S34*S37 + 2*S12*S14*S27*S46 + 
     -    S14*S23*S27*S46 + S12*S27*S34*S46 + S12*S14*S37*S46 + 
     -    S17*(S23*S26*S34 + S12*S34*S36 - S14*S23*(S36 + S46)) + 
     -    S14*S16*S23*S47 - 2*S12*S14*S26*S47 - S14*S23*S26*S47 - 
     -    S12*S26*S34*S47 - S12*S14*S36*S47 + 
     -    S13*(S17*(S23 + S24)*S26 - S16*(S23 + S24)*S27 - 
     -       S12*S27*S36 + S14*S27*S36 - S24*S27*S36 + S12*S26*S37 - 
     -       S14*S26*S37 + S24*S26*S37 + S12*S27*S46 + S14*S27*S46 - 
     -       S24*S27*S46 + S12*S37*S46 - S12*S26*S47 - S14*S26*S47 + 
     -       S24*S26*S47 - S12*S36*S47))/
     -  (2.*S13*S14*S23*(S13 + S14 + S34))

      BDPQB=BDPQB +
     -      (S12*(-(S14*S17*S26) - S17*S23*S26 - S17*S24*S26 + 
     -      S14*S16*S27 + S16*S23*S27 + S16*S24*S27 - S14*S17*S36 - 
     -      S17*S24*S36 + S14*S27*S36 + S24*S27*S36 + S14*S16*S37 + 
     -      S16*S24*S37 - S14*S26*S37 - S24*S26*S37 - S17*S23*S46 + 
     -      S23*S27*S46 + S16*S23*S47 - S23*S26*S47 + 
     -      S13*(S27*S46 - S17*(S26 + S46) - S26*S47 + 
     -         S16*(S27 + S47)) + 
     -      S12*(S27*S36 - S26*S37 + S27*S46 - 
     -         S17*(2*S26 + S36 + S46) - S26*S47 + 
     -         S16*(2*S27 + S37 + S47))))/(8.*S13*S14*S23*S24)

      BDPQB=BDPQB - 
     -     (S17*S26*S34 - S16*S27*S34 + 
     -     (S13 + S14)*(S27*(S36 + S46) - S26*(S37 + S47)))/
     -     (2.*S13*S14*(S13 + S14 + S34))

      BDPQB=BDPQB - 
     -       ((S12 + S14 + S24)*((S17*S26 - S16*S27)*
     -        (S13 + S14 + S23 + S24 + 2*S34) + 
     -       S12*(-(S27*S36) + S26*S37 - S27*S46 + 
     -          S17*(2*S26 + S36 + S46) + S26*S47 - 
     -          S16*(2*S27 + S37 + S47))))/
     -  (4.*S14*S24*(S13 + S14 + S34)*(S23 + S24 + S34))

      END
C-----------------------------------------------------------------------
      FUNCTION ERTGAXB(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTGAXB,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=2*DOT(P,6,7)

      ERTGAXB= -((S12**2*(2*S13*S23 + (S23 + 2*S24)*S34)*S36 - 
     -       S23*(S16*S24*(-S23 + S24)*S34 - S14**2*S24*S36 + 
     -          S14*S23*(S16*(S23 - 2*S24) + (S23 + S24)*S36)) - 
     -       S13**2*(S14*S23*S26 + 
     -          S24*(-(S24*S26) + S23*(2*S26 + S36))) + 
     -       S13*(S14*(S16*S23**2 + S24*S26*S34 + S23**2*(S26 + S36) - 
     -             S23*S24*(2*S26 + S46)) + 
     -          S24*(S16*S23*(2*S23 - S24) - S24**2*S36 + 
     -             S23**2*(S36 + S46) + S23*(-2*S26*S34 + S24*S46))) + 
     -       S12*(S16*S23*(S23 - 2*S24)*S34 + 
     -          (S23**2 + S24**2)*S34*S36 - 
     -          S14*(S23**2 + 2*S23*S24 + S24*S34)*S36 - 
     -          S13**2*S23*(2*S26 + S46) + 
     -          S13*(S16*S23*(2*S23 + S34) - 
     -             S24*(S26*S34 + 2*S24*S36) + S23**2*(2*S36 + S46) + 
     -             S23*(-2*S26*S34 + 2*S14*S36 + S24*S36 + S34*S36 + 
     -                3*S24*S46))))*S67)/
     -  (2.*S13*S14*S23*S24*(S13 + S14 + S34)*S36)

      ERTGAXB = ERTGAXB +
     -         ((S14**2*S16*S23 - S13**3*S26 + 
     -      S14*S34*(-(S12*S16) + 2*S16*S23 + (S23 + S24)*S36) - 
     -      S34**2*(S16*S24 + S12*(2*S16 + S36)) + 
     -      S13**2*(S16*S23 - 2*S14*S26 - 2*S26*S34 + S12*S36 + 
     -         2*S23*S36 + 2*S24*S36 + S12*S46) + 
     -      S13*(-(S14**2*S26) + 
     -         S14*(2*S16*S23 - 2*S26*S34 + S12*S36 + 2*S23*S36 + 
     -            2*S24*S36 + S12*S46) + 
     -         S34*(-(S12*S16) + 2*S16*S23 - S26*S34 + 2*S23*S36 + 
     -            2*S24*S36 + 2*S12*S46 + S23*S46 + S24*S46)))*S67)/
     -  (2.*S13*S14*(S13 + S14 + S34)**2*S36)

      ERTGAXB = ERTGAXB +
     -      ((-(S14*(S13 + S23)*S24*(-(S16*S23) + S13*S26)*
     -         (S13 + S14 + S23 + S24 + 2*S34)) + 
     -      S12**2*(S13 - S23)*(S13*S23 + S14*S24)*S36 + 
     -      S12*(-(S14*S23*S24*
     -            (-(S16*S23) + (S14 + S23 + S24 + 2*S34)*S36)) - 
     -         S13**3*S23*(S26 + S46) + 
     -         S13**2*(-(S23**2*S26) + S16*S23*(S23 + S34) + 
     -            S14*S24*(-S26 + S36) + 
     -            S23*(-3*S26*S34 + S14*S36 + S24*S36 - 2*S34*S46)) + 
     -         S13*(S16*S23*(S23**2 + S14*S24 + 3*S23*S34 + 2*S34**2) + 
     -            S14**2*S24*S36 + 
     -            S14*(-(S23*S24*S26) - S23**2*S36 + 
     -               S24*(S24 + 2*S34)*S36) + 
     -            S23*(-2*S26*S34**2 + S23**2*S46 - 
     -               S23*(S26*S34 + S24*S36 - 2*S34*S46)))))*S67)/
     -  (4.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34)*S36)

      ERTGAXB = ERTGAXB +
     -    ((S12 + S14 + S24)*(-(S13**2*S26) + S23*(S16*S23 - S12*S36) + 
     -      S13*(S16*S23 - S23*S26 + S12*S36))*S67)/
     -  (4.*S13*S14*S23*S24*S36)         

      ERTGAXB = -ERTGAXB

      END
C-----------------------------------------------------------------------
      FUNCTION ERTC(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT C FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTC,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      ERTC=-(5*S12*S34**2+2*S12*S24*S34+2*S12*S23*S34+2*S12*S14*S34+
     $     2*S12*S13*S34+4*S12**2*S34-S13*S24**2+S14*S23*S24+
     $     S13*S23*S24+S13*S14*S24-S12*S14*S24-S13**2*S24-
     $     3*S12*S13*S24-S14*S23**2-S14**2*S23+S13*S14*S23-
     $     3*S12*S14*S23-S12*S13*S23)/(4*S134*S234*S34**2)+
     $     (3*S12*S34**2-3*S13*S24*S34+3*S12*S24*S34+3*S14*S23*S34-
     $     S13*S24**2-S12*S23*S34+6*S12*S14*S34+2*S12*S13*S34-
     $     2*S12**2*S34+S14*S23*S24-3*S13*S23*S24-2*S13*S14*S24+
     $     4*S12*S14*S24+2*S12*S13*S24+3*S14*S23**2+2*S14**2*S23+
     $     2*S14**2*S12+2*S12**2*S14+6*S12*S14*S23-2*S12*S13**2-
     $     2*S12**2*S13)/(4*S13*S134*S234*S34)+
     $     (2*S12*S34**2-2*S13*S24*S34+S12*S24*S34+4*S13*S23*S34+
     $     4*S12*S14*S34+2*S12*S13*S34+2*S12**2*S34-S13*S24**2+
     $     3*S14*S23*S24+4*S13*S23*S24-2*S13*S14*S24+4*S12*S14*S24+
     $     2*S12*S13*S24+2*S14*S23**2+4*S13*S23**2+2*S13*S14*S23+
     $     2*S12*S14*S23+4*S12*S13*S23+2*S12*S14**2+4*S12**2*S13+
     $     4*S12*S13*S14+2*S12**2*S14)/(4*S13*S134*S24*S34)
      ERTC=ERTC-
     $     (S12*S34**2-2*S14*S24*S34-2*S13*S24*S34-S14*S23*S34+
     $     S13*S23*S34+S12*S14*S34+2*S12*S13*S34-2*S14**2*S24-
     $     4*S13*S14*S24-4*S13**2*S24-S14**2*S23-S13**2*S23+
     $     S12*S13*S14-S12*S13**2)/(2*S13*S34*S134**2)+
     $     (S12*S34**2-4*S14*S24*S34-2*S13*S24*S34-2*S14*S23*S34-
     $     4*S13*S23*S34-4*S12*S14*S34-4*S12*S13*S34-2*S13*S14*S24+
     $     2*S13**2*S24+2*S14**2*S23-2*S13*S14*S23-S12*S14**2-
     $     6*S12*S13*S14-S12*S13**2)/(4*S34**2*S134**2)
      END
C-----------------------------------------------------------------------
      FUNCTION BDPQC(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) PART OF ERT C FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQC,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=2*DOT(P,6,7)

      BDPQC=(-(S14**2*S17*S26) - S14*S17*S23*S26-2*S14*S17*S24*S26 + 
     -    S14**2*S16*S27 + S14*S16*S23*S27 + 2*S14*S16*S24*S27 - 
     -    S14*S17*S26*S34 + S17*S23*S26*S34 - S17*S24*S26*S34 + 
     -    S14*S16*S27*S34 - S16*S23*S27*S34 + S16*S24*S27*S34 - 
     -    2*S14*S17*S24*S36 + S14**2*S27*S36 + 2*S14*S23*S27*S36 + 
     -    S14*S24*S27*S36 - 2*S17*S24*S34*S36 + S14*S27*S34*S36 + 
     -    2*S14*S16*S24*S37 - S14**2*S26*S37 - 2*S14*S23*S26*S37 - 
     -    S14*S24*S26*S37 + 2*S16*S24*S34*S37 - S14*S26*S34*S37 + 
     -    S14*S17*S23*S46 + 2*S14*S23*S27*S46 + S17*S23*S34*S46 - 
     -    S14*S16*S23*S47 - 2*S14*S23*S26*S47 - S16*S23*S34*S47 + 
     -    S13*(2*S16*S23*S27 + S16*S24*S27 + S16*S27*S34 + 
     -       4*S23*S27*S36 + 2*S24*S27*S36 + 2*S27*S34*S36 - 
     -       4*S23*S26*S37 - 2*S24*S26*S37 - 2*S26*S34*S37 + 
     -       2*S23*S27*S46 - S24*S27*S46 - 2*S23*S37*S46 - 
     -       S17*(2*S23*S26 + S24*S26 + S26*S34 + S23*S46 - 
     -          2*S24*S46) + S16*S23*S47 - 2*S16*S24*S47 - 
     -       2*S23*S26*S47 + S24*S26*S47 + 2*S23*S36*S47 + 
     -       S14*(-2*S17*S26 + 2*S16*S27 + S27*S36 - S26*S37 - 
     -          2*S27*S46 + 2*S26*S47)) + 
     -    S12*(S14*(2*S16*S27 + S27*S36 - S26*S37 + 2*S27*S46 - 
     -          S17*(2*S26 + S46) + S16*S47 - 2*S26*S47) + 
     -       S13*(4*S16*S27 + 2*S27*S36 - 2*S26*S37 + S27*S46 - 
     -          S37*S46 - 2*S17*(2*S26 + S46) + 2*S16*S47 - S26*S47 + 
     -          S36*S47) + S34*
     -        (-2*S17*S26 + 2*S16*S27 + S27*S36 - S26*S37 - 3*S17*S46 - 
     -          2*S37*S46 + 3*S16*S47 + 2*S36*S47)))/
     -  (8.*S13*S24*S34*(S13 + S14 + S34))

      BDPQC=BDPQC+
     -   (S13*S14*(S17*S26 - S16*S27 + 4*S27*S46 - 4*S26*S47) + 
     -    S13**2*(-(S17*S26) + S16*S27 + S27*S36 - S26*S37 + 
     -       4*S27*S46 - 4*S26*S47) + 
     -    S13*S34*(2*S17*S26 - 2*S16*S27 - S27*S36 + S26*S37 + 
     -       2*S27*S46 - 2*S26*S47) - 
     -    (S14 + S34)*(-(S17*S26*S34) + S16*S27*S34 - 
     -       S14*S27*(S36 + 2*S46) + S14*S26*(S37 + 2*S47)))/
     -  (4.*S13*S34*(S13 + S14 + S34)**2)

      BDPQC=BDPQC+
     -    (-(S14**2*S17*S26)-3*S14*S17*S23*S26-2*S14*S17*S24*S26 + 
     -    S14**2*S16*S27 + 3*S14*S16*S23*S27 + 2*S14*S16*S24*S27 + 
     -    S13**2*(S17*S26 - S16*S27) - 4*S14*S17*S26*S34 - 
     -    S17*S24*S26*S34 + 4*S14*S16*S27*S34 + S16*S24*S27*S34 - 
     -    2*S17*S26*S34**2 + 2*S16*S27*S34**2 + 2*S14*S17*S24*S36 + 
     -    3*S14*S23*S27*S36 + S14*S24*S27*S36 + 2*S17*S24*S34*S36 + 
     -    S14*S27*S34*S36 - 2*S14*S16*S24*S37 - 3*S14*S23*S26*S37 - 
     -    S14*S24*S26*S37 - 2*S16*S24*S34*S37 - S14*S26*S34*S37 - 
     -    2*S14*S17*S23*S46 - 2*S14*S23*S37*S46 + 2*S14*S16*S23*S47 + 
     -    2*S14*S23*S36*S47 + 
     -    S13*(S17*(S23 - 2*S24)*S26 - S16*(S23 - 2*S24)*S27 - 
     -       3*S23*S27*S46 - S24*S27*S46 - 3*S27*S34*S46 + 
     -       2*S24*S37*S46 + 3*S23*S26*S47 + S24*S26*S47 + 
     -       3*S26*S34*S47 - 2*S24*S36*S47) + 
     -    S12*(S34*(-(S27*S36) + S26*S37 + 2*S27*S46 + 3*S37*S46 + 
     -          S17*(2*S26 - S36 + 3*S46) + 
     -          S16*(-2*S27 + S37 - 3*S47) - 2*S26*S47 - 3*S36*S47) + 
     -       S13*(S27*S36 - S26*S37 - S37*S46 + 
     -          S17*(2*S26 + S36 + S46) + S36*S47 - 
     -          S16*(2*S27 + S37 + S47)) + 
     -       S14*(3*S27*S36 - 3*S26*S37 + 2*S27*S46 - S37*S46 - 
     -          S17*(2*S26 + S36 + S46) - 2*S26*S47 + S36*S47 + 
     -          S16*(2*S27 + S37 + S47))))/
     -  (8.*S13*S34*(S13 + S14 + S34)*(S23 + S24 + S34))

      BDPQC=BDPQC+
     -    ((-(S17*S26) + S16*S27)*S34**2 + 
     -    S14**2*(S17*S26 - S16*S27 + 2*S27*S36 - 2*S26*S37) + 
     -    S13**2*(S17*S26 - S16*S27 + 2*S27*S46 - 2*S26*S47) + 
     -    2*S14*S34*(2*S17*S26 - 2*S16*S27 - S27*S36 + S26*S37 - 
     -       2*S27*S46 + 2*S26*S47) + 
     -    2*S13*(S14*(3*S17*S26 - 3*S16*S27 - S27*S36 + S26*S37 - 
     -          S27*S46 + S26*S47) + 
     -       S34*(2*S17*S26 - 2*S16*S27 - 2*S27*S36 + 2*S26*S37 - 
     -          S27*S46 + S26*S47)))/(8.*S34**2*(S13 + S14 + S34)**2)

      BDPQC=BDPQC+
     -      (S13*(-(S17*(S24*(3*S26 + S36) + S23*(S26 - S46))) - 
     -       S23*S27*S46 + S24*S27*S46 - S27*S34*S46 + 2*S24*S37*S46 + 
     -       S23*S26*S47 - S24*S26*S47 + S26*S34*S47 - 2*S24*S36*S47 + 
     -       S16*(S23*S27 + 3*S24*S27 + S24*S37 - S23*S47)) + 
     -    S14*(S23*S27*S36 - S24*S27*S36 - S27*S34*S36 - S23*S26*S37 + 
     -       S24*S26*S37 + S26*S34*S37 - 2*S23*S37*S46 - 
     -       S17*(S24*(S26 - S36) + S23*(3*S26 + S46)) + 
     -       2*S23*S36*S47 + 
     -       S16*(3*S23*S27 + S24*S27 - S24*S37 + S23*S47)) + 
     -    S34*(S17*S26*S34 - S16*S27*S34 + S17*S24*S36 - S16*S24*S37 + 
     -       S17*S23*S46 - S16*S23*S47 + 
     -       2*S12*(-(S27*S36) + S26*S37 - S27*S46 + 
     -          S17*(2*S26 + S36 + S46) + S26*S47 - 
     -          S16*(2*S27 + S37 + S47))))/
     -  (8.*S34**2*(S13 + S14 + S34)*(S23 + S24 + S34))

      END
C-----------------------------------------------------------------------
      FUNCTION ERTGAXC(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTGAXC,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=2*DOT(P,6,7)

      ERTGAXC= ((-(S13**4*S24*S26) + 
     -      S14*S23*(S14**2*S23*(S16 - S36) + 
     -         S34*(S12*S16*(2*S23 + S24 - S34) + 
     -            S16*S23*(S23 + 2*S24 + S34) + 2*S12**2*S36 + 
     -            S12*S23*S36) - 
     -         S14*(S16*(S23**2 + S23*S24 + S12*S34 - 2*S23*S34) + 
     -            S23*(3*S12 + 2*S23 + 2*S24 + S34)*S36)) + 
     -      S13**3*S24*(S16*S23 - 2*S14*S26 - 2*S23*S26 + S24*S26 - 
     -         2*S26*S34 + S12*S36 - S24*S36 + S23*S46) + 
     -      S13**2*(S14**2*S23*S26 + 
     -         S14*(2*S16*S23*S24 + 
     -            S24*(-2*S26*S34 + (2*S12 + S24)*S36) + 
     -            S23**2*(3*S26 + S46) - 
     -            S23*(-2*S26*S34 + S24*(S26 + 2*S36 - S46) + S12*S46))
     -          - S24*(S16*S23*(-2*S23 + S24 - S34) + 2*S23*S26*S34 + 
     -            S24*S26*S34 + S26*S34**2 + 2*S23*S24*S36 + 
     -            2*S24**2*S36 + S24*S34*S36 - 2*S23**2*S46 - 
     -            2*S23*S24*S46 - S23*S34*S46 + 
     -            S12*(S26*S34 + 2*S23*S36 + 4*S24*S36 - 3*S34*S36 - 
     -               5*S23*S46))) +S13*(-(S14**3*S23*S26) + 
     -         S14**2*S23*(-(S16*S23) + S24*S26 - 3*S26*S34 + 
     -            2*S24*S36 + S12*S46 + S23*(S26 - 3*S36 + S46)) + 
     -         S24*S34*(-(S16*S23*S24) + 4*S12**2*S36 + 
     -            S12*(-4*S16*S23 + S26*S34 + 2*S24*S36 + 2*S34*S36 - 
     -           S23*S46)) +S14*(-(S23**2*S26*S34) - 3*S23*S24*S26*S34 - 
     -            2*S23*S26*S34**2 - 2*S24*S26*S34**2 + 
     -            S16*S23*(-3*S23**2 + S23*(S24 - 3*S34) + 
     -               (S12 + S24)*S34) - 4*S23**3*S36 - 
     -            4*S23**2*S24*S36 - S23*S24**2*S36 - 
     -            4*S23**2*S34*S36 + 2*S24**2*S34*S36 + 
     -            S23**2*S24*S46 + S23**2*S34*S46 + 2*S23*S24*S34*S46 + 
     -            S12*(S24*S34*S36 + S23**2*(-5*S36 + 2*S46) - 
     -               S23*(4*S26*S34 + 7*S24*S36 + 2*S34*S36 + S24*S46 - 
     -                  S34*S46)))))*S67)/
     -  (8.*S13*S14*S23*S24*S34*(S13 + S14 + S34)*S36)

      ERTGAXC = ERTGAXC +
     -      ((3*S13**4*S26 + S14*
     -       (S14**2*S23*S36 - 
     -         S14*S34*(2*S16*S23 + (S12 + S23 + 2*S24)*S36) + 
     -         S34**2*(2*S16*S24 + S12*(2*S16 + S36))) - 
     -      S13**3*(3*S16*S23 - 3*S14*S26 - 6*S26*S34 + 3*S12*S36 + 
     -         4*S23*S36 + 4*S24*S36 + 3*S12*S46 + 2*S23*S46 + S24*S46)
     -       + S13*(S16*(S23 + 3*S24)*S34**2 + 
     -         S14**2*(4*S26*S34 + 5*S23*S36 - 4*S24*S36 - 2*S23*S46) - 
     -         S14*S34*(7*S16*S23 + 3*S16*S24 - 4*S26*S34 + 2*S23*S36 + 
     -            10*S24*S36 + 4*S23*S46 + 2*S24*S46) + 
     -         S12*(3*S14**2*S36 + S34**2*(5*S16 + S36) - 
     -            S14*S34*(S16 + S36 + 2*S46))) - 
     -      S13**2*(S14*(3*S16*S23 - 7*S26*S34 + 4*S12*S36 + 
     -            4*S23*S36 + 12*S24*S36 - S12*S46 - 3*S24*S46) + 
     -         S34*(4*S16*S23 - S16*S24 - 3*S26*S34 + 4*S23*S36 + 
     -            4*S24*S36 + 4*S23*S46 + 3*S24*S46 + 
     -            S12*(-3*S16 + 2*S36 + 5*S46))))*S67)/
     -  (8.*S13*S14*S34*(S13 + S14 + S34)**2*S36)

      ERTGAXC = ERTGAXC -
     -    ((2*S13**4*S26 + S14*
     -        (-(S14*S16*(S12*S34 + S23*(S24 + S34))) + 
     -          S14*(S12*(3*S23 - S34) + S23*(3*S23 + S24 + 3*S34))*
     -           S36 + S14**2*S23*(S16 + 2*S36) - 
     -          S34*(S16*S23*(S23 + 2*S24 + 3*S34) + 2*S12**2*S36 + 
     -             S12*(S16*S23 - S16*S24 + 2*S23*S36 - S24*S36))) + 
     -       S13**3*(-2*S16*S23 + S14*S26 + S23*S26 + 3*S24*S26 + 
     -          5*S26*S34 - 2*S12*S36 - 2*S24*S36 - 2*S12*S46) + 
     -       S13**2*(-2*S14**2*S26 - 4*S12*S26*S34 + S23*S26*S34 + 
     -          3*S24*S26*S34 + 3*S26*S34**2 - 
     -          S16*(S23**2 + 3*S23*S24 - 2*S12*S34 + 5*S23*S34) + 
     -          3*S12*S23*S36 - S12*S24*S36 - 2*S24**2*S36 + 
     -          S12*S34*S36 - 2*S24*S34*S36 + 3*S12*S23*S46 + 
     -          S12*S24*S46 + S23*S24*S46 + S24**2*S46 - 
     -          3*S12*S34*S46 + S24*S34*S46 + 
     -          S14*(-(S16*S23) - S23*S26 + 2*S24*S26 + 2*S26*S34 + 
     -             S12*S36 + 2*S23*S36 + S12*S46)) - 
     -       S13*(S14**3*S26 - 
     -          S14**2*(2*S16*S23 + S26*S34 + S24*(S26 - 2*S36) + 
     -             S12*S36 + S12*S46) + 
     -          S34*(S16*(S23**2 + S24*(S24 + S34) + 
     -                S23*(4*S24 + 3*S34)) - 2*S12**2*S36 - 
     -             S12*(-4*S26*S34 + S16*(S23 - S24 + 3*S34) + 
     -                3*S23*S36 + 3*S24*S36 + 3*S34*S36 + 4*S23*S46)) + 
     -          S14*(-2*S23*S26*S34 - 3*S24*S26*S34 - 4*S26*S34**2 + 
     -             S16*(-S23**2 + S12*S34 + 2*S23*(S24 + S34)) + 
     -             S23*S24*S36 + S24**2*S36 - 2*S23*S34*S36 + 
     -             3*S24*S34*S36 + S23**2*S46 + S23*S24*S46 + 
     -             S23*S34*S46 + 
     -             S12*(-2*S26*S34 + 2*S23*S36 - S24*S36 - 3*S34*S36 + 
     -                S23*S46 + S24*S46))))*S67)/
     -  (8.*S13*S14*S34*(S13 + S14 + S34)*(S23 + S24 + S34)*S36)        

      ERTGAXC = ERTGAXC +
     -  ((S16*(4*S12 + S23 + 3*S24)*S34**2 + 
     -      S14**2*(S26*S34 + 7*S23*S36 - S23*S46) - 
     -      S14*S34*(S16*(3*S23 + S24) - S26*S34 + 8*S12*S36 + 
     -         S23*S36 + 8*S24*S36 + S23*S46) + 
     -      S13**2*(3*S26*S34 + S24*(-7*S36 + S46)) - 
     -      S13*(S16*(3*S23 + S24)*S34 + 
     -         S14*(-4*S26*S34 - 7*S23*S36 + 7*S24*S36 + S23*S46 - 
     -            S24*S46) + 
     -         S34*(-3*S26*S34 - 4*S12*S36 - 4*S23*S36 + 3*S24*S36 + 
     -            4*S12*S46 + 4*S23*S46 + 3*S24*S46)))*S67)/
     -  (8.*S34**2*(S13 + S14 + S34)**2*S36)

      ERTGAXC = ERTGAXC -
     - ((S14*(-(S16*(3*S23 + S24)*S34) + S23*S26*S34 + 
     -          S34*(S24*S26 + 2*S26*S34 - 8*S12*S36) + 
     -          S23**2*(7*S36 - S46) + S23*(S24 + 2*S34)*(7*S36 - S46))
     -        + S14**2*(S26*S34 + 7*S23*S36 - S23*S46) + 
     -       S13*(-(S16*(3*S23 + S24)*S34) + 4*S14*S26*S34 + 
     -          3*S23*S26*S34 + 3*S24*S26*S34 + 6*S26*S34**2 - 
     -          7*S23*S24*S36 - 7*S24**2*S36 + 4*S12*S34*S36 - 
     -          14*S24*S34*S36 + S14*(S23 - S24)*(7*S36 - S46) + 
     -          S23*S24*S46 + S24**2*S46 - 4*S12*S34*S46 + 2*S24*S34*S46
     -          ) + S13**2*(3*S26*S34 + S24*(-7*S36 + S46)) - 
     -       S34*(S16*(3*S23**2 + 4*S23*S24 + S24**2 - 4*S12*S34 + 
     -             6*S23*S34 + 2*S24*S34) + 
     -          4*S12*(S26*S34 + S23*S36 - 2*S24*S36 - S23*S46)))*S67)/
     -  (16.*S34**2*(S13 + S14 + S34)*(S23 + S24 + S34)*S36) 

      ERTGAXC = -ERTGAXC

      END
C-----------------------------------------------------------------------
      FUNCTION BDPG(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE BDP FUNCTION FOR THE WHOLE GLUON CHANNEL
      INTEGER I,J,K,L,SCHEME,NF
      DOUBLE PRECISION BDPG,GCF,GCA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S124,S123,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S123=S12+S13-S23
      S124=S12+S14-S24
      S67=S12+S13+S14-S23-S24-S34

      GCF=(2*(S36 + S46))/(S13*S14) - 
     -  ((S34 - S36 + S37)*(2*S23 + S34 - S36 + S37)*(S24 + S34 + S67)*
     -     (2*S46 + S67))/(S123**2*S13*S14*S23) - 
     -  ((3*S123 - 3*S24 + 4*S26 - 6*S34)*S34 + 2*S124*S36 + 
     -     (S123 - 2*S24 + 2*S26 - 4*S34)*S67)/(S13*S14*S23) + 
     -  (S34*(S123**2 + S34*(-4*S26 + 3*S34 - 4*S36) + 
     -       S123*(2*S26 - 3*S34 + 2*S36 - S67) - 
     -       2*(S26 - S34 + S36)*S67))/(S13*S14*S23*S24) + 
     -  ((S34 + S67)*(S23**2*S24*(3*S34 + 2*S36 - 2*S46 - S67) + 
     -       (S24**3 + S24**2*(-2*S26 + 3*S34 - 2*S36) + 
     -          S34**2*(-2*S26 + S34 - S36 - S46) - 
     -          2*S24*S34*(2*S26 - 2*S34 + S36 + S46))*(S34 + S67) + 
     -       S23*S24*(S34*(-2*S26 + 3*S34 - 2*S36) + 
     -          (2*S26 + S34 + 2*S36)*S67)))/(S123*S124*S13*S14*S23*S24)
      GCF=GCF-((S24**3*(S34 + S67) + 
     -      S23**2*(3*S34**2 + 5*S34*S67 + 2*S67**2 + 
     -         S24*(3*S34 + S67)) + 
     -      S24**2*(3*S34**2 + 3*S34*S67 - 2*S26*(S34 + S67) - 
     -         2*S36*(S36 - S37 + S67)) + 
     -      S34*(S34 + S67)*(4*S34**2 - 2*S26*(3*S34 - S36 + S37) + 
     -         (S36 - S37)*(S36 - S37 + 2*S46 + S67) - 
     -         S34*(5*S36 - 3*S37 + 4*S46 + S67)) - 
     -      S24*(-6*S34**3 + S34**2*(3*S36 - 3*S37 + 8*S46 - 2*S67) + 
     -         2*S26*S34*(3*S34 - S36 + S37 + 2*S67) + 
     -         (S36 - S37)*(2*S36*S46 - 2*S37*S46 + 3*S36*S67 - 
     -            S37*S67) + 
     -         S34*(S36**2 - S37**2 - 8*S36*S46 + 8*S37*S46 - 
     -            3*S36*S67 + 3*S37*S67 + 2*S46*S67)) + 
     -      S23*(4*S24**2*S36 - 
     -         (S34 + S67)*(-6*S34**2 + 
     -            S34*(3*S36 - 3*S37 + 2*S46 - S67) + 
     -            2*S26*(2*S34 + S67) + S67*(S36 - S37 + 2*S46 + S67))
     -          - S24*(4*S26*S34 - 6*S34**2 - 4*S36*S46 + 4*S37*S46 - 
     -            7*S36*S67 + S37*S67 + 2*S46*S67 + S67**2 + 
     -            S34*(-S36 - 3*S37 + 8*S46 + S67))))/
     -    (S123*S13*S14*S23*S24))
      GCF=8*GCF
     -     -4*S34*(S36**2 - S37**2 + S46**2 - S47**2)/(S13*S14*S23*S24)

      GCA=8/S12 + 4/S13 + 20/S23 - (4*(2*S14 - S23 - 2*S34))/(S13*S24)- 
     -  (8*S34)/(S13*S14) - (4*S34*(S12 + S24 - 2*(S26 - 2*S34 + S36)))/
     -   (S13*S14*S23) - (4*(S14 - 2*S23 + S24 - 2*S46))/(S12*S13) - 
     -  (4*S34**2*(2*S12 - 2*S26 + S34 - S36 - S46))/
     -   (S13*S14*S23*S24) + 
     -  (8*(S12 + S13)*(S12 + 2*S13)*(S14 - S24 - S34 + 2*S46))/
     -   (S12*S123**2*S23) + 
     -  (4*(S12 - S14 - S24 + 3*S34 - 2*S36 + 2*S46))/(S13*S23) + 
     -  (8*(2*S13 + S14 - S24 - S26 + 3*S34 - S36 + 3*S46))/(S12*S23) - 
     -  (4*(S14**2 - S23**2 + 2*S26*S34 - 2*S14*(S26 + S46) + 
     -       2*S23*(S26 - S34 + S46)))/(S12*S13*S24) + 
     -  (2*(4*S12 + S13 + S14 - 2*(2*S26 - 6*S34 + S36 + S46)))/
     -   (S23*S24) - (4*(S13**2 + S14**2 - 
     -       S13*(2*S14 + 4*S34 + S36 - S46) - 
     -       S14*(4*S34 - S36 + S46) + 2*S34*(3*S26 - S34 + S36 + S46)))
     -    /(S12*S23*S24) - (4*
     -     (S14**2 + 2*S12*(S14 - 2*S34) + 2*S34*(2*S26 + S36 + S46) - 
     -       S14*(2*S26 + S34 + 2*S46)))/(S13*S23*S24)
      GCA=GCA+(-2*(S12**2 + 2*S13*S14 + S12*(S13 + S14))*
     -    (2*S12**2*((S13 - S14)**2 - 2*(S13 + S14)*S34 + 2*S34**2) - 
     -      2*S13*S14*(S13*(2*S34 + S36 - S46) + 
     -         S14*(2*S34 - S36 + S46) - 2*S34*(2*S26 - S34 + S36 + S46)
     -         ) + S12*(S13**3 + S14**3 - 
     -         S13**2*(S14 + 2*S26 + S34 + 2*S36) + 
     -         2*S14*S34*(2*S26 + S36 + S46) - 
     -         2*S34**2*(2*S26 - S34 + S36 + S46) - 
     -         S14**2*(2*S26 + S34 + 2*S46) + 
     -         S13*(-S14**2 + 2*S34*(2*S26 + S36 + S46) + 
     -            2*S14*(2*S26 - 7*S34 + S36 + S46)))))/
     -  (S12*S123*S124*S13*S14*S23*S24)
      GCA=GCA+(4*(S12**3*(2*S13**2 + 2*S14**2 - 4*S13*(S14 + S34) - 
     -         S14*(S24 + 4*S34) + S34*(S24 + 4*S34)) + 
     -      S12**2*(3*S13**3 + S14**3 - 
     -         S13**2*(5*S14 - 2*S24 + 2*S26 + 5*S34 + 2*S36) + 
     -         S14**2*(S24 - 2*S26 - S34 - 2*S46) + 
     -         S14*(S24**2 + S24*(-3*S34 + 2*S36 - 2*S46) + 
     -            2*S34*(2*S26 + S36 + S46)) + 
     -         S34*(S24**2 - 2*S24*(S26 - 2*S34 + S36) - 
     -            2*S34*(2*S26 - S34 + S36 + S46)) + 
     -         S13*(S14**2 - S24**2 - S24*S34 + 
     -            2*S14*(-3*S24 + 2*S26 - 9*S34 + S36 + S46) + 
     -            2*S34*(2*S26 + 2*S34 + S36 + S46))) + 
     -      S12*S13*(S13**3 + 3*S14**3 - S24**3 + 2*S24**2*S26 - 
     -         S24**2*S34 + 4*S24*S34**2 - 4*S26*S34**2 + 2*S34**3 + 
     -         S13**2*(S14 + 3*S24 - 2*S26 - S34 - 2*S36) + 
     -         2*S24**2*S36 - 2*S24*S34*S36 - 2*S34**2*S36 - 
     -         2*S34**2*S46 - 
     -         S14**2*(3*S24 + 2*S26 + 9*S34 - 2*S36 + 4*S46) - 
     -         S13*(5*S14**2 + S24**2 + 2*S24*(S26 + S34 + S36) + 
     -            S14*(7*S24 - 4*S26 + 22*S34 - 4*S46) - 
     -            2*S34*(2*S26 + S36 + S46)) + 
     -         S14*(5*S24**2 + S24*(2*S26 - 5*S34 + 6*S36 - 12*S46) + 
     -            2*S34*(6*S26 + 3*S36 + 2*S46))) + 
     -      S13**2*(S14**3 + S13**2*(S14 + S24) + 
     -         S24*(-S24**2 + 2*S26*S34 + 2*S24*(S26 - S34 + S36)) - 
     -         2*S13*(S14**2 + S24*(S26 + S36) + 
     -            S14*(S24 + S26 + 3*S34 + 2*S36 - S46)) - 
     -         2*S14**2*(3*S24 - S26 + 5*S34 - 2*S36 + S46) + 
     -         2*S14*(2*S24**2 + S24*(2*S26 - 3*S34 + 2*S36 - 6*S46) + 
     -            S34*(7*S26 + 3*(-S34 + S36 + S46))))))/
     -  (S12*S123*S13*S14*S23*S24)
      GCA=GCA-(2*(S14*S23 + S13*S24 - S12*S34)*
     -    (S36**2 - S37**2 + S46**2 - S47**2))/(S12*S13*S14*S23*S24)

      BDPG=(CA*GCA+CF*GCF)*S67
      END
C-----------------------------------------------------------------------
      FUNCTION BDPGAXP(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE BDP FUNCTION FOR THE WHOLE GLUON CHANNEL
C---This needs to be called in an atisimetric way between the quarks
C--- in order to get the correct full-contribution!!!
      INTEGER I,J,K,L,SCHEME,NF
      DOUBLE PRECISION BDPGAXP,GCF,GCA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S124,S123,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S123=S12+S13-S23
      S124=S12+S14-S24
      S67=S12+S13+S14-S23-S24-S34

C=============================CF PART===================================     
      GCF = (-2*S37*(S34 + S67))/(S23*S24) + 
     -  (2*S34*S37*(2*S26 + S34 + S37 + S47)*S67*(S34 + S67)**2)/
     -   (S123*S124*S13*S14*S23*S24) + 
     -  ((S34 + S67)**2*(2*(S14 - S34)*S37*S47 + 
     -       (S14**2 - S14*(2*S26 + 3*S34 + 4*S47) + 
     -          2*(S26**2 + S26*S34 + S34**2 - S26*S37 - S34*S37 - 
     -             S37**2 + (3*(S26 + S34) + S37)*S47 + S47**2))*S67))/
     -   (S123*S124*S13*S23*S24) - 
     -  ((S34 + S67)*(2*S13*S47*(S34 + S47) - 
     -       2*S34*S47*(S26 + S34 + S47) - S13*(S34 + 2*S47)*S67 + 
     -       (S124*S34 + 2*(S26 + S47)*(2*S34 + S47))*S67 + 
     -       (S124 + 2*(S26 + S37))*S67**2))/(S123*S14*S23*S24) + 
     -  (S34*S67*(-S123**2 + 2*S37*(2*S26 + S37 + S47 - S67) + 
     -       S123*(-2*S26 + S34 - 2*S47 + S67)))/(S13*S14*S23*S24)

      GCF = GCF + (2*S34*S37*(S26 + S34 + S37) - 
     -     (S123*S34 + 2*((S124 + 2*S26 - S34)*S34 + 
     -           (S124 + S26 + S34)*S37 + S37**2))*S67 - 
     -     (S123 + 2*(S26 - S34 - S37 + S47))*S67**2 + 
     -     S14*(-2*S37*(S34 + S37) + (S34 + 2*S37)*S67))/(S13*S23*S24)
     -   + (2*(S124 - S34 - S67)*(-S14 + S34 + S67)*
     -      (2*S47**2 - 2*S47*S67 + S67**2) + 
     -     S123*(2*S47*(S124*S34 + 2*S14*S34 - 2*S26*S34 - 4*S34**2 + 
     -           4*S14*S37 - 3*S34*S37 - S34*S47 + 2*S13*(S34 + S47)) + 
     -        (-3*S13*S34 - 6*S14*S34 + 8*S26*S34 + 6*S34**2 - 
     -           4*S14*S37 + 6*S34*S37 - 
     -           2*(2*(S13 + S14 - S26) + S37)*S47 + 2*S47**2 + 
     -           S124*(S34 + 2*S47))*S67 - 
     -        (S124 + S13 + 2*(S14 - 2*(S26 + 2*S34 + S37) - S47))*
     -         S67**2 + 4*S67**3))/(S123**2*S23*S24) + 
     -  ((S34 + S67)*(4*S37*S67*(2*S26 + 2*S34 + S37 + S47 + S67) + 
     -       S13*(S37*(4*S47 - 6*S67) + S67*(-3*S34 + 2*S47 + S67))))/
     -   (S123*S124*S23*S24) + 
     -  (S34*S67*(S34 + S67)*
     -     (-S124**2 + S124*(-2*S26 + S34 - 2*S37 + S67) + 
     -       2*(S47*(S34 + S47) + S37*(-S37 + S67) + 
     -          S26*(S34 - 2*S37 + 2*S47 + S67))))/
     -   (S123*S13*S14*S23*S24)

      GCF = GCF + ((-S124 + S34 + S67)**2*(-S14 + S34 + S67)*
     -     (2*S47**2 - 2*S47*S67 + S67**2) + 
     -    S123*(S124**2*(S34 + 2*S47 - S67)*S67 - 
     -       (S34 + S67)*(2*S14*(S34 + 3*S37)*S47 - 
     -          2*S34*(S26 + S34 + 3*S37)*S47 + S14**2*S67 - 
     -          S14*(2*S26 + 5*S34 + 2*S37 + 6*S47)*S67 + 
     -          2*(S26**2 + 2*S26*S34 + 2*S34**2 - S26*S37 + S34*S37 - 
     -             S37**2 + (4*S26 + 2*S34 + S37)*S47 + S47**2)*S67 + 
     -          (3*S34 + 2*S37 - 2*S47)*S67**2 + S67**3) + 
     -       S124*(2*S14*(S34 + 2*S37)*S47 - 
     -          2*S34*(S26 + S34 + 2*S37)*S47 - 
     -          2*S14*(S34 + S37 + S47)*S67 + 
     -          (S34*(S34 + 4*S37 - 4*S47) + 2*S26*(S34 + S47))*S67 + 
     -          (3*S34 + 2*S37 - 4*S47)*S67**2 + 2*S67**3)))/
     -  (S123**2*S13*S23*S24)

      GCF=8*(GCF + (S34*(S36**2 + S37**2)*S67)/(S13*S14*S23*S24))

c     (this last piece is the only simmetric part in 1-2, the gluons,
c      the rest is antisymmetric. Same thing happens in CA. We could 
c      use that to reduce the expression even more...)

C=============================CA PART===================================
      GCA = (2*S34*S37*(2*S26 + S34 + S37 + S47)*S67*(S34 + S67)**3)/
     -   (S12*S123*S124*S13*S14*S23*S24) + 
     -  (2*(S124 + S13 - S34 - S67)*(S124 + 2*S13 - S34 - S67)*
     -     (-S14 + S34 + S67)*(2*S47**2 - 2*S47*S67 + S67**2))/
     -   (S12*S123**2*S23*S24) + 
     -  (2*(S124*(-(S34*(2*(S26 + S34) + S37)) + S13*(S34 - 2*S47) + 
     -           S14*(2*S34 + S37 - 2*S47) + 2*(S26 + S37)*S47 + S47**2)
     -          + S37*(4*S14**2 + S34*(2*S26 + S34 + S37 + S47)) + 
     -        S13*(S34**2 + 4*S14*S37 - 2*S34*(S37 - 3*S47) + 
     -           S47*(4*S26 + 8*S37 + S47))) + 
     -     (6*S124**2 + 4*S13**2 + 
     -        S124*(8*S13 - 2*S14 + 6*S26 - 13*S34 + 6*S37 - 2*S47) + 
     -        S13*(6*S26 - 13*S34 - 10*S37 + 12*S47) + 
     -        2*S37*(10*S26 + 3*S34 + 5*(S37 + S47)))*S67 + 
     -     2*(-3*S124 + 2*(S14 + S47))*S67**2)/(S12*S23*S24) + 
     -  ((S34 + S67)**2*(-2*(S14 - S34)**2*S37*S47 + 
     -       (-S14**3 + 2*S14**2*(S26 + 2*(S34 + S47)) + 
     -          2*S34*(S26**2 + S34**2 - 2*S34*S37 - 2*S37**2 + 
     -             4*S34*S47 + 2*S47**2 + S26*(S34 - 3*S37 + 5*S47)) - 
     -          S14*(2*S26**2 + 5*S34**2 - 2*S34*S37 - 2*S37**2 + 
     -             10*S34*S47 + 2*S47**2 + S26*(4*S34 - 2*S37 + 6*S47)))
     -         *S67 + (S14**2 - S14*(2*S26 + 3*S34 + 4*S47) + 
     -          2*(S26**2 + S26*S34 + S34**2 - S26*S37 - S34*S37 - 
     -             S37**2 + (3*(S26 + S34) + S37)*S47 + S47**2))*S67**2)
     -     )/(S12*S123*S124*S13*S23*S24)

      GCA = GCA + ((S34 + S67)*(2*S47*(S34*(-S124 - S13 + S34)*
     -           (-S13 + S26 + S34) + (S13 - S34)*S34*S37 + 
     -          (S13 - S34)*(S124 + S13 - S34)*S47) + 
     -       2*(S34*(S124**2 + (S26 - S34)*(S13 + S26 - S34) + 
     -             (-2*S13 + 5*S26 + 3*S34)*S37 + 2*S37**2 + 
     -             S124*(3*S26 - S34 + S37)) + 
     -          (-S13**2 - 3*S34*(S26 + S34) + 
     -             S124*(-S13 + S26 + 2*S34) + S13*(S26 + 2*S34 + S37))*
     -           S47 + (S124 - 2*S34)*S47**2)*S67 + 
     -       (S124**2 + S13**2 + S124*(S13 + 2*S26 - 3*S34 + 2*S37) - 
     -          2*S13*(S34 + S37 - S47) + 
     -          2*(S26**2 - 3*S26*S34 + S34**2 + 3*S26*S37 + S34*S37 + 
     -             S37**2 - 2*S26*S47 - 3*S34*S47 + S37*S47 - 2*S47**2))
     -         *S67**2 - (S124 + 2*(S26 + S37))*S67**3))/
     -   (S12*S123*S14*S23*S24) + 
     -  (S34*S67*(S123**3 + S123**2*
     -        (S124 + 2*S26 - 3*S34 + 2*S47 - 3*S67) + 
     -       2*S37*(6*S26 + S34 + 3*(S37 + S47) - 2*S67)*(S34 + S67) - 
     -       S123*(2*S124*S37 + 
     -          2*(S34 + S37 - S47)*(2*S26 - S34 + S37 + S47) + 
     -          2*(2*S26 - 2*S34 - S37 + 3*S47)*S67 - 2*S67**2)))/
     -   (S12*S13*S14*S23*S24)

      GCA = GCA + (2*S47*(4*S14**2*(S34 - S37) - 4*S13**2*(S14-S34+S47)+ 
     -       S34**2*(4*S26 + 4*S34 + 3*S37 + S47) - 
     -       S14*S34*(4*S26 + 9*S34 - 2*S37 + 3*S47) + 
     -       S13*(-4*S14**2 - S34*(4*S26 + 7*S34 + 8*S37 - 4*S47) + 
     -          2*S14*(2*S26 + 6*S34 + 4*S37 + S47))) + 
     -    2*(-S14**3 + S13**2*(S14 - 2*S34 + 8*S47) + 
     -       S14**2*(2*S26 + 6*S34 - S37 + 11*S47) - 
     -       S13*(2*S14**2 + 3*S34*(S26 - 2*(S34 + S37)) + 
     -          S14*(-4*S26 + 3*(S34 + S37) - 11*S47) + 
     -          (4*S26 + 23*S34 + 6*S37)*S47 - 2*S47**2) - 
     -       S14*(2*S26**2 + 5*S26*S34 + 5*S34**2 - 6*S26*S37 - 
     -          5*S34*S37 - 4*S37**2 + 2*(7*S26 + 14*S34 + S37)*S47 + 
     -          7*S47**2) + S34*
     -        (-S34**2 - 2*(3*S37 - 5*S47)*(S37 + S47) + 
     -          S26*(5*S34 - 12*S37 + 20*S47) + S34*(-7*S37 + 23*S47)))*
     -     S67 + (-8*S13**2 + 2*S14**2 + 
     -       2*S13*(S14 - 3*S26 + 9*S34 + 3*S37 - 13*S47) - 
     -       S14*(6*S26 + 7*S34 - 2*S37 + 22*S47) + 
     -       2*(8*S26*S34 - 3*S34**2 - 8*S26*S37 - 4*S34*S37 - 
     -          4*S37**2 + (12*S26 + 22*S34 + S37)*S47 + 7*S47**2))*
     -     S67**2 + (8*S13 - S14 + 2*(3*S26 - 3*S34 + S37 + 5*S47))*
     -     S67**3 - 2*S67**4 + 
     -    S124**2*(-2*S47*(S14 - S34 + S47) + 
     -       (5*S14 - 7*S34 + 6*S47)*S67 - 8*S67**2) + 
     -    S124*(-2*S47*(4*S14**2 + S34*(4*S26 + 5*S34 + 3*S37) + 
     -          S13*(4*S14 - 4*S34 + 5*S47) - 
     -          S14*(4*S26 + 10*S34 + 3*(S37 + S47))) - 
     -       (2*S14**2 + S34*(10*S26 - 9*S34 + 6*S37) + 
     -          S13*(-7*S14 + 7*S34 - 20*S47) + 
     -          2*(4*S26 + 13*S34 + S37)*S47 + 4*S47**2 - 
     -          S14*(4*S26 - 5*S34 + 4*S37 + 10*S47))*S67 - 
     -       (14*S13 + 3*S14 + 6*S26 - 19*S34 + 6*S37 + 12*S47)*
     -        S67**2 + 10*S67**3))/(S12*S123*S23*S24)

      GCA = GCA + (S124**3*(S14 - 2*S34 - S67)*S67 + 
     -    2*S124**2*((-S14 + S34)*(S14 - S26 - S34 - S47)*S47 - 
     -       (S14*(S34 - S37) + S34*(S26 - 2*S34 + 2*S37) + S26*S47 + 
     -          S47**2)*S67 - (S14 - 3*S34 + S37 - S47)*S67**2 + S67**3)
     -      - 2*(S34 + S67)*((S14 - S34)*S47*
     -        (S14*(S34 - 2*S37) - S34*(S26 + S34 - 2*S37 + S47)) + 
     -       (-S14**3 + S14**2*(2*S26 + 4*S34 + 5*S47) + 
     -          S34*(2*S26**2 + 2*S34**2 - 2*S34*S37 - 4*S37**2 + 
     -             9*S34*S47 + 6*S47**2 + 3*S26*(S34 - 2*S37 + 4*S47))
     -           - S14*(2*S26**2 + 5*S34**2 - S34*S37 - 2*S37**2 + 
     -             12*S34*S47 + 3*S47**2 + S26*(4*S34 - 2*S37 + 7*S47)))
     -         *S67 + (S14**2 + 2*S26**2 + 3*S26*S34 + 2*S34**2 - 
     -          2*S26*S37 + S34*S37 - 2*S37**2 + 7*S26*S47 + 
     -          5*S34*S47 + 2*S37*S47 + 3*S47**2 - 
     -          S14*(2*S26 + 3*S34 + S37 + 4*S47))*S67**2 + 
     -       (S37 - S47)*S67**3) + 
     -    S124*(2*(S14 - S34)*
     -        (S14*(2*S34 - S37) + S34*(-2*S26 - 2*S34 + S37 - 2*S47))*
     -        S47 + (-S14**3 + 2*S14**2*(S26 + 2*S34 + 4*S47) + 
     -          2*S34*(S26**2 + 2*(S34 - S37)*S37 + 6*S34*S47 + 
     -             6*S47**2 + 3*S26*(S34 - S37 + 3*S47)) - 
     -          2*S14*(S26**2 + 2*S34**2 + S34*S37 - S37**2 + 
     -             9*S34*S47 + 3*S47**2 + S26*(2*S34 - S37 + 5*S47)))*
     -        S67 + (S14**2 + 2*S26**2 + 6*S26*S34 - 3*S34**2 - 
     -          2*S26*S37 + 10*S34*S37 - 2*S37**2 + 
     -          2*(5*S26 + S34 + S37)*S47 + 6*S47**2 - 
     -          S14*(2*S26 + S34 + 4*(S37 + S47)))*S67**2 + 
     -       (S14 - 4*(S34 - S37 + S47))*S67**3 - S67**4))/
     -  (S12*S123*S13*S23*S24)

      GCA = GCA + (-2*(S34*(S26+S34-S37) + S123*S37 + S14*(S37-S34))*
     -     (-(S14*(S34 + S37)) + S34*(S26 + S34 + S37)) - 
     -    2*S34*(-(S14*(S34 + 3*S37)) + S34*(S26 + S34 + 3*S37))*S47 + 
     -    S124**3*S67 + (3*S123**2*S34 - S14**2*(S34 + 4*S37) + 
     -       S14*(6*S26*S34 - S34**2 + 4*S26*S37 + 8*S34*S37 - 
     -          4*S34*S47 + 6*S37*S47) + 
     -       2*S34*(S26**2 + S34**2 - 3*S34*S37 - 3*S37**2 - 
     -          5*S26*(S34 + S37) + 7*S26*S47 + S34*S47 - 2*S37*S47 + 
     -          3*S47**2) + S123*
     -        (S14*(S34 - 2*S37) + 
     -          2*(4*S26*S34 - 3*S34**2 + S26*S37 + S34*S37 + S37**2 + 
     -             2*S34*S47)))*S67 + 
     -    (2*S123**2 + 2*S123*S14 + 2*S14**2 + 4*S123*S26 + 4*S26**2 - 
     -       8*S123*S34 - 5*S14*S34 - 10*S26*S34 + 6*S34**2 - 
     -       2*S123*S37 + 2*S14*S37 - 6*S26*S37 - 6*S37**2 + 
     -       2*(2*S123 - 2*(S14 - 3*S26 + S34) + S37)*S47 + 4*S47**2)*
     -     S67**2 - 2*(S123 + 2*S26 - S34 - 2*S37 + 3*S47)*S67**3 + 
     -    S124*(2*(S14 - S26 - S34)*
     -        (S14*(S34 + S37) - S34*(S26 + S34 + S37)) + 
     -       2*(-(S14*(S34 + 2*S37)) + S34*(S26 + S34 + 2*S37))*S47 - 
     -       (S123**2 + S14**2 + 2*S26**2 - 3*S14*S34 - 2*S26*S34 - 
     -          2*S26*S37 + 4*S34*S37 - 2*S37**2 - 
     -          2*(S14 - 3*S26 + 2*S34)*S47 + 2*S47**2 + 
     -          S123*(S14 + 2*(S26 - 2*S34 - S37 + S47)))*S67 + 
     -       (S123 + 2*S26 + S34 - 6*S37 + 6*S47)*S67**2 + S67**3) - 
     -    S124**2*S67*(3*S34 + 2*(-S37 + S47 + S67)))/(S12*S13*S23*S24)

      GCA = GCA + (S34*S67*(S34 + S67)*
     -     (S124**3 + 2*S124**2*(S26 - S34 + S37 - S67) + 
     -       S124*(S34**2 + 2*S37**2 - 2*S47**2 - 
     -          2*S34*(S37 + S47 - S67) - 4*S37*S67 + S67**2 - 
     -          4*S26*(S34 - S37 + S47 + S67)) + 
     -       2*(S34 + S67)*(-(S34*(S37 - 2*S47)) + 2*S47**2 + 
     -          S37*(-2*S37 + S67) + S26*(S34 - 4*S37 + 4*S47 + S67))))/
     -   (S12*S123*S13*S14*S23*S24) + 
     -  (2*S14**3*S67*(S34 + S67) + 
     -     2*S14**2*S37*(S34 + S67)*(4*S47 + S67) + 
     -     4*S37*S67*(S34 + S67)*
     -      (S34*(6*S26 + 5*S34 + 3*S37) + 
     -        (4*S26 + 5*S34 + 2*S37)*S67 + S67**2) + 
     -     S13*(10*S34**2*S37*S47 + 
     -        2*(4*S14*S37*(2*S26 + 3*S34 + S37) + 
     -           3*S14**2*(S34 + 2*S47) + 
     -           S34*S47*(6*S26 + 7*S34 + 6*S37 + 4*S47))*S67 + 
     -        2*(8*S14*S37 + S47*(6*S26 + 10*S34 + S37 + 4*S47))*
     -         S67**2 + 6*S47*S67**3 + S67**4) + 
     -     S14*S67*(10*S34**3 + 4*S26**2*(S34 + S67) + 
     -        15*S34**2*(2*S47 + S67) + 
     -        2*S26*(S34 + S67)*(3*S34 + 10*S47 + S67) + 
     -        2*S47*S67*(4*S47 + 7*S67) + 
     -        4*S34*(2*S47**2 + 11*S47*S67 + S67**2)) + 
     -     2*S13**2*((2*S26 + 5*S34 + 7*S37)*S67*(S34 + S67) + 
     -        S14*(4*S37*S47 + S67*(2*S47 + S67))))/
     -   (S12*S123*S124*S23*S24)

      GCA=4*(GCA + ((S14*S23 + S13*S24 - S12*S34)*(S36**2+S37**2)*S67)/
     -  (S12*S13*S14*S23*S24))

C============================ Final Sum ================================
      BDPGAXP=(CA*GCA+CF*GCF)
      END
C-----------------------------------------------------------------------
      FUNCTION ERTD(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT D FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTD,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S123,S134
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S123=S12+S13+S23
      S134=S13+S14+S34
      ERTD=(S13*S23*S34+S12*S23*S34-S12**2*S34+S13*S23*S24+
     $     2*S12*S23*S24-S14*S23**2+S12*S13*S24+S12*S14*S23+
     $     S12*S13*S14)/(S13**2*S123**2)-
     $     (S12*S34**2-S13*S24*S34+S12*S24*S34-S14*S23*S34-
     $     S12*S23*S34-S13*S24**2+S14*S23*S24-S13*S23*S24-
     $     S13**2*S24+S14*S23**2)/(S13**2*S123*S134)
      END
C-----------------------------------------------------------------------
      FUNCTION BDPQD1(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED ERT D FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L     
      DOUBLE PRECISION BDPQD1,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47

      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)


      BDPQD1=(S17*(-(S26*S34**2) - S24**2*S36 + S23*S24*(S26 + S46) + 
     -       S23*S34*(S36 + S46) + S24*S34*(S26 + S36 + 2*S46)) - 
     -    S16*(-(S27*S34**2) - S24**2*S37 + S23*S24*(S27 + S47) + 
     -       S23*S34*(S37 + S47) + S24*S34*(S27 + S37 + 2*S47)))/
     -  (2.*S23**2*(S23 + S24 + S34)**2)

      BDPQD1=BDPQD1+(S12*S23*(S17*S46 + S27*S46 - (S16 + S26)*S47) + 
     -    S12*S13*(2*S17*S46 + S27*S46 + S37*S46 - 2*S16*S47 - 
     -       S26*S47 - S36*S47) + S12**2*(-(S37*S46) + S36*S47) + 
     -    S13*(S17*S23*S46 - S13*S27*S46 + S23*S37*S46 - S16*S23*S47 + 
     -       S13*S26*S47 - S23*S36*S47))/
     -  (2.*S23**2*(S12 + S13 + S23)**2)

      BDPQD1=BDPQD1+
     -  (S13*S17*S26*S34 + S17*S23*S26*S34 - S13*S16*S27*S34 - 
     -    S16*S23*S27*S34 - S13*S17*S24*S36 + S17*S23*S24*S36 - 
     -    2*S13*S24*S27*S36 + S13*S16*S24*S37 - S16*S23*S24*S37 + 
     -    2*S13*S24*S26*S37 - 2*S13*S17*S24*S46 + S13*S23*S27*S46 - 
     -    S13*S24*S27*S46 + S13*S27*S34*S46 + 2*S13*S16*S24*S47 - 
     -    S13*S23*S26*S47 + S13*S24*S26*S47 - S13*S26*S34*S47 + 
     -    S14*S23*(S27*S46 + S37*S46 + S17*(S26 + S36 + 2*S46) - 
     -       S26*S47 - S36*S47 - S16*(S27 + S37 + 2*S47)) + 
     -    S12*(2*S27*S34*S36 - 2*S26*S34*S37 + S23*S37*S46 + 
     -       S24*S37*S46 - S34*S37*S46 + 
     -       S17*(-(S26*S34) + S24*S36 - 2*S34*S46) - S23*S36*S47 - 
     -       S24*S36*S47 + S34*S36*S47 + 
     -       S16*(S27*S34 - S24*S37 + 2*S34*S47)))/
     -  (2.*S23**2*(S12 + S13 + S23)*(S23 + S24 + S34))

      END
C-----------------------------------------------------------------------
      FUNCTION BDPQD2(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED ERT D FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L     
      DOUBLE PRECISION BDPQD2,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47

      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)


      BDPQD2=(-((S14 + S34)*(-(S17*S26*S34)+S16*S27*S34-S14*S27*S36 + 
     -         S14*S26*S37)) + 
     -    S13*(S27*S34*(S36 + S46) - S26*S34*(S37 + S47) + 
     -       S14*(S17*S26 - S16*S27 - S27*S46 + S26*S47)))/
     -  (2.*S13**2*(S13 + S14 + S34)**2)

      BDPQD2=BDPQD2+(S12*S13*(S17*S46 + S27*S46 - (S16 + S26)*S47)+ 
     -    S12**2*(-(S37*S46) + S36*S47) + 
     -    S12*S23*(S17*S46 - S37*S46 + (-S16 + S36)*S47) + 
     -    S23*(S17*S23*S46 - S16*S23*S47 + 
     -       S13*(-(S27*S46) - S37*S46 + (S26 + S36)*S47)))/
     -  (2.*S13**2*(S12 + S13 + S23)**2)

      BDPQD2=BDPQD2+(S12*(-(S17*S26*S34) + S16*S27*S34-S14*S27*S36 + 
     -       S14*S26*S37 + S14*S37*S46 + S34*S37*S46 - S14*S36*S47 - 
     -       S34*S36*S47) + S23*
     -     (-(S14*S27*S36) + S14*S26*S37 - 
     -       S17*(S26*S34 + (S14 + S34)*S46) + 
     -       S16*(S27*S34 + (S14 + S34)*S47)) + 
     -    S13*(-(S14*S27*S36) + S24*S27*S36 + S14*S26*S37 - 
     -       S24*S26*S37 + S12*S37*S46 - S24*S37*S46 + 
     -       S17*(-(S26*S34) - S23*S46 + S24*(S26 + 2*S36 + S46)) - 
     -       S12*S36*S47 + S24*S36*S47 + 
     -       S16*(S27*S34 + S23*S47 - S24*(S27 + 2*S37 + S47))))/
     -  (2.*S13**2*(S12 + S13 + S23)*(S13 + S14 + S34))

      END
C-----------------------------------------------------------------------
      FUNCTION BDPQD2AX(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED ERT D FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L     
      DOUBLE PRECISION BDPQD2AX,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67

      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S67=2*DOT(P,6,7)

      BDPQD2AX = (((S14 + S34)*(-(S17*S26*S34) 
     -   - S16*S27*S34 + S14*S27*S36 + 
     -          S14*S26*S37) + 
     -       S13*(S27*S34*(S36 + S46) + S26*S34*(S37 + S47) - 
     -          S14*(S17*S26 + S16*S27 + S27*S46 + S26*S47)))/
     -     (S13 + S14 + S34)**2 + 
     -    (S12*S13*(S17*S46 + S27*S46 + (S16 + S26)*S47) + 
     -       S12*S23*(S17*S46 - S37*S46 + (S16 - S36)*S47) - 
     -       S12**2*(S37*S46 + S36*S47) + 
     -       S23*(S17*S23*S46 + S16*S23*S47 - 
     -          S13*(S27*S46 + S37*S46 + (S26 + S36)*S47)))/
     -     (S12 + S13 + S23)**2 + 
     -    (S23*(S17*S26*S34 + S16*S27*S34 - S14*(S27*S36 + S26*S37) - 
     -          S17*(S14 + S34)*S46 - S16*(S14 + S34)*S47) + 
     -       S12*(S17*S26*S34 + S16*S27*S34 - S14*S27*S36 - 
     -          S14*S26*S37 + S14*S37*S46 + S34*S37*S46 + S14*S36*S47 + 
     -          S34*S36*S47) + 
     -       S13*(2*S14*S26*S27 - 2*S26*S27*S34 + S14*S27*S36 + 
     -          S24*S27*S36 + S14*S26*S37 + S24*S26*S37 - S12*S37*S46 - 
     -          S24*S37*S46 + 
     -          S17*(-(S26*S34) + S23*S46 + S24*(-S26 + S46)) - 
     -          S12*S36*S47 - S24*S36*S47 - 2*S12*S46*S47 + 
     -          2*S23*S46*S47 + 
     -          S16*(-(S24*S27) - S27*S34 + S23*S47 + S24*S47) - 
     -          2*S14*S23*S67 + 2*S12*S34*S67))/
     -     ((S12 + S13 + S23)*(S13 + S14 + S34)))/(2.*S13**2)

      END
C-----------------------------------------------------------------------
      FUNCTION ERTE(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED ERT E FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTE,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S123,S124,S134,S234

      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S123=S12+S13+S23
      S124=S12+S14+S24
      S134=S13+S14+S34
      S234=S23+S24+S34
      ERTE=(S12*S23*S34-S12*S24*S34+S12*S14*S34+S12*S13*S34+S13*S24**2-
     $     S14*S23*S24+S13*S23*S24+S13*S14*S24+S13**2*S24-S14*S23**2-
     $     S14**2*S23-S13*S14*S23)/(S13*S23*S123*S134)-
     $     S12*(S12*S34-S23*S24-S13*S24-S14*S23-
     $     S14*S13)/(S13*S23*S123**2)-
     $     (S14+S13)*(S24+S23)*S34/(S13*S23*S134*S234)
      END
C-----------------------------------------------------------------------
      FUNCTION BDPQE(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT E FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQE,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)

      BDPQE=(-(S16*S23*S27*S34) + S14*S23*S27*S36 - S14*S16*S23*S37 - 
     -    S14*S23*S26*S37 + S12*S16*S34*S37 + S14*S23*S27*S46 + 
     -    S12*S27*S34*S46 - S12*S14*S37*S46 + 
     -    S17*(S23*S26*S34 - S12*S34*S36 + S14*S23*(S36 + S46)) - 
     -    S14*S16*S23*S47 - S14*S23*S26*S47 - S12*S26*S34*S47 + 
     -    S12*S14*S36*S47 + S13*
     -     (S17*S26*S34 - S16*S27*S34 - S24*S27*S36 + S24*S26*S37 - 
     -       S24*S27*S46 - S12*S37*S46 - S17*S24*(S36 + S46) + 
     -       S24*S26*S47 + S12*S36*S47 + S16*S24*(S37 + S47)))/
     -  (2.*S13*S14*(S12 + S13 + S23)*(S13 + S14 + S34))

      BDPQE=BDPQE-
     -    (S16*S23*S27*S34 - S14*S23*S27*S36 + S12*S17*S34*S36 + 
     -    S14*S16*S23*S37 + S14*S23*S26*S37 - S12*S16*S34*S37 - 
     -    S14*S23*S27*S46 - S12*S27*S34*S46 + S12*S14*S37*S46 - 
     -    S17*S23*(S26*S34 + S14*(S36 + S46)) + S14*S16*S23*S47 + 
     -    S14*S23*S26*S47 + S12*S26*S34*S47 - S12*S14*S36*S47 + 
     -    S13*(-(S17*S26*S34) + S16*S27*S34 + S24*S27*S36 - 
     -       S24*S26*S37 + S24*S27*S46 + S12*S37*S46 + 
     -       S17*S24*(S36 + S46) - S24*S26*S47 - S12*S36*S47 - 
     -       S16*S24*(S37 + S47)))/
     -  (2.*S13*S23*(S12 + S13 + S23)*(S13 + S14 + S34))

      BDPQE=BDPQE+(S34*(-(S17*S26*S34) + S16*S27*S34 - 
     -      (S13 + S14)*(S27*(S36 + S46) - S26*(S37 + S47))))/
     -  (2.*S13*S14*(S13 + S14 + S34)**2)

      BDPQE=BDPQE-
     -    (S12*(S17*S23*S46 + S23*S27*S46 - S12*S37*S46 - S16*S23*S47 - 
     -      S23*S26*S47 + S12*S36*S47 + 
     -      S13*(S17*S46 + S27*S46 - (S16 + S26)*S47)))/
     -  (2.*S13*S23*(S12 + S13 + S23)**2)

      BDPQE=BDPQE+
     -     (S12*(S14*S17*S36 + S17*S24*S36 + S14*S27*S36 + S24*S27*S36 - 
     -      S16*S24*S37 - S24*S26*S37 - S14*(S16 + S26)*S37 + 
     -      S17*(S13 + S23)*S46 + S13*S27*S46 + S23*S27*S46 - 
     -      S13*S16*S47 - S16*S23*S47 - S13*S26*S47 - S23*S26*S47))/
     -  (4.*S13*S14*(S12 + S13 + S23)*(S12 + S14 + S24))

      BDPQE=BDPQE-
     -      (S34*(-(S14*S27*S36) + S16*S23*S37 + S16*S24*S37 + 
     -      S14*S26*S37 - S14*S27*S46 - S17*(S23 + S24)*(S36 + S46) - 
     -      S13*S27*(S36 + S46) + S16*S23*S47 + S16*S24*S47 + 
     -      S14*S26*S47 + S13*S26*(S37 + S47)))/
     -  (4.*S13*S23*(S13 + S14 + S34)*(S23 + S24 + S34))


      END
C-----------------------------------------------------------------------
      FUNCTION BDPQFVV(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE UN-POLARIZED ERT F FUNCTION (FOR VECTOR-CURRENT ONLY)
C     ONLY NON-ZERO IF THE CHARGE OF THE FINAL STATE IS OBSERVED 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQFVV,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S67=2*DOT(P,6,7)

      BDPQFVV=(-2*S14*S23*S26*S27+2*S13*S26*S27*S34-S14*S23*S27*S36 - 
     -    S13*S24*S27*S36 - S12*S27*S34*S36 - S14*S23*S26*S37 - 
     -    S13*S24*S26*S37 - S12*S26*S34*S37 - 2*S14*S23*S36*S37 + 
     -    2*S12*S24*S36*S37 + S12*S23*S27*S46 + S13*S23*S27*S46 + 
     -    S14*S23*S27*S46 - S13*S24*S27*S46 + S13*S23*S37*S46 + 
     -    S12*S24*S37*S46 + S17*
     -     (-2*S23**2*S46 - S12*S34*S46 + 
     -       S23*(-2*S16*S34 + S14*S36 + S34*S36 + S24*(S26 + S36) + 
     -          S14*S46) + S13*(S26*S34 - S24*(S36 + S46))) + 
     -    S12*S23*S26*S47 + S13*S23*S26*S47 + S14*S23*S26*S47 - 
     -    S13*S24*S26*S47 + S13*S23*S36*S47 + S12*S24*S36*S47 - 
     -    2*S12*S23*S46*S47 + 
     -    S16*(S13*S27*S34 - S13*S24*S37 - 2*S23**2*S47 - S13*S24*S47 - 
     -       S12*S34*S47 + S23*
     -        (S34*S37 + S24*(S27 + S37) + S14*(S37 + S47))) + 
     -    2*S14*S23**2*S67 + 2*S12*S23*S34*S67)/
     -  (2.*S13*(S12 + S13 + S23)*S24*(S23 + S24 + S34))

      BDPQFVV=BDPQFVV
     -      -(-(S13*S17*S24*S26) - S12*S17*S26*S34 - 2*S12*S26*S27*S34 + 
     -     S12*S17*S24*S36 - S13*S17*S24*S36 + S12*S24*S27*S36 + 
     -     S12*S17*S34*S36 + S12*S24*S26*S37 + S12*S13*S17*S46 + 
     -     S17*S23*S24*S46 + S12*S13*S27*S46 + S12*S23*S27*S46 - 
     -     S13*S24*S27*S46 + S12*S27*S34*S46 - 2*S12**2*S37*S46 - 
     -     S13*S24*S37*S46 + S12*S34*S37*S46 + S12*S13*S26*S47 + 
     -     S12*S23*S26*S47 - S13*S24*S26*S47 + S12*S26*S34*S47 - 
     -     2*S12**2*S36*S47 - S13*S24*S36*S47 + S12*S34*S36*S47 - 
     -     2*S12*S23*S46*S47 + 
     -     S16*(2*S17*S23*S24 - S14*S23*S27 - S13*S24*S27 - 
     -        2*S12*S17*S34 - S12*S27*S34 + S12*S14*S37 + S12*S24*S37 - 
     -        S13*S24*S37 + S12*S34*S37 + S12*S13*S47 + S23*S24*S47) + 
     -     2*S12**2*S34*S67 + 
     -     S14*(-(S17*S23*S26) + 2*S13*S26*S27 + S12*S17*S36 + 
     -        S13*S27*S36 + S13*S26*S37 - 2*S12*S36*S37 - S23*S37*S46 - 
     -        S23*S36*S47 + 2*S12*S23*S67))/
     -  (2.*S13*(S12 + S13 + S23)*S24*(S12 + S14 + S24))

      BDPQFVV=BDPQFVV
     -    +(S13*(-(S27*S34*S36) + S17*S24*(S26 + S36) - S26*S34*S37 + 
     -       S16*S24*(S27 + S37) - S17*S23*S46 + S24*S27*S46 - 
     -       S27*S34*S46 + S24*S37*S46 - S16*S23*S47 + S24*S26*S47 - 
     -       S26*S34*S47 + S24*S36*S47 - 2*S23*S46*S47) + 
     -    S14*(S17*S23*S26 + S16*S23*S27 + 2*S26*S27*S34 - 
     -       S24*S27*S36 - S24*S26*S37 - 2*S24*S36*S37 - S27*S34*S46 + 
     -       S23*S37*S46 - S26*S34*S47 + S23*S36*S47 - 2*S23*S34*S67) - 
     -    S34*(S17*(-2*S26*S34 + S23*S36 + S24*(S36 + S46)) + 
     -       S16*(-2*S17*S23 - 2*S27*S34 + S23*S37 + S24*S37 + 
     -          S12*(S27 + S37) + S24*S47) + 
     -       S12*(S17*(S26 + S36) + S27*S46 - S37*S46 + S26*S47 - 
     -          2*S46*S47 - S36*(2*S37 + S47) + 2*S34*S67)))/
     -  (2.*S13*S24*(S13 + S14 + S34)*(S23 + S24 + S34))
     
      BDPQFVV=BDPQFVV
     -    -(-(S17*S24*S26*S34) + S13*S17*S24*S36 + S13*S24*S27*S36 + 
     -     S12*S27*S34*S36 + S13*S24*S26*S37 + S12*S26*S34*S37 + 
     -     S13*S17*S24*S46 + S13*S24*S27*S46 + S12*S17*S34*S46 - 
     -     S12*S13*S37*S46 + S13*S24*S26*S47 - S12*S13*S36*S47 - 
     -     2*S12*S13*S46*S47 + 
     -     S16*(-2*S17*S24*S34 - S24*S27*S34 + S13*S24*S37 + 
     -        S13*S24*S47 + S12*S34*S47) + 
     -     2*S14**2*(S27*S36 + S26*S37 - S23*S67) - 
     -     S14*(S13*S17*S26 - 2*S26*S27*S34 + S12*S17*S36 + 
     -        S17*S23*S36 + S17*S24*S36 + S23*S27*S36 + S23*S26*S37 - 
     -        2*S12*S36*S37 - S17*S23*S46 + S13*S27*S46 + S23*S27*S46 + 
     -        S27*S34*S46 + S24*S37*S46 + S13*S26*S47 + S23*S26*S47 + 
     -        S26*S34*S47 + S24*S36*S47 - 2*S23*S46*S47 + 
     -        S16*(-2*S17*S23 + S13*S27 + S12*S37 + S23*S37 + S24*S37 - 
     -           S23*S47) + 2*S12*S34*S67))/
     -  (2.*S13*S24*(S12 + S14 + S24)*(S13 + S14 + S34))
          
      END      
C-----------------------------------------------------------------------
      FUNCTION BDPQFAA(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE UN-POLARIZED ERT F FUNCTION (FOR VECTOR-CURRENT ONLY)
C     ONLY NON-ZERO IF THE CHARGE OF THE FINAL STATE IS OBSERVED 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQFAA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S67=2*DOT(P,6,7)

C---- AAAA-like Couplings
      BDPQFAA=-(-2*S14*S23*S26*S27-2*S13*S26*S27*S34-3*S14*S23*S27*S36+ 
     -     S13*S24*S27*S36 + S12*S27*S34*S36 - 3*S14*S23*S26*S37 + 
     -     S13*S24*S26*S37 + S12*S26*S34*S37 - 2*S14*S23*S36*S37 - 
     -     2*S12*S24*S36*S37 + S12*S23*S27*S46 + S13*S23*S27*S46 - 
     -     S14*S23*S27*S46 + S13*S24*S27*S46 + 2*S12*S23*S37*S46 + 
     -     S13*S23*S37*S46 - S12*S24*S37*S46 + 
     -     S17*(-2*S23**2*S46 + S12*S34*S46 + 
     -        S23*(2*S16*S34 + 2*S26*S34 - S14*S36 + S34*S36 + 
     -           S24*(S26 + S36) - S14*S46) + 
     -        S13*(-(S26*S34) + S24*(S36 + S46))) + S12*S23*S26*S47 + 
     -     S13*S23*S26*S47 - S14*S23*S26*S47 + S13*S24*S26*S47 + 
     -     2*S12*S23*S36*S47 + S13*S23*S36*S47 - S12*S24*S36*S47 + 
     -     2*S12*S23*S46*S47 + 
     -     S16*(-2*S23**2*S47 + S12*S34*S47 + 
     -        S23*(2*S27*S34 - S14*S37 + S34*S37 + S24*(S27 + S37) - 
     -           S14*S47) + S13*(-(S27*S34) + S24*(S37 + S47))) + 
     -     2*S14*S23**2*S67 - 2*S12*S23*S34*S67)/
     -  (2.*S13*(S12 + S13 + S23)*S24*(S23 + S24 + S34))

      BDPQFAA=BDPQFAA
     -      -(S13*S17*S24*S26 - 3*S12*S17*S26*S34 - 2*S12*S26*S27*S34 + 
     -     S12*S17*S24*S36 + S13*S17*S24*S36 + S12*S24*S27*S36 - 
     -     S12*S17*S34*S36 + S12*S24*S26*S37 + S12*S13*S17*S46 + 
     -     2*S12*S17*S23*S46 - S17*S23*S24*S46 + S12*S13*S27*S46 + 
     -     S12*S23*S27*S46 + S13*S24*S27*S46 - S12*S27*S34*S46 - 
     -     2*S12**2*S37*S46 + S13*S24*S37*S46 - S12*S34*S37*S46 + 
     -     S12*S13*S26*S47 + S12*S23*S26*S47 + S13*S24*S26*S47 - 
     -     S12*S26*S34*S47 - 2*S12**2*S36*S47 + S13*S24*S36*S47 - 
     -     S12*S34*S36*S47 + 2*S12*S23*S46*S47 + 
     -     S16*(S14*S23*S27 + S13*S24*S27 - 3*S12*S27*S34 - 
     -        2*S17*(S23*S24 + S12*S34) + S12*S14*S37 + S12*S24*S37 + 
     -        S13*S24*S37 - S12*S34*S37 + S12*S13*S47 + 2*S12*S23*S47 - 
     -        S23*S24*S47) + 2*S12**2*S34*S67 + 
     -     S14*(2*S12*S27*S36 + S17*(S23*S26 + S12*S36) + 
     -        2*S12*S26*S37 + 2*S12*S36*S37 - 
     -        S13*(2*S26*S27 + S27*S36 + S26*S37) + S23*S37*S46 + 
     -        S23*S36*S47 - 2*S12*S23*S67))/
     -  (2.*S13*(S12 + S13 + S23)*S24*(S12 + S14 + S24))

      BDPQFAA=BDPQFAA
     -    -(S13*(S27*S34*S36 + S17*S24*(S26 + S36) + S26*S34*S37 + 
     -        S16*S24*(S27 + S37) - S17*S23*S46 + S24*S27*S46 + 
     -        S27*S34*S46 + S24*S37*S46 - S16*S23*S47 + S24*S26*S47 + 
     -        S26*S34*S47 + S24*S36*S47 - 2*S23*S46*S47) + 
     -     S14*(S17*S23*S26 + S16*S23*S27 + 2*S26*S27*S34 - 
     -        S24*S27*S36 + 2*S27*S34*S36 - S24*S26*S37 + 
     -        2*S26*S34*S37 - 2*S24*S36*S37 + S27*S34*S46 + 
     -        S23*S37*S46 + S26*S34*S47 + S23*S36*S47 - 2*S23*S34*S67)
     -      - S34*(S17*(2*S26*S34 - S24*(S36 + S46) - 
     -           S23*(S36 + 2*S46)) - 
     -        S16*(2*S17*S23 - 2*S27*S34 + S23*S37 + S24*S37 - 
     -           S12*(S27 + S37) + 2*S23*S47 + S24*S47) + 
     -        S12*(S17*(S26 + S36) + 2*S36*S37 + S27*S46 + 3*S37*S46 + 
     -           S26*S47 + 3*S36*S47 + 2*S46*S47 - 2*S34*S67)))/
     -  (2.*S13*S24*(S13 + S14 + S34)*(S23 + S24 + S34))
     
      BDPQFAA=BDPQFAA
     -    -(-(S17*S24*S26*S34) + S13*S17*S24*S36 + S13*S24*S27*S36 + 
     -     S12*S27*S34*S36 + S13*S24*S26*S37 + S12*S26*S34*S37 + 
     -     S13*S17*S24*S46 + S13*S24*S27*S46 + S12*S17*S34*S46 - 
     -     S12*S13*S37*S46 + S13*S24*S26*S47 - S12*S13*S36*S47 - 
     -     2*S12*S13*S46*S47 + 
     -     S16*(-2*S17*S24*S34 - S24*S27*S34 + S13*S24*S37 + 
     -        S13*S24*S47 + S12*S34*S47) - 
     -     2*S14**2*(S27*S36 + S26*S37 - S23*S67) + 
     -     S14*(S13*S17*S26 + 2*S17*S26*S34 + 2*S26*S27*S34 + 
     -        S12*S17*S36 - S17*S23*S36 + S17*S24*S36 - S23*S27*S36 - 
     -        S23*S26*S37 + 2*S12*S36*S37 - 3*S17*S23*S46 + 
     -        S13*S27*S46 - S23*S27*S46 + S27*S34*S46 + 2*S12*S37*S46 + 
     -        S24*S37*S46 + S13*S26*S47 - S23*S26*S47 + S26*S34*S47 + 
     -        2*S12*S36*S47 + S24*S36*S47 - 2*S23*S46*S47 + 
     -        S16*(-2*S17*S23 + S13*S27 + 2*S27*S34 + S12*S37 - 
     -           S23*S37 + S24*S37 - 3*S23*S47) - 2*S12*S34*S67))/
     -  (2.*S13*S24*(S12 + S14 + S24)*(S13 + S14 + S34))
          
      END      
C-----------------------------------------------------------------------
      FUNCTION BDPQFAV(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE UN-POLARIZED ERT F FUNCTION (FOR VECTOR-CURRENT ONLY)
C     ONLY NON-ZERO IF THE CHARGE OF THE FINAL STATE IS OBSERVED 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQFAV,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S67=2*DOT(P,6,7)

C---- AVAV-like Couplings
      BDPQFAV=(-(S14*S23*S27*S36) + S13*S24*S27*S36 - S12*S27*S34*S36 + 
     -    S14*S23*S26*S37 - S13*S24*S26*S37 + S12*S26*S34*S37 + 
     -    S12*S23*S27*S46 + S13*S23*S27*S46 - S14*S23*S27*S46 + 
     -    S13*S24*S27*S46 + S13*S23*S37*S46 - S12*S24*S37*S46 + 
     -    S17*(S12*S34*S46 - 
     -       S23*(2*S26*S34 + S14*S36 + S34*S36 + S24*(S26 + S36) + 
     -          S14*S46) + S13*(-(S26*S34) + S24*(S36 + S46))) - 
     -    S12*S23*S26*S47 - S13*S23*S26*S47 + S14*S23*S26*S47 - 
     -    S13*S24*S26*S47 - S13*S23*S36*S47 + S12*S24*S36*S47 + 
     -    S16*(-(S12*S34*S47) + 
     -       S23*(2*S27*S34 + S14*S37 + S34*S37 + S24*(S27 + S37) + 
     -          S14*S47) + S13*(S27*S34 - S24*(S37 + S47))))/
     -  (2.*S13*(S12 + S13 + S23)*S24*(S23 + S24 + S34))

      BDPQFAV=BDPQFAV
     -      +(S12*S17*S26*S34 - S12*S16*S27*S34 + S12*S17*S24*S36 + 
     -    S12*S24*S27*S36 + S12*S17*S34*S36 - S12*S16*S24*S37 - 
     -    S12*S24*S26*S37 - S12*S16*S34*S37 - S17*S23*S24*S46 + 
     -    S12*S23*S27*S46 - S12*S27*S34*S46 - S12*S34*S37*S46 + 
     -    S16*S23*S24*S47 - S12*S23*S26*S47 + S12*S26*S34*S47 + 
     -    S12*S34*S36*S47 + S14*
     -     (2*S12*S27*S36 + S13*S27*S36 + S17*(S23*S26 + S12*S36) - 
     -       2*S12*S26*S37 - S13*S26*S37 - S16*(S23*S27 + S12*S37) + 
     -       S23*S37*S46 - S23*S36*S47) + 
     -    S13*(-(S17*S24*(S26 + S36)) + S16*S24*(S27 + S37) + 
     -       S12*S17*S46 + S12*S27*S46 + S24*S27*S46 + S24*S37*S46 - 
     -       S12*S16*S47 - S12*S26*S47 - S24*S26*S47 - S24*S36*S47))/
     -  (2.*S13*(S12 + S13 + S23)*S24*(S12 + S14 + S24))

      BDPQFAV=BDPQFAV
     -    +(S14*(-(S17*S23*S26) + S16*S23*S27 - S24*S27*S36 + 
     -       S24*S26*S37 + S27*S34*S46 - S23*S37*S46 - S26*S34*S47 + 
     -       S23*S36*S47) + S13*
     -     (S27*S34*S36 - S17*S24*(S26 + S36) - S26*S34*S37 + 
     -       S16*S24*(S27 + S37) + S17*S23*S46 + S24*S27*S46 + 
     -       S27*S34*S46 + S24*S37*S46 - S16*S23*S47 - S24*S26*S47 - 
     -       S26*S34*S47 - S24*S36*S47) + 
     -    S34*(S17*(S24*(S36 + S46) + S23*(S36 + 2*S46)) + 
     -       S12*(S17*(S26 + S36) - S16*(S27 + S37) - S27*S46 - 
     -          S37*S46 + S26*S47 + S36*S47) - 
     -       S16*(S24*(S37 + S47) + S23*(S37 + 2*S47))))/
     -  (2.*S13*S24*(S13 + S14 + S34)*(S23 + S24 + S34))
     
      BDPQFAV=BDPQFAV
     -    +(-(S16*S24*S27*S34) - S14*S23*S27*S36 + S12*S27*S34*S36 + 
     -    S12*S14*S16*S37 + S14*S16*S23*S37 + S14*S16*S24*S37 + 
     -    S14*S23*S26*S37 - S12*S26*S34*S37 - S14*S23*S27*S46 + 
     -    S14*S27*S34*S46 + 2*S12*S14*S37*S46 + S14*S24*S37*S46 + 
     -    S17*(S24*(S26*S34 - S14*S36) - S14*S23*(S36 + S46) - 
     -       S12*(S14*S36 + S34*S46)) + S14*S16*S23*S47 + 
     -    S14*S23*S26*S47 + S12*S16*S34*S47 - S14*S26*S34*S47 - 
     -    2*S12*S14*S36*S47 - S14*S24*S36*S47 + 
     -    S13*(S24*S27*S36 - S16*S24*S37 - S24*S26*S37 + S24*S27*S46 + 
     -       S12*S37*S46 + S17*S24*(S36 + S46) - S16*S24*S47 - 
     -       S24*S26*S47 - S12*S36*S47 + 
     -       S14*(-(S17*S26) + S16*S27 + S27*S46 - S26*S47)))/
     -  (2.*S13*S24*(S12 + S14 + S24)*(S13 + S14 + S34))
          
      END      
C-----------------------------------------------------------------------
      FUNCTION BDPQFVA(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE UN-POLARIZED ERT F FUNCTION (FOR VECTOR-CURRENT ONLY)
C     ONLY NON-ZERO IF THE CHARGE OF THE FINAL STATE IS OBSERVED 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQFVA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S67=2*DOT(P,6,7)

C---- AVAV-like Couplings
      BDPQFVA=(-(S14*S23*S27*S36) + S13*S24*S27*S36 - S12*S27*S34*S36 + 
     -    S14*S23*S26*S37 - S13*S24*S26*S37 + S12*S26*S34*S37 - 
     -    S12*S23*S27*S46 - S13*S23*S27*S46 - S14*S23*S27*S46 + 
     -    S13*S24*S27*S46 - 2*S12*S23*S37*S46 - S13*S23*S37*S46 - 
     -    S12*S24*S37*S46 + S17*
     -     (-(S13*S26*S34) + S13*S24*S36 + S13*S24*S46 + S12*S34*S46 + 
     -       S23*(S34*S36 + S24*(S26 + S36) - S14*(S36 + S46))) + 
     -    S12*S23*S26*S47 + S13*S23*S26*S47 + S14*S23*S26*S47 - 
     -    S13*S24*S26*S47 + 2*S12*S23*S36*S47 + S13*S23*S36*S47 + 
     -    S12*S24*S36*S47 - S16*
     -     (-(S13*S27*S34) + S13*S24*S37 + S13*S24*S47 + S12*S34*S47 + 
     -       S23*(S34*S37 + S24*(S27 + S37) - S14*(S37 + S47))))/
     -  (2.*S13*(S12 + S13 + S23)*S24*(S23 + S24 + S34))

      BDPQFVA=BDPQFVA
     -      +(-(S12*S17*S26*S34) + S12*S16*S27*S34 + S12*S17*S24*S36 + 
     -    S12*S24*S27*S36 - S12*S17*S34*S36 - S12*S16*S24*S37 - 
     -    S12*S24*S26*S37 + S12*S16*S34*S37 + 2*S12*S17*S23*S46 + 
     -    S17*S23*S24*S46 + S12*S23*S27*S46 + S12*S27*S34*S46 + 
     -    S12*S34*S37*S46 - 2*S12*S16*S23*S47 - S16*S23*S24*S47 - 
     -    S12*S23*S26*S47 - S12*S26*S34*S47 - S12*S34*S36*S47 + 
     -    S14*(-(S17*S23*S26) + S16*S23*S27 + S12*S17*S36 - 
     -       S13*S27*S36 - S12*S16*S37 + S13*S26*S37 - S23*S37*S46 + 
     -       S23*S36*S47) + S13*
     -     (S12*S27*S46 - S24*S27*S46 - S24*S37*S46 + 
     -       S17*(S24*(S26 + S36) + S12*S46) - S12*S26*S47 + 
     -       S24*S26*S47 + S24*S36*S47 - S16*(S24*(S27 + S37) + S12*S47)
     -       ))/(2.*S13*(S12 + S13 + S23)*S24*(S12 + S14 + S24))

      BDPQFVA=BDPQFVA
     -    +(S14*(S17*S23*S26 - S16*S23*S27 + S24*S27*S36+2*S27*S34*S36 - 
     -       S24*S26*S37 - 2*S26*S34*S37 + S27*S34*S46 + S23*S37*S46 - 
     -       S26*S34*S47 - S23*S36*S47) + 
     -    S13*(S27*S34*S36 + S17*S24*(S26 + S36) - S26*S34*S37 - 
     -       S16*S24*(S27 + S37) - S17*S23*S46 - S24*S27*S46 + 
     -       S27*S34*S46 - S24*S37*S46 + S16*S23*S47 + S24*S26*S47 - 
     -       S26*S34*S47 + S24*S36*S47) + 
     -    S34*(S17*(S23*S36 + S24*(S36 + S46)) + 
     -       S12*(-(S17*(S26 + S36)) + S16*(S27 + S37) + S27*S46 + 
     -          S37*S46 - S26*S47 - S36*S47) - 
     -       S16*(S23*S37 + S24*(S37 + S47))))/
     -  (2.*S13*S24*(S13 + S14 + S34)*(S23 + S24 + S34))
     
      BDPQFVA=BDPQFVA
     -    +(S34*(S17*S24*S26 - S16*S24*S27 + S12*S27*S36 - S12*S26*S37 - 
     -       S12*S17*S46 + S12*S16*S47) + 
     -    S14*(-(S23*S27*S36) + S23*S26*S37 - S23*S27*S46 - 
     -       S27*S34*S46 - S24*S37*S46 + 
     -       S17*(2*S26*S34 + S12*S36 - S23*S36 + S24*S36 - S23*S46) + 
     -       S23*S26*S47 + S26*S34*S47 + S24*S36*S47 + 
     -       S16*(-2*S27*S34 - S12*S37 + S23*S37 - S24*S37 + S23*S47))
     -     + S13*(S24*S27*S36 - S16*S24*S37 - S24*S26*S37 + 
     -       S24*S27*S46 + S12*S37*S46 + S17*S24*(S36 + S46) - 
     -       S16*S24*S47 - S24*S26*S47 - S12*S36*S47 + 
     -       S14*(S17*S26 - S16*S27 - S27*S46 + S26*S47)))/
     -  (2.*S13*S24*(S12 + S14 + S24)*(S13 + S14 + S34))
          
      END      
C-----------------------------------------------------------------------
      FUNCTION LEIA(P,Q,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIA,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,Q2,
     $     P1K,P2K,P3K,P4K,KK
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIA = 0
      LEIA = LEIA + P1K*P2K * (  - 2*Q2*P12*P25 - Q2*P14*
     +     P25 - Q2*P23*P25 + P12*P13*P25 + P12*P24*P25 + 2*P12*P25*
     +     P34 - P13*P23*P25 - 2*P13*P24*P25 - P13*P25*P34 + 2*P14*P15*
     +     P24 - 2*P14*P23*P25 - P14*P24*P25 - 2*P15*P23*P25 - 2*P15*P24
     +     *P25 - 2*P15*P25*P34 - P24*P25*P34 )
      LEIA = LEIA + P1K*P3K * (  - Q2*P12 + Q2*P24 + P12*
     +     P13 + P14*P24 + 2*P15*P24 - P24*P34 )*P25
      LEIA = LEIA + P1K*P4K * (  - Q2*P12 - Q2*P23 + P12*
     +     P13 + P12*P34 - P14*P23 )*P25
      LEIA = LEIA + P1K**2 * ( Q2*P23 + Q2*P24 - P13*P23
     +     - P13*P24 - P24*P34 )*P25
      LEIA = LEIA + P2K*P3K * (  - Q2*P12*P25 - Q2*P14*
     +     P25 - 2*P12*P15*P25 + P12*P24*P25 + P12*P25*P34 + 2*P14*P15*
     +     P24 - P14*P23*P25 - 4*P15*P23*P25 - 2*P15*P24*P25 - 2*P15*P25
     +     *P34 )
      LEIA = LEIA + P2K*P4K * (  - Q2*P12*P25 + Q2*P13*
     +     P25 + P12*P24*P25 + 2*P13*P15*P25 + P13*P23*P25 - P13*P25*P34
     +     + 2*P14*P15*P24 - 2*P15*P23*P25 - 2*P15*P24*P25 )
      LEIA = LEIA + P2K**2 * ( Q2*P13 + Q2*P14 + 2*P13*P15
     +     - P13*P24 - P13*P34 - P14*P24 + 2*P15*P34 )*P25
      LEIA = LEIA + P3K*P4K * (  - 2*P12*P15 + P12*P34 - P13*
     +     P24 - P14*P23 - 2*P15*P23 - P15*P25 )*P25
      LEIA = LEIA + P3K**2 * ( P14 + 2*P15 )*P24*P25
      LEIA = LEIA + P4K**2 * ( P13*P23*P25 )
c$$$      LEIA = LEIA + KK * ( 2*Q2*P12*P14*P25 + 2*Q2*P12*P23
c$$$     +     *P25 + 2*Q2*P12**2*P25 - 2*Q2*P14*P15*P24 + 2*Q2*P14*
c$$$     +     P23*P25 + P12*P13*P25*P34 + 4*P12*P15*P23*P25 + 2*P12*P15*P25
c$$$     +     *P34 + P12*P24*P25*P34 - P13*P14*P23*P25 - 2*P13*P14*P24*P25
c$$$     +     - 2*P13*P15*P24*P25 - 2*P13*P23*P24*P25 - P13*P24**2*P25 - 
c$$$     +     P13**2*P24*P25 + 2*P14*P15*P23*P25 - P14*P23*P24*P25 + 4*P15*
c$$$     +     P23*P24*P25 + 4*P15*P23*P25*P34 + 4*P15*P23**2*P25 + 2*P15*
c$$$     +     P24*P25*P35 + P15*P25**2*P34 )/4
      LEIA = LEIA/(P15*P25**2*P13*P24)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIGAXA(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIGAXA,P(4,7),DOT,A6,A6COEF,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=-2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=-2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=-2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=-2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=-2*DOT(P,6,7)
      A6=A6COEF(P,I,J,K,L)

      LEIGAXA=(S14**3*S23**2*S26 + 
     -     S13*S14*S24*(2*S13*S26**2 - 2*S26**2*S34 - 
     -        S16*S23*(S24 + 2*S26 + S34) - 2*S12*S26*S36 + 
     -        2*S24*S26*S36 + S12*S23*S46 + S13*S23*S46 + S23**2*S46 - 
     -        2*S23*S36*S46) - 
     -     S14**2*S23*(-(S13*S24*(S26 + S36)) + S16*S23*(S24 - 2*S46) + 
     -        S13*S23*S46 + S23*
     -         (-(S26*S34) + S24*S36 + (S12 + S23 - 2*S36)*S46)) + 
     -     2*S13*S24*(S16**2*S23*S24 + 
     -        S46*(-2*S13*S24*S26 + S12**2*S36 + 
     -           S12*(-(S13*S26) + S26*S34 + S24*S36 - S23*S46)) - 
     -        S16*(S13*S24*S26 + S24*(-(S26*S34) + S24*S36 - S23*S46) + 
     -           S12*(-2*S26*S34 + S24*S36 + S23*S46))))/
     -   (4.*S13*S14**2*S23**2*S24) + 
     -  (A6*(S14**3*S23**2 - 
     -       S14*S23*S24*(-2*S16*S23 + S23*(S24 + S34) + 
     -          2*(S13*S26 - S26*S34 + S12*S36 - S24*S36)) + 
     -       S14**2*S23**2*(S34 - 2*S46) + 
     -       2*S24*(S13**2*S24*S26 - 
     -          S13*(S16*S23*S24 + 
     -             S12*(S26*S34 + S24*S36 - 2*S23*S46) + 
     -             S24*(-(S26*S34) + S24*S36 - S23*S46)) + 
     -          S34*(-2*S16*S23*S24 + S12**2*S36 + 
     -             S12*(-(S16*S23) - S26*S34 + S24*S36 + S23*S46)))))/
     -   (4.*S14**2*S23**2*S24) 

      LEIGAXA = LEIGAXA +
     -    (A6*(2*S14**2*S23*(S23 - S24) + 2*S12**2*(S14*S23 - S24*S34) + 
     -       S14*S23*(2*S23**2 + 
     -          S24*(-5*S13 + 4*S16 + 2*S24 - 4*S26 + S34 + 6*S36) + 
     -          S23*(3*S24 - 4*S26 + 2*S34 - 4*S46)) + 
     -       S12*(2*S14**2*S23 + S13*S24*(2*S24 - S34) - 
     -          S24*S34*(3*S23 + 2*S24 - 4*S26 + S34 - 2*S36) + 
     -          2*S14*S23*(2*S23 + 2*S24 - 2*S26 + S34 - 2*S46)) + 
     -       S24*(S13**2*S24 + S13*S23*(3*S24 + 2*S34) + 
     -          S13*S24*(2*S24 - 4*S26 + S34 - 2*S36) - 
     -          2*S23*S34*(2*S16 + S23 - 2*S26 - S34 + 2*S46))))/
     -   (4.*S14*S23*S24*(S13 + S14 + S34)) - 
     -  (S13**2*S24*(S16*S24 + S14*S26 + 2*S24*S26 - S12*S46 + 
     -        2*S23*S46) + S13*
     -      (2*S14**2*S23*S26 + 
     -        S14*(-4*S16*S23*S24 + 2*S12*S23*S26 + 2*S23**2*S26 - 
     -           3*S23*S24*S26 + 4*S23*S26**2 + 2*S23*S26*S34 + 
     -           S24*S26*S34 - 4*S23*S24*S36 - 2*S24**2*S36 - 
     -           2*S24*S26*S36 + 4*S23*S24*S46 + 4*S23*S26*S46) + 
     -        S24*(-2*S12*S26*S34 - 2*S23*S26*S34 + 4*S26**2*S34 - 
     -           4*S24*S26*S36 + S16*S24*(S34 + 2*S36) + 
     -           2*S23*S34*S46 + 4*S24*S36*S46 + 
     -           S12*(S23 - S34 + 2*S36)*S46 - 4*S23*S46**2 - 
     -           S16*S23*(S24 + 4*S46))) - 
     -     2*S14*S23*(-2*S16**2*S24 - S14**2*S26 + S24*S26*S34 - 
     -        2*S26**2*S34 - S26*S34**2 - S24**2*S36 + 4*S12*S26*S36 + 
     -        4*S23*S26*S36 + 2*S24*S26*S36 + S24*S34*S36 + 
     -        2*S26*S34*S36 - 2*S24*S36**2 + S12*S24*S46 + 
     -        S23*S24*S46 + 2*S12*S26*S46 + 2*S23*S26*S46 + 
     -        S12*S34*S46 + S23*S34*S46 + 2*S12*S36*S46 + 
     -        2*S23*S36*S46 + 
     -        S14*(-2*S26**2 + S16*(S24 + 2*S26) - 2*S26*S34 + 
     -           2*S26*S36 + S24*(S26 + S36) + S12*S46 + S23*S46) + 
     -        S16*(-S24**2 + S24*(2*S26 + S34 - 4*S36) + 
     -           2*(2*S23*S26 + S26*S34 + S23*S46 + S12*(2*S26 + S46))))
     -     )/(4.*S13*S14*S23*S24*(S13 + S14 + S34))

      LEIGAXA = LEIGAXA -
     -     (A6*(S12 + S23 + S24 - 2*S26)*(S14 - S34))/
     -   (2.*S14*(S13 + S14 + S34)) - 
     -  (S13**2*S26*(S14 + S34) + 
     -     S14*(S14 + S34)*(-(S16*S24) + S14*S26 + S26*S34 - S24*S36 + 
     -        S12*S46 + S23*S46 - 2*S26*S46) + 
     -     S13*(2*S14**2*S26 + 
     -        S14*(-(S16*S24) + 3*S26*S34 - S24*S36 + S12*S46 + 
     -           S23*S46 - 2*S26*S46) + 
     -        S26*S34*(-2*S16 + S34 - 2*(S36 + S46))))/
     -   (2.*S13*S14*(S13 + S14 + S34)**2)      


      LEIGAXA = LEIGAXA +
     -  (A6*(2*S13**4*S23*S24 + 
     -       S14**2*S23**2*(2*S14**2 + 2*S23*S26 + S23*S34 + 
     -          2*S26*S34 + S34**2 - 2*S16*(S23 + 2*S24 + S34) + 
     -          S14*(-4*S16 + S23 + 4*S26 + 3*S34 - 4*S36) - 
     -          2*S23*S36 - 2*S24*S36 - 2*S34*S36) + 
     -       2*S12**2*(S13**2*S23*S24 - S13*S14*S23*(S23 + 2*S34) - 
     -          S14*S23*S34*(S14 - 2*S16 + 2*S26 + S34 - 2*S36) + 
     -          S13*S24*S34*(2*S16 + 2*S23 + S24 - 2*S26 + S34 + 2*S36))
     -         + S13**3*S24*(3*S14*S23 + S23*S24 - S24**2 - 
     -          2*S16*(2*S23 + S24) + 4*S23*S26 + 10*S24*S26 + 
     -          4*S23*S34 - S24*S34 + 2*S24*S36 - 4*S23*S46) - 
     -       S13*S23*(S14**3*(-3*S23 + 2*S24) + 
     -          4*S16*S24*(S23 + 4*S24)*S34 + 
     -          S14**2*(S23**2 + 2*S16*(S23 - 2*S24) + 
     -             S24*(S34 - 6*S36) + S23*(5*S24 - 6*S26 - 4*S46)) + 
     -          S14*(2*S23**3 + 4*S16*(S23**2 + S23*S34 + S24*S34) - 
     -             S24*(2*S24**2 + 12*S26*S34 + 
     -                S24*(-4*S26 + S34 + 2*S36)) - 
     -             S23*(S24**2 - 7*S24*S34 + 4*S26*S34 - 2*S34**2 - 
     -                2*S24*S36 + 4*S34*S46) + 
     -             S23**2*(3*S24 - 4*(S26 - S34 + S46)))) + 
     -       S13**2*(S14**2*S23*(S23 - S24) + 
     -          S14*(4*S23**2*S26 + 2*S24**2*(2*S26 + S36) + 
     -             S23*(5*S24**2 + 4*S26*S34 + 
     -                S24*(7*S34 + 2*S36 - 8*S46))) - 
     -          S24*(S23**2*S24 + 
     -             2*S16*(2*S23**2 + S24*(2*S24 + S34) + 
     -                S23*(7*S24 + 2*S34)) + 
     -             S24*(2*S24**2 + S34*(-10*S26 + S34 - 2*S36) + 
     -                S24*(-4*S26 + 3*S34 + 4*S36)) + 
     -             S23*(3*S24**2 - 2*S34*(2*S26 + S34 - 2*S46) - 
     -                2*S24*(S26 + 2*S46))))))/
     -   (8.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34)) 

      LEIGAXA = LEIGAXA + 
     -       (A6*(S12*(2*S14**3*S23*(S23 - S34) + 
     -          S14**2*S23*(-4*S16*(S23 - S34) + 
     -             S23*(4*S26 + S34 - 4*S36) - 
     -             S34*(5*S13 + 4*S26 + 3*S34 - 4*S36)) + 
     -          S13*S24*(4*S13**2*S23 + 
     -             S13*(-2*S24**2 + 4*S23*S26 + 7*S23*S34 - 
     -                10*S26*S34 + S34**2 + 
     -                S16*(-8*S23 - 4*S24 + 2*S34) + 
     -                S24*(4*S26 - S34 - 4*S36) - 2*S34*S36 + 8*S23*S46)
     -               + S34*(S23**2 + 2*S24**2 - 4*S24*S26 + 3*S24*S34 - 
     -                10*S26*S34 + S34**2 + 
     -                2*S16*(-5*S23 + 2*S24 + S34) + 4*S24*S36 - 
     -                2*S34*S36 + S23*(5*S24 + 2*S26 + 4*(S34 + S46))))
     -           + S14*(S13**2*S23*(6*S24 - S34) + 
     -             S23*S34*(-2*S26*S34 - S34**2 + 
     -                2*S16*(S23 + 2*S24 + S34) - 
     -                S23*(2*S26 + S34 - 2*S36) + 2*S24*S36 + 2*S34*S36)
     -               - S13*(4*S23**3 + 
     -                2*S16*S23*(2*S23 + 2*S24 + S34) + 
     -                2*S24*S34*(2*S26 + S36) + 
     -                S23**2*(6*S24 - 8*S26 + 7*S34 - 8*S46) - 
     -                2*S23*(2*S24*(S26 - 2*S36) + 
     -                   S34*(S26 - 2*S34 + 2*S46)))))))/
     -   (8.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34))

      LEIGAXA = LEIGAXA +  
     -  (-2*S13**3*S24*(S16*S23 + 2*S26**2 + S23*S46) + 
     -     S14*S23*(-(S14**2*
     -           (2*S16*S23 + S23*S26 + S26*S34 - 2*S24*S36)) + 
     -        S16*(2*S16*(2*S23**2 + 3*S23*S24 + 2*S23*S34 + S24*S34) - 
     -           S24*(-2*S26*S34 + S34**2 + 2*S24*S36 - 2*S34*S36 + 
     -              S23*(-2*S26 + S34 + 2*S36))) + 
     -        4*S12**2*S36*(-S16 + S26 + S46) + 
     -        S14*(4*S16**2*S23 - 2*S26**2*S34 - S26*S34**2 + 
     -           4*S12*S26*S36 + 2*S24*S26*S36 + 2*S24*S34*S36 + 
     -           2*S26*S34*S36 - 4*S24*S36**2 - 
     -           S16*(-2*S12*S34 + S24*S34 - 2*S26*S34 + 
     -              S23*(S24 + 2*(S26 + S34)) + 4*S12*S36 + 4*S24*S36)
     -            + S12*S34*S46 + 4*S12*S36*S46 + 
     -           S23*(-2*S26**2 - S26*S34 + 6*S26*S36 + S12*S46 + 
     -              4*S36*S46)) + 
     -        S12*(4*S16**2*S23 + 
     -           (2*S26*S34 + S34**2 + 2*S24*S36 - 2*S34*S36)*S46 + 
     -           S23*((S34 + 2*S36)*S46 + 2*S26*(2*S36 + S46)) + 
     -           2*S16*(2*S26*S34 + S34**2 - 2*S23*S36 - 4*S24*S36 - 
     -              S23*S46 - S34*(2*S36 + S46)))) + 
     -     S13**2*(S14**2*S23*S26 + 
     -        S14*S24*(-(S16*(S23 + 2*S26)) + 
     -           S26*(S24 + 2*S26 + S34 + 2*S36) + 
     -           S23*(3*S26 + 2*S36 - 4*S46)) - 
     -        S12*(S14*S23*S46 + S24**2*S46 + 
     -           2*S16*S24*(S23 + 2*S36 + S46) + 
     -           S24*(4*S26**2 - 4*S26*S36 + 3*S23*S46 + 6*S26*S46 + 
     -              S34*S46 + 2*S36*S46)) + 
     -        S24*(2*S16**2*(2*S23 + S24) + 
     -           S16*(S24*(S24 - 6*S26 + S34 - 6*S36) + 
     -              S23*(S24 + 4*S26 - 2*S34 + 8*S46)) + 
     -           2*(S23*S24*S26 + S24**2*S26 - 2*S26**2*S34 + 
     -              S23*S46*(-S34 + 2*S46) + 
     -              S24*(-2*S26**2 + S26*(S34 - 8*S46) - 2*S36*S46)))))/
     -   (8.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34))

      LEIGAXA = LEIGAXA +
     -       (S13*(S14**3*S23*S26 + 
     -        S14**2*(S16*S23*(-2*S23 + S24 - 2*S26) - S23**2*S26 + 
     -           2*S24*S26*S36 + 
     -           S23*(2*S26**2 + 2*S24*(S26 + 2*S36 - S46) - S12*S46))
     -         + S14*(-2*S16**2*S23*S24 + 2*S23**3*S26 + 
     -           S23**2*S24*S26 - S23*S24**2*S26 - 4*S23**2*S26**2 + 
     -           2*S23*S24*S26**2 + 2*S23**2*S26*S34 + 
     -           2*S23*S24*S26*S34 + S24**2*S26*S34 - 
     -           10*S24*S26**2*S34 + S24*S26*S34**2 - 
     -           2*S23**2*S24*S36 - 4*S23*S24**2*S36 - 2*S24**3*S36 + 
     -           12*S24**2*S26*S36 - 2*S24**2*S34*S36 - 
     -           2*S24*S26*S34*S36 + 4*S24**2*S36**2 - 
     -           S16*(S23**2*(3*S24 + 4*S26) + 
     -              2*S24*(S26*S34 + S24*S36) + 
     -              2*S23*(S24**2 + 2*S26*S34 + 
     -                 S24*(4*S26 + S34 - 4*S46))) + 2*S23**3*S46 + 
     -           4*S23**2*S24*S46 + 2*S23*S24**2*S46 - 
     -           8*S23**2*S26*S46 + 2*S23**2*S34*S46 - 
     -           4*S23**2*S46**2 + 
     -           S12*(-2*S24*S36*S46 + 2*S16*S23*(-S24 + S34 + S46) + 
     -              S23**2*(2*S26 + 3*S46) + 
     -              2*S23*(S24*S26 - 2*S26**2 + S26*(S34 - 5*S46) + 
     -                 (S34 - 2*S46)*S46))) + 
     -        S24*(4*S12**2*S36*(-S16 + S26 + S46) + 
     -           S16*(-(S23**2*S24) + 2*S16*S24*(3*S23 + S34) + 
     -              S24*(S24*(S34 - 8*S36) + 
     -                 S34*(6*S26 + S34 - 2*S36)) + 
     -              S23*(-S24**2 + 4*S26*S34 + 2*S24*(S26 + 4*S46))) + 
     -           S12*(4*S16**2*S23 - 2*S24*S26*S34 - 2*S26*S34**2 + 
     -              4*S24*S26*S36 + 4*S26*S34*S36 + S23**2*S46 - 
     -              S24*S34*S46 + 10*S26*S34*S46 - S34**2*S46 + 
     -              4*S24*S36*S46 + 2*S34*S36*S46 - 
     -              2*S16*(-6*S26*S34 + 6*S24*S36 + S23*(S34 - S46) + 
     -                 S34*S46) - 
     -              S23*(2*S26*(S34 + S46) + S46*(-S24 + 2*S34 + 4*S46))
     -              ))))/
     -   (8.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34))

      END
C-----------------------------------------------------------------------
      FUNCTION A6COEF(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION A6COEF,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=2*DOT(P,6,7)

      A6COEF=(S14**2*S23*S26 + (S13*S24 - S12*S34)*(S16*S24 - S12*S46) - 
     -    S14*(S16*S23*S24 + S13*S24*S26 + S12*S26*S34 - 
     -       2*S12*S24*S36 + S12*S23*S46))/
     -  (S14**2*S23**2 + (S13*S24 - S12*S34)**2 - 
     -    2*S14*S23*(S13*S24 + S12*S34))

      END
C-----------------------------------------------------------------------
      FUNCTION LEIB(P,Q,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIB,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIB = 0
      LEIB = LEIB + P1K*P2K * (  - 4*Q2*P14*P24*P45 - 2*P12
     +     *P15*P25*P34 - P12*P15*P25*P35 - P12*P15*P25*P45 + 8*P13*P15*
     +     P23*P24 + 4*P13*P15*P24*P34 + 4*P13*P15*P24**2 + 4*P14*P15*
     +     P23*P24 + 2*P14*P15*P24*P45 + 2*P14*P24*P25*P45 - 4*P15*P23*
     +     P24*P25 + 4*P15*P23*P24*P34 )
      LEIB = LEIB + P1K*P3K * (  - P12*P14*P15*P25 - 2*P12*
     +     P14*P24*P45 + P12*P15*P24*P25 + 4*P12*P15*P24*P34 - P12**2*
     +     P15*P25 - 4*P13*P15*P24**2 + 4*P14*P15*P23*P24 )
      LEIB = LEIB + P1K*P4K * (  - P12*P13*P15*P25 - 2*P12*
     +     P14*P24*P45 + P12*P15*P23*P25 - P12**2*P15*P25 + 4*P13*P15*
     +     P23*P24 + 4*P14*P15*P23*P24 )
      LEIB = LEIB + P1K**2 * ( 2*Q2*P14*P24*P45 - 2*P12*P14*
     +     P24*P45 + P12*P15**2*P25 - 2*P14*P24*P25*P45 - 4*P15*P23*P24*
     +     P34 )
      LEIB = LEIB + P2K*P3K * ( 4*P12*P13*P15*P24 + P12*P14*
     +     P15*P25 - 2*P12*P14*P24*P45 - P12*P15*P24*P25 - P12**2*P15*
     +     P25 + 4*P13*P15*P23*P24 + 4*P13*P15*P24**2 )
      LEIB = LEIB + P2K*P4K * (  - 4*P12*P13*P15*P24 + P12*
     +     P13*P15*P25 - 8*P12*P14*P15*P24 - 2*P12*P14*P24*P45 - P12*P15
     +     *P23*P25 - 4*P12*P15*P24*P34 - P12**2*P15*P25 - 4*P13*P14*P15
     +     *P24 + 4*P13*P15*P23*P24 + 4*P13*P15*P24**2 - 4*P13**2*P15*
     +     P24 )
      LEIB = LEIB + P2K**2 * ( 2*Q2*P14*P24*P45 - 2*P12*P14*
     +     P24*P45 + P12*P15*P25**2 - 4*P13*P14*P15*P24 - 4*P13*P15*P24*
     +     P34 - 4*P13**2*P15*P24 - 2*P14*P15*P24*P45 )
      LEIB = LEIB + P3K*P4K * (  - 2*P12*P25 - 4*P14*P24 )
     +     *P12*P15
      LEIB = LEIB + P3K**2 * (  - 4*P12*P14*P15*P24 )
c$$$      LEIB = LEIB + KK * ( 2*Q2*P12*P14*P24*P45 + Q2*
c$$$     +     P12**2*P15*P25 - 2*Q2*P13*P15*P23*P24 - 2*Q2*P14*P15*P23*
c$$$     +     P24 + 2*Q2*P15*P23*P24*P25 + P12*P13*P14*P15*P25 - 4*P12*
c$$$     +     P13*P15*P23*P24 - 2*P12*P13*P15*P24*P34 + 2*P12*P14*P15*P24*
c$$$     +     P34 + 4*P12*P14*P15*P24**2 + P12*P15*P23*P24*P25 - 2*P12*P15*
c$$$     +     P23*P24*P34 + 2*P12*P15*P24**2*P34 - 2*P13*P14*P15*P23*P24 + 
c$$$     +     2*P13*P14*P15*P24**2 - 2*P13*P15*P23*P24**2 - 2*P13*P15*
c$$$     +     P24**3 + 2*P13**2*P15*P24**2 + 2*P14*P15*P23*P24**2 + 2*P14*
c$$$     +     P15*P23**2*P24 - 2*P14**2*P15*P23*P24 - 2*P15**2*P23*P24*P25
c$$$     +     )/2
      LEIB = LEIB/(2*P13*P14*P15*P23*P24*P25)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIGAXB(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIGAXB,P(4,7),DOT,A6,A6COEF,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=-2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=-2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=-2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=-2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=-2*DOT(P,6,7)
      A6=A6COEF(P,I,J,K,L)

      LEIGAXB= (-(S14*S16*S23**2*S24) + S14**2*S23**2*S26 - 
     -     4*S14*S16*S23**2*S26 + 
     -     S13**2*(S16*S24*(2*S23 + S24) + 
     -        S26*(-6*S23*S24 - S14*(2*S23 + S24) + 4*S24*S26)) - 
     -     6*S16**2*S23*S24*S34 + S16*S23**2*S24*S34 - 
     -     4*S16**2*S24**2*S34 + S16*S23*S24**2*S34 + 
     -     2*S14*S16*S23*S26*S34 + S14*S23**2*S26*S34 + 
     -     4*S14*S16*S24*S26*S34 - S14*S23*S24*S26*S34 - 
     -     2*S16*S23*S24*S26*S34 + 2*S14*S23*S26**2*S34 + 
     -     6*S14*S16*S23*S24*S36 + 4*S14*S16*S24**2*S36 - 
     -     2*S14**2*S23*S26*S36 - 4*S14*S23**2*S26*S36 - 
     -     4*S14**2*S24*S26*S36 + 
     -     2*S12**2*S24*(S13 + 2*S14 + S34)*S36 - 
     -     4*S14*S16*S23**2*S46 - 4*S14*S16*S23*S24*S46 - 
     -     4*S14*S23**2*S26*S46 + 
     -     S13*(-4*S16**2*S23*S24 - S14*S16*S23*(S24 - 4*S26) - 
     -        S14**2*S23*S26 + 4*S24*S26*(S26*S34 - S24*(S36 + S46) + 
     -           S23*(-S34 + S36 + S46)) + 
     -     S16*S24*(2*S23**2 + S23*(S24 + 4*S26 + S34 - 4*S36 + 4*S46) + 
     -           S24*(-8*S26 + S34 - 2*S36 + 4*S46)) + 
     -        S14*S26*(2*S23**2 + S24*(8*S26 - S34 + 2*S36 + 4*S46) - 
     -           S23*(5*S24 + S34 - 4*(S36 + S46)))) - 
     -     S12*(-2*S14**2*S24*S36 + 
     -        S13**2*(4*S23*S26 + 2*S24*S26 + 2*S23*S46 - S24*S46) + 
     -        S13*(4*S23*S26*S34 + 4*S24*S26*S34 + 2*S23*S24*S36 - 
     -           2*S24**2*S36 - 8*S23*S26*S36 + 4*S24*S26*S36 + 
     -           2*S16*(S24*S34 + S23*(S24 - 4*S26 - 2*S46)) + 
     -           3*S23*S24*S46 - 4*S23*S26*S46 + S23*S34*S46 - 
     -           S24*S34*S46 - 4*S23*S36*S46 + 2*S24*S36*S46 + 
     -           4*S24*S46**2 + 
     -           S14*(2*S23*S26 + 4*S24*S26 - 4*S24*S36 + S23*S46 - 
     -              2*S24*S46)) + S14*(2*S16*S24*(2*S23 + S34) + 
     -           2*S24*(S36*(-S24 - S34 + 2*S36) + S26*(S34 + 4*S36)) - 
     -           S23**2*S46 + 2*S23*(S26*S34 + (S24 - S36)*S46)) + 
     -        S34*(2*S16*S23*(S24 - 2*S26 - S46) + 
     -           2*S23*S26*(S34 - 2*S36 - S46) - S23**2*S46 + 
     -           S23*S24*(2*S36 + 3*S46) - 
     -           2*S24*(-(S26*S34) + S24*S36 + 2*S26*S46) + 
     -           2*S16*S24*(-4*S26 + S34 - 2*(S36 + S46)))))/
     -   (4.*S13*S14*S23*S24*(S13 + S14 + S34))

      LEIGAXB = LEIGAXB + 
     -  (A6*(-(S13**3*S24*(2*S23 + S24)) + 
     -       2*S12**2*(2*S13**2*S23 + S13*S14*S23 + S14*S23*S34 + 
     -          S13*(2*S23 + S24)*S34 + (S23 + S24)*S34**2) - 
     -       S14*S23*(-2*S16*(S23 + 2*S24)*S34 + 
     -          S23*(S23 + S24 - 2*S26)*S34 + 
     -          S14*(S23**2 - 4*S23*S26 + 2*S23*S36 + 4*S24*S36)) + 
     -       S13**2*(2*S14*S23*(S23 + S24 - 2*S26) + 
     -          S24*(4*S16*S23 + 4*S23**2 + 
     -             S23*(5*S24 - 8*S26 - S34 - 8*S46) + 
     -             S24*(4*S26 - S34 + 2*S36 - 4*S46))) + 
     -       S13*(S14**2*S23*(S23 + 2*S24 - 4*S26) + 
     -          S24*(2*S16*(5*S23 + 2*S24) + 3*S23*(S23 + S24 - 2*S26))*
     -           S34 + S14*(-2*S23**3 + 
     -             S23*(4*S16*S24 + 4*S24**2 - 12*S24*S26 + S24*S34 - 
     -                8*S26*S34) + 4*S24*(-2*S26*S34 + S24*S36) + 
     -             S23**2*(2*S24 + 4*S26 + S34 + 4*S46))) + 
     -       S12*(-2*S14**2*S23**2 + 
     -          (-2*S16*(S23 + 2*S24) + S23*(S23 + S24 - 2*S26))*
     -           S34**2 - S14*S34*
     -           (S23**2 + 2*S23*(S24 + 2*S26 - S36) - 4*S24*S36) + 
     -          S13**2*(4*S14*S23 + 4*S23**2 + S24*(-2*S24 + S34) + 
     -             2*S23*(4*S24 - 4*S26 + S34 - 2*S46)) + 
     -          S13*(2*S14**2*S23 + S14*S23*(6*S24 - 4*S26 + 3*S34) + 
     -             S34*(4*S23**2 + S23*(5*S24 - 8*S26 + S34) + 
     -                S24*(-2*S24 - 4*S26 + S34 - 2*S36 + 4*S46))))))/
     -   (4.*S13*S14*S23*S24*(S13 + S14 + S34)) 

      LEIGAXB = LEIGAXB -    
     -  (A6*(S12 + S23 + S24 - 2*S26))/(2.*S14) + 
     -  (2*S13**3*S26 + S13**2*
     -      (4*S14*S26 - S16*(S23 + S24 + 2*S26) + 4*S26*S34 + 
     -        S12*S36 - 4*S26*S36 + S12*S46 - 4*S26*S46) + 
     -     S34*(-(S14*S16*(S23 + S24)) + S14**2*S26 - 
     -        S16*(S23 + S24 - 2*S26)*S34 + S12*S14*(S36 + S46) + 
     -        S12*S34*(S36 + S46) + S14*S26*(S34 - 2*(S36 + S46))) + 
     -     S13*(2*S14**2*S26 - 
     -        S14*(S16*(S23 + S24 + 2*S26) - 5*S26*S34 - 
     -           S12*(S36 + S46) + 4*S26*(S36 + S46)) + 
     -        2*S34*(-(S16*(S23 + S24)) + S12*(S36 + S46) + 
     -           S26*(S34 - 2*(S36 + S46)))))/
     -   (2.*S13*S14*(S13 + S14 + S34)**2)

      LEIGAXB = LEIGAXB +
     -  (A6*(-(S13**4*S23*S24)+S13**3*(S14*(S23**2 - S23*S24 + S24**2) + 
     -          S23*S24*(2*S16 + 2*S23 + S24 + 2*S26)) + 
     -       S14**2*S23*S24*(-(S23*(S23 + S24 - 2*S26)) + 
     -          2*S16*(S23 + 2*S34) + S14*(S23 - 4*S36)) - 
     -       2*S12**2*(S13**2*S23*(S24 + S34) + 
     -          S14*S24*S34*(S14 - 2*S16 - S23 - S24 + 2*S26 - 2*S36) - 
     -          S13*(S14*(S23**2 - S24*S34) + 
     -             S23*S34*(2*S16 + S23 - 2*S26 + 2*S36))) + 
     -       S13**2*(S14**2*(3*S23**2 + S23*S24 + S24**2) - 
     -          S23*S24*(6*S16*S23 + S23**2 + 
     -             2*(S24**2 - 2*S24*S26 + S24*S34 - 2*S26*S34 + 
     -                2*S24*S36) + S23*(3*S24 - 2*S26 - 4*S46)) - 
     -          S14*(2*S23**3 + 2*S16*(S23**2 + S24**2) + 
     -             2*S23**2*S26 + 2*S23*S24*(S24 - 2*S26 + S34) + 
     -             S24**2*(S24 - 6*S26 + 4*S46))) + 
     -       S13*(S14**3*S23*(2*S23 + S24) + 
     -          4*S16*S23*S24*(-2*S23 + S24)*S34 + 
     -          S14*(S16*(6*S23**3 - 4*S23**2*S24 - 6*S23*S24**2 + 
     -                4*S23*S24*S34 + 4*S24**2*S34) + 
     -             S23*(S23**3 - S24*(S24**2 - 2*S24*S26 + 4*S26*S34) + 
     -                S23*(-S24**2 + 2*S24*S34 + 4*S26*S34 + 
     -                   4*S24*S36) + S23**2*(S24 - 2*S26 - 4*S46))) - 
     -          S14**2*(S23**3 + 2*S16*S23*(2*S23 + S24) - 
     -             2*S23**2*(S24 + S34) - 4*S24*(-2*S26*S34 + S24*S36) + 
     -             2*S23*(S24*S26 + 2*S26*S34 - 2*S24*S46))) + 
     -       S12*(S13**3*S23*(2*S23 + S34) + 
     -          S14*S24*(2*S14**2*S23 + S23*(S23 + S24 - 2*S26)*S34 - 
     -             2*S16*S34*(S23 + 2*S34) - 
     -             S14*(4*S16*S23 + 2*S23**2 - 4*S34*S36 + 
     -                S23*(2*S24 - 4*S26 + S34 + 4*S36))) + 
     -          S13**2*(S14*(4*S23**2 + 2*S23*S24 + 2*S24**2 - 
     -                S23*S34 - S24*S34) - 
     -             S23*(2*S23**2 + 4*S24**2 - 4*S24*S26 + 5*S24*S34 + 
     -                2*S26*S34 + 2*S16*(2*S23 + S34) + 4*S24*S36 + 
     -                4*S23*(S24 - S26 - 2*S46))) + 
     -          S13*(S14**2*
     -              (4*S23**2 + 6*S23*S24 + S24*(2*S24 - S34)) + 
     -             S23*S34*(-6*S16*S23 - S23**2 + 8*S16*S24 + S23*S24 + 
     -                2*S23*S26 - 4*S24*S26 - 4*S26*S34 + 4*S24*S36 + 
     -                4*S23*S46) + 
     -             S14*(S16*(-4*S23**2 + 2*S24*(-2*S24 + S34) + 
     -                   S23*(-8*S24 + 4*S34)) + 
     -                S23**2*(-2*S24 + 5*S34 - 4*S36) + 
     -                S24*(-2*S24**2 - 6*S26*S34 + 
     -                   S24*(4*S26 + S34 - 4*S36) + 4*S34*S46) + 
     -                S23*(-6*S24**2 - 8*S26*S34 + 8*S24*(S26 + S46)))))
     -       ))/(8.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34))

      LEIGAXB = LEIGAXB +       
     -     (S13**3*S23*(S16*S24 + S14*S26) + 
     -     S13**2*(S14**2*(S23 + S24)*S26 - 
     -        S14*(S16*(2*S23**2 + S24**2 + S23*(S24 + 2*S26)) + 
     -           S23*S26*(-2*S23 + S24 - 2*(S26 + S34 - 2*S46))) - 
     -        S23*S24*(2*S16**2 + S16*(2*S23 + S24 + 2*S26) - 
     -           2*S26*(S23 + S24 - 2*S26 + S34 - 6*S46))) + 
     -     2*S12**2*S36*(S14*S24*
     -         (S14 - 2*S16 - S23 - S24 + 2*S26 + 2*S46) + 
     -        S13*(-2*S16*S23 + S14*S24 + 2*S23*(S26 + S46))) + 
     -     S13*(2*S14*S16**2*(2*S23**2 + 2*S23*S24 + S24**2) + 
     -        S14**3*S24*S26 - 
     -        S14**2*S16*(2*S23**2 + 2*S23*S24 + S24*(S24 + 2*S26)) + 
     -        S16*S23*S24*(6*S16*S23 - S23**2 + 4*S26*S34 - 8*S24*S36 - 
     -           S23*(S24 - 2*S26 + 2*S34 - 8*S46)) + 
     -        S14*S23*S26*(-S23**2 + 
     -           2*(S24**2 - 2*S24*S26 - 2*S26*S34 + 4*S24*S36) + 
     -           S23*(S24 + 2*S26 + 4*S46)) + 
     -        S14**2*S26*(S23**2 - 4*S23*(S24 - S36) + 
     -           S24*(-S24 + 6*S26 + 4*S46)) + 
     -        S14*S16*(4*S23*S24*(S24 - S36) + 
     -           S24**2*(S24 - 6*S26 + 4*S46) + 
     -           S23**2*(S24 - 6*S26 - 2*S34 + 4*S46))) + 
     -     S14*S24*(S16*S24*(-(S23*(S23 + S24 - 2*S26)) + 
     -           2*S16*(S23 - 2*S34)) - S14**2*S26*(S23 + 4*S36) + 
     -        S14*(S23*(S23 + S24 - 2*S26)*S26 + 
     -           4*S16*(S26*S34 + S24*S36) + 
     -           S16*S23*(S24 - 2*(S26 + 4*S46)))) + 
     -     S12*(-(S13**3*S23*S46) - 
     -        S13**2*(-2*S23**2*S26 + 2*S16*S23*(S23 - S46) + 
     -           S14*S24*(2*S26 + S46) + 
     -           S23*(4*S26**2 - 2*S26*S34 + S14*S46 + 10*S26*S46 - 
     -              S24*(2*S26 + S46))) + 
     -        S14*S24*(4*S16**2*S23 + 
     -           2*S16*(S23**2 + 4*S26*S34 - 2*S24*S36 + 
     -              S23*(S24 - S46) + 2*S34*S46) + 
     -           S14*(-2*S16*S23 + 4*S36*(-S26 + S46) + 
     -              S23*(2*S36 + S46)) + 
     -           S23*((S23 + S24)*S46 + 2*S34*(S36 + S46) + 
     -              2*S26*(2*S36 + S46))) - 
     -        S13*(S14**2*S24*(2*S26 + S46) + 
     -           S14*(-2*S23*S26*(S24 + 2*S36) + S23**2*S46 + 
     -              S24*(-2*S24*S26 + 4*S26**2 + 2*S24*S36 + 
     -                 2*S34*S36 + S24*S46 + 6*S26*S46 + 2*S34*S46 + 
     -                 4*S46**2) + 
     -              2*S16*(S23**2 + S23*(S24 + 2*S36) + 
     -                 S24*(2*S36 + S46))) - 
     -           S23*(4*S16**2*S23 + 
     -              S24*(4*S26*S36 + (S23 + 4*S36)*S46) + 
     -              2*S16*(4*S26*S34 - 6*S24*S36 + S23*(-S34 + S46)) + 
     -              S46*(S23**2 + 4*S26*S34 - 2*S23*(S26 + 2*S46))))))/
     -   (8.*S13*S14*S23*S24*(S13 + S14 + S34)*(S23 + S24 + S34)) 

      LEIGAXB = LEIGAXB + 
     -  (-(S14*S16*S23*S24) + 6*S16**2*S23*S24 + S16*S23**2*S24 + 
     -     S16*S23*S24**2 + S14**2*S23*S26 - 6*S14*S16*S23*S26 - 
     -     S14*S23**2*S26 - S14*S23*S24*S26 - 2*S16*S23*S24*S26 + 
     -     2*S14*S23*S26**2 + S13**2*(S16*S24 - S14*S26) - 
     -     4*S16**2*S24*S34 + 4*S14*S16*S26*S34 + 4*S14*S16*S24*S36 - 
     -     4*S14**2*S26*S36 - 8*S14*S16*S23*S46 + 
     -     2*S12**2*S36*(S14 - 2*S16 - S24 + 2*S26 + 2*S46) + 
     -     S13*(-(S14**2*S26) + S14*S16*(S24 + 2*S26) - 
     -        S16*S24*(2*S16 + 2*S23 + S24 + 2*S26 - 4*S46) + 
     -        S14*S26*(2*S23 + S24 + 2*S26 + 4*S46)) + 
     -     S12*(4*S16**2*S23 + 2*S16*S23*S24 + 4*S16**2*S34 + 
     -        2*S16*S24*S34 + 8*S16*S26*S34 + 2*S23*S26*S34 + 
     -        2*S24*S26*S34 - 4*S26**2*S34 + 2*S14**2*S36 - 
     -        8*S16*S24*S36 - 2*S24**2*S36 + 4*S24*S26*S36 + 
     -        S13**2*S46 - 2*S16*S23*S46 - S23**2*S46 + S23*S24*S46 + 
     -        2*S23*S26*S46 + 4*S16*S34*S46 + 
     -        S14*(-2*S26*S34 - 2*S16*(S23 + S34 + 2*S36) + S23*S46 + 
     -           4*S36*S46) - 
     -        S13*(2*S14*S26 - 2*S23*S26 - 2*S24*S26 + 4*S26**2 + 
     -           S14*S46 + S24*S46 + 6*S26*S46 + 4*S46**2 + 
     -           2*S16*(S23 + S34 + S46))))/(8.*S13*S14*S23*S24) 

      LEIGAXB = LEIGAXB + 
     -  (A6*(-(S13**3*S24) + 
     -       2*S12**2*S34*(S13 + S14 - 2*S16 - S23 - S24 + 2*S26 + 
     -          2*S36) + S14*S23*
     -        (S23*(S23 + S24 - 2*S26) + S16*(6*S23 + 4*S34) - 
     -          S14*(S23 + 4*S36)) + 
     -       S13**2*(S14*(S23 - S24) + 
     -          S24*(2*S16 + 2*S23 + S24 + 2*S26 - 4*S46)) + 
     -       S13*(S14**2*S23 - 
     -          S24*(S23*(S23 + S24 - 2*S26) + S16*(6*S23 - 4*S34)) - 
     -          2*S14*(S16*S23 + S23**2 + S23*S26 + 4*S26*S34 - 
     -             2*S24*S36 - 2*S23*S46)) + 
     -       S12*(-2*S14**2*S23 + S13**2*(2*S23 - 2*S24 + S34) - 
     -          S34*(S23*(S23 + S24 - 2*S26) + S16*(6*S23 + 4*S34)) + 
     -          S14*(4*S16*S23 + 2*S23**2 + 
     -             S23*(2*S24 - 4*S26 + S34 - 4*S36) + 4*S34*S36) - 
     -          S13*(2*S23**2 + 2*S14*S24 - 2*S24**2 + 4*S24*S26 - 
     -             S14*S34 + S24*S34 + 2*S26*S34 + 
     -             2*S16*(2*S23 - 2*S24 + S34) + 4*S24*S36 - 
     -             4*S34*S46 - 4*S23*(S26 + 2*S46)))))/
     -   (8.*S13*S14*S23*S24)  

      END
C-----------------------------------------------------------------------
      FUNCTION LEIC(P,Q,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIC,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIC = 0
      LEIC = LEIC + P1K*P2K * (  - 4*P12*P13*P15*P25*P34 + 6*
     +     P12*P13*P24*P25*P34 - 2*P12*P14*P15*P25*P34 - 2*P12*P14*P24*
     +     P25*P34 - 2*P12*P15*P25*P34**2 + 2*P12*P24*P25*P34**2 + 4*P13
     +     *P14*P15*P24*P34 - 2*P13*P14*P15*P25*P34 - 3*P13*P14*P23*P24*
     +     P25 - P13*P14*P24**2*P25 + P13*P14**2*P15*P24 - 2*P13*P15*P23
     +     *P25*P34 + P13*P15*P24*P25*P34 - P13*P15*P24*P34**2 - P13*P15
     +     *P25*P34**2 + P13*P23*P24*P25*P34 - 2*P13*P24**2*P25*P34 + 6*
     +     P13**2*P14*P15*P24 - P13**2*P23*P24*P25 + P13**2*P24*P25*P34
     +     - 3*P13**2*P24**2*P25 + P13**3*P15*P24 - P14*P15*P23*P25*P34
     +     - 2*P14*P15*P24*P25*P34 - 2*P14*P15*P25*P34**2 - 3*P14*P23*
     +     P24*P25*P34 + P14*P24*P25*P34**2 - 2*P14*P24**2*P25*P34 - 
     +     P14**2*P15*P25*P34 - P14**2*P24*P25*P34 + P15*P23*P25*P34**2
     +     + P15*P24*P25*P34**2 + P15*P25*P34**3 + P24*P25*P34**3 - 
     +     P24**2*P25*P34**2 )
      LEIC = LEIC + P1K*P3K * ( 3*P12*P13*P34 - P12*P14*P34
     +     - P12*P34**2 + P13*P14*P24 + 2*P13*P15*P34 - P13**2*P24 - 
     +     P14*P15*P34 + 3*P14*P24*P34 + P15*P34**2 - 3*P24*P34**2 )
     +     *P24*P25
      LEIC = LEIC + P1K*P4K * (  - 2*P12*P13*P15*P34 + 3*P12*
     +     P13*P24*P34 - P12*P14*P15*P34 - P12*P14*P24*P34 - 3*P12*P15*
     +     P34**2 + 3*P12*P24*P34**2 - P13*P14*P23*P24 - P13*P15*P23*P34
     +     - 2*P13*P15*P24*P34 - 2*P13*P23*P24*P34 + P13**2*P23*P24 + 2
     +     *P14*P15*P23*P34 - 3*P14*P23*P24*P34 + P23*P24*P34**2 )*P25
      LEIC = LEIC + P1K**2 * ( 2*P13*P15 - 3*P13*P23 - 3*P13*
     +     P24 + P14*P15 + P14*P23 + P14*P24 + 3*P15*P34 + P23*P34 - 3*
     +     P24*P34 )*P24*P25*P34
      LEIC = LEIC + P2K*P3K * (  - 2*P12*P13*P15*P25*P34 + 
     +     P12*P13*P24*P25*P34 - P12*P14*P15*P25*P34 - 3*P12*P14*P24*P25
     +     *P34 - P12*P15*P25*P34**2 + P12*P24*P25*P34**2 + 2*P13*P14*
     +     P15*P24*P34 - P13*P14*P15*P25*P34 - P13*P14*P23*P24*P25 - 2*
     +     P13*P14*P24*P25*P34 + P13*P14*P24**2*P25 - 2*P13*P14**2*P15*
     +     P24 - 4*P13*P15*P23*P25*P34 - 2*P13*P15*P25*P34**2 + 2*P13**2
     +     *P14*P15*P24 - 2*P14*P15*P23*P25*P34 - 3*P14*P15*P24*P25*P34
     +     - 2*P14*P15*P25*P34**2 - 3*P14*P23*P24*P25*P34 - P14*P24**2*
     +     P25*P34 - P14**2*P24*P25*P34 )
      LEIC = LEIC + P2K*P4K * (  - P12*P13*P15*P25*P34 + 2*
     +     P12*P13*P24*P25*P34 - 2*P12*P14*P15*P25*P34 - 2*P12*P14*P24*
     +     P25*P34 - 2*P12*P24*P25*P34**2 + 4*P13*P14*P15*P24*P34 - P13*
     +     P14*P15*P25*P34 + P13*P14*P24*P25*P34 - 2*P13*P15*P23*P25*P34
     +     - 3*P13*P15*P24*P25*P34 - P13*P15*P25*P34**2 + 3*P13*P23*P24
     +     *P25*P34 - 2*P13*P24*P25*P34**2 + P13*P24**2*P25*P34 + 2*
     +     P13**2*P14*P15*P24 - 2*P13**2*P15*P24*P34 + 2*P13**2*P15*P25*
     +     P34 + P13**2*P23*P24*P25 - P13**2*P24**2*P25 - 2*P13**3*P15*
     +     P24 - 2*P14*P15*P23*P25*P34 - 4*P14*P15*P24*P25*P34 )
      LEIC = LEIC + P2K**2 * ( 2*P13*P14*P15 + P13*P14*P24 + 
     +     P13*P15*P34 - P13*P24*P34 + 2*P13**2*P15 - P13**2*P24 + 2*P14
     +     *P15*P34 + 2*P14*P24*P34 + 2*P14**2*P15 + 2*P14**2*P24 )
     +     *P25*P34
      LEIC = LEIC + P3K*P4K * (  - P12*P13*P15*P34 + 2*P12*
     +     P13*P24*P34 - P12*P14*P15*P34 - P12*P15*P34**2 + 2*P12*P24*
     +     P34**2 - 2*P13*P14*P23*P24 - 2*P13*P15*P23*P34 - 2*P13*P24**2
     +     *P34 - 2*P13**2*P24**2 - 2*P14*P23*P24*P34 )*P25
      LEIC = LEIC + P3K**2 * ( 2*P13*P14*P24 + 2*P13*P15*P34
     +     + 2*P14*P24*P34 )*P24*P25
      LEIC = LEIC + P4K**2 * ( 2*P12*P15*P34 + 2*P13*P23*P24
     +     + 2*P23*P24*P34 )*P13*P25
c$$$      LEIC = LEIC + KK * (  - 4*P12*P13*P14*P15*P24*P34 + 4*
c$$$     +     P12*P13*P14*P15*P25*P34 + 3*P12*P13*P14*P23*P24*P25 - 2*P12*
c$$$     +     P13*P14*P24*P25*P34 + P12*P13*P14*P24**2*P25 - P12*P13*P14**2
c$$$     +     *P15*P24 + 4*P12*P13*P15*P23*P25*P34 + P12*P13*P15*P24*P34**2
c$$$     +     + 2*P12*P13*P15*P25*P34**2 - 2*P12*P13*P23*P24*P25*P34 - 2*
c$$$     +     P12*P13*P24*P25*P34**2 - 6*P12*P13**2*P14*P15*P24 + P12*
c$$$     +     P13**2*P23*P24*P25 - 4*P12*P13**2*P24*P25*P34 + 3*P12*P13**2*
c$$$     +     P24**2*P25 - P12*P13**3*P15*P24 + 2*P12*P14*P15*P23*P25*P34
c$$$     +     + 4*P12*P14*P15*P24*P25*P34 + 5*P12*P14*P15*P25*P34**2 + 6*
c$$$     +     P12*P14*P23*P24*P25*P34 + P12*P14*P24*P25*P34**2 + 4*P12*P14*
c$$$     +     P24**2*P25*P34 + 2*P12*P14**2*P15*P25*P34 + 2*P12*P14**2*P24*
c$$$     +     P25*P34 - P12*P15*P24*P25*P34**2 + P12*P15*P25*P34**3 - P12*
c$$$     +     P23*P24*P25*P34**2 + 3*P12*P24**2*P25*P34**2 + 4*P12**2*P13*
c$$$     +     P15*P25*P34 - 6*P12**2*P13*P24*P25*P34 + 2*P12**2*P14*P15*P25
c$$$     +     *P34 + 2*P12**2*P14*P24*P25*P34 + 2*P12**2*P15*P25*P34**2 - 2
c$$$     +     *P12**2*P24*P25*P34**2 )/4
c$$$      LEIC = LEIC + KK * (  - 2*P13*P14*P15*P23*P24*P34 + 2*
c$$$     +     P13*P14*P15*P23*P25*P34 + P13*P14*P15*P24*P25*P34 - 4*P13*P14
c$$$     +     *P15*P24**2*P34 + 3*P13*P14*P23*P24*P25*P34 - P13*P14*P23*
c$$$     +     P24**2*P25 + P13*P14*P23**2*P24*P25 - 3*P13*P14*P24**2*P25*
c$$$     +     P34 + 2*P13*P14**2*P15*P23*P24 + P13*P14**2*P23*P24*P25 + 2*
c$$$     +     P13*P15*P23*P24*P25*P34 + 4*P13*P15*P23*P25*P34**2 + 4*P13*
c$$$     +     P15*P23**2*P25*P34 - P13*P15*P24*P25*P34**2 + 3*P13*P15*
c$$$     +     P24**2*P25*P34 - 3*P13*P23*P24**2*P25*P34 + 2*P13*P24**2*P25*
c$$$     +     P34**2 - P13*P24**3*P25*P34 - 2*P13**2*P14*P15*P23*P24 - 2*
c$$$     +     P13**2*P14*P15*P24**2 - P13**2*P14*P23*P24*P25 - P13**2*P14*
c$$$     +     P24**2*P25 - 2*P13**2*P15*P24*P25*P34 + 2*P13**2*P15*P24**2*
c$$$     +     P34 - P13**2*P23*P24**2*P25 + P13**2*P24**2*P25*P34 + P13**2*
c$$$     +     P24**3*P25 + 2*P13**3*P15*P24**2 + P13**3*P24**2*P25 + 5*P14*
c$$$     +     P15*P23*P24*P25*P34 + P14*P15*P23*P25*P34**2 + 2*P14*P15*
c$$$     +     P23**2*P25*P34 + 4*P14*P15*P24**2*P25*P34 + 2*P14*P23*P24*P25
c$$$     +     *P34**2 + )/4
c$$$      LEIC = LEIC + KK * ( P14*P23*P24**2*P25*P34 + 3*P14*
c$$$     +     P23**2*P24*P25*P34 - P14**2*P15*P23*P25*P34 + 3*P14**2*P23*
c$$$     +     P24*P25*P34 )/4
      LEIC = LEIC/(2*P15*P25**2*P13*P34**2*P24)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIGAXC(P,I,J,K,L)
      IMPLICIT NONE
C     EVALUATE THE POLARIZED (QUARK) ERT A FUNCTION 
C     WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIGAXC,P(4,7),DOT,A6,A6COEF,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=-2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=-2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=-2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=-2*DOT(P,ABS(L),7)*SIGN(1,L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S67=-2*DOT(P,6,7)
      A6=A6COEF(P,I,J,K,L)

      LEIGAXC=  (A6*(-3*S13**3*S24 - 2*S12**2*(S13 + S14 - 3*S34)*S34 + 
     -       S14*S23*(-2*S14**2 + 
     -          S14*(4*S16 + 9*S23 + 4*S24 - 8*S26 - 3*S34 - 10*S36) + 
     -          S34*(16*S16 + S23 + 2*S24 - 4*S26 - S34 + 2*S36)) + 
     -       S13**2*(6*S23**2 + 8*S24**2 + S14*(S23 + 3*S24) - 
     -          8*S24*S26 + 2*S24*S34 - 4*S26*S34 + 6*S24*S36 + 
     -          S23*(11*S24 - 12*S26 + 6*S34 - 8*S46) - 8*S24*S46) + 
     -       S12*(2*S14**2*(S23 + S34) + 
     -          S13*S14*(4*S23 + 2*S24 + 5*S34) + 
     -          S13**2*(6*S23 + 4*S24 + 5*S34) + 
     -          S34**2*(-16*S16 - S23 - 2*S24 + 4*S26 + S34 - 2*S36) + 
     -          S14*S34*(-4*S16 - 15*S23 - 4*S24 + 8*S26 + 3*S34 + 
     -             10*S36) + 
     -          S13*S34*(-4*S16 - S23 - 14*S24 + 8*S26 + 4*S34 - 
     -             2*S36 + 16*S46)) + 
     -       S13*(S14**2*(-3*S23 + 4*S24) + 
     -          S34*(2*S23**2 + 4*S16*S24 + 3*S24*S34 - 4*S26*S34 + 
     -             S23*(3*S24 - 4*S26 + 2*S34) - 2*S24*S36) + 
     -          S14*(5*S23**2 + 4*S16*(S23 - 2*S24) - 2*S24**2 + 
     -             4*S24*S26 + 7*S24*S34 - 28*S26*S34 + 10*S24*S36 + 
     -             S23*(-5*S24 - 12*S26 + 2*S36 + 8*S46)))))/
     -   (16.*S13*S24*S34*(S13 + S14 + S34)) 

      LEIGAXC= LEIGAXC + 
     -  (7*S14*S16*S23*S24 - 9*S14**2*S23*S26 + 8*S14*S16*S23*S26 + 
     -     S14*S16*S24*S34 - 12*S16**2*S24*S34 - S16*S23*S24*S34 - 
     -     2*S16*S24**2*S34 - S14**2*S26*S34 + 12*S14*S16*S26*S34 - 
     -     S14*S23*S26*S34 + 2*S14*S24*S26*S34 + 4*S16*S24*S26*S34 - 
     -     4*S14*S26**2*S34 + S16*S24*S34**2 - S14*S26*S34**2 + 
     -     10*S14*S16*S24*S36 - 10*S14**2*S26*S36 + 8*S14*S23*S26*S36 - 
     -     2*S16*S24*S34*S36 + 2*S14*S26*S34*S36 + 2*S14**2*S23*S46 - 
     -     24*S14*S16*S23*S46 + 8*S14*S23*S26*S46 + 2*S14*S23*S34*S46 - 
     -     2*S12**2*(S13*S36 + S14*S36 + S34*(S36 + 4*S46)) + 
     -     S13**2*(3*S16*S24 + S14*S26 - 2*S23*(3*S26 + S46) - 
     -        4*(S26**2 + S26*(S34 + S46) + S24*(-2*S36 + S46))) + 
     -     S13*(-(S14**2*S26) - 6*S24**2*S36 - 
     -        4*S26*S34*(S26 + S34 - 2*S36 - S46) + 
     -        S14*(S16*S24 - 5*S23*S26 + 4*S24*S26 - 4*S26**2 - 
     -           6*S26*S34 + 4*S24*S36 + 2*S26*S36 + 4*S23*S46 - 
     -           4*S24*S46 + 12*S26*S46) - 
     -        4*S24*(2*S36**2 + S26*(S34 - 2*S36 + S46) + 
     -           S34*(-S36 + S46)) + 
     -        2*S23*((3*S24 + S34 + 4*S36)*S46 + 
     -           S26*(-S34 + 8*S36 + 4*S46)) + 
     -        S16*(-3*S23*S24 - 14*S24**2 + 8*S26*S34 + 
     -           8*S23*(2*S26 + S46) + 
     -           2*S24*(4*S26 + S34 - 5*S36 + 8*S46))) + 
     -     S12*(-2*S14**2*S36 + S13**2*(2*S26 + S46) + 
     -        S14*(10*S26*S34 + 2*S16*(S23 + S34) - 14*S24*S36 + 
     -           4*S26*S36 - 2*S34*S36 + 7*S23*S46 - S34*S46 + 
     -           14*S36*S46) + 
     -        S13*(4*S26*S34 + 2*S16*(S23 + S34) + 2*S24*S36 + 
     -           4*S26*S36 + S23*S46 + 14*S24*S46 + 2*S34*S46 + 
     -           2*S36*S46 - 16*S46**2 + S14*(2*S26 - 2*S36 + 3*S46)) + 
     -        S34*(2*S24*S36 - S23*S46 + 2*S24*S46 - S34*S46 + 
     -           2*S36*S46 + 
     -           2*S16*(S23 + 4*S24 - 4*S26 + S34 + 6*S46) + 
     -           2*S26*(S34 - 2*(S36 + S46)))))/
     -   (16.*S13*S24*S34*(S13 + S14 + S34))

      LEIGAXC= LEIGAXC +
     - (-8*S14**2*S16*S23 + 8*S14*S16**2*S23 + 4*S14*S16*S23**2 + 
     -     2*S14*S16*S23*S24 + 2*S14**2*S23*S26 - 10*S14*S16*S23*S34 + 
     -     12*S16**2*S23*S34 - 2*S14*S16*S24*S34 - 3*S16*S23*S24*S34 - 
     -     2*S14**2*S26*S34 + 7*S14*S23*S26*S34 - 8*S14*S26**2*S34 - 
     -     4*S16*S23*S34**2 - S16*S24*S34**2 + S14*S26*S34**2 + 
     -     4*S14**2*S23*S36 - 8*S14*S16*S23*S36 + 4*S14**2*S24*S36 - 
     -     4*S14*S16*S24*S36 + 4*S14**2*S26*S36 + 8*S14*S24*S26*S36 + 
     -     4*S14*S23*S34*S36 + 2*S16*S24*S34*S36 - 2*S14*S26*S34*S36 - 
     -     8*S14**2*S23*S46 + 8*S14*S16*S23*S46 - 4*S14*S23*S34*S46 + 
     -     4*S16*S23*S34*S46 - 8*S14*S24*S36*S46 + 8*S14*S23*S46**2 + 
     -     S13**2*(-(S16*(2*S23 + S24)) + 3*S14*S26 + 4*S23*S26 - 
     -        8*S26**2 + 2*S23*S46 + 2*S24*S46) - 
     -     2*S12**2*(S13*S36 + 3*S14*S36 + S34*(S36 + 2*S46)) + 
     -     S13*(4*S16**2*S23 - 2*S14**2*S26 - 8*S26**2*S34 - 
     -        2*S16*S24*(2*S24 - 8*S26 + S34 - 3*S36) + 8*S24*S26*S36 + 
     -        8*S24*S26*S46 + 2*S24*S34*S46 + 2*S23*S34*(2*S26 + S46) - 
     -        S16*S23*(5*S24 + 6*S34 + 4*S46) + 
     -        S14*(-4*S16*S23 + 11*S23*S26 + 2*S24*S26 - 20*S26**2 + 
     -           4*S26*S34 + 4*S23*S36 - 6*S26*S36 - 2*S23*S46 + 
     -           2*S24*S46)) + 
     -     S12*(4*S14**2*S36 + S13**2*(2*S26 - S46) + 
     -        S34*(2*S16*(S23 + 2*S24 - 8*S26 + S34 - 6*S36) - 
     -           2*S24*S36 + 2*S26*(S34 - 4*S46) + 5*S23*S46 - 
     -           S34*S46 - 6*S36*S46) + 
     -        2*S14*(4*S26*S34 + S16*(3*S23 + S34 - 4*S36) - 
     -           2*S23*S36 - 5*S24*S36 + 10*S26*S36 + S34*S36 + 
     -           4*S36**2 + 4*S23*S46 - 2*S36*S46) + 
     -        S13*(6*S14*S26 + 4*S26*S34 + 2*S16*(S23 + S34 - 2*S36) - 
     -           2*S24*S36 + 8*S26*S36 + 7*S23*S46 + 4*S24*S46 - 
     -           2*S34*S46 - 2*S36*S46 - 2*S14*(S36 + S46))))/
     -   (16.*S14*S23*S34*(S13 + S14 + S34))
     
      LEIGAXC= LEIGAXC - 
     -  (A6*(-(S13**3*(2*S23 + S24)) + 2*S12**2*(S13 - S34)*S34 + 
     -       S13**2*(-5*S14*S23 + 4*S16*S23 + 4*S23**2 + S23*S24 - 
     -          2*S24**2 - 8*S23*S26 + 4*S24*S26 - 4*S23*S34 + 2*S24*S36
     -          ) + S12*(6*S13*S14*S23 + 2*S14*S23*S34 + 
     -          S13**2*(2*S23 - 2*S24 + S34) - 
     -          S34**2*(S23 + 2*S24 - 4*S26 + S34 - 2*S36) + 
     -          S13*S34*(3*S23 + 4*S24 - 2*(2*S26 + S36))) + 
     -       S23*(-4*S16*S34**2 + 2*S14**2*(S23 - 2*(S34 + S36)) + 
     -          S14*S34*(4*S16 + 3*S23 - 8*S26 - 3*S34 + 2*S36 + 8*S46))
     -         + S13*(-4*S14**2*S23 + 
     -          S14*S23*(8*S16 + 11*S23 + 8*S24 - 20*S26 - 12*S34 - 
     -             2*S36 + 4*S46) + 
     -          S34*(8*S16*S23 + 4*S23**2 + 
     -             S24*(2*S24 - 4*S26 + S34 - 2*S36) + 
     -             S23*(3*S24 - 8*S26 - 2*S34 + 8*S46)))))/
     -   (16.*S14*S23*S34*(S13 + S14 + S34))

      LEIGAXC= LEIGAXC +
     -       (A6*(S12 + S23 + S24 - 2*S26)*
     -     (3*S13**2 + 3*S13*S34 + 4*S14*S34))/
     -   (8.*S14*S34*(S13 + S14 + S34)) + 
     -  (-4*S13**4*S26 + S13**3*
     -      (-12*S14*S26 + S16*(S23 + S24 + 6*S26) - 8*S26*S34 - 
     -        S12*S36 + 8*S26*S36 - S12*S46 + 8*S26*S46) + 
     -     S13*(-(S14**2*(S16*(2*S23 - 7*S24 + 6*S26) + 8*S26*S34 - 
     -             2*S12*S36 - 9*S24*S36 + 10*S26*S36 + 7*S12*S46 + 
     -             9*S23*S46 - 8*S26*S46)) + 
     -        S14*S34*(S16*(S23 + 9*S24 + 2*S26) - 6*S26*S34 - 
     -           S12*S36 + 8*S24*S36 + 4*S26*S36 - 9*S12*S46 - 
     -           8*S23*S46 + 20*S26*S46) + 
     -        S34**2*(S16*(S23 + S24 - 2*S26) - S12*(S36 + S46))) + 
     -     2*S13**2*(-4*S14**2*S26 + 
     -        S14*(4*S16*(S24 + S26) - 9*S26*S34 + 4*S24*S36 + 
     -           4*S26*S36 - 4*S12*S46 - 4*S23*S46 + 12*S26*S46) + 
     -        S34*(S16*(S23 + S24 + 2*S26) - 2*S26*S34 - 
     -           S12*(S36 + S46) + 4*S26*(S36 + S46))) + 
     -     S14*(-(S14**2*(S16*S23 - S12*S36 - S24*S36 + 
     -             2*S26*(S34 + S36) + S23*S46)) + 
     -        S34**2*(S16*(S23 + S24 - 2*S26) - S12*(S36 + S46)) + 
     -        S14*S34*(S16*(S24 + 2*S26) + S24*S36 - S12*S46 - 
     -           S23*S46 + S26*(-2*S34 + 2*S36 + 4*S46))))/
     -   (8.*S13*S14*S34*(S13 + S14 + S34)**2)  

      LEIGAXC= LEIGAXC +
     - (A6*(-(S13**3*(4*S23 + 5*S24)) + 
     -       2*S12**2*S34*(S13 + S14 + S34) + 
     -       S14*S23*(2*S14**2 - 
     -          S14*(4*S16 + 5*S23 + 4*S24 - 8*S26 + S34 - 2*S36) - 
     -          S34*(3*S23 + 2*S24 - 4*S26 + S34 - 2*S36)) + 
     -       S12*(-2*S14**2*(S23 + S34) + S13**2*(-2*S23 + 7*S34) + 
     -          S14*S34*(4*S16 + 3*S23 + 4*S24 - 8*S26 + S34 - 2*S36) + 
     -          S34**2*(3*S23 + 2*S24 - 4*S26 + S34 - 2*S36) - 
     -          S13*(S14*(6*S23 + 2*S24 - 3*S34) + 
     -             S34*(16*S16 + 5*S23 + 4*S24 - 12*S26 - 2*S34 + 2*S36)
     -             )) + S13**2*
     -        (2*S23**2 + 7*S23*S24 + 4*S24**2 + 8*S16*(S23 + S24) - 
     -          S14*(5*S23 + 9*S24) - 4*S23*S26 - 8*S24*S26 + 
     -          2*S23*S34 + 4*S26*S34 + 2*S24*S36 - 4*S23*S46) + 
     -       S13*(S14**2*(S23 - 4*S24) + 
     -          S34*(-2*S23**2 - 4*S16*(S23 + 2*S24) - 
     -             S24*(2*S24 - 4*S26 + S34 - 2*S36) + 
     -             S23*(-3*S24 + 4*S26 + 2*S34 - 4*S46)) + 
     -          S14*(-5*S23**2 + 2*S24**2 + 4*S16*(S23 + 2*S24) - 
     -             4*S24*S26 - S24*S34 + 8*S26*S34 + 2*S24*S36 - 
     -             S23*(S24 - 8*S26 - 2*S36 + 4*S46)))))/
     -   (16.*S13*S34*(S13 + S14 + S34)*(S23 + S24 + S34))
     
      LEIGAXC= LEIGAXC + 
     -  (-2*S14**2*S16*S23 + 4*S14*S16**2*S23 + 6*S14*S16*S23**2 + 
     -     S14*S16*S23*S24 + S14**2*S23*S26 + 4*S14*S16*S23*S34 - 
     -     4*S16**2*S23*S34 + 2*S16*S23**2*S34 - 3*S14*S16*S24*S34 + 
     -     12*S16**2*S24*S34 - 3*S16*S23*S24*S34 - S14**2*S26*S34 - 
     -     4*S14*S16*S26*S34 + 3*S14*S23*S26*S34 - 4*S16*S23*S26*S34 + 
     -     2*S14*S24*S26*S34 - 4*S14*S26**2*S34 + 2*S16*S23*S34**2 - 
     -     3*S16*S24*S34**2 + S14*S26*S34**2 - 4*S16*S26*S34**2 + 
     -     4*S14**2*S24*S36 - 14*S14*S16*S24*S36 - 2*S14*S23*S24*S36 - 
     -     2*S14*S24**2*S36 + 6*S14**2*S26*S36 + 12*S14*S23*S26*S36 + 
     -     4*S14*S24*S26*S36 + 2*S14*S24*S34*S36 + 10*S16*S24*S34*S36 + 
     -     2*S14*S26*S34*S36 - 4*S12**2*(S13 + S14 + S34)*S36 - 
     -     8*S14*S24*S36**2 + 16*S14*S16*S23*S46 + 2*S14*S23**2*S46 + 
     -     2*S14*S23*S24*S46 + 2*S14*S23*S34*S46 + 8*S14*S23*S36*S46 + 
     -     S13**2*(4*S16*S23 + 5*S16*S24 - S14*S26 + 4*S23*S26 - 
     -        2*S24*S26 - 4*S24*S46) - 
     -     S13*(8*S16**2*(S23 + S24) + S14**2*S26 - 2*S23*S26*S34 + 
     -        4*S24*S26*S34 - 4*S26**2*S34 + 2*S26*S34**2 + 
     -        6*S23*S24*S36 + 6*S24**2*S36 + 6*S24*S34*S36 + 
     -        S16*(-2*S23**2 + 8*S24**2 + 4*S26*S34 + 
     -           2*S24*(-6*S26 + 3*S34 + S36) + 
     -           S23*(11*S24 + 4*S26 - 6*S34 - 4*S46)) - 6*S23**2*S46 - 
     -        6*S23*S24*S46 + 12*S23*S26*S46 + 4*S24*S26*S46 - 
     -        6*S23*S34*S46 + 4*S24*S34*S46 - 8*S26*S34*S46 - 
     -        8*S24*S36*S46 + 8*S23*S46**2 + 
     -        S14*(-2*S16*S23 - 5*S16*S24 - 3*S23*S26 + 4*S24*S26 + 
     -           4*S26**2 + 4*S26*S34 - 4*S24*S36 - 2*S26*S36 + 
     -           4*S24*S46 + 4*S26*S46)) + 
     -     S12*(S13**2*(4*S26 + S46) + 
     -        S34*(4*S24*S36 - 4*S26*S36 + 2*S34*S36 + 
     -           2*S16*(2*S23 - 4*S26 + S34 + 2*S36 - 6*S46) + 
     -           4*S26*S46 + S34*S46 - 10*S36*S46 - S23*(2*S36 + S46))
     -         + S14*(2*S16*(2*S23 + S34 - 2*S36) + 12*S26*S36 - 
     -           2*S34*S36 + S34*S46 - 2*S36*S46 - S23*(6*S36 + S46)) + 
     -        S13*(4*S14*S26 + 4*S26*S34 - 2*S23*S36 + 4*S24*S36 + 
     -           2*S34*S36 + S14*S46 + 7*S23*S46 + 8*S24*S46 - 
     -           8*S26*S46 + 10*S34*S46 - 2*S36*S46 + 
     -           2*S16*(2*S23 + S34 + 4*(S36 + S46)))))/
     -   (16.*S13*S34*(S13 + S14 + S34)*(S23 + S24 + S34))

      LEIGAXC= LEIGAXC + 
     - (A6*(-(S14**3*S23) + 2*S14**2*S16*S23 + 3*S14**2*S23**2 + 
     -       S13**3*(2*S23 + S24) - S14**2*S23*S34 + 
     -       6*S14*S16*S23*S34 - 3*S14*S23**2*S34 - 2*S14*S16*S24*S34 - 
     -       3*S14*S23*S24*S34 - S14*S24**2*S34 + 2*S14**2*S26*S34 + 
     -       4*S14*S23*S26*S34 + 2*S14*S24*S26*S34 + 4*S12**2*S34**2 - 
     -       2*S14*S23*S34**2 + 4*S16*S23*S34**2 - S14*S24*S34**2 + 
     -       4*S16*S24*S34**2 - 2*S14*S26*S34**2 - 6*S14**2*S23*S36 + 
     -       2*S14*S23*S34*S36 + 
     -       S13**2*(-S23**2 - 3*S23*S24 + S24**2 - 4*S16*(S23 + S24) + 
     -          S14*(4*S23 + 3*S24) + 2*S23*S26 + 6*S24*S26 + 
     -          3*S23*S34 - 2*S24*S34 - 2*S26*S34 + 2*S24*S36 - 
     -          4*S23*S46 - 8*S24*S46) + 
     -       S12*(-(S14**2*S23) + S13**2*(S23 - S24 - 2*S34) + 
     -          S13*S14*(S23 - S34) - 
     -          S14*S34*(2*S16 + 9*S23 + S24 - 2*S26 + S34 - 6*S36) + 
     -          S34**2*(-2*S16 + 2*S23 + S24 - 6*S26 + S34 - 2*S36) + 
     -          S13*S34*(8*S16 + 4*S23 - 7*S24 - 8*S26 - 2*S36 + 8*S46))
     -         + S13*(S14**2*(S23 + 2*S24) + 
     -          S34*(-S23**2 + 10*S16*S24 - 2*S24**2 - 2*S26*S34 + 
     -             2*S24*(2*S26 - S34 + S36) + 
     -             S23*(-2*S24 + 2*S26 + S34 - 4*S46)) - 
     -          S14*(S23**2 + 9*S23*S24 + 2*S16*(S23 + 2*S24) + 
     -             2*(S24**2 - 2*S24*S26 + 9*S26*S34 - 5*S24*S36) - 
     -             2*S23*(2*S34 + S36 + 2*S46)))))/
     -   (8.*S14*S34*(S13 + S14 + S34)*(S23 + S24 + S34))

      LEIGAXC= LEIGAXC +
     -  (S14**2*(S16*S23 - 3*S23*S26 + 2*S23*S36 - 8*S26*S36) - 
     -     S13**2*(S16*(2*S23 + S24) + S14*S26 - 2*S23*S26 - 
     -        2*S24*S26 + 4*S26**2 - 2*S26*S34 - S12*S46 + 2*S23*S46) + 
     -     S14*(-2*S16**2*S23 + 4*S23*S26*S34 + 3*S24*S26*S34 - 
     -        6*S26**2*S34 + 3*S26*S34**2 - 2*S23*S24*S36 - 
     -        2*S24**2*S36 + 4*S24*S26*S36 + 2*S23*S34*S36 - 
     -        2*S24*S34*S36 - 2*S26*S34*S36 + 4*S24*S36**2 + 
     -        S16*(6*S26*S34 + 2*S12*S36 + 8*S24*S36 + 
     -           S23*(5*S24 + S34 - 14*S46)) + 2*S23**2*S46 + 
     -        2*S23*S24*S46 + 2*S23*S34*S46 - 4*S23*S36*S46 + 
     -        S12*(4*S26*S34 - 8*S24*S36 - 2*S26*S36 + 3*S23*S46 + 
     -           6*S36*S46)) + 
     -     S34*(-6*S16**2*S24 + 
     -        S12*(S24*S36 + 4*S26*S36 + S34*S36 + S23*(S36 - S46) - 
     -           4*S12*S46 + 6*S26*S46 + 2*S36*S46) - 
     -        S16*(S23**2 + 2*S26*S34 + S23*(-2*S26 + S34) + 
     -           2*S24*S36 - 2*S12*(2*S24 + 2*S26 + 3*S46))) + 
     -     S13*(4*S16**2*(S23 + S24) - S14**2*S26 + 2*S23*S26*S34 + 
     -        2*S24*S26*S34 - 4*S26**2*S34 + 2*S26*S34**2 + 
     -        S12*S23*S36 + S12*S24*S36 + 4*S12*S26*S36 + S12*S34*S36 - 
     -        2*S12*S23*S46 + 3*S12*S24*S46 + 4*S12*S26*S46 - 
     -        4*S24*S26*S46 - 2*S23*S34*S46 - 2*S12*S36*S46 - 
     -        4*S24*S36*S46 - 8*S12*S46**2 + 4*S23*S46**2 + 
     -        S14*(-(S16*(S23 + S24)) + 7*S24*S26 - 2*S26**2 + 
     -           2*S26*S34 + 2*S26*S36 + 2*S23*(2*S26 + S36 - S46) + 
     -           S12*S46 + 10*S26*S46) - 
     -        S16*(S23**2 + 3*S24**2 + 2*S26*S34 + 4*S12*S36 + 
     -           S24*(6*S26 + 6*S36 - 8*S46) + 4*S12*S46 - 
     -           S23*(S24 + 2*S26 - 3*S34 + 8*S46))))/
     -   (8.*S14*S34*(S13 + S14 + S34)*(S23 + S24 + S34))

      LEIGAXC= LEIGAXC +
     - (A6*(3*S13 + S14)*(S12 + S23 + S24 - 2*S26))/
     -   (8.*S34*(S13 + S14 + S34)) + 
     -  (-(S14**2*(7*S16*S23 + 8*S26*S34 - 7*S12*S36 - 7*S24*S36 + 
     -          14*S26*S36 + 7*S23*S46)) + 
     -     S13*S34*(7*S16*S24 - 8*S16*S26 + 4*S26*S34 + 7*S24*S36 - 
     -        8*S26*S36 - 7*S12*S46 - 7*S23*S46 + 6*S26*S46) + 
     -     S13**2*(7*S16*S24 + 4*S26*S34 + 7*S24*S36 - 7*S12*S46 - 
     -        7*S23*S46 + 14*S26*S46) + 
     -     S13*S14*(-7*S16*(S23 - S24) - 
     -        2*S26*(2*S34 + 7*S36 - 7*S46) + 
     -        7*(S12*S36 + 2*S24*S36 - S12*S46 - 2*S23*S46)) + 
     -     S14*S34*(S16*(-7*S23 + 16*S26) + 
     -        2*S26*(-4*S34 + S36 + 8*S46) + 
     -        7*(S12*S36 + S24*S36 - S23*S46)))/
     -   (8.*S34**2*(S13 + S14 + S34)**2)
     
      LEIGAXC= LEIGAXC +
     - (A6*((8*S12**2 + 2*S16*(S23 + 3*S24) + 
     -          S12*(-8*S16 + 3*S23 + S24 - 8*S26))*S34**2 - 
     -       S14*S34*(S23**2 + 2*S23*S24 + S24**2 + 
     -          2*S16*(-5*S23 + S24) - 2*S23*S26 - 2*S24*S26 - 
     -          2*S23*S34 + 6*S26*S34 + S12*(14*S23 + S34 - 16*S36)) + 
     -       S14**2*(8*S23**2 + 3*S23*S34 + S24*S34 - 16*S23*S36) + 
     -       S13**2*(8*S24**2 + 3*S23*S34 + S24*(S34 - 16*S46)) - 
     -       S13*S34*(6*S16*S23 + 3*S23**2 + 18*S12*S24 - 14*S16*S24 + 
     -          6*S23*S24 + 3*S24**2 - 6*S23*S26 - 6*S24*S26 + 
     -          3*S12*S34 + 2*S24*S34 + 2*S26*S34 - 16*S12*S46) + 
     -       2*S13*S14*(-16*S26*S34 + S24*(S34 + 8*S36) + 
     -          S23*(-8*S24 + 3*S34 + 8*S46))))/
     -   (16.*S34**2*(S13 + S14 + S34)*(S23 + S24 + S34)) + 
     -  (S14**2*(-(S23*S26) + 7*S23*S36 + 7*S24*S36 - 16*S26*S36) - 
     -     7*S13**2*(S23*S46 + S24*(S26 + S46)) - 
     -     S13*(11*S24*S26*S34 - 8*S26**2*S34 + 4*S26*S34**2 - 
     -        7*S12*S24*S36 + 7*S24**2*S36 + 7*S24*S34*S36 + 
     -        S16*S24*(15*S24 + 15*S34 + 14*S36 - 16*S46) + 
     -        S16*S23*(7*S24 - 4*S34 - 14*S46) - 7*S23**2*S46 - 
     -        15*S12*S24*S46 + 14*S24*S26*S46 - 7*S12*S34*S46 + 
     -        7*S24*S34*S46 - 10*S26*S34*S46 + 16*S12*S46**2 + 
     -        S23*(4*S26*S34 + 7*S24*S36 - 7*S24*S46 + 14*S26*S46)) + 
     -     S14*(8*S12*S26*S34 + 15*S23*S26*S34 + 8*S24*S26*S34 - 
     -        16*S26**2*S34 + 8*S26*S34**2 - 7*S12*S23*S36 - 
     -        16*S12*S24*S36 - 7*S23*S24*S36 - 7*S24**2*S36 + 
     -        14*S23*S26*S36 + 14*S24*S26*S36 - 7*S12*S34*S36 + 
     -        7*S23*S34*S36 - 10*S26*S34*S36 + 
     -        S16*(7*S23**2 - 8*S24*S34 + 16*S26*S34 + 2*S24*S36 + 
     -           S23*(15*S24 + 11*S34 - 18*S46)) + S12*S23*S46 + 
     -        7*S23**2*S46 + 7*S23*S24*S46 + 7*S23*S34*S46 + 
     -        16*S12*S36*S46 + 
     -        S13*(S24*(S26 + 7*S36 - 7*S46) + 
     -           7*S23*(S26 + S36 - S46) + 16*S26*S46)) + 
     -     S34*(-8*S16**2*S23 + 
     -        2*S16*(2*S23*S34 - 4*S24*S34 + 5*S24*S36 + 
     -           4*S12*(S24 + S36) - 5*S23*S46) + 
     -        S12*(7*S24*S36 - 8*S26*S36 - 8*S12*S46 - 7*S23*S46 + 
     -           16*S26*S46)))/
     -   (16.*S34**2*(S13 + S14 + S34)*(S23 + S24 + S34))                      

      END
C-----------------------------------------------------------------------
      FUNCTION LEID(P,Q,II,JJ,KKK,LL)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L,II,JJ,KKK,LL
      DOUBLE PRECISION LEID,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
C---CONVERT FROM ERT NOTATION (Q'QQ'BARQBAR) TO LEIDEN (QQBARQ'BARQ')
      I=JJ
      J=LL
      K=KKK
      L=II
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEID = 0
      LEID = LEID + P1K*P2K * (  - 2*P12*P25*P34 - 4*P13*P14*
     +     P15 - 2*P13*P15*P34 + 2*P13*P24*P25 - 2*P14*P15*P34 + 2*P14*
     +     P23*P25 )/P25
      LEID = LEID + P1K*P3K * (  - P12*P34 - 2*P14*P24 + P24*P25 )
      LEID = LEID + P1K*P4K * (  - P12*P34 - 2*P13*P23 + P23*P25 )
      LEID = LEID + P1K**2 * ( P23 + P24 )*P34
      LEID = LEID + P2K*P3K * (  - P12*P25*P34 - 4*P13*P14*P15
     +     - 2*P13*P15*P34 + 3*P14*P15*P25 - 2*P14*P15*P34 - 2*P14*P24*
     +     P25 )/P25
      LEID = LEID + P2K*P4K * (  - P12*P25*P34 - 4*P13*P14*P15
     +     + 3*P13*P15*P25 - 2*P13*P15*P34 - 2*P13*P23*P25 - 2*P14*P15*
     +     P34 )/P25
      LEID = LEID + P2K**2 * ( P13 + P14 )*P34
      LEID = LEID + P3K*P4K * (  - 2*P12*P34 + 2*P13*P24 + 2*
     +     P14*P23 )
      LEID = LEID + P3K**2 * (  - 2*P14*P24 )
      LEID = LEID + P4K**2 * (  - 2*P13*P23 )
c$$$      LEID = LEID + KK * ( Q2*P12*P25*P34 - Q2*P13*P24*P25
c$$$     +     - Q2*P14*P23*P25 + 2*P12*P13*P15*P25 - P12*P13*P24*P25 - 2
c$$$     +     *P12*P13**2*P15 + 2*P12*P14*P15*P25 - P12*P14*P23*P25 - 2*P12
c$$$     +     *P14**2*P15 + P12*P25*P34**2 + P12**2*P25*P34 + 2*P13*P14*P23
c$$$     +     *P25 + 2*P13*P14*P24*P25 + 2*P13*P15*P23*P25 + 2*P13*P23*P24*
c$$$     +     P25 - P13*P24*P25*P34 - 2*P13**2*P15*P23 - 2*P13**2*P15*P24
c$$$     +     + 2*P14*P15*P24*P25 + 2*P14*P23*P24*P25 - P14*P23*P25*P34 - 
c$$$     +     2*P14**2*P15*P23 - 2*P14**2*P15*P24 )/(4*P25)
      LEID = LEID/(P15*P25*P34**2)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIQAXD(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT E FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIQAXD,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)

      LEIQAXD = (S13*S14*(S17*S26 - S16*S27 - S27*S36 + S26*S37) + 
     -       S13**2*(S27*S46 - S26*S47) + 
     -       S13*S34*(S17*S26 - S16*S27 - 2*S27*S36 + 2*S26*S37 - 
     -          S27*S46 + S26*S47) + 
     -       S34*(-(S17*S26*S34) + S16*S27*S34 - S14*S27*(S36 + S46) + 
     -          S14*S26*(S37 + S47)))/(S14**2*(S13 + S14 + S34)**2) + 
     -    (S14*(-2*S23*S27*S36 + 
     -          S17*(S26*S34 + S23*(S26 - S36) - S24*S36) + 
     -          2*S23*S26*S37 + 
     -          S16*(-(S23*S27) - S27*S34 + S23*S37 + S24*S37) - 
     -          S13*S27*S46 - S23*S27*S46 + S12*S37*S46 + S23*S37*S46 + 
     -          S13*S26*S47 + S23*S26*S47 - S12*S36*S47 - S23*S36*S47)
     -        + S12*(2*S27*S34*S36 - 2*S26*S34*S37 - S13*S27*S46 + 
     -          S13*S37*S46 - S34*S37*S46 - S17*S34*(S26 + 2*S46) + 
     -          S13*S26*S47 - S13*S36*S47 + S34*S36*S47 + 
     -          S16*S34*(S27 + 2*S47)) + 
     -       S24*(S13*S27*(2*S36 + S46) + 
     -          S17*(S26*S34 + S13*S36 - S34*S36 + 2*S13*S46) - 
     -          S13*S26*(2*S37 + S47) - 
     -          S16*(S27*S34 + S13*S37 - S34*S37 + 2*S13*S47)))/
     -     (S14**2*(S12 + S14 + S24)*(S13 + S14 + S34)) 

      LEIQAXD = LEIQAXD + 
     -    ((-(S17*S26) + S16*S27)*S34**2 + S14**2*(S27*S36 - S26*S37) + 
     -       S14*S34*(S17*S26 - S16*S27 - S27*S36 + S26*S37 - 
     -          2*S27*S46 + 2*S26*S47) + 
     -       S13*(-(S27*S34*(S36 + S46)) + S26*S34*(S37 + S47) + 
     -          S14*(S17*S26 - S16*S27 - S27*S46 + S26*S47)))/
     -     (S13**2*(S13 + S14 + S34)**2) + 
     -    (S12*(-(S14*S27*S36) - S17*S34*(S26 + 2*S36) + S14*S26*S37 + 
     -          S16*S34*(S27 + 2*S37) + 2*S27*S34*S46 - S14*S37*S46 + 
     -          S34*S37*S46 - 2*S26*S34*S47 + S14*S36*S47 - S34*S36*S47)
     -         + S13*(-(S14*S27*S36) - S24*S27*S36 + S14*S26*S37 + 
     -          S24*S26*S37 - 2*S24*S27*S46 - S12*S37*S46 - 
     -          S24*S37*S46 + 
     -          S17*(S26*S34 + S24*(S26 - S46) - S23*S46) + 
     -          2*S24*S26*S47 + S12*S36*S47 + S24*S36*S47 + 
     -          S16*(-(S24*S27) - S27*S34 + S23*S47 + S24*S47)) + 
     -       S23*(S14*S27*(S36 + 2*S46) + 
     -          S17*(S26*S34 + 2*S14*S36 + S14*S46 - S34*S46) - 
     -          S14*S26*(S37 + 2*S47) - 
     -          S16*(S27*S34 + 2*S14*S37 + S14*S47 - S34*S47)))/
     -     (S13**2*(S12 + S13 + S23)*(S13 + S14 + S34)) 

      LEIQAXD = LEIQAXD + 
     -    (S12*S14*(-(S17*S36) - S27*S36 + (S16 + S26)*S37) + 
     -       S12*S24*(-(S17*S36) - 2*S27*S36 + S16*S37 + 2*S26*S37 + 
     -          S37*S46 - S36*S47) + S12**2*(-(S37*S46) + S36*S47) + 
     -       S24*(S17*S24*S36 - S16*S24*S37 + 
     -          S14*(-(S27*S36) + S26*S37 + S37*S46 - S36*S47)))/
     -     (S14**2*(S12 + S14 + S24)**2) - 
     -    (S12*S13*(S17*S46 + S27*S46 - (S16 + S26)*S47) + 
     -       S12*S23*(S17*S46 + 2*S27*S46 + S37*S46 - S16*S47 - 
     -          2*S26*S47 - S36*S47) + S12**2*(-(S37*S46) + S36*S47) + 
     -       S23*(-(S17*S23*S46) + S16*S23*S47 + 
     -          S13*(S27*S46 + S37*S46 - (S26 + S36)*S47)))/
     -     (S13**2*(S12 + S13 + S23)**2)

      LEIQAXD = LEIQAXD/2.D0

      END
C-----------------------------------------------------------------------
      FUNCTION LEIE(P,Q,II,JJ,KKK,LL)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L,II,JJ,KKK,LL
      DOUBLE PRECISION LEIE,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
C---CONVERT FROM ERT NOTATION (Q'QQ'BARQBAR) TO LEIDEN (QQBARQ'BARQ')
      I=JJ
      J=LL
      K=KKK
      L=II
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIE = 0
      LEIE = LEIE + P1K*P2K * (  - 2*Q2*P13*P15*P25 + Q2*
     +     P23*P25**2 + 2*P13*P14*P15*P25 - 4*P14*P15*P25*P45 - P14*P23*
     +     P25**2 + 2*P14**2*P15*P45 - 2*P15*P23*P25*P34 + 2*P15*P25**2*
     +     P45 - P23*P25**2*P45 )/(P13*P25*P45)
      LEIE = LEIE + P1K*P3K * (  - 2*Q2*P13*P15 + Q2*P23*
     +     P25 - 2*P12*P14*P15 - P14*P23*P25 - 2*P15*P23*P34 + 2*P15*P25
     +     *P45 - P23*P25*P45 )/(P13*P45)
      LEIE = LEIE + P1K*P4K * ( 2*Q2*P25 + 2*P13*P15 - 2*P14*
     +     P25 - 2*P15*P25 - 2*P15*P34 - 2*P25*P45 )*P23/(P13*P45)
      LEIE = LEIE + P1K**2 * (  - 2*P15*P23*P34 )/(P13*P45)
      LEIE = LEIE + P2K*P3K * (  - 2*P12*P15*P25 + 2*P13*P15*
     +     P25 + 2*P14*P15*P45 - 2*P23*P25**2 )*P14/(P13*P25*P45)
      LEIE = LEIE + P2K*P4K * ( Q2*P23*P25**2 + 2*P12*P15*P25
     +     *P34 - 2*P13*P15*P24*P25 + 2*P14*P15*P23*P25 - 2*P14*P15*P25*
     +     P45 - P14*P23*P25**2 + 2*P14**2*P15*P45 - P15*P23*P25**2 )
     +     /(P13*P25*P45)
      LEIE = LEIE + P2K**2 * ( 2*P13*P15 - P23*P25 )*P14/(P13*P45)
      LEIE = LEIE + P3K*P4K * (  - 2*Q2*P13*P15 + Q2*P23*
     +     P25 - 2*P12*P14*P15 + 2*P13*P15*P23 - P14*P23*P25 - 3*P15*P23
     +     *P25 + 2*P15*P25*P45 )/(P13*P45)
      LEIE = LEIE + P3K**2 * (  - 2*P12*P15 - P23*P25 )*P14/(P13*P45)
      LEIE = LEIE + P4K**2 * ( 2*P15*P23 )/P45
c$$$      LEIE = LEIE + KK * ( Q2*P13*P15*P24*P25 + Q2*P14*P23*
c$$$     +     P25**2 - Q2*P14**2*P15*P45 + Q2*P15*P23*P25**2 + Q2*P23
c$$$     +     *P25**2*P45 - Q2**2*P23*P25**2 - P12*P13*P15*P25*P34 + P12
c$$$     +     *P14*P15*P23*P25 - P12*P15*P24*P25*P34 - P12*P15*P25*P34**2 - 
c$$$     +     P12**2*P15*P25*P34 - P13*P14*P15*P23*P25 - P13*P14*P15*P24*
c$$$     +     P25 - P13*P15*P23*P24*P25 - P14*P15*P23*P24*P25 + P14*P15*P23
c$$$     +     *P25*P34 + P14*P15*P25*P35*P45 - P15*P23*P25**2*P45 )
c$$$     +     /(2*P13*P25*P45)
      LEIE = LEIE/(P15*P25*P34)
      END
C-----------------------------------------------------------------------
      FUNCTION BDPQEW(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE E FUNCTION POLARIZED CONTRIBUTIONS TO W-EXCHANGE
      INTEGER I,J,K,L
      DOUBLE PRECISION BDPQEW,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)

      BDPQEW=(-(S16*S23*S27*S34) + S14*S23*S27*S36 - S14*S16*S23*S37 -
     -    S14*S23*S26*S37 + S12*S16*S34*S37 + S14*S23*S27*S46 +
     -    S12*S27*S34*S46 - S12*S14*S37*S46 +
     -    S17*(S23*S26*S34 - S12*S34*S36 + S14*S23*(S36 + S46)) -
     -    S14*S16*S23*S47 - S14*S23*S26*S47 - S12*S26*S34*S47 +
     -    S12*S14*S36*S47 + S13*
     -     (S17*S26*S34 - S16*S27*S34 - S24*S27*S36 + S24*S26*S37 -
     -       S24*S27*S46 - S12*S37*S46 - S17*S24*(S36 + S46) +
     -       S24*S26*S47 + S12*S36*S47 + S16*S24*(S37 + S47)))/
     -  (2.*S13*S14*(S12 + S13 + S23)*(S13 + S14 + S34))


      BDPQEW=BDPQEW+(S34*(-(S17*S26*S34) + S16*S27*S34 -
     -      (S13 + S14)*(S27*(S36 + S46) - S26*(S37 + S47))))/
     -  (S13*S14*(S13 + S14 + S34)**2)


      BDPQEW=BDPQEW+
     -     (S12*(S14*S17*S36 + S17*S24*S36 + S14*S27*S36 + S24*S27*S36 -
     -      S16*S24*S37 - S24*S26*S37 - S14*(S16 + S26)*S37 +
     -      S17*(S13 + S23)*S46 + S13*S27*S46 + S23*S27*S46 -
     -      S13*S16*S47 - S16*S23*S47 - S13*S26*S47 - S23*S26*S47))/
     -  (2.*S13*S14*(S12 + S13 + S23)*(S12 + S14 + S24))

      BDPQEW=BDPQEW+
     -     (-(S16*S24*S27*S34) + S13*S24*S27*S36 + S12*S27*S34*S36 -
     -    S13*S16*S24*S37 - S13*S24*S26*S37 - S12*S26*S34*S37 +
     -    S13*S24*S27*S46 + S12*S13*S37*S46 +
     -    S17*(S24*S26*S34 - S12*S34*S46 + S13*S24*(S36 + S46)) -
     -    S13*S16*S24*S47 - S13*S24*S26*S47 + S12*S16*S34*S47 -
     -    S12*S13*S36*S47 + S14*
     -     (S17*S26*S34 - S16*S27*S34 - S23*S27*S36 + S23*S26*S37 -
     -       S23*S27*S46 + S12*S37*S46 - S17*S23*(S36 + S46) +
     -       S23*S26*S47 - S12*S36*S47 + S16*S23*(S37 + S47)))/
     -  (2.*S13*S14*(S12 + S14 + S24)*(S13 + S14 + S34))

      END
C-----------------------------------------------------------------------
      FUNCTION ERTEW(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE E FUNCTION UNPOLARIZED CONTRIBUTIONS TO W-EXCHANGE
C   (ERT+LEI).
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTEW,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S,
     $     S16,S17,S26,S27,S36,S37,
     $     S46,S47,S67
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S16=2*DOT(P,ABS(I),6)*SIGN(1,I)
      S17=2*DOT(P,ABS(I),7)*SIGN(1,I)
      S26=2*DOT(P,ABS(J),6)*SIGN(1,J)
      S27=2*DOT(P,ABS(J),7)*SIGN(1,J)
      S36=2*DOT(P,ABS(K),6)*SIGN(1,K)
      S37=2*DOT(P,ABS(K),7)*SIGN(1,K)
      S46=2*DOT(P,ABS(L),6)*SIGN(1,L)
      S47=2*DOT(P,ABS(L),7)*SIGN(1,L)
      S67=2*DOT(P,6,7)

      ERTEW=(-(S17*S23*S26*S34) + S14*S17*S23*S36 + S14*S23*S27*S36 +
     -    S12*S17*S34*S36 + S14*S23*S26*S37 - 2*S12*S14*S36*S37 +
     -    S14*S17*S23*S46 + S14*S23*S27*S46 + S12*S27*S34*S46 -
     -    S12*S14*S37*S46 + S14*S23*S26*S47 + S12*S26*S34*S47 -
     -    S12*S14*S36*S47 + S16*
     -     (-2*S17*S23*S34 - S23*S27*S34 + S14*S23*S37 + S12*S34*S37 +
     -       S14*S23*S47 + S13*(S27*S34 - S24*S37 - S24*S47)) +
     -    S13*(S17*S26*S34 + 2*S26*S27*S34 - S24*S27*S36 -
     -       S24*S27*S46 + S12*S37*S46 - S17*S24*(S36 + S46) +
     -       S12*S36*S47 + 2*S12*S46*S47 - S24*S26*(S37 + S47) -
     -       2*S12*S34*S67))/
     -  (2.*S13*S14*(S12 + S13 + S23)*(S13 + S14 + S34))

      ERTEW=ERTEW+(S34*(S17*S26*S34 + S16*S27*S34 -
     -      (S13 + S14)*(S27*(S36 + S46) + S26*(S37 + S47))))/
     -  (S13*S14*(S13 + S14 + S34)**2)


      ERTEW=ERTEW+
     -     (S12*(-2*S26*S27*S34 + S14*S27*S36 + S24*S27*S36 +
     -      S14*S26*S37 + S24*S26*S37 + S13*S27*S46 + S23*S27*S46 -
     -      2*S12*S37*S46 + S17*
     -       (-2*S26*S34 + S14*S36 + S24*S36 + S13*S46 + S23*S46) +
     -      S13*S26*S47 + S23*S26*S47 - 2*S12*S36*S47 +
     -      S16*(-2*S17*S34 - 2*S27*S34 + S14*S37 + S24*S37 + S13*S47 +
     -         S23*S47) + 2*S12*S34*S67))/
     -  (2.*S13*S14*(S12 + S13 + S23)*(S12 + S14 + S24))

      ERTEW=ERTEW+
     -     (-(S17*S24*S26*S34) + S13*S17*S24*S36 + S13*S24*S27*S36 +
     -    S12*S27*S34*S36 + S13*S24*S26*S37 + S12*S26*S34*S37 +
     -    S13*S17*S24*S46 + S13*S24*S27*S46 + S12*S17*S34*S46 -
     -    S12*S13*S37*S46 + S13*S24*S26*S47 - S12*S13*S36*S47 -
     -    2*S12*S13*S46*S47 +
     -    S16*(-2*S17*S24*S34 - S24*S27*S34 + S13*S24*S37 +
     -       S13*S24*S47 + S12*S34*S47 +
     -       S14*(S27*S34 - S23*S37 - S23*S47)) +
     -    S14*(S17*S26*S34 + 2*S26*S27*S34 - S23*S27*S36 +
     -       2*S12*S36*S37 - S23*S27*S46 + S12*S37*S46 -
     -       S17*S23*(S36 + S46) + S12*S36*S47 - S23*S26*(S37 + S47) -
     -       2*S12*S34*S67))/
     -  (2.*S13*S14*(S12 + S14 + S24)*(S13 + S14 + S34))

      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION DILOG(X)
      IMPLICIT NONE
c$$$C---RETURN THE DILOGARITHM (MAXIMUM ERROR = 0.0058, MAX FRAC ERROR = 1%)
c$$$      DOUBLE PRECISION DILOG,X,PISQO6,OX,XX,L2
c$$$      DATA PISQO6,L2/2*0/
c$$$      IF (PISQO6.EQ.0) PISQO6=(ATAN(1D0)*4)**2/6
c$$$      IF (L2.EQ.0) L2=LOG(2D0)
c$$$      IF (X.LT.-1) THEN
c$$$        XX=1/X
c$$$      ELSE
c$$$        XX=X
c$$$      ENDIF
c$$$      IF (XX.LT.-0.5) THEN
c$$$        OX=1+XX
c$$$        DILOG=-PISQO6/2-L2*LOG(-XX)
c$$$     $       -OX**2/4-5*OX**3/24-OX**4/6-131*OX**5/960
c$$$      ELSEIF (XX.LT.0.5) THEN
c$$$        DILOG=XX+XX**2/4+XX**3/9
c$$$      ELSEIF (XX.LT.1) THEN
c$$$        OX=1-XX
c$$$        DILOG=PISQO6-LOG(OX)*LOG(XX)-OX-OX**2/4-OX**3/9
c$$$      ELSEIF (XX.EQ.1) THEN
c$$$        DILOG=PISQO6
c$$$      ELSE
c$$$        WRITE (*,*) 'DILOG CALLED FOR X=',X
c$$$        DILOG=0
c$$$      ENDIF
c$$$      IF (X.LT.-1) DILOG=-DILOG-PISQO6-LOG(-X)**2/2
      double precision dilog,rsp,x,pisqo6
      data pisqo6/0/
      if (pisqo6.eq.0) pisqo6=(atan(1d0)*4)**2/6
      if (x.lt.1) then
        dilog=rsp(x)
      elseif (x.gt.1) then
c---if dilog is complex, return its real part
        dilog=pisqo6-log(x)*log(x-1)-rsp(1-x)
      else
        dilog=pisqo6
      endif
      END
      
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c
c spence function
c
      block data splint
      implicit double precision (a-h,o-z)
      common/spint/a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,zeta2
      data a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,zeta2/
     1 -0.250000000000000D0,
     2 -0.111111111111111D0,
     3 -0.010000000000000D0,
     4 -0.017006802721088D0,
     5 -0.019444444444444D0,
     6 -0.020661157024793D0,
     7 -0.021417300648069D0,
     8 -0.021948866377231D0,
     9 -0.022349233811171D0,
     1 -0.022663689135191D0,
     2  1.644934066848226D0/
      end
c
c spence function taking only real arguments
c
      function rsp(x)
      implicit double precision(a-h,o-z)
      common/spint/a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,zeta2
      x2=x*x
      if(x.gt.1.D0)then
        write(*,*)' argument greater than 1 passed to spence function'
        rsp=0.D0
        return
      endif
      if(x2.gt.1.D0.and.x.gt.0.5D0)then
        y=(x-1.D0)/x
        z=-log(1.D0-y)
        z2=z*z
        rsp=z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
     3 +zeta2-log(x)*log(1.D0-x)+0.5D0*log(x)**2
        return
      elseif(x2.gt.1.D0.and.x.le.0.5D0)then
        y=1.D0/x
        z=-log(1.D0-y)
        z2=z*z
        rsp=-z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
     3 -zeta2-0.5D0*log(-x)**2
        return
      elseif(x2.eq.1.D0)then
        rsp=zeta2
        return
      elseif(x2.le.1.D0.and.x.gt.0.5D0)then
        y=1.D0-x
        z=-log(1.D0-y)
        z2=z*z
        rsp=-z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
     3 +zeta2+z*log(1.D0-x)
       return
      elseif(x2.le.1.D0.and.x.le.0.5D0)then
        y=x
        z=-log(1.D0-y)
        z2=z*z
        rsp=z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
        return
      else
        write(*,*)' illegal x value in spence function'
        rsp=0.D0
      endif
      return
      end
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE GFLIN1(ID,X,W)
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NSIZE=200000,NMAX=2000)
      COMMON /GBOOK/ A(NSIZE)
      IF (ID.GT.NMAX) RETURN
      IS=INT(A(ID+2)+0.5)
      A(IS+9)=A(IS+9)+1.
      IOX=2
      IF(X.LT.A(IS+2)) IOX=1
      IF(X.GE.A(IS+3)) IOX=3
      A(IS+12+IOX)=A(IS+12+IOX)+W
      IF(IOX.NE.2) RETURN
      IX=INT((X-A(IS+2))/A(IS+4)+0.5)
      DX=(X-A(IS+2))/A(IS+4)+0.5-IX
      IF (IX.EQ.0) THEN
        A(IS+19+IX)=A(IS+19+IX)+W
      ELSEIF (IX.EQ.A(IS+1)) THEN
        A(IS+18+IX)=A(IS+18+IX)+W
      ELSE
        A(IS+18+IX)=A(IS+18+IX)+(1-DX)*W
        A(IS+19+IX)=A(IS+19+IX)+DX*W
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------      
      
      
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RANGEN(N,R)
      IMPLICIT NONE
C---RANDOM NUMBER GENERATOR
C   USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C   RETURNS A VECTOR OF N RANDOM VALUES
C   IF (N.EQ.0) THE FIRST TWO VALUES IN R SET THE SEEDS
C   IF (N.LT.0) PRINT THE CURRENT VALUES OF THE SEEDS
      DOUBLE PRECISION R(*)
      INTEGER N,ISEED(3)
      if (N.GT.0) then
         call RM48(R,N)
      else if (N.LT.0) then
         call RM48UT(iseed(1),iseed(2),iseed(3))
         WRITE (6,'(I10,A,I10,I11,I11)') -N-1,', ISEED=',
     $        iseed(1),iseed(2),iseed(3)
      else ! N=0
        IF(NINT(R(1)) .eq. 0) then
           !-- retrieve seed --
           call RM48UT(iseed(1),iseed(2),iseed(3))
           R(1) = ISEED(1)
           R(2) = ISEED(2)
           R(3) = ISEED(3)
        else
           !-- set seed -------
           ISEED(1)=NINT(R(1))
           ISEED(2)=NINT(R(2))
           ISEED(3)=NINT(R(3))
           call RM48IN(iseed(1),iseed(2),iseed(3))
        end if
      end if
      !DATA ISEED/12345,678900/
      !DATA ISEED/1277158507, 1826842337/
      !DATA ISEED/1542788427,  474274578/
c$$$      IF (N.LT.0) WRITE (*,'(I10,A,I10,I11)') -N-1,', ISEED=',ISEED
c$$$      IF (N.GT.0) THEN
c$$$        DO I=1,N
c$$$          K=ISEED(1)/53668
c$$$          ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
c$$$          IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
c$$$          K=ISEED(2)/52774
c$$$          ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
c$$$          IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
c$$$          IZ=ISEED(1)-ISEED(2)
c$$$          IF (IZ.LT.1) IZ=IZ+2147483562
c$$$          R(I)=DBLE(IZ)*4.656613D-10
c$$$        ENDDO
c$$$      ELSEIF (N.EQ.0) THEN
c$$$        IF(NINT(R(1)) .eq. 0) then
c$$$           !-- retrieve seed --
c$$$           R(1) = ISEED(1)
c$$$           R(2) = ISEED(2)
c$$$        else
c$$$           !-- set seed -------
c$$$           ISEED(1)=NINT(R(1))
c$$$           ISEED(2)=NINT(R(2))
c$$$        end if
c$$$      ENDIF
      END


      SUBROUTINE RM48(RVEC,LENV)
C     Double-precision version of
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        based on RANMAR, modified by F. James, to generate vectors
C        of pseudorandom numbers RVEC of length LENV, where the numbers
C        in RVEC are numbers with at least 48-bit mantissas.
C   Input and output entry points: RM48IN, RM48UT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RM48:                                    ++
C!!!      CALL RM48 (RVEC, LEN)     returns a vector RVEC of LEN     ++
C!!!                   64-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RM48IN(I1,N1,N2)   initializes the generator from one ++
C!!!                   64-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++ 
C!!!                    output by RM48UT)                            ++ 
C!!!      CALL RM48UT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++  
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C for 32-bit machines, use IMPLICIT DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RVEC(*)
      COMMON/R48ST1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL,TWOM49, ONE, ZERO
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RM48 without RM48IN.
      IJKL = 54217137
c      IJKL = 12345678
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RM48IN(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RM48, may be called before
C         generating pseudorandom numbers with RM48.   The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      WRITE(6,'(A,I10,2X,2I10)') ' RM48 INITIALIZED:',IJKL,NTOT,NTOT2
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      ONE = 1.
      HALF = 0.5
      ZERO = 0.
      DO 2 II= 1, 97
      S = 0.
      T = HALF
      DO 3 JJ= 1, 48
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = HALF*T
    2 U(II) = S
      TWOM49 = T
      TWOM24 = ONE
      DO 4 I24= 1, 24
    4 TWOM24 = HALF*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
      WRITE(6,'(A,I15)') ' RM48IN SKIPPING OVER ',NOW
          DO 40 IDUM = 1, NTOT
          UNI = U(I97)-U(J97)
          IF (UNI .LT. ZERO)  UNI=UNI+ONE
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. ZERO)  C=C+CM
   40     CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. ZERO)  UNI=UNI+ONE
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. ZERO)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. ZERO) UNI=UNI+ONE
      RVEC(IVEC) = UNI
C             Replace exact zeros by 2**-49
         IF (UNI .EQ. ZERO)  THEN
            RVEC(IVEC) = TWOM49
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY RM48UT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RANGEN2(N,R)
      IMPLICIT NONE
C---RANDOM NUMBER GENERATOR
C   USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C   RETURNS A VECTOR OF N RANDOM VALUES
C   IF (N.EQ.0) THE FIRST TWO VALUES IN R SET THE SEEDS
C   IF (N.LT.0) PRINT THE CURRENT VALUES OF THE SEEDS
      DOUBLE PRECISION R(*)
      INTEGER N,I,ISEED(2),K,IZ
      DATA ISEED/12345,678900/
      IF (N.LT.0) WRITE (*,'(I10,A,I10,I11)') -N-1,', ISEED=',ISEED
      IF (N.GT.0) THEN
        DO I=1,N
          K=ISEED(1)/53668
          ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
          IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
          K=ISEED(2)/52774
          ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
          IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
          IZ=ISEED(1)-ISEED(2)
          IF (IZ.LT.1) IZ=IZ+2147483562
          R(I)=DBLE(IZ)*4.656613D-10
        ENDDO
      ELSEIF (N.EQ.0) THEN
        ISEED(1)=NINT(R(1))
        ISEED(2)=NINT(R(2))
      ENDIF
      END       

CC-----------------------------------------------------------------------
      SUBROUTINE BORNPROJ(S,P,WWW,*)
      IMPLICIT NONE
C---GENERATE A TWO-PARTON CONFIGURATION
C   WITHIN THE PHASE-SPACE LIMITS GIVEN IN CUTS()
C   THE WEIGHT GIVES THE TOTAL VOLUME OF PHASE-SPACE
      DOUBLE PRECISION S,P(4,7),W,R(1),
     $     Q2,QJAC,Y,YJAC,COMBPDFVEC,COMBPDFAX
      INTEGER SCHEME,NF,I,J,NEV,N,LPOL,VOGT,ICH,INNLO,IFL,EWFLAG
      INTEGER IND(-6:6)
      DOUBLE PRECISION X,Z,ZJAC,U,NRM
      DOUBLE PRECISION F(-6:6),WWW(0:5)
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGMU,CLqns1,CLg1,C2qns_reg1,
     $     C2qns_sing1,C2qns_delta1,C2g1,
     $     F1(-6:6),SIGMA,DELTA,SIGMA1,DELTA1,GL,GL1, 
     $     F2NLO,FLNLO,F2LO,FLLO,
     $     F2NLOq,FLNLOq,F2NLOg,FLNLOg, ZMAX, SUMC,
     $     CLqs1,C2qs_reg1,C2qs_sing1,C2qs_delta1
      DOUBLE PRECISION F3LO,CF3,CF3min,CF31,C3qns_reg1,C3qns_sing1,
     $     CF2min,C3qns_delta1,F3NLOq,SUMCEW,DOT,C4qns_reg1,C4qns_sing1,
     $     C4qns_delta1,C4qs_reg1,C4qs_sing1,C4qs_delta1,GLNLOq,G4NLOq,
     $     G4LO,GLLO,SUMCEWAX,GLNLOg
      DOUBLE PRECISION SCALEREN,SCALEFACT,QQ2
      DOUBLE PRECISION MUR,MUF,SCREN,SCFACT,ALPHAS,ALPHASPDF,A2PI,
     $     NLOPREF
      DOUBLE PRECISION G2LO,G2NLOg,G2NLOQ,G2NLO,Ee,Ep,ALPHAE

      DOUBLE PRECISION LOGZ,DIMZ,TRIMZ,SMZ,LOG1MZ,DI1MZ,DIZ,TRIZ,
     $       TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,ZETA2,ZETA3   
      DOUBLE PRECISION C2QNSMIN_REG,DLOG,C2QPS_REG,C2G_REG
      DOUBLE PRECISION C2QNSPLUS_SING,C2QNSPLUS_REG,C2QNSPLUS_DELTA
      DOUBLE PRECISION  L3,L2,L1,L0,L3X,L2X,L1X,L0X   
      DOUBLE PRECISION C2QNS_SING2,C2QNS_DELTA2,C2QNS_REG2P,C2QNS_REG2M,
     $        C2QS_SING2,C2QS_DELTA2,C2QS_REG2,C2G_REG2,
     $     G2NNLOq,G2NNLOg,G2NNLO,C2QNS_SING, C2QNS_REG,C2QNS_DELTA, 
     $     C2G_DELTA2,
     $     F2NNLOq,F2NNLOg,F2NNLO,CLQNSP_delta2,CLQNSM_delta2,CKM2(3,3),
     $     C3qns_sing2,C3qns_delta2,C3qns_reg2p,C3qns_reg2m,F3NNLOq,
     $     G4NNLOq,GLNNLOq,
     $     C4QNS_SING2,C4QNS_DELTA2,C4QNS_REG2P,C4QNS_REG2M,C4QPS_REG,
     $     C4QS_SING2,C4QS_DELTA2,C4QS_REG2
      DOUBLE PRECISION FLNNLOq,FLNNLOG,FLNNLO,CLqns2,CLqs2,CLG2     
      COMPLEX*16 WGPLG
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,GZ,MZ,LEPCHARGE,
     #                 GAMMAZ,PROPZ,INTGZ,SIGMAW,DELTAW,SIGMAW1,DELTAW1
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG 

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/NEV/NEV
      COMMON/JETSCALES/SCREN,SCFACT,QQ2
      COMMON/FLAGPOL/LPOL
      COMMON/EeEp/Ee,Ep
      COMMON/PBORNPAR/X,Y,Q2,YJAC,QJAC
      COMMON/CHANNEL/ICH ,IFL    
      COMMON/INNLO/INNLO
      COMMON/ALPHAQED/ALPHAE
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

C     ========================== START UP ==============================

C---  MASSIVE BOSON PROPAGATOR AND WIDTH
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

C---DEFINE RIEMANN ZETAS
      ZETA2=PI**2/6.0D0
      ZETA3=1.2020569031D0

      if ((1-X).lt.CUTOFF) then
        Return 1
      endif    
C  
      CALL RANGEN(1,R)
C---BRING JACOBIAN WEIGHT FROM LO
      W=QJAC*YJAC    
C
C--- CALL JET FUNCTION TO GET SCALES           
      CALL JETALG(P)
C SCALES DEFINED IN JETALG AND ALPHAS
      MUR=SQRT(SCREN)
      MUF=SQRT(SCFACT)
c      ALPHAS= alphasPDF(MUR)
c      A2PI=ALPHAS/2/PI

C TOGGLE VOGT CALCULATION
      VOGT=0
C X-SECTION IN PB Gev to nb to pb
      NRM=389385*1000
C CALL LHAPDF AND COMPUTE PDFS
c      CALL EVOLVEPDF(X,MUF,F1)

C GET PDF AND ALPHAS VALUES ALREADY CALCULATED FOR MATTWO
      CALL GETPDF(F1,ALPHAS)
      A2PI=ALPHAS/2/PI
C      IF(ABS(Q2-100.D0).LT.5.D0) THEN
C        WRITE(6,*) 'POLDIS PDF at x=',X,' mu2=',MUF**2,' Q2=',Q2
C        WRITE(6,*) '  xf(2)=',F1(2),' xf(1)=',F1(1)
C        WRITE(6,*) '  xf(-2)=',F1(-2),' xf(-1)=',F1(-1)
C        WRITE(6,*) '  xf(3)=',F1(3),' xf(-3)=',F1(-3)
C        WRITE(6,*) '  xf(4)=',F1(4),' xf(-4)=',F1(-4)
C        WRITE(6,*) '  xf(5)=',F1(5),' xf(-5)=',F1(-5)
C      ENDIF



C     ======================= PDF COMBINATIONS =========================

C     CREATE NEW MOMENTUM FRACTION FOR HIGHER ORDERS AND BUILD ITS PDF
      N=1                
      ZMAX=1-CUTOFF
      Z=X*(ZMAX/X)**(R(1)**N)
      ZJAC=LOG(ZMAX/X)*N*(R(1))**(N-1)*Z
                     
      U=X/Z

      CALL EVOLVEPDF(U,MUF,F) 
C SELECT INITIAL STATE CHANNELS      
      IF (ICH.EQ.0) THEN
      ELSEIF (ICH.EQ.1) THEN ! GLUONS TO 0, ONLY QUARKS
      F(0)=0D0
            IF (IFL.NE.0) THEN
            DO I=1,NF
            IF(I.NE.IFL) THEN
               F(I)=0D0
               F(-I)=0D0
           ENDIF
           ENDDO
           ENDIF
      ELSEIF (ICH.EQ.2) THEN ! QUARKS TO 0, ONLY GLUONS
           DO I=1,NF
               F(I)=0D0
               F(-I)=0D0
c               F(0)=ABS(F(0))  !debug Gluon Sign
           ENDDO
      ELSE
         WRITE(6,*)"WRONG CHANNEL= ",ICH
         STOP
      ENDIF         

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
        SUMCEW = SUMCEW + PROPZ*(CVZE**2+CAZE**2)*(CVZ(I)**2+CAZ(I)**2)
     #             +2*EQ(I)*PROPZ*INTGZ*(CVZ(I)*CVZE)
        SUMCEWAX = SUMCEWAX + PROPZ*(4*CVZ(I)*CAZ(I)*CVZE*CAZE)
     #             +2*EQ(I)*PROPZ*INTGZ*(CAZ(I)*CAZE)  
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
C     INDEXES OF FLAVOR TO USE (this was a mistake)
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
      DELTA = DELTA     + COMBPDFVEC(P,F,0)/U - SIGMA
      DELTA1 = DELTA1   + COMBPDFVEC(P,F1,0)/X - SIGMA1
      DELTAW = DELTAW   + COMBPDFAX(P,F,0)/U - SIGMAW
      DELTAW1 = DELTAW1 + COMBPDFAX(P,F1,0)/X - SIGMAW1

C     COMBINATIONS FOR MINUS COEFFICIENTS (ONLY RELEVANT FOR W)      
      CF2min = COMBPDFVEC(P,F,1)/U
      CF3min = COMBPDFAX(P,F,1)/U

C     POR FAVOR CAMBIAR A FUTURO      
      CF3 = DELTAW + SIGMAW
      CF31 = DELTAW1 + SIGMAW1

C     -- GLUON -- (SUMCEWAX DOESNT CONTRIBUTE HERE!)
      GL1 = GL1+ F1(0)/X * (SUMC+SUMCEW)/NF
      GL  = GL +  F(0)/U * (SUMC+SUMCEW)/NF 


C     ====================== STRUCTURE FNCTIONS ========================

C COMPUTE STRUCTURE FUNCTIONS
      IF (LPOL.EQ.0) THEN
        F2LO=(DELTA1+SIGMA1)*X
        FLLO=0D0
        F3LO=(DELTAW1+SIGMAW1)*X

        WWW(0)=W*2*PI/Y/Q2**2*ALPHAE**2 *NRM/NEV*   
     $                   ( (1+(1-Y)**2)*F2LO - Y**2*FLLO
     $                     + (1-(1-Y)**2)*F3LO )
C      IF (ABS(Q2-100D0).LT.5D0) THEN
C        WRITE(6,*) 'DIAG: Q2=',Q2,' Y=',Y,' X=',X,
C     $    ' F2LO=',F2LO,' dSig=',
C     $    2*PI/Y/Q2**2*ALPHAE**2*NRM*
C     $    ((1+(1-Y)**2)*F2LO - Y**2*FLLO
C     $        +(1-(1-Y)**2)*F3LO), ' pb/GeV2'
C        WRITE(6,*) 'PDF: xfu=',F1(2),' xfd=',F1(1),
C     $  ' at x=',X,' mu2=',MUF**2
C      ENDIF

        
C -- POLARIZED CASE
      ELSEIF (LPOL.Eq.1) THEN !wip
        G2LO=(DELTA1+SIGMA1)*X ! G2 MEANS 2*X*G1
        GLLO=0D0
        G4LO=(DELTAW1+SIGMAW1)*X

        WWW(0)=W*2*PI/Y/Q2**2*ALPHAE**2 *NRM/NEV*   
     $                           ( (1-(1-Y)**2)*G2LO 
     $                            - Y**2*GLLO + (1+(1-Y)**2)*G4LO )              
      ENDIF

C START HIGHER ORDER PART (NLO-NNLO)

C     BUILD THE COEFFICIENTS           
c     (in units of alphas/4/Pi)
      LOGMU=Log(Q2/MUF**2)  

C     ADD WEIGHT FROM Z FRACTION
      W=W*ZJAC 


C ----- UNPOLARIZED INCLUSIVE CROSS SECTION -----
       
      IF (LPOL.EQ.0) THEN
C LONGITUDINAL       
       CLqns1=CF*4*z
       CLqs1=CLqns1
       CLg1=NF*TR*16*z*(1-z)
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
     #                -2*(1+z)*LOGMU) !debug(chequear termino LOGMU)        
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
       FLNLOq=X*DELTA*CLqns1/Z*A2PI/2 *W !ONLY REGULAR PART NS 
     $        + X*SIGMA*CLqs1/Z*A2PI/2 *W !ONLY REGULAR PART S
c     
       F2NLOq=X*A2PI/2*  W*  (C2qns_reg1*SIGMA/Z  !REGULAR PART NS
     $   +    C2qns_sing1*(SIGMA/Z-SIGMA1)! PLUS DISTRIBUTION NS
     $   + SIGMA1* C2qns_delta1/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*A2PI/2*  W*  (C2qs_reg1*DELTA/Z  !REGULAR PART  S
     $   +    C2qs_sing1*(DELTA/Z-DELTA1)! PLUS DISTRIBUTION  S
     $   + DELTA1* C2qs_delta1/ZJAC )   ! DELTA TERM   S

       F3NLOq=X*A2PI/2*  W*  (C3qns_reg1*CF3/Z  !REGULAR PART NS
     $   +    C3qns_sing1*(CF3/Z-CF31)! PLUS DISTRIBUTION NS
     $   + CF31* C3qns_delta1/ZJAC )   ! DELTA TERM  NS 


C GLUON
       F2NLOg=X*(GL)*C2g1/Z*A2PI/2 *W !ONLY REGULAR PART  G
       FLNLOg=X*(GL)*CLg1/Z*A2PI/2 *W !ONLY REGULAR PART  G
C             
       FLNLO=FLNLOq+FLNLOg
       F2NLO=F2NLOq+F2NLOg      
              
       WWW(1)=WWW(0)+ 2*PI/Y/Q2**2*ALPHAE**2 *NRM/NEV*   
     $     ( (1+(1-Y)**2)*F2NLO -Y**2*FLNLO+ (1-(1-Y)**2)*F3NLOq)    


C -- SECOND ORDER IN ALPHAS --

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
C
       F2NNLOq=X*(A2PI/2)**2 * (C2qns_reg2p*DELTA/Z  !REGULAR PART NS
     $   + C2qns_reg2m*CF2min/Z  !REGULAR PART NS   
     $   + C2qns_sing2*(DELTA/Z-DELTA1)! PLUS DISTRIBUTION NS
     $   + DELTA1* C2qns_delta2/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*(A2PI/2)**2 *  (C2qs_reg2*SIGMA/Z  !REGULAR PART  S
     $   +    C2qs_sing2*(SIGMA/Z-SIGMA1)! PLUS DISTRIBUTION  S
     $   + SIGMA1* C2qs_delta2/ZJAC )   ! DELTA TERM   S

       F2NNLOg=X*((GL)*C2G_REG2/Z  !ONLY REGULAR PART  G
     $   + GL1* C2G_delta2/ZJAC ) *(A2PI/2)**2  ! DELTA TERM  GLU VOGT

       F2NNLO=(F2NNLOq+F2NNLOg)*W


C LONGITUDINAL STRUCTURE FUNCTION    
      CALL CFL(Z,VOGT,CLQS2,CLQNS2,CLQNSP_DELTA2,CLG2)

C QUARK
         FLNNLOq=X*(DELTA*CLqns2/Z+ !REGULAR PART NS
     $     DELTA1*(CLQNSP_delta2)/ZJAC )*(A2PI/2)**2 !DELTA PART 
     $        + X*(SIGMA*(CLqs2+CLqns2)/Z+ !REGULAR PART S
     $     SIGMA1*(CLQNSP_delta2)/ZJAC )*(A2PI/2)**2 !DELTA PART S

C GLUON
       FLNNLOg=X*((GL)*CLg2/Z) *(A2PI/2)**2   !ONLY REGULAR PART  G

C             
       FLNNLO=(FLNNLOq+FLNNLOg)*W

C F3 STRUCTURE FUNCTION
      CALL C2POLNS(Z,C3qns_sing2,C3qns_delta2,C3qns_reg2p,C3qns_reg2m)     

       F3NNLOq=X*(A2PI/2)**2 *  W*  (C3qns_reg2p*CF3/Z  !REGULAR PART NS
     $   + C3qns_reg2m*CF3min/Z  !REGULAR PART NS   
     $   + C3qns_sing2*(CF3/Z-CF31)! PLUS DISTRIBUTION NS
     $   + CF31* C3qns_delta2/ZJAC )   ! DELTA TERM  NS 


C INCLUSIVE NNLO CROSS SECTION

       WWW(2)=WWW(1)+ 2*PI/Y/Q2**2*ALPHAE**2 *NRM/NEV*   
     $     ( (1+(1-Y)**2)*F2NNLO -Y**2*FLNNLO+ (1-(1-Y)**2)*F3NNLOq)

C Scale Dependence
              WWW(2)= WWW(2)+ (WWW(1)-WWW(0)) * 
     $         LOG(MUR**2/MUF**2)*(11*CA-2*NF)/(6)*(A2PI)           

 
      IF (INNLO.NE.1) WWW(2)=WWW(1)  



C ----- POLARIZED INCLUSIVE CROSS SECTION -----

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

c     
       G2NLOq=X*A2PI/2*  W*  (C2qns_reg1*delta/Z  !REGULAR PART NS
     $   +    C2qns_sing1*(delta/Z-delta1)! PLUS DISTRIBUTION NS
     $   + delta1* C2qns_delta1/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*A2PI/2*  W*  (C2qs_reg1*sigma/Z  !REGULAR PART  S
     $   +    C2qs_sing1*(sigma/Z-sigma1)! PLUS DISTRIBUTION  S
     $   + sigma1* C2qs_delta1/ZJAC )   ! DELTA TERM   S


C EW PART (NO GLUONS NOR TRIANGLES) 
       GLNLOq=X*DELTAW*CLqns1/Z*A2PI/2 *W !ONLY REGULAR PART NS 
     $        + X*SIGMAW*CLqs1/Z*A2PI/2 *W !ONLY REGULAR PART S
c
       G4NLOq=X*A2PI/2*  W*  (C4qns_reg1*deltaw/Z  !REGULAR PART NS
     $   +    C4qns_sing1*(deltaw/Z-deltaw1)! PLUS DISTRIBUTION NS
     $   + deltaw1* C4qns_delta1/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*A2PI/2*  W*  (C4qs_reg1*sigmaw/Z  !REGULAR PART  S
     $   +    C4qs_sing1*(sigmaw/Z-sigmaw1)! PLUS DISTRIBUTION  S
     $   + sigmaw1* C4qs_delta1/ZJAC )   ! DELTA TERM   S


C GLUON
       G2NLOg=X*(GL)*C2g1/Z*A2PI/2 *W !ONLY REGULAR PART  G

       G2NLO=G2NLOq+G2NLOg
              
       NLOPREF=2*PI/Y/Q2**2*ALPHAE**2 *NRM/NEV
       WWW(1)=WWW(0)+ NLOPREF *
     $     ( (1-(1-Y)**2)*G2NLO
     $       + (1+(1-Y)**2)*G4NLOq - Y**2 * GLNLOq )
       CALL ACCUMULATE_NLO_BLOCKS(NLOPREF,G2NLOq,G2NLOg,G4NLOq,
     $  GLNLOq,Y)


C -- SECOND ORDER IN ALPHAS --

      LOGZ=DLOG(Z)  
C      
      DIMZ=WGPLG(1,1,-Z)
      TRIMZ=WGPLG(2,1,-Z)
      SMZ=WGPLG(1,2,-Z)   
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
      

c NON SINGLET G1 (WE KEEP CALLING THE COEFFICIENT C2)
       CALL C2POLNS(Z,C2QNS_SING2,C2QNS_DELTA2,C2QNS_REG2P,C2QNS_REG2M)

c PURE SINGLET G1
       CALL C2POLPS(Z,C2QPS_REG)
        
       C2QS_SING2 = C2QNS_SING2
       C2QS_DELTA2 = C2QNS_DELTA2
       C2QS_REG2 = C2QNS_REG2P + C2QPS_REG
        
C GLUON G1
       CALL C2POLG(Z,C2G_REG2)


C G1 NNLO STRUCTURE FUNCTION (WE CALL IT G2 FOR "REASONS")

       G2NNLOq=X*(A2PI/2)**2 * (C2qns_reg2p*DELTA/Z  !REGULAR PART NS
     $   + C2qns_reg2m*CF2min/Z  !REGULAR PART MINUS
     $   + C2qns_sing2*(delta/Z-delta1)! PLUS DISTRIBUTION NS
     $   + delta1* C2qns_delta2/ZJAC )   ! DELTA TERM  NS
c     
     $   +  X*(A2PI/2)**2 *  (C2qs_reg2*sigma/Z  !REGULAR PART  S
     $   +    C2qs_sing2*(sigma/Z-sigma1)! PLUS DISTRIBUTION  S
     $   + sigma1* C2qs_delta2/ZJAC )   ! DELTA TERM   S

       G2NNLOg=X*(GL)*C2G_REG2/Z*(A2PI/2)**2  !ONLY REGULAR PART  G

       G2NNLO=(G2NNLOq+G2NNLOg)*W


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

       GLNNLOq=X*DELTAW*CLqns2/Z*(A2PI/2)**2 *W !ONLY REGULAR PART NS 
     $        + X*SIGMAW*CLqns2/Z*(A2PI/2)**2 *W !ONLY REGULAR PART S (equal to ns)
c
       G4NNLOq=X*(A2PI/2)**2*  W*  (C4qns_reg2p*deltaw/Z  !REGULAR PART NS
     $   + C4qns_reg2m*CF3min/Z           !REGULAR PART MINUS  
     $   + C4qns_sing2*(deltaw/Z-deltaw1) !PLUS DISTRIBUTION NS
     $   + deltaw1* C4qns_delta2/ZJAC )   !DELTA TERM  NS
c     
     $   +  X*(A2PI/2)**2*  W*  (C4qs_reg2*sigmaw/Z  !REGULAR PART  S
     $   +    C4qs_sing2*(sigmaw/Z-sigmaw1)! PLUS DISTRIBUTION  S
     $   + sigmaw1* C4qs_delta2/ZJAC )   ! DELTA TERM   S



C INCLUSIVE NNLO CROSS SECTION

              WWW(2)= WWW(1)+ 2*PI/Y/Q2**2*ALPHAE**2 *NRM/NEV*   
     $     ( (1-(1-Y)**2)*G2NNLO
     $       + (1+(1-Y)**2)*G4NNLOq - Y**2 * GLNNLOq )   
C Scale Dependence
              WWW(2)= WWW(2) + (WWW(1)-WWW(0)) * 
     $         LOG(MUR**2/MUF**2)*(11*CA-2*NF)/(6)*(A2PI)  
     
       IF (INNLO.NE.1) WWW(2)=WWW(1)    
              
       
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
      INTEGER VOGT,SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
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
       C2QNSPLUS_REG =   ! ACTUALLY PLUS+MINUS
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
      INTEGER VOGT,SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
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
      INTEGER VOGT,SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
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
      INTEGER VOGT,SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3


      IF (vogt.eq.1) then
C LONGITUDINAL STRUCTURE FUNCTION FROM VOGT + SCALE FROM VN
c      pure singlet
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
c      pure singlet
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
      INTEGER SCHEME,NF
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON/ZETLOG/LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3


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
      C2QNS_REG2M = -C2QNSMIN_REG

      END
C-----------------------------------------------------------------------
      SUBROUTINE C2POLPS(Z,C2QPS_REG)
      IMPLICIT NONE
C---CALCULATE THE C2PS COEFFICIENTS
      DOUBLE PRECISION Z,C2QPS_REG
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      INTEGER SCHEME,NF
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
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
      INTEGER SCHEME,NF
      DOUBLE PRECISION LOGZ,LOGMU,DIMZ,TRIMZ,SMZ,DIZ,TRIZ,LOG1MZ,DI1MZ,
     $    TRI1MZ,S1MZ,LOG1PZ,TRIPCO,TRIMCO,L0,L1,L2,L3,L0X,L1X,L2X,L3X,
     $    ZETA2,ZETA3

      COMMON/COLFAC/CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
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
      DOUBLEPRECISION FUNCTION COMBPDFVEC(P,F,IMINUS)
C---COMBINE THE QUARK PDFS FOR THE F2 LIKE STRUCTURE FUNCTIONS.  
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF,IMINUS
      DOUBLE PRECISION P(4,7),F(-6:6),DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3),CHARGE
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG


      IF (IMINUS.EQ.1) THEN
      CHARGE = - LEPCHARGE
      ELSE
      CHARGE = LEPCHARGE
      ENDIF 

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

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
      DOUBLEPRECISION FUNCTION COMBPDFAX(P,F,IMINUS)
C---COMBINE THE QUARK PDFS FOR THE F3 LIKE STRUCTURE FUNCTIONS. 
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF,IMINUS
      DOUBLE PRECISION P(4,7),F(-6:6),DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3),CHARGE
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

      IF (IMINUS.EQ.1) THEN
      CHARGE = - LEPCHARGE
      ELSE
      CHARGE = LEPCHARGE
      ENDIF 

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

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
      SUBROUTINE COMBINEGLUON(P,M,CS,CSAX)
C---COMBINE THE GLUON CHANNEL CROSS SECTION
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF
      DOUBLE PRECISION P(4,7),M(-6:6),CS,CSAX,DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3)
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

      IF (EWFLAG.LT.3) THEN
        DO I=1,NF
C         Photon Part: 
          M(0)=M(0)+EQ(I)**2*CS
          IF (EWFLAG.GE.1) THEN
C           Z Part: VVVV-like coupligs:        
            M(0)=M(0)+CS*PROPZ*(CVZ(I)**2+CAZ(I)**2)*(CVZE**2+CAZE**2)     
C           Z Part: AVAV-like coupligs:
            M(0)=M(0)+CSAX*PROPZ*(4*CVZ(I)*CAZ(I)*CVZE*CAZE)
C           Photon/Z interference terms (null for W)
            M(0)=M(0)+2*EQ(I)*CS*PROPZ*INTGZ*(CVZ(I)*CVZE)
     $             +2*EQ(I)*CSAX*PROPZ*INTGZ*(CAZ(I)*CAZE)
          ENDIF
        ENDDO
      ELSEIF (LEPCHARGE.LT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: VVVV-like coupligs:        
            M(0)=M(0)+CS*PROPZ*(CVZ(2*I)**2+CAZ(2*I)**2)
     $                  *(CVZE**2+CAZE**2)*CKM2(I,J)         
C           W Part: AVAV-like coupligs:
            M(0)=M(0)+CSAX*PROPZ*(4*CVZ(2*I)*CAZ(2*I)*CVZE*CAZE)
     $                    *CKM2(I,J)  
          ENDDO
        ENDDO
      ELSEIF (LEPCHARGE.GT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: VVVV-like coupligs:        
            M(0)=M(0)+CS*PROPZ*(CVZ(2*I-1)**2+CAZ(2*I-1)**2)
     $                  *(CVZE**2+CAZE**2)*CKM2(J,I)            
C           W Part: AVAV-like coupligs:
            M(0)=M(0)+CSAX*PROPZ*(4*CVZ(2*I-1)*CAZ(2*I-1)*CVZE*CAZE)
     $                    *CKM2(J,I)
          ENDDO
        ENDDO
      ENDIF  
      END
C-----------------------------------------------------------------------
      SUBROUTINE COMBINEGTRIANG(P,M,CS,CSAX)
C---COMBINE THE GLUON CHANNEL CROSS SECTION
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF
      DOUBLE PRECISION P(4,7),M(-6:6),CS,CSAX,DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3)
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

      DO I=1,NF
C       Z Part: VVVV-like coupligs:        
        M(0)=M(0)+CS*PROPZ*(CAZ(I)**2)*(CVZE**2+CAZE**2)
C       Z Part: AVAV-like coupligs:
        M(0)=M(0)+CSAX*PROPZ*(2*CVZ(I)*CAZ(I)*CVZE*CAZE)
C       Photon/Z interference terms
        M(0)=M(0)+1*EQ(I)*CSAX*PROPZ*INTGZ*(CAZ(I)*CAZE)
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE COMBINEQI(P,M,CS,CSAX)
C---COMBINE THE QUARK CROSS SECTION WITH THE BOSON COUPLED TO THE
C   INITIAL QUARK LINE.  
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF
      DOUBLE PRECISION P(4,7),M(-6:6),CS,CSAX,DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3)
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

      IF (EWFLAG.LT.3) THEN
        DO I=-6,6
C         Photon Part: 
          M(I)=M(I)+EQ(I)**2*CS
          IF (EWFLAG.GE.1) THEN
C           Z Part: VVVV-like coupligs:        
            M(I)=M(I)+CS*PROPZ*(CVZ(I)**2+CAZ(I)**2)
     $                *(CVZE**2+CAZE**2)
C           Z Part: AVAV-like coupligs:
            M(I)=M(I)+CSAX*PROPZ*(4*CVZ(I)*CAZ(I)*CVZE*CAZE)
C           Photon/Z interference terms (null for W)
            M(I)=M(I)+2*EQ(I)*CS*PROPZ*INTGZ*(CVZ(I)*CVZE)
     $               +2*EQ(I)*CSAX*PROPZ*INTGZ*(CAZ(I)*CAZE)
          ENDIF
        ENDDO
      ELSEIF (LEPCHARGE.LT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: VVVV-like coupligs: 
            M(2*I)=M(2*I)+CS*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(2*I)**2+CAZ(2*I)**2)*CKM2(I,J)
            M(-2*I+1)=M(-2*I+1)+CS*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(-2*I+1)**2+CAZ(-2*I+1)**2)*CKM2(J,I)
C           W Part: AVAV-like coupligs:
            M(2*I)=M(2*I)+CSAX*PROPZ*(4*CVZE*CAZE)
     $                    *(CVZ(2*I)*CAZ(2*I)*CKM2(I,J))
            M(-2*I+1)=M(-2*I+1)+CSAX*PROPZ*(4*CVZE*CAZE)
     $                     *(CVZ(-2*I+1)*CAZ(-2*I+1)*CKM2(J,I))
          ENDDO
        ENDDO
      ELSEIF (LEPCHARGE.GT.0) THEN
        DO I=1,3
          DO J=1,3
C           W Part: VVVV-like coupligs: 
            M(2*I-1)=M(2*I-1)+CS*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(2*I-1)**2+CAZ(2*I-1)**2)*CKM2(J,I)
            M(-2*I)=M(-2*I)+CS*PROPZ*(CVZE**2+CAZE**2)
     $                  *(CVZ(-2*I)**2+CAZ(-2*I)**2)*CKM2(I,J)
C           W Part: AVAV-like coupligs:
            M(2*I-1)=M(2*I-1)+CSAX*PROPZ*(4*CVZE*CAZE)
     $                    *(CVZ(2*I-1)*CAZ(2*I-1)*CKM2(J,I))
            M(-2*I)=M(-2*I)+CSAX*PROPZ*(4*CVZE*CAZE)
     $                     *(CVZ(-2*I)*CAZ(-2*I)*CKM2(I,J))
          ENDDO
        ENDDO
      ENDIF
      END
C-----------------------------------------------------------------------
      SUBROUTINE COMBINEQITRIANG(P,M,CS,CSAX)
C---COMBINE THE TRIANGLE QUARK CONTRIBUTIONS TO VIRTUAL PART  
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF
      DOUBLE PRECISION P(4,7),M(-6:6),CS,CSAX,DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3)
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

      DO I=-6,6
C       Z Part: VVVV-like coupligs:        
        M(I)=M(I)+CS*PROPZ*(CAZ(I)**2)*(CVZE**2+CAZE**2)
C       Z Part: AVAV-like coupligs:
        M(I)=M(I)+CSAX*PROPZ*(2*CVZ(I)*CAZ(I)*CVZE*CAZE)
C       Photon/Z interference terms
        M(I)=M(I)+1*EQ(I)*CSAX*PROPZ*INTGZ*(CAZ(I)*CAZE)
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE COMBINEQF(P,M,CS,CSAX)
C---COMBINE THE QUARK CROSS SECTION WITH THE BOSON COUPLED TO A FINAL
C   QUARK LINE
      IMPLICIT NONE
      INTEGER I,J,K,EWFLAG,SCHEME,NF
      DOUBLE PRECISION P(4,7),M(-6:6),CS,CSAX,DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3)
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)


      IF (EWFLAG.LT.3) THEN
        DO I=-6,6
          IF (I.NE.0) THEN
          DO J=1,NF
C           Photon Part: 
            M(I)=M(I)+EQ(J)**2*CS
            IF (EWFLAG.GE.1) THEN
C             Z Part: VVVV-like coupligs:        
              M(I)=M(I)+CS*PROPZ*(CVZ(J)**2+CAZ(J)**2)
     #                 *(CVZE**2+CAZE**2)
C             Z Part: AVAV-like coupligs:
              M(I)=M(I)+CSAX*PROPZ*(4*CVZ(J)*CAZ(J)*CVZE*CAZE)
C             Photon/Z interference terms (null for W)
              M(I)=M(I)+2*EQ(J)*CS*PROPZ*INTGZ*(CVZ(J)*CVZE)
     $                 +2*EQ(J)*CSAX*PROPZ*INTGZ*(CAZ(J)*CAZE)
            ENDIF
          ENDDO
          ENDIF
        ENDDO


      ELSEIF (LEPCHARGE.LT.0) THEN
        DO I=-6,6
           IF (I.NE.0) THEN          
           DO J=1,3
             DO K=1,3
C                W Part: VVVV-like coupligs:        
                 M(I)=M(I)+CS*PROPZ*(CVZ(2*J)**2+CAZ(2*J)**2)*
     #            (CVZE**2+CAZE**2)*CKM2(J,K)
C                W Part: AVAV-like coupligs:
                 M(I)=M(I)+CSAX*PROPZ*(4*CVZ(2*J)*CAZ(2*J)
     #                      *CVZE*CAZE)*CKM2(J,K)           
             ENDDO
           ENDDO
           ENDIF        
        ENDDO    

      ELSEIF (LEPCHARGE.GT.0) THEN
        DO I=-6,6
           IF (I.NE.0) THEN          
           DO J=1,3
             DO K=1,3
C                W Part: VVVV-like coupligs:        
                 M(I)=M(I)+CS*PROPZ*(CVZ(2*J-1)**2+CAZ(2*J-1)**2)*
     #            (CVZE**2+CAZE**2)*CKM2(K,J)
C                W Part: AVAV-like coupligs:
                 M(I)=M(I)+CSAX*PROPZ*(4*CVZ(2*J-1)
     #                      *CAZ(2*J-1)*CVZE*CAZE)*CKM2(K,J)

             ENDDO
           ENDDO
           ENDIF           
        ENDDO    
      ENDIF

      END
C-----------------------------------------------------------------------
      SUBROUTINE COMBINEQIF(P,M,CS,AX1,AX2)
C---COMBINE THE QUARK CROSS SECTION WITH THE BOSON COUPLED TO BOTH THE 
C   FINAL AND INITIAL QUARK LINES. AX1 (I) AND AX2 (F) INDICATE THE
C   TYPE OF COUPLING TO THE BOSON (VECTOR: AX=0, AXIAL: AX=1)
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF,AX1,AX2
      DOUBLE PRECISION P(4,7),M(-6:6),CS,DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3)
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

      DO I=-6,6
        IF (I.NE.0) THEN
        DO J=1,NF
C         Photon Part: 
          M(I)=M(I)+EQ(J)*EQ(I)*CS*(1-AX1)*(1-AX2)
          IF (EWFLAG.GE.1) THEN
C           Z Part: VVVV-like coupligs:        
            M(I)=M(I)+CS*PROPZ*(CVZE**2+CAZE**2)
     $         * (CVZ(J)*CVZ(I)*(1-AX1)*(1-AX2) + CAZ(J)*CAZ(I)*AX1*AX2)       
C           Z Part: AVAV-like coupligs:
            M(I)=M(I)+CS*PROPZ*(CVZE*CAZE)*2
     $         * (CVZ(J)*CAZ(I)*(1-AX2)*AX1 + CAZ(J)*CVZ(I)*(1-AX1)*AX2)      
C           Photon/Z interference terms
            M(I)=M(I)+CS*PROPZ*INTGZ*( EQ(J)*CVZ(I)*CVZE*(1-AX1)*(1-AX2)
     $                                +EQ(I)*CVZ(J)*CVZE*(1-AX1)*(1-AX2)
     $                                +EQ(I)*CAZ(J)*CAZE*(1-AX1)*AX2
     $                                +EQ(J)*CAZ(I)*CAZE*AX1*(1-AX2) )
          ENDIF
        ENDDO
        ENDIF
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE COMBINEQIW(P,M,CS,CSAX)
C---COMBINE THE QUARK CROSS SECTION WITH THE BOSON COUPLED TO THE
C   INITIAL QUARK LINE. SAME AS COMBINEQI, BUT USES THE OPOSITE 
C   FLAVOURS FOR W!  
      IMPLICIT NONE
      INTEGER I,J,EWFLAG,SCHEME,NF
      DOUBLE PRECISION P(4,7),M(-6:6),CS,CSAX,DOT
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      DOUBLE PRECISION CVZ(-6:6),CAZ(-6:6),CVZE,CAZE,QZ,QQZ,GZ,MZ,
     $                 GAMMAZ,PROPZ,INTGZ,LEPCHARGE,CKM2(3,3) 
      COMMON /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      COMMON /EWCOUPLINGS/ CVZ,CAZ,CVZE,CAZE,MZ,GAMMAZ,LEPCHARGE,CKM2,
     $                     EWFLAG

C     MASSIVE BOSON SQUARE AND INTERFERENCE PROPAGATORS MODIFIERS     
      PROPZ=DOT(P,5,5)**2/((DOT(P,5,5)-MZ**2)**2+(MZ*GAMMAZ)**2)
      INTGZ=LEPCHARGE*(DOT(P,5,5)-MZ**2)/DOT(P,5,5)

C     ELECTRON (LEPTON)
      IF (LEPCHARGE.GT.0) THEN
      DO I=1,3
        DO J=1,3
C         W Part: VVVV-like coupligs:        
          M(2*I)=M(2*I)+CS*PROPZ*(CVZ(2*I)**2+CAZ(2*I)**2)
     $              *(CVZE**2+CAZE**2)*CKM2(I,J)
          M(-2*I+1)=M(-2*I+1)+CS*PROPZ*(CVZ(-2*I+1)**2+CAZ(-2*I+1)**2)
     $              *(CVZE**2+CAZE**2)*CKM2(J,I)
C         W Part: AVAV-like coupligs:
          M(2*I)=M(2*I)+CSAX*PROPZ
     #                 *(4*CVZ(2*I)*CAZ(2*I)*CVZE*CAZE)*CKM2(I,J)
          M(-2*I+1)=M(-2*I+1)+CSAX*PROPZ
     #                 *(4*CVZ(-2*I+1)*CAZ(-2*I+1)*CVZE*CAZE)*CKM2(J,I)     
        ENDDO
      ENDDO

C     POSITRON (ANTI-LEPTON)
      ELSEIF (LEPCHARGE.LT.0) THEN
      DO I=1,3
        DO J=1,3
C         W Part: VVVV-like coupligs:        
          M(2*I-1)=M(2*I-1)+CS*PROPZ*(CVZ(2*I-1)**2+CAZ(2*I-1)**2)
     $              *(CVZE**2+CAZE**2)*CKM2(J,I)
          M(-2*I)=M(-2*I)+CS*PROPZ*(CVZ(-2*I)**2+CAZ(-2*I)**2)
     $              *(CVZE**2+CAZE**2)*CKM2(I,J)
C         W Part: AVAV-like coupligs:
          M(2*I-1)=M(2*I-1)+CSAX*PROPZ
     #                 *(4*CVZ(2*I-1)*CAZ(2*I-1)*CVZE*CAZE)*CKM2(J,I)
          M(-2*I)=M(-2*I)+CSAX*PROPZ
     #                 *(4*CVZ(-2*I)*CAZ(-2*I)*CVZE*CAZE)*CKM2(I,J)     
        ENDDO
      ENDDO        
      ENDIF

      END
C-----------------------------------------------------------------------      
      SUBROUTINE EVOPDF(P,S,NF,SCREN,SCFACT,F,ALPHAS,ISAVE)
C     THIS ROUTINE EVOLVES THE PDFS AND GETS ALPHAS. IT CAN SAVE THE
C     RESULTS TO AVOID REDOING THE PDF EVOLUTION WHEN IT'S NOT NECESARY   
      IMPLICIT NONE
      DOUBLE PRECISION P(4,7),S,SCREN,SCFACT,Q2   
      DOUBLE PRECISION ETA,MUF,MUR,DOT,ALPHASPDF
      DOUBLE PRECISION F(-6:6),FNEW(-6:6),FSAVE(-6:6),ALPHAS,ALPHASAVE 
      INTEGER ISAVE,ICH,IFL,I,NF
      COMMON/CHANNEL/ICH,IFL
 
      SAVE FSAVE,ALPHASAVE

C--- SCALES DEFINED IN JETALG
      MUR=SQRT(SCREN)
      MUF=SQRT(SCFACT)
C--- CALCULATE MOMENTUM FRACTION AND KINEMATICS
      ETA=2*DOT(P,1,6)/S
C CALL LHAPDF TO COMPUTE PDFS AND GET ALPHAS
      CALL EVOLVEPDF(ETA,MUF,F)
      ALPHAS= alphasPDF(MUR)
C      ALPHAS=0.118d0 ! AP: FIXED THIS FOR COMPARISONS TO HERWIG
C SELECT INITIAL STATE CHANNELS  
      IF (ICH.EQ.0) THEN
      ELSEIF (ICH.EQ.1) THEN ! GLUONS TO 0, ONLY QUARKS
      F(0)=0D0
            IF (IFL.NE.0) THEN
            DO I=1,NF
            IF(I.NE.IFL) THEN
               F(I)=0D0
               F(-I)=0D0
           ENDIF
           ENDDO
           ENDIF
      ELSEIF (ICH.EQ.2) THEN ! QUARKS TO 0, ONLY GLUONS
           DO I=1,NF
               F(I)=0D0
               F(-I)=0D0
           ENDDO
      ELSE
         WRITE(6,*)"WRONG CHANNEL= ",ICH
         STOP
      ENDIF
c     SAVE VARIABLES FOR LATER RETRIVAL 
      IF (ISAVE.EQ.1) THEN
      FSAVE=F
      ALPHASAVE=ALPHAS
      ENDIF

      RETURN

      ENTRY GETPDF(F,ALPHAS) 
      F=FSAVE
      ALPHAS=ALPHASAVE
      END
C-----------------------------------------------------------------------
      SUBROUTINE GETWEIGHT(N,NA,IT,P,S,WEIGHT,SCREN,SCFACT,Q2,NF,W,IW)
C     THE INTEGER FLAG IW CONTROLS WHAT TO DO WITH THE PDF CALL:
C       IW=0: CALLS GETPDF TO GET PREVIOUSLY SAVED PDFS VALUES
C       IW=1: CALLS EVOPDF AND SAVES THE PDFS VALUES
C       IW=2: CALLS EVOPDF AND DOES NOT SAVE THE PDFS VALUES     
      IMPLICIT NONE
      DOUBLE PRECISION P(4,7),S,SCREN,SCFACT,Q2,W(0:5),WEIGHT(-6:6)      
      DOUBLE PRECISION ETA,AS,A2PI,MUF,MUR,TOTWGT
      DOUBLE PRECISION CA,PI,DOT,F(-6:6),ALPHAE  
      INTEGER NF,I,N,NA,IT,INNLO,IW    
      COMMON/INNLO/INNLO
      COMMON/ALPHAQED/ALPHAE    
      CA=3
      PI=ATAN(1D0)*4
C
      W(0)=0d0
      W(1)=0d0
      W(2)=0d0

C SCALES DEFINED IN JETALG
      MUR=SQRT(SCREN)
      MUF=SQRT(SCFACT)
C---CALCULATE MOMENTUM FRACTION AND KINEMATICS
      ETA=2*DOT(P,1,6)/S
C---CALL THE PDFS AND ALPHAS      
      IF (IW.GE.1) THEN
         CALL EVOPDF(P,S,NF,SCREN,SCFACT,F,AS,2-IW)
      ELSE
         CALL GETPDF(F,AS)
      ENDIF
C GET THE COUPLING CONSTANT
c      ALPHAS= alphasPDF(MUR)
      A2PI=AS/2/PI       
C COMPUTE WEIGHTS      
        TOTWGT=0
C      IF(N.NE.5) THEN  
        DO I=-NF,NF
          TOTWGT=TOTWGT+WEIGHT(I)*f(I)
        ENDDO

C DEFINE WEIGHTS ACCORDING TO POWER OF ALPHAS^(0,1,2) (NA)
C
C W(0) WEIGHT FOR ORDER ALPHAS^0 CONTRIBUTION 
C W(1) WEIGHT FOR ORDER ALPHAS^0+ALPHAS^1 CONTRIBUTION 
C W(2) WEIGHT FOR ORDER ALPHAS^0+ALPHAS^1+ALPHAS^2 CONTRIBUTION 
C
C AND ACCOUNT FOR BETA FUNCTION CONTRIBUTION AT ALPHAS^2 ~ BETA_0*ORDER(1)        
      IF (NA.EQ.0) then  
              W(0)= TOTWGT*ALPHAE**2     
              W(1)=  W(0)   
              W(2)=  W(0) 
      ELSEIF (NA.EQ.1) then       
              W(0)= 0D0
              W(1)= TOTWGT*ALPHAE**2  * A2PI      
              W(2)= W(1) + TOTWGT*ALPHAE**2  * A2PI  
     $                   * LOG(MUR**2/Q2)*(11*CA-2*NF)/6 * A2PI    
      IF (INNLO.NE.1) W(2)=W(1)    
      ELSEIF (NA.EQ.2) then  
              W(0)= 0D0                     
              W(1)= 0D0                     
              W(2)=TOTWGT*ALPHAE**2  *(A2PI)**2
      ENDIF    
      END
      
      

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
         
      
             
       
