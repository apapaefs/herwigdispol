      PROGRAM NATIVE_HELICITY_CHECK
      IMPLICIT NONE
      INCLUDE 'coupl.inc'
      INCLUDE '../../Source/MODEL/input.inc'
      INCLUDE '../../Source/MODEL/mp_coupl.inc'
      INCLUDE '../../Source/MODEL/mp_input.inc'
      INCLUDE 'MadLoopParams.inc'
      REAL*8 P(0:3,4),ANS(0:3,0:1),MU2,AIN,AEIN,PREC(0:1)
      COMPLEX*16 JAMPL(3,2,1),W(20,4)
      COMMON/ML5_0_0_1_JAMPL/JAMPL
      COMMON/ML5_0_0_1_W/W
      REAL*8 PS(0:3,4)
      COMMON/ML5_0_0_1_PSPOINT/PS
      INTEGER HELC(4,16),H,I,K,N,RC,SAVEDMODE
      CHARACTER*512 RESOURCEPATH
      LOGICAL BYPASS_CHECK,ALWAYS_TEST_STABILITY
      COMMON/ML5_0_0_1_BYPASS_CHECK/BYPASS_CHECK,ALWAYS_TEST_STABILITY
      COMMON/ML5_0_0_1_HELCONFIGS/HELC
      CALL SETPARA('param_card.dat')
      RESOURCEPATH='../MadLoop5_resources'
      CALL SETMADLOOPPATH(RESOURCEPATH)
C     Positive-energy incoming then outgoing momenta, in GeV.
      READ(*,*) MU2,AIN,AEIN
      DO I=1,4
        READ(*,*) P(:,I)
      ENDDO
      AS=AIN
      AEWM1=1D0/AEIN
      MU_R=DSQRT(MU2)
      G=DSQRT(4D0*3.141592653589793D0*AIN)
      INCLUDE '../../Source/MODEL/actualize_mp_ext_params.inc'
      CALL COUP()
      WRITE(*,*) 'COUPLING',G,ABS(GAL(1))
      DO N=1,4
        CALL ML5_0_0_1_SLOOPMATRIX_THRES(P,ANS,-1D0,PREC,RC)
      ENDDO
      WRITE(*,*) 'TOTAL',ANS(1,0),PREC(0),RC
      DO I=1,4
        WRITE(*,*) 'MOMENTUM',I,P(:,I)
      ENDDO
      DO H=1,16
        CALL ML5_0_0_1_SLOOPMATRIXHEL_THRES(P,H,ANS,
     $       -1D0,PREC,RC)
        WRITE(*,*) 'HEL',H,HELC(:,H),ANS(1,0),PREC(0),RC
C       Refresh the complex colour amplitude in the original frame: the
C       stability call may leave its final rotated-point amplitude cached.
        BYPASS_CHECK=.TRUE.
        SAVEDMODE=CTMODERUN
        CTMODERUN=1
        CALL ML5_0_0_1_SLOOPMATRIXHEL(P,H,ANS)
        CTMODERUN=SAVEDMODE
        BYPASS_CHECK=.FALSE.
        WRITE(*,*) 'AMP',H,JAMPL(1,1,1),JAMPL(1,2,1)
        DO I=1,4
C         Export the actual external wavefunctions used above. MadLoop's
C         phase-space repair can zero tiny beam transverse components;
C         rebuilding VXXXXX from the uncorrected input gives wrong phases.
          WRITE(*,*) 'POL',I,(W(K,I),K=5,8)
          WRITE(*,*) 'MOM_USED',I,PS(:,I)
        ENDDO
      ENDDO
      END
