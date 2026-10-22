      PROGRAM GAMMAJETLOOP_CHECK
      IMPLICIT NONE
      INCLUDE 'coupl.inc'
      INCLUDE 'MadLoopParams.inc'
      REAL*8 P(0:3,4),ANS(0:3,0:1),E,C,ST,PHI,PREC(0:1)
      COMPLEX*16 JAMPL(3,2,1),W(20)
      COMMON/ML5_0_0_1_JAMPL/JAMPL
      INTEGER HELC(4,16),H,I,K,N,RC,NS,SAVEDMODE
      CHARACTER*512 RESOURCEPATH
      LOGICAL BYPASS_CHECK,ALWAYS_TEST_STABILITY
      COMMON/ML5_0_0_1_BYPASS_CHECK/BYPASS_CHECK,ALWAYS_TEST_STABILITY
      COMMON/ML5_0_0_1_HELCONFIGS/HELC
      CALL SETPARA('param_card.dat')
      RESOURCEPATH='../MadLoop5_resources'
      CALL SETMADLOOPPATH(RESOURCEPATH)
      READ(*,*) E,C,PHI
      ST=DSQRT(1D0-C*C)
      P=0D0
      DO I=1,4
        P(0,I)=E
      ENDDO
      P(3,1)=E
      P(3,2)=-E
      P(1,3)=E*ST*DCOS(PHI)
      P(2,3)=E*ST*DSIN(PHI)
      P(3,3)=E*C
      DO I=1,3
        P(I,4)=-P(I,3)
      ENDDO
      MU_R=2D0*E
      CALL UPDATE_AS_PARAM()
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
          NS=1
          IF(I.LE.2) NS=-1
          CALL VXXXXX(P(0,I),0D0,HELC(I,H),NS,W)
          WRITE(*,*) 'POL',I,(W(K),K=5,8)
        ENDDO
      ENDDO
      END
