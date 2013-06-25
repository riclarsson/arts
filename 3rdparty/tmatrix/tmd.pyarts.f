C   The code was originally  developed by Michael Mishchenko at the NASA      
C   Goddard Institute for Space Studies, New York. This research
C   was funded by the NASA Radiation Science Program.
                                                                       
C   The code can be used without limitations in any not-for-       
C   profit scientific research.  The only request is that in any      
C   publication using the code the source of the code be acknowledged  
C   and relevant references be made. 

C   This version code has been modified by Cory Davis (University 
C   of Edinburgh) for inclusion in the the PyARTS atmospheric radiative 
C   transfer package                                  
                                                                       
C CPD:8/12/03-removed STOP statements in TMD, added ERRMSG output variable.
      SUBROUTINE TMD(RAT,NDISTR,AXMAX,NPNAX,B,GAM,NKMAX,EPS,NP,LAM,MRR,
     &     MRI,DDELT,NPNA,NDGS,R1RAT,R2RAT,QUIET,REFF,VEFF,CEXT,CSCA,
     &     WALB,ASYMM,F11,F22,F33,F44,F12,F34,ERRMSG)
                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'tmd.par.f'
      REAL*8  LAM,MRR,MRI,X(NPNG2),W(NPNG2),S(NPNG2),SS(NPNG2),
     *        AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8  XG(1000),WG(1000),TR1(NPN2,NPN2),TI1(NPN2,NPN2),
     &        ALPH1(NPL),ALPH2(NPL),ALPH3(NPL),ALPH4(NPL),BET1(NPL),
     &        BET2(NPL),XG1(2000),WG1(2000),
     &        AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BE1(NPL),BE2(NPL),
     &     F11(NPNA),F22(NPNA),F33(NPNA),F44(NPNA),F12(NPNA),F34(NPNA)
      REAL*4
     &     RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),
     &     RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),
     &     IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),
     &     IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
      INTEGER QUIET
      CHARACTER ERRMSG*100
Cf2py intent(out) REFF,VEFF,CEXT,CSCA,W,ASYMM,F11,F22,F33,F44,F12,F34,ERRMSG
Cf2py depend(npna) F11,F22,F33,F44,F12,F34 
      COMMON /CT/ TR1,TI1
      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
      COMMON /CHOICE/ ICHOICE
      P=DACOS(-1D0)
 
C  OPEN FILES *******************************************************
 
C      OPEN (6,FILE='test')
C      OPEN (10,FILE='tmatr.write')
 
C  INPUT DATA ********************************************************
 
C      RAT=1.0 D0 
C      NDISTR=4
C      AXMAX=10D0
C      NPNAX=1
C      B=0.1D0
C      GAM=0.5D0
C      NKMAX=-1
C      EPS=0.33D0
C      NP=-1
C      LAM=1208.840556451613D0
C      MRR=1.78145992756 d0
C      MRI=0.00512840878218D0 
C      DDELT=0.001D0 
C      NPNA=19 
C      NDGS=2

      ICHOICE=1
      NCHECK=0
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1
C      WRITE (6,5454) ICHOICE,NCHECK
 5454 FORMAT ('ICHOICE=',I1,'  NCHECK=',I1)
      DAX=AXMAX/NPNAX
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      if (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
C      PRINT*, 8000, RAT
 8000 FORMAT ('RAT=',F8.6)
      IF(QUIET.EQ.0) THEN      
         IF(NP.EQ.-1.AND.EPS.GE.1D0) PRINT 7000,EPS
         IF(NP.EQ.-1.AND.EPS.LT.1D0) PRINT 7001,EPS
         IF(NP.GE.0) PRINT 7100,NP,EPS
         IF(NP.EQ.-2.AND.EPS.GE.1D0) PRINT 7150,EPS
         IF(NP.EQ.-2.AND.EPS.LT.1D0) PRINT 7151,EPS
      
         PRINT 7400, LAM,MRR,MRI
         PRINT 7200, DDELT
      ENDIF
 7000 FORMAT('RANDOMLY ORIENTED OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('RANDOMLY ORIENTED PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('RANDOMLY ORIENTED CHEBYSHEV PARTICLES, T',
     &       I1,'(',F5.2,')')
 7150 FORMAT('RANDOMLY ORIENTED OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('RANDOMLY ORIENTED PROLATE CYLINDERS, D/L=',F11.7)
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 7400 FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
      DDELT=0.1D0*DDELT
      DO 600 IAX=1,NPNAX
         AXI=AXMAX-DAX*DFLOAT(IAX-1)
         R1=R1RAT*AXI
         R2=R2RAT*AXI
         NK=MAX(IDINT(AXI*NKMAX/AXMAX+2),1) !MAX call added by Cory to avoid occasional        
         IF (NK.GT.1000) THEN               !problems with nkmax=-1
            WRITE (ERRMSG,8001) NK
            PRINT *,ERRMSG
            RETURN
         ENDIF
C used to be a STOP here
C         IF (NK.GT.1000) RETURN
         IF (NDISTR.EQ.3) CALL POWER (AXI,B,R1,R2)
 8001    FORMAT ('NK=',I4,' I.E., IS GREATER THAN 1000. ',
     &           'EXECUTION TERMINATED.')
         CALL GAUSS (NK,0,0,XG,WG)
         Z1=(R2-R1)*0.5D0
         Z2=(R1+R2)*0.5D0
         Z3=R1*0.5D0
         IF (NDISTR.EQ.5) GO TO 3
         DO I=1,NK
            XG1(I)=Z1*XG(I)+Z2
            WG1(I)=WG(I)*Z1
         ENDDO
         GO TO 4
    3    DO I=1,NK
            XG1(I)=Z3*XG(I)+Z3
            WG1(I)=WG(I)*Z3
         ENDDO
         DO I=NK+1,2*NK
            II=I-NK
            XG1(I)=Z1*XG(II)+Z2
            WG1(I)=WG(II)*Z1
         ENDDO
         NK=NK*2
    4    CALL DISTRB (NK,XG1,WG1,NDISTR,AXI,B,GAM,R1,R2,QUIET,
     &                REFF,VEFF,P)
         IF (QUIET.EQ.0) THEN
            PRINT 8002,R1,R2
 8002       FORMAT('R1=',F10.6,'   R2=',F10.6)
            IF (DABS(RAT-1D0).LE.1D-6) PRINT 8003, REFF,VEFF
            IF (DABS(RAT-1D0).GT.1D-6) PRINT 8004, REFF,VEFF
 8003       FORMAT('EQUAL-VOLUME-SPHERE REFF=',F8.4,'   VEFF=',F7.4)
 8004       FORMAT('EQUAL-SURFACE-AREA-SPHERE REFF=',F8.4,
     &           '   VEFF=',F7.4)
            PRINT 7250,NK
 7250       FORMAT('NUMBER OF GAUSSIAN QUADRATURE POINTS ',
     &           'IN SIZE AVERAGING =',I4)
         ENDIF
         DO I=1,NPL
            ALPH1(I)=0D0
            ALPH2(I)=0D0
            ALPH3(I)=0D0
            ALPH4(I)=0D0
            BET1(I)=0D0
            BET2(I)=0D0
         ENDDO      
         CSCAT=0D0
         CEXTIN=0D0
         L1MAX=0
         DO 500 INK=1,NK
            I=NK-INK+1
            A=RAT*XG1(I)
            XEV=2D0*P*A/LAM
            IXXX=XEV+4.05D0*XEV**0.333333D0
            INM1=MAX0(4,IXXX)
            IF (INM1.GE.NPN1) THEN
               WRITE (ERRMSG,7333) NPN1
               PRINT *,ERRMSG
               RETURN
            END IF
C            IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  
     &       '.  EXECUTION TERMINATED')
            QEXT1=0D0
            QSCA1=0D0
            DO 50 NMA=INM1,NPN1
               NMAX=NMA
               MMAX=1
               NGAUSS=NMAX*NDGS
               IF (NGAUSS.GT.NPNG1) THEN
                  WRITE (ERRMSG,7340) NGAUSS
                  PRINT *,ERRMSG
                  RETURN
               END IF
c               IF (NGAUSS.GT.NPNG1) STOP
 7340          FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',
     &                '  EXECUTION TERMINATED')
 7334          FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
 7335 FORMAT('                              NMAX1 =', I3,'  DC2=',D8.2,
     &      '  DC1=',D8.2)
               CALL CONST1(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,
     &                   DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                      DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                          +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
               ENDDO
               DSCA=DABS((QSCA1-QSCA)/QSCA)
               DEXT=DABS((QEXT1-QEXT)/QEXT)
C              PRINT 7334, NMAX,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               NMIN=DFLOAT(NMAX)/2D0+1D0
               DO 10 N=NMIN,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  DQSCA=DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                      +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  DQEXT=(TR1NN+TR1NN1)*DN1
                  DQSCA=DABS(DQSCA/QSCA)
                  DQEXT=DABS(DQEXT/QEXT)
                  NMAX1=N
                  IF (DQSCA.LE.DDELT.AND.DQEXT.LE.DDELT) GO TO 12
   10          CONTINUE
   12          CONTINUE
c              PRINT 7335, NMAX1,DQSCA,DQEXT
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 55
               IF (NMA.EQ.NPN1) THEN
                  WRITE (ERRMSG,7333) NPN1
                  PRINT *,ERRMSG
                  RETURN
               END IF
c               IF (NMA.EQ.NPN1) STOP      
   50       CONTINUE
   55       NNNGGG=NGAUSS+1
            IF (NGAUSS.EQ.NPNG1) PRINT 7336
            MMAX=NMAX1
            DO 150 NGAUS=NNNGGG,NPNG1
               NGAUSS=NGAUS
               NGGG=2*NGAUSS
 7336          FORMAT('WARNING: NGAUSS=NPNG1')
 7337          FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
               CALL CONST1(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,
     &                   DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                      DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO 104 N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                          +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104          CONTINUE
               DSCA=DABS((QSCA1-QSCA)/QSCA)
               DEXT=DABS((QEXT1-QEXT)/QEXT)
c              PRINT 7337, NGGG,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 155
               IF (NGAUS.EQ.NPNG1) PRINT 7336
  150       CONTINUE
  155       CONTINUE
            QSCA=0D0
            QEXT=0D0
            NNM=NMAX*2
            DO 204 N=1,NNM
               QEXT=QEXT+TR1(N,N)
  204       CONTINUE
            IF (NMAX1.GT.NPN4) THEN
               WRITE (ERRMSG,7550) NMAX1
               PRINT *,ERRMSG
            END IF
 7550       FORMAT ('NMAX1 = ',I3, ', i.e. greater than NPN4.',
     &              ' Execution terminated')
C            IF (NMAX1.GT.NPN4) STOP              
            DO 213 N2=1,NMAX1
               NN2=N2+NMAX
               DO 213 N1=1,NMAX1
                  NN1=N1+NMAX
                  ZZ1=TR1(N1,N2)
                  RT11(1,N1,N2)=ZZ1
                  ZZ2=TI1(N1,N2)
                  IT11(1,N1,N2)=ZZ2
                  ZZ3=TR1(N1,NN2)
                  RT12(1,N1,N2)=ZZ3
                  ZZ4=TI1(N1,NN2)
                  IT12(1,N1,N2)=ZZ4
                  ZZ5=TR1(NN1,N2)
                  RT21(1,N1,N2)=ZZ5
                  ZZ6=TI1(NN1,N2)
                  IT21(1,N1,N2)=ZZ6
                  ZZ7=TR1(NN1,NN2)
                  RT22(1,N1,N2)=ZZ7
                  ZZ8=TI1(NN1,NN2)
                  IT22(1,N1,N2)=ZZ8
                  QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                 +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
  213       CONTINUE
C           PRINT 7800,0,DABS(QEXT),QSCA,NMAX
            DO 220 M=1,NMAX1
               CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &                    DDR,DRR,DRI,NMAX,NCHECK)
               NM=NMAX-M+1
               NM1=NMAX1-M+1
               M1=M+1
               QSC=0D0
               DO 214 N2=1,NM1
                  NN2=N2+M-1
                  N22=N2+NM
                  DO 214 N1=1,NM1
                     NN1=N1+M-1
                     N11=N1+NM
                     ZZ1=TR1(N1,N2)
                     RT11(M1,NN1,NN2)=ZZ1
                     ZZ2=TI1(N1,N2)
                     IT11(M1,NN1,NN2)=ZZ2
                     ZZ3=TR1(N1,N22)
                     RT12(M1,NN1,NN2)=ZZ3
                     ZZ4=TI1(N1,N22)
                     IT12(M1,NN1,NN2)=ZZ4
                     ZZ5=TR1(N11,N2)
                     RT21(M1,NN1,NN2)=ZZ5
                     ZZ6=TI1(N11,N2)
                     IT21(M1,NN1,NN2)=ZZ6
                     ZZ7=TR1(N11,N22)
                     RT22(M1,NN1,NN2)=ZZ7
                     ZZ8=TI1(N11,N22)
                     IT22(M1,NN1,NN2)=ZZ8
                     QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                       +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2D0
  214          CONTINUE
               NNM=2*NM
               QXT=0D0
               DO 215 N=1,NNM
                  QXT=QXT+TR1(N,N)*2D0
  215          CONTINUE
               QSCA=QSCA+QSC
               QEXT=QEXT+QXT
C              PRINT 7800,M,DABS(QXT),QSC,NMAX
 7800          FORMAT(' m=',I3,'  qxt=',d12.6,'  qsc=',d12.6,
     &                '  nmax=',I3)
  220       CONTINUE
            COEFF1=LAM*LAM*0.5D0/P
            CSCA=QSCA*COEFF1
            CEXT=-QEXT*COEFF1
c           PRINT 7880, NMAX,NMAX1
 7880       FORMAT ('nmax=',I3,'   nmax1=',I3)
            CALL GSP (NMAX1,CSCA,LAM,AL1,AL2,AL3,AL4,BE1,BE2,LMAX)
            L1M=LMAX+1
            L1MAX=MAX(L1MAX,L1M)
            WGII=WG1(I)
            WGI=WGII*CSCA
            DO 250 L1=1,L1M
               ALPH1(L1)=ALPH1(L1)+AL1(L1)*WGI
               ALPH2(L1)=ALPH2(L1)+AL2(L1)*WGI
               ALPH3(L1)=ALPH3(L1)+AL3(L1)*WGI
               ALPH4(L1)=ALPH4(L1)+AL4(L1)*WGI
               BET1(L1)=BET1(L1)+BE1(L1)*WGI
               BET2(L1)=BET2(L1)+BE2(L1)*WGI
  250       CONTINUE
            CSCAT=CSCAT+WGI
            CEXTIN=CEXTIN+CEXT*WGII
C           PRINT 6070, I,NMAX,NMAX1,NGAUSS
 6070       FORMAT(4I6)
  500    CONTINUE
         DO 510 L1=1,L1MAX
            ALPH1(L1)=ALPH1(L1)/CSCAT
            ALPH2(L1)=ALPH2(L1)/CSCAT
            ALPH3(L1)=ALPH3(L1)/CSCAT
            ALPH4(L1)=ALPH4(L1)/CSCAT
            BET1(L1)=BET1(L1)/CSCAT
            BET2(L1)=BET2(L1)/CSCAT
  510    CONTINUE
         WALB=CSCAT/CEXTIN
         CALL HOVENR(L1MAX,ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,QUIET)
         ASYMM=ALPH1(2)/3D0
         IF (QUIET.EQ.0) THEN
            PRINT 9100,CEXTIN,CSCAT,WALB,ASYMM
 9100       FORMAT('CEXT=',D12.6,2X,'CSCA=',D12.6,2X,
     &           2X,'W=',D12.6,2X,'<COS>=',D12.6)
            IF (WALB.GT.1D0) PRINT 9111
 9111       FORMAT ('WARNING: W IS GREATER THAN 1')
c         WRITE (10,580) WALB,L1MAX
C         DO L=1,L1MAX
c            WRITE (10,575) ALPH1(L),ALPH2(L),ALPH3(L),ALPH4(L),
C     &                     BET1(L),BET2(L)           
C         ENDDO   
 575        FORMAT(6D14.7)
 580        FORMAT(D14.8,I8)
         ENDIF
         LMAX=L1MAX-1
         CALL MATR (ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,LMAX,NPNA,QUIET,
     &        F11,F22,F33,F44,F12,F34)
  600 CONTINUE
c      ITIME=MCLOCK()
c      TIME=DFLOAT(ITIME)/6000D0
c      IF (QUIET.EQ.0) PRINT 1001,TIME
c 1001 FORMAT (' time =',F8.2,' min')
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE CONST1 (NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'tmd.par.f'
      REAL*8 X(NPNG2),W(NPNG2),X1(NPNG1),W1(NPNG1),
     *        X2(NPNG1),W2(NPNG1),
     *        S(NPNG2),SS(NPNG2),
     *        AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)
 
      DO 10 N=1,NMAX
           NN=N*(N+1)
           AN(N)=DFLOAT(NN)
           D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
           DD(N)=D
           DO 10 N1=1,N
                DDD=D*DD(N1)*0.5D0
                ANN(N,N1)=DDD
                ANN(N1,N)=DDD
   10 CONTINUE
      NG=2*NGAUSS
      IF (NP.EQ.-2) GO  TO 11
      CALL GAUSS(NG,0,0,X,W)
      GO TO 19
   11 NG1=DFLOAT(NGAUSS)/2D0
      NG2=NGAUSS-NG1
      XX=-DCOS(DATAN(EPS))
      CALL GAUSS(NG1,0,0,X1,W1)
      CALL GAUSS(NG2,0,0,X2,W2)
      DO 12 I=1,NG1
         W(I)=0.5D0*(XX+1D0)*W1(I)
         X(I)=0.5D0*(XX+1D0)*X1(I)+0.5D0*(XX-1D0)
   12 CONTINUE
      DO 14 I=1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
   14 CONTINUE
      DO 16 I=1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
   16 CONTINUE
   19 DO 20 I=1,NGAUSS
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y
           SS(NG-I+1)=Y
           Y=DSQRT(Y)
           S(I)=Y
           S(NG-I+1)=Y
   20 CONTINUE
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE VARY (LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,
     *                 R,DR,DDR,DRR,DRI,NMAX)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),R(NPNG2),DR(NPNG2),MRR,MRI,LAM,
     *        Z(NPNG2),ZR(NPNG2),ZI(NPNG2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),DDR(NPNG2),
     *        DRR(NPNG2),DRI(NPNG2),
     *        DY(NPNG2,NPN1)
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
      CHARACTER*60 ERRMSG
      NG=NGAUSS*2
      IF (NP.EQ.-1) CALL RSP1(X,NG,NGAUSS,A,EPS,NP,R,DR)
      IF (NP.GE.0) CALL RSP2(X,NG,A,EPS,NP,R,DR)
      IF (NP.EQ.-2) CALL RSP3(X,NG,NGAUSS,A,EPS,R,DR)
      PI=P*2D0/LAM
      PPI=PI*PI
      PIR=PPI*MRR
      PII=PPI*MRI
      V=1D0/(MRR*MRR+MRI*MRI)
      PRR=MRR*V
      PRI=-MRI*V
      TA=0D0
      DO 10 I=1,NG
           VV=DSQRT(R(I))
           V=VV*PI
           TA=MAX(TA,V)
           VV=1D0/V
           DDR(I)=VV
           DRR(I)=PRR*VV
           DRI(I)=PRI*VV
           V1=V*MRR
           V2=V*MRI
           Z(I)=V
           ZR(I)=V1
           ZI(I)=V2
   10 CONTINUE
      IF (NMAX.GT.NPN1) THEN
         WRITE (ERRMSG,9000) NMAX,NPN1
         PRINT *, ERRMSG
      END IF
C     IF (NMAX.GT.NPN1) STOP
 9000 FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
      TB=TA*DSQRT(MRR*MRR+MRI*MRI)
      TB=DMAX1(TB,DFLOAT(NMAX))
      NNMAX1=1.2D0*DSQRT(DMAX1(TA,DFLOAT(NMAX)))+3D0
      NNMAX2=(TB+4D0*(TB**0.33333D0)+1.2D0*DSQRT(TB))
      NNMAX2=NNMAX2-NMAX+5
      CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RSP1 (X,NG,NGAUSS,REV,EPS,NP,R,DR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      A=REV*EPS**(1D0/3D0)
      AA=A*A
      EE=EPS*EPS
      EE1=EE-1D0
      DO 50 I=1,NGAUSS
          C=X(I)
          CC=C*C
          SS=1D0-CC
          S=DSQRT(SS)
          RR=1D0/(SS+EE*CC)
          R(I)=AA*RR
          R(NG-I+1)=R(I)
          DR(I)=RR*C*S*EE1
          DR(NG-I+1)=-DR(I)
   50 CONTINUE
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RSP2 (X,NG,REV,EPS,N,R,DR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      DNP=DFLOAT(N)
      DN=DNP*DNP
      DN4=DN*4D0
      EP=EPS*EPS
      A=1D0+1.5D0*EP*(DN4-2D0)/(DN4-1D0)
      I=(DNP+0.1D0)*0.5D0
      I=2*I
      IF (I.EQ.N) A=A-3D0*EPS*(1D0+0.25D0*EP)/
     *              (DN-1D0)-0.25D0*EP*EPS/(9D0*DN-1D0)
      R0=REV*A**(-1D0/3D0)
      DO 50 I=1,NG
         XI=DACOS(X(I))*DNP
         RI=R0*(1D0+EPS*DCOS(XI))
         R(I)=RI*RI
         DR(I)=-R0*EPS*DNP*DSIN(XI)/RI
   50 CONTINUE
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RSP3 (X,NG,NGAUSS,REV,EPS,R,DR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )
      A=H*EPS
      DO 50 I=1,NGAUSS
         CO=-X(I)
         SI=DSQRT(1D0-CO*CO)
         IF (SI/CO.GT.A/H) GO TO 20
         RAD=H/CO
         RTHET=H*SI/(CO*CO)
         GO TO 30
   20    RAD=A/SI
         RTHET=-A*CO/(SI*SI)
   30    R(I)=RAD*RAD
         R(NG-I+1)=R(I)
         DR(I)=-RTHET/RAD
         DR(NG-I+1)=-DR(I)
   50 CONTINUE
      RETURN
      END
 
C************************************************************************
 
      SUBROUTINE BESS (X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),XR(NG),XI(NG),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),
     *        AJ(NPN1),AY(NPN1),AJR(NPN1),AJI(NPN1),
     *        ADJ(NPN1),ADY(NPN1),ADJR(NPN1),
     *        ADJI(NPN1)
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
 
      DO 10 I=1,NG
           XX=X(I)
           CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
           CALL RYB(XX,AY,ADY,NMAX)
           YR=XR(I)
           YI=XI(I)
           CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,NNMAX2)
           DO 10 N=1,NMAX
                J(I,N)=AJ(N)
                Y(I,N)=AY(N)
                JR(I,N)=AJR(N)
                JI(I,N)=AJI(N)
                DJ(I,N)=ADJ(N)
                DY(I,N)=ADY(N)
                DJR(I,N)=ADJR(N)
                DJI(I,N)=ADJI(N)
   10 CONTINUE
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RJB(X,Y,U,NMAX,NNMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),U(NMAX),Z(800)
      L=NMAX+NNMAX
      XX=1D0/X
      Z(L)=1D0/(DFLOAT(2*L+1)*XX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         Z(I1)=1D0/(DFLOAT(2*I1+1)*XX-Z(I1+1))
    5 CONTINUE
      Z0=1D0/(XX-Z(1))
      Y0=Z0*DCOS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1
      DO 10 I=2,NMAX
         YI1=Y(I-1)
         YI=YI1*Z(I)
         U(I)=YI1-DFLOAT(I)*YI*XX
         Y(I)=YI
   10 CONTINUE
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RYB(X,Y,V,NMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),V(NMAX)
      C=DCOS(X)
      S=DSIN(X)
      X1=1D0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3D0*X3+X1)*C-3D0*X2*S
      NMAX1=NMAX-1
      DO 5 I=2,NMAX1
    5     Y(I+1)=DFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
      V(1)=-X1*(C+Y1)
      DO 10 I=2,NMAX
  10       V(I)=Y(I-1)-DFLOAT(I)*X1*Y(I)
      RETURN
      END
 
C**********************************************************************
C                                                                     *
C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
C   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE CJB (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 YR(NMAX),YI(NMAX),UR(NMAX),UI(NMAX)
      REAL*8 CYR(NPN1),CYI(NPN1),CZR(1200),CZI(1200),
     *       CUR(NPN1),CUI(NPN1)
      L=NMAX+NNMAX
      XRXI=1D0/(XR*XR+XI*XI)
      CXXR=XR*XRXI
      CXXI=-XI*XRXI 
      QF=1D0/DFLOAT(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1
      DO I=1,L1
         I1=L-I
         QF=DFLOAT(2*I1+1)
         AR=QF*CXXR-CZR(I1+1)
         AI=QF*CXXI-CZI(I1+1)
         ARI=1D0/(AR*AR+AI*AI)
         CZR(I1)=AR*ARI
         CZI(I1)=-AI*ARI
      ENDDO   
      AR=CXXR-CZR(1)
      AI=CXXI-CZI(1)
      ARI=1D0/(AR*AR+AI*AI)
      CZ0R=AR*ARI
      CZ0I=-AI*ARI
      CR=DCOS(XR)*DCOSH(XI)
      CI=-DSIN(XR)*DSINH(XI)
      AR=CZ0R*CR-CZ0I*CI
      AI=CZ0I*CR+CZ0R*CI
      CY0R=AR*CXXR-AI*CXXI
      CY0I=AI*CXXR+AR*CXXI
      CY1R=CY0R*CZR(1)-CY0I*CZI(1)
      CY1I=CY0I*CZR(1)+CY0R*CZI(1)
      AR=CY1R*CXXR-CY1I*CXXI
      AI=CY1I*CXXR+CY1R*CXXI
      CU1R=CY0R-AR
      CU1I=CY0I-AI
      CYR(1)=CY1R
      CYI(1)=CY1I
      CUR(1)=CU1R
      CUI(1)=CU1I
      YR(1)=CY1R
      YI(1)=CY1I
      UR(1)=CU1R
      UI(1)=CU1I
      DO I=2,NMAX
         QI=DFLOAT(I)
         CYI1R=CYR(I-1)
         CYI1I=CYI(I-1)
         CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
         CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
         AR=CYIR*CXXR-CYII*CXXI
         AI=CYII*CXXR+CYIR*CXXI
         CUIR=CYI1R-QI*AR
         CUII=CYI1I-QI*AI
         CYR(I)=CYIR
         CYI(I)=CYII
         CUR(I)=CUIR
         CUI(I)=CUII
         YR(I)=CYIR
         YI(I)=CYII
         UR(I)=CUIR
         UI(I)=CUII
      ENDDO   
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
 
      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      REAL*4 PLUS(NPN6*NPN4*NPN4*8)
      COMMON /TMAT/ PLUS,
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1D0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG ( X(I1), NMAX, 0, DV1, DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR12=0D0
                AR21=0D0
                AI12=0D0
                AI21=0D0
                GR12=0D0
                GR21=0D0
                GI12=0D0
                GI21=0D0
                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0D0) GO TO 205

                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
 
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    URI=DR(I)
                    RRI=RR(I)
 
                    F1=RRI*A22
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
  200           CONTINUE
 
  205           AN12=ANN(N1,N2)*FACTOR
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
  300 CONTINUE
 
      TPIR=PIR
      TPII=PII
      TPPI=PPI
 
      NM=NMAX
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=0D0
                TQI(K1,KK2)=0D0
                TRGQR(K1,KK2)=0D0
                TRGQI(K1,KK2)=0D0
 
                TQR(KK1,K2)=0D0
                TQI(KK1,K2)=0D0
                TRGQR(KK1,K2)=0D0
                TRGQI(KK1,K2)=0D0
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
      CALL TT(NMAX,NCHECK)
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE TMATR (M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
 
      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      REAL*4 PLUS(NPN6*NPN4*NPN4*8)
      COMMON /TMAT/ PLUS,
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1D0
      NM=NMAX+NMAX
      DO 5 N=1,NM
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG (X(I1),NMAX,M,DV1,DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           WR=W(I)*R(I)
           DS(I)=S(I)*QM*WR
           DSS(I)=SS(I)*QMM
           RR(I)=WR
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0D0
                AR12=0D0
                AR21=0D0
                AR22=0D0
                AI11=0D0
                AI12=0D0
                AI21=0D0
                AI22=0D0
                GR11=0D0
                GR12=0D0
                GR21=0D0
                GR22=0D0
                GI11=0D0
                GI12=0D0
                GI21=0D0
                GI22=0D0
                SI=SIG(N1+N2)
 
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
                    AA2=A11*DSS(I)+A22
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1
                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1
 
                    C7R=C4R*DDRI
                    C7I=C4I*DDRI
                    B7R=B4R*DDRI
                    B7I=B4I*DDRI
 
                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII
                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII
 
                    URI=DR(I)
                    DSI=DS(I)
                    DSSI=DSS(I)
                    RRI=RR(I)
 
                    IF (NCHECK.EQ.1.AND.SI.GT.0D0) GO TO 150
 
                    E1=DSI*AA1
                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I
                    IF (NCHECK.EQ.1) GO TO 160
 
  150               F1=RRI*AA2
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
                    IF (NCHECK.EQ.1) GO TO 200
 
  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1
                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I
  200           CONTINUE
                AN12=ANN(N1,N2)*FACTOR
                R11(N1,N2)=AR11*AN12
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                R22(N1,N2)=AR22*AN12
                I11(N1,N2)=AI11*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                I22(N1,N2)=AI22*AN12
                RG11(N1,N2)=GR11*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                RG22(N1,N2)=GR22*AN12
                IG11(N1,N2)=GI11*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                IG22(N1,N2)=GI22*AN12
 
  300 CONTINUE
      TPIR=PIR
      TPII=PII
      TPPI=PPI
      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22
 
                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
 
      CALL TT(NM,NCHECK)
 
      RETURN
      END
 
C*****************************************************************
 
      SUBROUTINE VIG (X, NMAX, M, DV1, DV2)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN1),DV2(NPN1)
 
      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DO N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
      ENDDO   
      IF (M.NE.0) GO TO 20
      D1=1D0
      D2=X  
      DO N=1,NMAX  
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1 
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      RETURN
   20 QMM=DFLOAT(M*M)
      DO I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
      ENDDO   
      D1=0D0
      D2=A 
      DO N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      RETURN
      END 
 
C**********************************************************************
C                                                                     *
C   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              *
C                                                                     *
C   INPUT INFORTMATION IS IN COMMON /CTT/                             *
C   OUTPUT INFORMATION IS IN COMMON /CT/                              *
C                                                                     *
C**********************************************************************
 
      SUBROUTINE TT(NMAX,NCHECK)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(NPN2,NPN2),B(NPN2),WORK(NPN2),
     *       QR(NPN2,NPN2),QI(NPN2,NPN2),
     *       RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *       A(NPN2,NPN2),C(NPN2,NPN2),D(NPN2,NPN2),E(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      COMPLEX*16 ZQ(NPN2,NPN2),ZW(NPN2)
      COMPLEX*16 ZQR(NPN2,NPN2),ZAFAC(NPN2,NPN2),ZT(NPN2,NPN2),
     &           ZTHETA(NPN2,NPN2)
      INTEGER IPIV(NPN2),IPVT(NPN2)
      COMMON /CHOICE/ ICHOICE
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
      NDIM=NPN2
      NNMAX=2*NMAX
      IF (ICHOICE.EQ.2) GO TO 5
 
C	Inversion from NAG-LIB or Waterman's method
 
	DO I=1,NNMAX
	   DO J=1,NNMAX
	      ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
	      ZAFAC(I,J)=ZQ(I,J)
	   ENDDO
	ENDDO
	IF (ICHOICE.EQ.1) THEN
	   INFO=0
C           CALL F07ARF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
           CALL ZGETRF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
           IF (INFO.NE.0) WRITE (6,1100) INFO
C           CALL F07AWF(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
           CALL ZGETRI(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
C           IF (INFO.NE.0) WRITE (6,1100) INFO
 1100      FORMAT ('WARNING:  info=', i2)
	   DO I=1,NNMAX
	      DO J=1,NNMAX
	         TR=0D0
	         TI=0D0
	         DO K=1,NNMAX
                    ARR=RGQR(I,K)
                    ARI=RGQI(I,K)
                    AR=ZQ(K,J)
                    AI=DIMAG(ZQ(K,J))
                    TR=TR-ARR*AR+ARI*AI
                    TI=TI-ARR*AI-ARI*AR
                 ENDDO
	         TR1(I,J)=TR
	         TI1(I,J)=TI
	      ENDDO
	   ENDDO
 
	   ELSE
	   IFAIL=0
C          CALL F01RCF(NNMAX,NNMAX,ZAFAC,NPN2,ZTHETA,IFAIL)
C          CALL F01REF('S',NNMAX,NNMAX,NNMAX,ZAFAC,
C    &                 NPN2,ZTHETA,ZW,IFAIL)
	   DO I=1,NNMAX
	      DO J=1,NNMAX
	         ZQ(I,J)=DCMPLX(DREAL(ZAFAC(I,J)),
     &                  -DIMAG(ZAFAC(I,J)))
	      ENDDO
	   ENDDO
	   DO I=1,NNMAX
	      DO J=1,NNMAX
                 IF (I.LE.NNMAX/2.AND.I.EQ.J) THEN
	            D(I,J)=-1D0
                    ELSE IF (I.GT.NNMAX/2.AND.I.EQ.J) THEN
	               D(I,J)=1D0
	               ELSE
	               D(I,J)=0D0
	         ENDIF
              ENDDO
           ENDDO
	   DO I=1,NNMAX
	      DO J=1,NNMAX
	         ZT(I,J)=DCMPLX(0D0,0D0)
	         DO K=1,NNMAX
	            ZT(I,J)=ZT(I,J)+D(I,I)
     &                     *ZQ(I,K)*D(K,K)*ZQ(J,K)
	         ENDDO
	         ZT(I,J)=0.5D0*(ZT(I,J)-D(I,J)**2)
	         TR1(I,J)=DREAL(ZT(I,j))
	         TI1(I,J)=DIMAG(ZT(i,j))
              ENDDO
	   ENDDO
	ENDIF
		
	GOTO 70
 
C  Gaussian elimination

    5 DO 10 N1=1,NNMAX
         DO 10 N2=1,NNMAX
            F(N1,N2)=QI(N1,N2)
   10 CONTINUE
      IF (NCHECK.EQ.1) THEN
          CALL INV1(NMAX,F,A)
        ELSE
          CALL INVERT(NDIM,NNMAX,F,A,COND,IPVT,WORK,B) 
      ENDIF
      CALL PROD(QR,A,C,NDIM,NNMAX)
      CALL PROD(C,QR,D,NDIM,NNMAX)
      DO 20 N1=1,NNMAX
           DO 20 N2=1,NNMAX
                C(N1,N2)=D(N1,N2)+QI(N1,N2)
   20 CONTINUE
      IF (NCHECK.EQ.1) THEN
          CALL INV1(NMAX,C,QI)
        ELSE
          CALL INVERT(NDIM,NNMAX,C,QI,COND,IPVT,WORK,B) 
      ENDIF
      CALL PROD(A,QR,D,NDIM,NNMAX)
      CALL PROD(D,QI,QR,NDIM,NNMAX)
 
      CALL PROD(RGQR,QR,A,NDIM,NNMAX)
      CALL PROD(RGQI,QI,C,NDIM,NNMAX)
      CALL PROD(RGQR,QI,D,NDIM,NNMAX)
      CALL PROD(RGQI,QR,E,NDIM,NNMAX)
      DO 30 N1=1,NNMAX
           DO 30 N2=1,NNMAX
                TR1(N1,N2)=-A(N1,N2)-C(N1,N2)
                TI1(N1,N2)= D(N1,N2)-E(N1,N2)
   30 CONTINUE
   70 RETURN
      END
 
C********************************************************************
 
      SUBROUTINE PROD(A,B,C,NDIM,N)
      REAL*8 A(NDIM,N),B(NDIM,N),C(NDIM,N),cij
      DO 10 I=1,N
           DO 10 J=1,N
                CIJ=0d0
                DO 5 K=1,N
                     CIJ=CIJ+A(I,K)*B(K,J)
    5           CONTINUE
                C(I,J)=CIJ
   10 CONTINUE
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE INV1 (NMAX,F,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'tmd.par.f'
      REAL*8  A(NPN2,NPN2),F(NPN2,NPN2),B(NPN1),
     *        WORK(NPN1),Q1(NPN1,NPN1),Q2(NPN1,NPN1),
     &        P1(NPN1,NPN1),P2(NPN1,NPN1)
      INTEGER IPVT(NPN1),IND1(NPN1),IND2(NPN1)
      NDIM=NPN1
      NN1=(DFLOAT(NMAX)-0.1D0)*0.5D0+1D0 
      NN2=NMAX-NN1
      DO 5 I=1,NMAX
         IND1(I)=2*I-1
         IF(I.GT.NN1) IND1(I)=NMAX+2*(I-NN1)
         IND2(I)=2*I
         IF(I.GT.NN2) IND2(I)=NMAX+2*(I-NN2)-1
    5 CONTINUE
      NNMAX=2*NMAX
      DO 15 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 15 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            Q1(J,I)=F(J1,I1)
            Q2(J,I)=F(J2,I2)
   15 CONTINUE
      CALL INVERT(NDIM,NMAX,Q1,P1,COND,IPVT,WORK,B)
      CALL INVERT(NDIM,NMAX,Q2,P2,COND,IPVT,WORK,B)
      DO 30 I=1,NNMAX
         DO 30 J=1,NNMAX
            A(J,I)=0D0
   30 CONTINUE
      DO 40 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 40 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            A(J1,I1)=P1(J,I)
            A(J2,I2)=P2(J,I)
   40 CONTINUE
      RETURN
      END
 
C*********************************************************************
 
      SUBROUTINE INVERT (NDIM,N,A,X,COND,IPVT,WORK,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),X(NDIM,N),WORK(N),B(N)
      INTEGER IPVT(N)
      CALL DECOMP (NDIM,N,A,COND,IPVT,WORK)
      IF (COND+1D0.EQ.COND) PRINT 5,COND
C     IF (COND+1D0.EQ.COND) STOP
    5 FORMAT(' THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY '
     *      ,'COND = ',D12.6)
      DO 30 I=1,N
           DO 10 J=1,N
                B(J)=0D0
                IF (J.EQ.I) B(J)=1D0
  10       CONTINUE
           CALL SOLVE (NDIM,N,A,B,IPVT)
           DO 30 J=1,N
                X(J,I)=B(J)
   30 CONTINUE
      RETURN
      END
 
C********************************************************************
 
      SUBROUTINE DECOMP (NDIM,N,A,COND,IPVT,WORK)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
      IPVT(N)=1
      IF(N.EQ.1) GO TO 80
      NM1=N-1
      ANORM=0D0
      DO 10 J=1,N
          T=0D0
          DO 5 I=1,N
              T=T+DABS(A(I,J))
    5     CONTINUE
          IF (T.GT.ANORM) ANORM=T
   10 CONTINUE
      DO 35 K=1,NM1
          KP1=K+1
          M=K
          DO 15 I=KP1,N
              IF (DABS(A(I,K)).GT.DABS(A(M,K))) M=I
   15     CONTINUE
          IPVT(K)=M
          IF (M.NE.K) IPVT(N)=-IPVT(N)
          T=A(M,K)
          A(M,K)=A(K,K)
          A(K,K)=T
          IF (T.EQ.0d0) GO TO 35
          DO 20 I=KP1,N
              A(I,K)=-A(I,K)/T
   20     CONTINUE
          DO 30 J=KP1,N
              T=A(M,J)
              A(M,J)=A(K,J)
              A(K,J)=T
              IF (T.EQ.0D0) GO TO 30
              DO 25 I=KP1,N
                  A(I,J)=A(I,J)+A(I,K)*T
   25         CONTINUE
   30     CONTINUE
   35 CONTINUE
      DO 50 K=1,N
          T=0D0
          IF (K.EQ.1) GO TO 45
          KM1=K-1
          DO 40 I=1,KM1
              T=T+A(I,K)*WORK(I)
   40     CONTINUE
   45     EK=1D0
          IF (T.LT.0D0) EK=-1D0
          IF (A(K,K).EQ.0D0) GO TO 90
          WORK(K)=-(EK+T)/A(K,K)
   50 CONTINUE
      DO 60 KB=1,NM1
          K=N-KB
          T=0D0
          KP1=K+1
          DO 55 I=KP1,N
              T=T+A(I,K)*WORK(K)
   55     CONTINUE
          WORK(K)=T
          M=IPVT(K)
          IF (M.EQ.K) GO TO 60
          T=WORK(M)
          WORK(M)=WORK(K)
          WORK(K)=T
   60 CONTINUE
      YNORM=0D0
      DO 65 I=1,N
          YNORM=YNORM+DABS(WORK(I))
   65 CONTINUE
      CALL SOLVE (NDIM,N,A,WORK,IPVT)
      ZNORM=0D0
      DO 70 I=1,N
          ZNORM=ZNORM+DABS(WORK(I))
   70 CONTINUE
      COND=ANORM*ZNORM/YNORM
      IF (COND.LT.1d0) COND=1D0
      RETURN
   80 COND=1D0
      IF (A(1,1).NE.0D0) RETURN
   90 COND=1D52
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE SOLVE (NDIM,N,A,B,IPVT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),B(N)
      INTEGER IPVT(N)
      IF (N.EQ.1) GO TO 50
      NM1=N-1
      DO 20 K=1,NM1
          KP1=K+1
          M=IPVT(K)
          T=B(M)
          B(M)=B(K)
          B(K)=T
          DO 10 I=KP1,N
              B(I)=B(I)+A(I,K)*T
   10     CONTINUE
   20 CONTINUE
      DO 40 KB=1,NM1
          KM1=N-KB
          K=KM1+1
          B(K)=B(K)/A(K,K)
          T=-B(K)
          DO 30 I=1,KM1
              B(I)=B(I)+A(I,K)*T
   30     CONTINUE
   40 CONTINUE
   50 B(1)=B(1)/A(1,1)
      RETURN
      END
 
C********************************************************************
C                                                                   *
C   CALCULATION OF THE EXPANSION COEFFICIENTS FOR (I,Q,U,V) -       *
C   REPRESENTATION.                                                 *
C                                                                   *
C   INPUT PARAMETERS:                                               *
C                                                                   *
C      LAM - WAVELENGTH OF LIGHT                                    *
C      CSCA - SCATTERING CROSS SECTION                              *
C      TR AND TI - ELEMENTS OF THE T-MATRIX. TRANSFERRED THROUGH    *
C                  COMMON /CTM/                                     *
C      NMAX - DIMENSION OF T(M)-MATRICES                            *
C                                                                   *
C   OUTPUT INFORTMATION:                                            *
C                                                                   *
C      ALF1,...,ALF4,BET1,BET2 - EXPANSION COEFFICIENTS             *
C      LMAX - NUMBER OF COEFFICIENTS MINUS 1                        *
C                                                                   *
C********************************************************************
 
      SUBROUTINE GSP(NMAX,CSCA,LAM,ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-B,D-H,O-Z),COMPLEX*16 (C)
      REAL*8 LAM,SSIGN(900)
      REAL*8  CSCA,SSI(NPL),SSJ(NPN1),
     &        ALF1(NPL),ALF2(NPL),ALF3(NPL),
     &        ALF4(NPL),BET1(NPL),BET2(NPL),
     &        TR1(NPL1,NPN4),TR2(NPL1,NPN4),
     &        TI1(NPL1,NPN4),TI2(NPL1,NPN4),
     &        G1(NPL1,NPN6),G2(NPL1,NPN6),
     &        AR1(NPN4),AR2(NPN4),AI1(NPN4),AI2(NPN4),
     &        FR(NPN4,NPN4),FI(NPN4,NPN4),FF(NPN4,NPN4)
      REAL*4 B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4),
     &       B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4),
     &       D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4),
     &       D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4),
     &       D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4),
     &       PLUS1(NPN6*NPN4*NPN4*8)         
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CIM(NPN1)
 
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
      COMMON /CBESS/ B1R,B1I,B2R,B2I    
      COMMON /SS/ SSIGN
      EQUIVALENCE ( PLUS1(1),TR11(1,1,1) )
      EQUIVALENCE (D1(1,1,1),PLUS1(1)),        
     &            (D2(1,1,1),PLUS1(NPL1*NPN4*NPN4+1)),
     &            (D3(1,1,1),PLUS1(NPL1*NPN4*NPN4*2+1)),
     &            (D4(1,1,1),PLUS1(NPL1*NPN4*NPN4*3+1)), 
     &            (D5R(1,1,1),PLUS1(NPL1*NPN4*NPN4*4+1)) 
 
      CALL FACT
      CALL SIGNUM
      LMAX=2*NMAX
      L1MAX=LMAX+1
      CI=(0D0,1D0)
      CIM(1)=CI
      DO 2 I=2,NMAX
         CIM(I)=CIM(I-1)*CI
    2 CONTINUE
      SSI(1)=1D0
      DO 3 I=1,LMAX
         I1=I+1
         SI=DFLOAT(2*I+1)
         SSI(I1)=SI
         IF(I.LE.NMAX) SSJ(I)=DSQRT(SI)
    3 CONTINUE
      CI=-CI
      DO 5 I=1,NMAX
         SI=SSJ(I)
         CCI=CIM(I)
         DO 4 J=1,NMAX
            SJ=1D0/SSJ(J)
            CCJ=CIM(J)*SJ/CCI
            FR(J,I)=CCJ
            FI(J,I)=CCJ*CI
            FF(J,I)=SI*SJ
    4    CONTINUE
    5 CONTINUE
      NMAX1=NMAX+1
 
C *****  CALCULATION OF THE ARRAYS B1 AND B2  *****
 
      K1=1
      K2=0
      K3=0
      K4=1
      K5=1
      K6=2
 
C     PRINT 3300, B1,B2
 3300 FORMAT (' B1 AND B2')
      DO 100 N=1,NMAX
 
C *****  CALCULATION OF THE ARRAYS T1 AND T2  *****
 
 
         DO 10 NN=1,NMAX
            M1MAX=MIN0(N,NN)+1
            DO 6 M1=1,M1MAX
               M=M1-1
               L1=NPN6+M
               TT1=TR11(M1,N,NN)
               TT2=TR12(M1,N,NN)
               TT3=TR21(M1,N,NN)
               TT4=TR22(M1,N,NN)
               TT5=TI11(M1,N,NN)
               TT6=TI12(M1,N,NN)
               TT7=TI21(M1,N,NN)
               TT8=TI22(M1,N,NN)
               T1=TT1+TT2
               T2=TT3+TT4
               T3=TT5+TT6
               T4=TT7+TT8
               TR1(L1,NN)=T1+T2
               TR2(L1,NN)=T1-T2
               TI1(L1,NN)=T3+T4
               TI2(L1,NN)=T3-T4
               IF(M.EQ.0) GO TO 6
               L1=NPN6-M
               T1=TT1-TT2
               T2=TT3-TT4
               T3=TT5-TT6
               T4=TT7-TT8
               TR1(L1,NN)=T1-T2
               TR2(L1,NN)=T1+T2
               TI1(L1,NN)=T3-T4
               TI2(L1,NN)=T3+T4
    6       CONTINUE
   10    CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****
 
         NN1MAX=NMAX1+N
         DO 40 NN1=1,NN1MAX
            N1=NN1-1
 
C  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****
 
            CALL CCG(N,N1,NMAX,K1,K2,G1)
            NNMAX=MIN0(NMAX,N1+N)
            NNMIN=MAX0(1,IABS(N-N1))
            KN=N+NN1
            DO 15 NN=NNMIN,NNMAX
               NNN=NN+1
               SIG=SSIGN(KN+NN)
               M1MAX=MIN0(N,NN)+NPN6
               AAR1=0D0
               AAR2=0D0
               AAI1=0D0
               AAI2=0D0
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=G1(M1,NNN)
                  RR1=TR1(M1,NN)
                  RI1=TI1(M1,NN)
                  RR2=TR2(M1,NN)
                  RI2=TI2(M1,NN)
                  IF(M.EQ.0) GO TO 12
                  M2=NPN6-M
                  RR1=RR1+TR1(M2,NN)*SIG
                  RI1=RI1+TI1(M2,NN)*SIG
                  RR2=RR2+TR2(M2,NN)*SIG
                  RI2=RI2+TI2(M2,NN)*SIG
   12             AAR1=AAR1+SSS*RR1
                  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2
                  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=FR(NN,N)
               XI=FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI
               AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI
               AI2(NN)=AAR2*XI+AAI2*XR
   15       CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****
 
            CALL CCG(N,N1,NMAX,K3,K4,G2)
            M1=MAX0(-N1+1,-N)
            M2=MIN0(N1+1,N)
            M1MAX=M2+NPN6
            M1MIN=M1+NPN6
            DO 30 M1=M1MIN,M1MAX
               BBR1=0D0
               BBI1=0D0
               BBR2=0D0
               BBI2=0D0
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)
                  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)
                  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1
               B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2
               B2I(NN1,M1,N)=BBI2
   30       CONTINUE
   40    CONTINUE
  100 CONTINUE
 
C  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****
 
C  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****
 
C     PRINT 3301
 3301 FORMAT(' D1, D2, ...')
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)
            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1=0D0
               DD2=0D0
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)
                  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)
                  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)
                  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
  150          CONTINUE
               D1(M1,NN,N)=DD1
               D2(M1,NN,N)=DD2
  180       CONTINUE
            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0D0
               DD4=0D0
               DD5R=0D0
               DD5I=0D0
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)
                  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)
                  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)
                  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
  183          CONTINUE
               D3(M1,NN,N)=DD3
               D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R
               D5I(M1,NN,N)=DD5I
  186       CONTINUE
  190    CONTINUE
  200 CONTINUE
 
C  *****  END OF THE CALCULATION OF THE D-ARRAYS *****
 
C  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****
 
C     PRINT 3303
 3303 FORMAT (' G1, G2, ...')
 
      DK=LAM*LAM/(4D0*CSCA*DACOS(-1D0))
      DO 300 L1=1,L1MAX
         G1L=0D0
         G2L=0D0
         G3L=0D0
         G4L=0D0
         G5LR=0D0
         G5LI=0D0
         L=L1-1
         SL=SSI(L1)*DK
         DO 290 N=1,NMAX
            NNMIN=MAX0(1,IABS(N-L))
            NNMAX=MIN0(NMAX,N+L)
            IF(NNMAX.LT.NNMIN) GO TO 290
            CALL CCG(N,L,NMAX,K1,K2,G1)
            IF(L.GE.2) CALL CCG(N,L,NMAX,K5,K6,G2)
            NL=N+L
            DO 280  NN=NNMIN,NNMAX
               NNN=NN+1
               MMAX=MIN0(N,NN)
               M1MIN=NPN6-MMAX
               M1MAX=NPN6+MMAX
               SI=SSIGN(NL+NNN)
               DM1=0D0
               DM2=0D0
               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=G1(M1,NNN)
                  IF(M.LT.0) SSS1=G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
  270          CONTINUE
               FFN=FF(NN,N)
               SSS=G1(NPN6+1,NNN)*FFN
               G1L=G1L+SSS*DM1
               G2L=G2L+SSS*DM2*SI
               IF(L.LT.2) GO TO 280
               DM3=0D0
               DM4=0D0
               DM5R=0D0
               DM5I=0D0
               MMAX=MIN0(N,NN+2)
               MMIN=MAX0(-N,-NN+2)
               M1MAX=NPN6+MMAX
               M1MIN=NPN6+MMIN
               DO 275 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  SSS1=G2(NPN6-M,NNN)
                  DM3=DM3+SSS1*D3(M1,NN,N)
                  DM4=DM4+SSS1*D4(M1,NN,N)
                  DM5R=DM5R+SSS1*D5R(M1,NN,N)
                  DM5I=DM5I+SSS1*D5I(M1,NN,N)
  275          CONTINUE
               G5LR=G5LR-SSS*DM5R
               G5LI=G5LI-SSS*DM5I
               SSS=G2(NPN4,NNN)*FFN
               G3L=G3L+SSS*DM3
               G4L=G4L+SSS*DM4*SI
  280       CONTINUE
  290    CONTINUE
         G1L=G1L*SL
         G2L=G2L*SL
         G3L=G3L*SL
         G4L=G4L*SL
         G5LR=G5LR*SL
         G5LI=G5LI*SL
         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2D0
         BET2(L1)=G5LI*2D0
         LMAX=L
         IF(DABS(G1L).LT.1D-6) GO TO 500
  300 CONTINUE
  500 CONTINUE
      RETURN
      END
 
C****************************************************************
 
C   CALCULATION OF THE QUANTITIES F(N+1)=0.5*LN(N!)
C   0.LE.N.LE.899
 
      SUBROUTINE FACT
      REAL*8 F(900)
      COMMON /FAC/ F
      F(1)=0D0
      F(2)=0D0
      DO 2 I=3,900
         I1=I-1
         F(I)=F(I1)+0.5D0*DLOG(DFLOAT(I1))
    2 CONTINUE
      RETURN
      END
 
C************************************************************
 
C   CALCULATION OF THE ARRAY SSIGN(N+1)=SIGN(N)
C   0.LE.N.LE.899
 
      SUBROUTINE SIGNUM
      REAL*8 SSIGN(900)
      COMMON /SS/ SSIGN
      SSIGN(1)=1D0
      DO 2 N=2,899 
         SSIGN(N)=-SSIGN(N-1)
    2 CONTINUE
      RETURN
      END
 
C******************************************************************
C
C   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
C   (N,M:N1,M1/NN,MM)
C   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
C   MM=M*K1+K2
C
C   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
C                               N.LE.NMAX
C                               N.GE.1
C                               N1.GE.0
C                               N1.LE.N+NMAX
C   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
C                                       COEFFICIENTS
C                               /M/.LE.N
C                               /M1/=/M*(K1-1)+K2/.LE.N1
C                               NN.LE.MIN(N+N1,NMAX)
C                               NN.GE.MAX(/MM/,/N-N1/)
C   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N
 
 
      SUBROUTINE CCG(N,N1,NMAX,K1,K2,GG)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 GG(NPL1,NPN6),CD(0:NPN5),CU(0:NPN5)
      IF(NMAX.LE.NPN4.
     &   AND.0.LE.N1.
     &   AND.N1.LE.NMAX+N.
     &   AND.N.GE.1.
     &   AND.N.LE.NMAX) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCG')
    1 NNF=MIN0(N+N1,NMAX)
      MIN=NPN6-N
      MF=NPN6+N
      IF(K1.EQ.1.AND.K2.EQ.0) MIN=NPN6
      DO 100 MIND=MIN,MF
         M=MIND-NPN6
         MM=M*K1+K2
         M1=MM-M
         IF(IABS(M1).GT.N1) GO TO 90
         NNL=MAX0(IABS(MM),IABS(N-N1))
         IF(NNL.GT.NNF) GO TO 90
         NNU=N+N1
         NNM=(NNU+NNL)*0.5D0
         IF (NNU.EQ.NNL) NNM=NNL
         CALL CCGIN(N,N1,M,MM,C)
         CU(NNL)=C  
         IF (NNL.EQ.NNF) GO TO 50
         C2=0D0
         C1=C
         DO 7 NN=NNL+1,MIN0(NNM,NNF)
            A=DFLOAT((NN+MM)*(NN-MM)*(N1-N+NN))
            A=A*DFLOAT((N-N1+NN)*(N+N1-NN+1)*(N+N1+NN+1))
            A=DFLOAT(4*NN*NN)/A
            A=A*DFLOAT((2*NN+1)*(2*NN-1))
            A=DSQRT(A)
            B=0.5D0*DFLOAT(M-M1)
            D=0D0
            IF(NN.EQ.1) GO TO 5
            B=DFLOAT(2*NN*(NN-1))
            B=DFLOAT((2*M-MM)*NN*(NN-1)-MM*N*(N+1)+
     &               MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN-1)*(NN-1))
            D=D*DFLOAT((2*NN-3)*(2*NN-1))
            D=DFLOAT((NN-MM-1)*(NN+MM-1)*(N1-N+NN-1))/D
            D=D*DFLOAT((N-N1+NN-1)*(N+N1-NN+2)*(N+N1+NN))
            D=DSQRT(D)
    5       C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CU(NN)=C
    7    CONTINUE
         IF (NNF.LE.NNM) GO TO 50
         CALL DIRECT(N,M,N1,M1,NNU,MM,C)
         CD(NNU)=C
         IF (NNU.EQ.NNM+1) GO TO 50
         C2=0D0
         C1=C
         DO 12 NN=NNU-1,NNM+1,-1
            A=DFLOAT((NN-MM+1)*(NN+MM+1)*(N1-N+NN+1))
            A=A*DFLOAT((N-N1+NN+1)*(N+N1-NN)*(N+N1+NN+2))
            A=DFLOAT(4*(NN+1)*(NN+1))/A
            A=A*DFLOAT((2*NN+1)*(2*NN+3))
            A=DSQRT(A)
            B=DFLOAT(2*(NN+2)*(NN+1))
            B=DFLOAT((2*M-MM)*(NN+2)*(NN+1)-MM*N*(N+1)
     &               +MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN+2)*(NN+2))
            D=D*DFLOAT((2*NN+5)*(2*NN+3))
            D=DFLOAT((NN+MM+2)*(NN-MM+2)*(N1-N+NN+2))/D
            D=D*DFLOAT((N-N1+NN+2)*(N+N1-NN-1)*(N+N1+NN+3))
            D=DSQRT(D)
            C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CD(NN)=C
   12    CONTINUE
   50    DO 9 NN=NNL,NNF
            IF (NN.LE.NNM) GG(MIND,NN+1)=CU(NN)
            IF (NN.GT.NNM) GG(MIND,NN+1)=CD(NN)
c           WRITE (6,*) N,M,N1,M1,NN,MM,GG(MIND,NN+1)
    9    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
 
C*********************************************************************
 
      SUBROUTINE DIRECT (N,M,N1,M1,NN,MM,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(900)
      COMMON /FAC/ F
      C=F(2*N+1)+F(2*N1+1)+F(N+N1+M+M1+1)+F(N+N1-M-M1+1)    
      C=C-F(2*(N+N1)+1)-F(N+M+1)-F(N-M+1)-F(N1+M1+1)-F(N1-M1+1)
      C=DEXP(C)
      RETURN
      END
 
C*********************************************************************
C
C   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
C   G=(N,M:N1,MM-M/NN,MM)
C   FOR GIVEN N,N1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
C                               /M/.LE.N
C                               /MM-M/.LE.N1
C                               /MM/.LE.N+N1
 
      SUBROUTINE CCGIN(N,N1,M,MM,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(900),SSIGN(900)
      COMMON /SS/ SSIGN
      COMMON /FAC/ F
      M1=MM-M
      IF(N.GE.IABS(M).
     &   AND.N1.GE.IABS(M1).
     &   AND.IABS(MM).LE.(N+N1)) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCGIN')
    1 IF (IABS(MM).GT.IABS(N-N1)) GO TO 100
      L1=N
      L2=N1
      L3=M
      IF(N1.LE.N) GO TO 50
      K=N
      N=N1
      N1=K
      K=M
      M=M1
      M1=K
   50 N2=N*2
      M2=M*2
      N12=N1*2
      M12=M1*2
      G=SSIGN(N1+M1+1)
     & *DEXP(F(N+M+1)+F(N-M+1)+F(N12+1)+F(N2-N12+2)-F(N2+2)
     &       -F(N1+M1+1)-F(N1-M1+1)-F(N-N1+MM+1)-F(N-N1-MM+1))
      N=L1
      N1=L2
      M=L3
      RETURN
  100 A=1D0
      L1=M
      L2=MM
      IF(MM.GE.0) GO TO 150
      MM=-MM
      M=-M
      M1=-M1
      A=SSIGN(MM+N+N1+1)
  150 G=A*SSIGN(N+M+1)
     &   *DEXP(F(2*MM+2)+F(N+N1-MM+1)+F(N+M+1)+F(N1+M1+1)
     &        -F(N+N1+MM+2)-F(N-N1+MM+1)-F(-N+N1+MM+1)-F(N-M+1)
     &        -F(N1-M1+1))
      M=L1
      MM=L2
      RETURN
      END
 
C*****************************************************************
 
      SUBROUTINE SAREA (D,RAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (D.GE.1) GO TO 10
      E=DSQRT(1D0-D*D)
      R=0.5D0*(D**(2D0/3D0) + D**(-1D0/3D0)*DASIN(E)/E)
      R=DSQRT(R)
      RAT=1D0/R
      RETURN
   10 E=DSQRT(1D0-1D0/(D*D))
      R=0.25D0*(2D0*D**(2D0/3D0) + D**(-4D0/3D0)*DLOG((1D0+E)/(1D0-E))
     &   /E)
      R=DSQRT(R)
      RAT=1D0/R
      RETURN
      END
 
c****************************************************************
 
      SUBROUTINE SURFCH (N,E,RAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(60),W(60)
      DN=DFLOAT(N)
      E2=E*E
      EN=E*DN
      NG=60
      CALL GAUSS (NG,0,0,X,W)
      S=0D0
      V=0D0
      DO 10 I=1,NG
         XI=X(I)
         DX=DACOS(XI)
         DXN=DN*DX
         DS=DSIN(DX)
         DSN=DSIN(DXN)
         DCN=DCOS(DXN)
         A=1D0+E*DCN
         A2=A*A
         ENS=EN*DSN
         S=S+W(I)*A*DSQRT(A2+ENS*ENS)
         V=V+W(I)*(DS*A+XI*ENS)*DS*A2
   10 CONTINUE
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0/4D0)**(1D0/3D0)
      RAT=RV/RS
      RETURN
      END
 
C********************************************************************
 
      SUBROUTINE SAREAC (EPS,RAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      RAT=(1.5D0/EPS)**(1D0/3D0)
      RAT=RAT/DSQRT( (EPS+2D0)/(2D0*EPS) )
      RETURN
      END
 
C********************************************************************
 
C  COMPUTATION OF R1 AND R2 FOR A POWER LAW SIZE DISTRIBUTION WITH
C  EFFECTIVE RADIUS A AND EFFECTIVE VARIANCE B
 
      SUBROUTINE POWER (A,B,R1,R2)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F
      COMMON AA,BB
      AA=A
      BB=B
      AX=1D-5
      BX=A-1D-5
      R1=ZEROIN (AX,BX,F,0D0)
      R2=(1D0+B)*2D0*A-R1
      RETURN
      END
 
C***********************************************************************
 
      DOUBLE PRECISION FUNCTION ZEROIN (AX,BX,F,TOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      EPS=1D0
   10 EPS=0.5D0*EPS
      TOL1=1D0+EPS
      IF (TOL1.GT.1D0) GO TO 10
   15 A=AX
      B=BX
      FA=F(A)
      FB=F(B)
   20 C=A
      FC=FA
      D=B-A
      E=D
   30 IF (DABS(FC).GE.DABS(FB)) GO TO 40
   35 A=B
      B=C
      C=A
      FA=FB
      FB=FC
      FC=FA
   40 TOL1=2D0*EPS*DABS(B)+0.5D0*TOL
      XM=0.5D0*(C-B)
      IF (DABS(XM).LE.TOL1) GO TO 90
   44 IF (FB.EQ.0D0) GO TO 90
   45 IF (DABS(E).LT.TOL1) GO TO 70
   46 IF (DABS(FA).LE.DABS(FB)) GO TO 70
   47 IF (A.NE.C) GO TO 50
   48 S=FB/FA
      P=2D0*XM*S
      Q=1D0-S
      GO TO 60
   50 Q=FA/FC
      R=FB/FC
      S=FB/FA
      P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))
      Q=(Q-1D0)*(R-1D0)*(S-1D0)
   60 IF (P.GT.0D0) Q=-Q
      P=DABS(P)
      IF ((2D0*P).GE.(3D0*XM*Q-DABS(TOL1*Q))) GO TO 70
   64 IF (P.GE.DABS(0.5D0*E*Q)) GO TO 70
   65 E=D
      D=P/Q
      GO TO 80
   70 D=XM
      E=D
   80 A=B
      FA=FB
      IF (DABS(D).GT.TOL1) B=B+D
      IF (DABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)
      FB=F(B)
      IF ((FB*(FC/DABS(FC))).GT.0D0) GO TO 20
   85 GO TO 30
   90 ZEROIN=B
      RETURN
      END
 
C***********************************************************************
 
      DOUBLE PRECISION FUNCTION F(R1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON A,B
      R2=(1D0+B)*2D0*A-R1
      F=(R2-R1)/DLOG(R2/R1)-A
      RETURN
      END
 
C**********************************************************************
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
C    N - NUMBER OF POINTS                                             *
C    Z - DIVISION POINTS                                              *
C    W - WEIGHTS                                                      *
C**********************************************************************
 
      SUBROUTINE GAUSS ( N,IND1,IND2,Z,W )
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      DATA A,B,C /1D0,2D0,3D0/
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.check*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END
 
C****************************************************************
 
      SUBROUTINE DISTRB (NNK,YY,WY,NDISTR,AA,BB,GAM,R1,R2,QUIET,REFF,                 
     &                   VEFF,PI)                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      INTEGER*8 QUIET
      REAL*8 YY(NNK),WY(NNK)                                                    
      IF (NDISTR.EQ.2) GO TO 100                                                
      IF (NDISTR.EQ.3) GO TO 200                                                
      IF (NDISTR.EQ.4) GO TO 300                                                
      IF (NDISTR.EQ.5) GO TO 360
      IF (QUIET.EQ.0) PRINT 1001,AA,BB,GAM                                                      
 1001 FORMAT('MODIFIED GAMMA DISTRIBUTION, alpha=',F6.4,'  r_c=',               
     &  F6.4,'  gamma=',F6.4)                                                   
      A2=AA/GAM                                                                 
      DB=1D0/BB
      DO 50 I=1,NNK                                                             
         X=YY(I)                                                             
         Y=X**AA                                                                
         X=X*DB
         Y=Y*DEXP(-A2*(X**GAM))                                                 
         WY(I)=WY(I)*Y                                                       
   50 CONTINUE                                                                  
      GO TO 400                                                                 
  100 IF (QUIET.EQ.0) PRINT 1002,AA,BB                                                          
 1002 FORMAT('LOG-NORMAL DISTRIBUTION, r_g=',F8.4,                
     &       '  [ln(sigma_g)]**2=', F6.4)          
      DA=1D0/AA                                                                 
      DO 150 I=1,NNK                                                            
         X=YY(I)                                                                
         Y=DLOG(X*DA)                                                          
         Y=DEXP(-Y*Y*0.5D0/BB)/X                                             
         WY(I)=WY(I)*Y                                                          
  150 CONTINUE                                                                  
      GO TO 400                                                                 
  200 IF (QUIET.EQ.0) PRINT 1003                                                                
 1003 FORMAT('POWER LAW DISTRIBUTION OF HANSEN & TRAVIS 1974')                 
      DO 250 I=1,NNK                                                            
         X=YY(I)                                                                
         WY(I)=WY(I)/(X*X*X)                                                 
  250 CONTINUE                                                                  
      GO TO 400                                                                 
  300 IF (QUIET.EQ.0) PRINT 1004,AA,BB                                                          
 1004 FORMAT ('GAMMA DISTRIBUTION,  a=',F6.3,'  b=',F6.4)
      B2=(1D0-3D0*BB)/BB                                                        
      DAB=1D0/(AA*BB)                                                          
      DO 350 I=1,NNK                                                            
         X=YY(I)                                                                
         X=(X**B2)*DEXP(-X*DAB)                                                 
         WY(I)=WY(I)*X                                                       
  350 CONTINUE                                                                  
      GO TO 400                                                                 
  360 IF (QUIET.EQ.0) PRINT 1005,BB
 1005 FORMAT ('MODIFIED POWER LAW DISTRIBUTION,  alpha=',D10.4)
      DO 370 I=1,NNK
         X=YY(I)
         IF (X.LE.R1) WY(I)=WY(I)
         IF (X.GT.R1) WY(I)=WY(I)*(X/R1)**BB
  370 CONTINUE
  400 CONTINUE                                                                  
      SUM=0D0
      DO 450 I=1,NNK
         SUM=SUM+WY(I)
  450 CONTINUE
      SUM=1D0/SUM
      DO 500 I=1,NNK
         WY(I)=WY(I)*SUM
  500 CONTINUE
      G=0D0
      DO 550 I=1,NNK
         X=YY(I)
         G=G+X*X*WY(I)
  550 CONTINUE
      REFF=0D0
      DO 600 I=1,NNK
         X=YY(I)
         REFF=REFF+X*X*X*WY(I)
  600 CONTINUE
      REFF=REFF/G
      VEFF=0D0
      DO 650 I=1,NNK
         X=YY(I)
         XI=X-REFF
         VEFF=VEFF+XI*XI*X*X*WY(I)
  650 CONTINUE
      VEFF=VEFF/(G*REFF*REFF)
      RETURN                                                                    
      END                                                                       
 
C*************************************************************
 
      SUBROUTINE HOVENR(L1,A1,A2,A3,A4,B1,B2,QUIET)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 QUIET
      REAL*8 A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
      DO 100 L=1,L1
         KONTR=1
         LL=L-1
         DL=DFLOAT(LL)*2D0+1D0
         DDL=DL*0.48D0
         AA1=A1(L)
         AA2=A2(L)
         AA3=A3(L)
         AA4=A4(L)
         BB1=B1(L)
         BB2=B2(L)
         IF(LL.GE.1.AND.DABS(AA1).GE.DL) KONTR=2
         IF(DABS(AA2).GE.DL) KONTR=2
         IF(DABS(AA3).GE.DL) KONTR=2
         IF(DABS(AA4).GE.DL) KONTR=2
         IF(DABS(BB1).GE.DDL) KONTR=2
         IF(DABS(BB2).GE.DDL) KONTR=2
         IF(KONTR.EQ.2) PRINT 3000,LL
         C=-0.1D0
         DO 50 I=1,11
            C=C+0.1D0
            CC=C*C
            C1=CC*BB2*BB2
            C2=C*AA4
            C3=C*AA3
            IF((DL-C*AA1)*(DL-C*AA2)-CC*BB1*BB1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL-C3)+C1.LE.-1D-4) KONTR=2
            IF((DL+C2)*(DL-C3)-C1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL+C3)-C1.LE.-1D-4) KONTR=2
            IF(KONTR.EQ.2) PRINT 4000,LL,C
   50    CONTINUE
  100 CONTINUE
      IF((KONTR.EQ.1).AND.(QUIET.EQ.0)) PRINT 2000
 2000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS SATISFIED')
 3000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3)
 4000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3,
     & '   A=',D9.2)
      RETURN
      END
 
C****************************************************************
 
C    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
C    COEFFICIENTS
 
C    A1,...,B2 - EXPANSION COEFFICIENTS
C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
C    N - NUMBER OF SCATTERING ANGLES
C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
C        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
      SUBROUTINE MATR(A1,A2,A3,A4,B1,B2,LMAX,NPNA,QUIET,F11,F22,F33,F44,
     &     F12,F34)
      INCLUDE 'tmd.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
C     The 201's in the following line were originally NPL, they have //CPD
C     been changed to get around an f2py bug.                        //CPD
      REAL*8 A1(201),A2(201),A3(201),A4(201),B1(201),B2(201),F11(NPNA),
     &     F22(NPNA),F33(NPNA),F44(NPNA),F12(NPNA),F34(NPNA)
      INTEGER QUIET
      N=NPNA
      DN=1D0/DFLOAT(N-1)
      DA=DACOS(-1D0)*DN
      DB=180D0*DN
      L1MAX=LMAX+1
      TB=-DB
      TAA=-DA
      IF (QUIET.EQ.0) THEN
         PRINT 1000
 1000    FORMAT(' ')
         PRINT 1001
 1001    FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',
     &        6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
         DO 10 L1=1,L1MAX
            L=L1-1
            PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
 10      CONTINUE
 1002    FORMAT(' ',I3,6F12.5)
         PRINT 1000
         PRINT 1003
 1003    FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',
     &        8X,'F44',8X,'F12',8X,'F34')
      ENDIF
      D6=DSQRT(6D0)*0.25D0
      DO 500 I1=1,N
         TAA=TAA+DA
         TB=TB+DB
         U=DCOS(TAA)
         F11(I1)=0D0
         F2=0D0
         F3=0D0
         F44(I1)=0D0
         F12(I1)=0D0
         F34(I1)=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11(I1)=F11(I1)+A1(L1)*PP1
            F44(I1)=F44(I1)+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12(I1)=F12(I1)+B1(L1)*PP4
            F34(I1)=F34(I1)+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1*(L*L-4))
            PL4=1D0/DFLOAT(L*(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400    CONTINUE
         F22(I1)=(F2+F3)*0.5D0
         F33(I1)=(F2-F3)*0.5D0
C        F22=F22/F11
C        F33=F33/F11
C        F44=F44/F11
C        F12=-F12/F11
C        F34=F34/F11
         IF (QUIET.EQ.0) PRINT 1004,TB,F11(I1),F22(I1),F33(I1),F44(I1),
     &        F12(I1),F34(I1)
  500 CONTINUE
      IF (QUIET.EQ.0) THEN
         PRINT 1000 
 1004    FORMAT(' ',F6.2,6F11.4)
      ENDIF
      RETURN
      END
