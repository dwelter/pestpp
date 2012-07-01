        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 30 08:48:55 2012
        MODULE DLANSVD__genmod
          INTERFACE 
            SUBROUTINE DLANSVD(JOBU,JOBV,M,N,K,KMAX,APROD,U,LDU,SIGMA,  &
     &BND,V,LDV,TOLIN,WORK,LWORK,IWORK,LIWORK,DOPTION,IOPTION,INFO,DPARM&
     &,IPARM)
              INTEGER(KIND=4) :: LIWORK
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDU
              CHARACTER(LEN=1) :: JOBU
              CHARACTER(LEN=1) :: JOBV
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: KMAX
              EXTERNAL APROD
              REAL(KIND=8) :: U(LDU,*)
              REAL(KIND=8) :: SIGMA(*)
              REAL(KIND=8) :: BND(*)
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: TOLIN
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: IWORK(LIWORK)
              REAL(KIND=8) :: DOPTION(*)
              INTEGER(KIND=4) :: IOPTION(*)
              INTEGER(KIND=4) :: INFO
              REAL(KIND=8) :: DPARM(*)
              INTEGER(KIND=4) :: IPARM(*)
            END SUBROUTINE DLANSVD
          END INTERFACE 
        END MODULE DLANSVD__genmod
