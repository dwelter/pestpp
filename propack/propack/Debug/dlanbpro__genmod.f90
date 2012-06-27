        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 27 16:46:28 2012
        MODULE DLANBPRO__genmod
          INTERFACE 
            SUBROUTINE DLANBPRO(M,N,K0,K,APROD,U,LDU,V,LDV,B,LDB,RNORM, &
     &DOPTION,IOPTION,WORK,IWORK,DPARM,IPARM,IERR)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDU
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K0
              INTEGER(KIND=4) :: K
              EXTERNAL APROD
              REAL(KIND=8) :: U(LDU,*)
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: RNORM
              REAL(KIND=8) :: DOPTION(*)
              INTEGER(KIND=4) :: IOPTION(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: IWORK(*)
              REAL(KIND=8) :: DPARM(*)
              INTEGER(KIND=4) :: IPARM(*)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE DLANBPRO
          END INTERFACE 
        END MODULE DLANBPRO__genmod
