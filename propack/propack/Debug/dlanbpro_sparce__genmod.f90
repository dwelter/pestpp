        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 27 15:37:04 2012
        MODULE DLANBPRO_SPARCE__genmod
          INTERFACE 
            SUBROUTINE DLANBPRO_SPARCE(M,N,K0,K,U,LDU,V,LDV,B,LDB,RNORM,&
     &DOPTION,IOPTION,WORK,IWORK,DPARM,IPARM,IERR)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDU
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K0
              INTEGER(KIND=4) :: K
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
            END SUBROUTINE DLANBPRO_SPARCE
          END INTERFACE 
        END MODULE DLANBPRO_SPARCE__genmod
