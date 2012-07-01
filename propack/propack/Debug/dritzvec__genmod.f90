        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 30 08:48:55 2012
        MODULE DRITZVEC__genmod
          INTERFACE 
            SUBROUTINE DRITZVEC(WHICH,JOBU,JOBV,M,N,K,DIM,D,E,S,U,LDU,V,&
     &LDV,WORK,IN_LWRK,IWORK)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDU
              CHARACTER(LEN=1) :: WHICH
              CHARACTER(LEN=1) :: JOBU
              CHARACTER(LEN=1) :: JOBV
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: DIM
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: U(LDU,*)
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: IN_LWRK
              INTEGER(KIND=4) :: IWORK(*)
            END SUBROUTINE DRITZVEC
          END INTERFACE 
        END MODULE DRITZVEC__genmod
