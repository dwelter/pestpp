        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 27 15:37:04 2012
        MODULE DGETU0__genmod
          INTERFACE 
            SUBROUTINE DGETU0(TRANSA,M,N,J,NTRY,U0,U0NORM,U,LDU,APROD,  &
     &DPARM,IPARM,IERR,ICGS,ANORMEST,WORK)
              CHARACTER(LEN=1) :: TRANSA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: NTRY
              REAL(KIND=8) :: U0(*)
              REAL(KIND=8) :: U0NORM
              REAL(KIND=8) :: U(*)
              INTEGER(KIND=4) :: LDU
              EXTERNAL APROD
              REAL(KIND=8) :: DPARM(*)
              INTEGER(KIND=4) :: IPARM(*)
              INTEGER(KIND=4) :: IERR
              INTEGER(KIND=4) :: ICGS
              REAL(KIND=8) :: ANORMEST
              REAL(KIND=8) :: WORK(*)
            END SUBROUTINE DGETU0
          END INTERFACE 
        END MODULE DGETU0__genmod
