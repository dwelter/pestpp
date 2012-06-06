        !COMPILER-GENERATED INTERFACE MODULE: Wed Apr 25 15:06:19 2012
        MODULE DBDQR__genmod
          INTERFACE 
            SUBROUTINE DBDQR(IGNORELAST,JOBQ,N,D,E,C1,C2,QT,LDQ)
              INTEGER(KIND=4) :: LDQ
              LOGICAL(KIND=4) :: IGNORELAST
              CHARACTER(LEN=1) :: JOBQ
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: C1
              REAL(KIND=8) :: C2
              REAL(KIND=8) :: QT(LDQ,*)
            END SUBROUTINE DBDQR
          END INTERFACE 
        END MODULE DBDQR__genmod
