        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 27 16:46:28 2012
        MODULE DGEMM_OVWR__genmod
          INTERFACE 
            SUBROUTINE DGEMM_OVWR(TRANSA,M,N,K,ALPHA,A,LDA,BETA,B,LDB,  &
     &DWORK,LDWORK)
              INTEGER(KIND=4) :: LDWORK
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANSA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: DWORK(LDWORK)
            END SUBROUTINE DGEMM_OVWR
          END INTERFACE 
        END MODULE DGEMM_OVWR__genmod
