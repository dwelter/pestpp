        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 04 08:43:31 2012
        MODULE DGEMM_OVWR_LEFT__genmod
          INTERFACE 
            SUBROUTINE DGEMM_OVWR_LEFT(TRANSB,M,N,K,ALPHA,A,LDA,BETA,B, &
     &LDB,DWORK,LDWORK)
              INTEGER(KIND=4) :: LDWORK
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANSB
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: DWORK(LDWORK)
            END SUBROUTINE DGEMM_OVWR_LEFT
          END INTERFACE 
        END MODULE DGEMM_OVWR_LEFT__genmod
