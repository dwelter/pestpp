        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 30 08:48:55 2012
        MODULE DREORTH__genmod
          INTERFACE 
            SUBROUTINE DREORTH(N,K,V,LDV,VNEW,NORMVNEW,INDEX,ALPHA,WORK,&
     &IFLAG)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: VNEW(*)
              REAL(KIND=8) :: NORMVNEW
              INTEGER(KIND=4) :: INDEX(*)
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE DREORTH
          END INTERFACE 
        END MODULE DREORTH__genmod
