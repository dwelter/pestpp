        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 30 08:48:55 2012
        MODULE DCGS__genmod
          INTERFACE 
            SUBROUTINE DCGS(N,K,V,LDV,VNEW,INDEX,WORK)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: VNEW(*)
              INTEGER(KIND=4) :: INDEX(*)
              REAL(KIND=8) :: WORK(*)
            END SUBROUTINE DCGS
          END INTERFACE 
        END MODULE DCGS__genmod
