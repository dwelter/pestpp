        !COMPILER-GENERATED INTERFACE MODULE: Mon Jul 02 20:34:13 2012
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
