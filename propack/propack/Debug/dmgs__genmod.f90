        !COMPILER-GENERATED INTERFACE MODULE: Mon Jul 02 20:34:13 2012
        MODULE DMGS__genmod
          INTERFACE 
            SUBROUTINE DMGS(N,K,V,LDV,VNEW,INDEX)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: VNEW(*)
              INTEGER(KIND=4) :: INDEX(*)
            END SUBROUTINE DMGS
          END INTERFACE 
        END MODULE DMGS__genmod
