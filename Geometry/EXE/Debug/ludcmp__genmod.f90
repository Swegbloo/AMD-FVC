        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 05 11:34:26 2024
        MODULE LUDCMP__genmod
          INTERFACE 
            SUBROUTINE LUDCMP(A,N,NP,INDX)
              INTEGER(KIND=4) :: NP
              COMPLEX(KIND=4) :: A(NP,NP)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: INDX(4000)
            END SUBROUTINE LUDCMP
          END INTERFACE 
        END MODULE LUDCMP__genmod
