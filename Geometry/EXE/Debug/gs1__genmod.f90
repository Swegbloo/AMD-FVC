        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 05 11:34:26 2024
        MODULE GS1__genmod
          INTERFACE 
            SUBROUTINE GS1(NROOT,ANU,AK,DP,XX,YY,XP,CGJ,AN,AR,R,GF,GNF, &
     &GXF,GYF,GZF)
              INTEGER(KIND=4) :: NROOT
              REAL(KIND=8) :: ANU
              REAL(KIND=8) :: AK
              REAL(KIND=8) :: DP
              REAL(KIND=8) :: XX(4)
              REAL(KIND=8) :: YY(4)
              REAL(KIND=8) :: XP(3)
              REAL(KIND=8) :: CGJ(3)
              REAL(KIND=8) :: AN(3)
              REAL(KIND=8) :: AR
              REAL(KIND=8) :: R(3,3)
              COMPLEX(KIND=4) :: GF
              COMPLEX(KIND=4) :: GNF
              COMPLEX(KIND=4) :: GXF
              COMPLEX(KIND=4) :: GYF
              COMPLEX(KIND=4) :: GZF
            END SUBROUTINE GS1
          END INTERFACE 
        END MODULE GS1__genmod
