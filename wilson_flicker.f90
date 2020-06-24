!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Wilson :           
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION E1,E2,H1,H2,I1,I2,X,Y,X1,Z1,Z2,PFBP,FA,epsi
       E1=U(1)
       H1=U(2)
       I1=U(3) 
       E2=U(4)
       H2=U(5)
       I2=U(6)
       X=U(7)
       Y=U(8)
       PFBP  = par(8)
       FA  = par(9)
       epsi  = par(10)

       X1=1/(1.0d0+dexp(-10*(X)))

       if ((PAR(7)*X1-PAR(2)*I2)>0.0d0) then
       Z1=((PAR(7)*X1-PAR(2)*I2)**2.0d0)
       else
       Z1=0
       end if
       if ((PAR(7)*X1-PAR(2)*I1)>0.0d0) then
       Z2=((PAR(7)*X1-PAR(2)*I1)**2.0d0)
       else
       Z2=0
       end if

       F(1)=( 1.0d0/PAR(5) ) * (-E1+((100.0d0*Z1*FA  )/((((10.0d0+H1+epsi)**2.0d0 )+(Z1*FA )))))
       F(2)=(1.0d0/PAR(3))*(-H1+PAR(1)*E1)     
       F(3)=(1.0d0/PAR(4))*(-I1+E1)
       F(4)=( 1.0d0/PAR(5) ) * (-E2+((100.0d0*Z2*FA  )/((((10.0d0+H2)**2.0d0 )+(Z2*FA )))))
       F(5)=(1.0d0/PAR(3))*(-H2+PAR(1)*E2)
       F(6)=(1.0d0/PAR(4))*(-I2+E2)

       F(7)=PFBP*(X+PAR(6)*Y-X*((X**2)+(Y**2)))-X*(1-PFBP)
       F(8)=PFBP*(-PAR(6)*X+Y-Y*((X**2)+(Y**2)))-Y*(1-PFBP)
      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values
!      T      :	  Not used here

      IMPLICIT NONE
      INTEGER NDIM
      DOUBLE PRECISION U(NDIM), PAR(*), T

! Initialize the equation parameters      
       PAR(1)=0
       PAR(2)=1.5
       PAR(3)=900.0d0
       PAR(4)=11.0d0
       PAR(5)=20.0d0
       PAR(6)=18*(2.0*3.14)*0.001
       PAR(7)=10.0d0
       PAR(8)=0.0d0                ! PFBP
       par(9)  = 0.0               ! A forcing amplitude
       par(10)  = 0.0              ! epsi
! Initialize the solution
       U(1)=0.
       U(2)=0.
       U(3)=0.
       U(4)=0.
       U(5)=0.
       U(6)=0.
       U(7)=0.
       U(8)=0.

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      

