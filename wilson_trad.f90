!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Wilson_trad :           
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
      DOUBLE PRECISION E1,E2,H1,H2,I1,I2,X1,X2,PFBP

       E1=U(1)
       H1=U(2)
       I1=U(3) 
       E2=U(4)
       H2=U(5)
       I2=U(6)
       PFBP  = par(7)

       if ((PAR(6)-PAR(2)*I2)>0.0d0) then
       X1=((PAR(6)-PAR(2)*I2)**2.0d0)
       else
       X1=0
       end if

       if ((PAR(6)-PAR(2)*I1)>0.0d0) then
       X2=((PAR(6)-PAR(2)*I1)**2.0d0)
       else
       X2=0
       end if

       F(1)=PFBP*((1.0d0/PAR(5))*(-E1+((100.0d0*X1)/(((10.0d0+H1)**2.0d0)+(X1)))))-E1*(1-PFBP)

       F(2)=PFBP*((1.0d0/PAR(3))*(-H1+PAR(1)*E1))-H1*(1-PFBP)  

       F(3)=PFBP*((1.0d0/PAR(4))*(-I1+E1))-I1*(1-PFBP)

       F(4)=PFBP*((1.0d0/PAR(5))*(-E2+((100.0d0*X2)/(((10.0d0+H2)**2.0d0)+(X2)))))-E2*(1-PFBP)

       F(5)=PFBP*((1.0d0/PAR(3))*(-H2+PAR(1)*E2))-H2*(1-PFBP)

       F(6)=PFBP*((1.0d0/PAR(4))*(-I2+E2))-I2*(1-PFBP)

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
       PAR(1)=0.0
       PAR(2)=1.5
       PAR(3)=900.0d0
       PAR(4)=11.0d0
       PAR(5)=20.0d0
       PAR(6)=10.0d0
       PAR(7)=0.0d0

! Initialize the solution
#       U(1)=0.
#       U(2)=0.
#       U(3)=0.
#       U(4)=0.
#       U(5)=0.
#       U(6)=0.
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
