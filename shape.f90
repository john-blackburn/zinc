! Part of Zinc FE package. Author: John Blackburn

module shape
  implicit none

  double precision, private, save :: xic(8),etac(8),muc(8),xic2d(4),etac2d(4)
  
  data  xic/-1, 1, 1,-1,-1, 1, 1,-1/
  data etac/-1,-1, 1, 1,-1,-1, 1, 1/
  data  muc/-1,-1,-1,-1, 1, 1, 1, 1/

  data xic2d /-1,+1,+1,-1/
  data etac2d/-1,-1,+1,+1/
  
contains

! ######################################################################

function dNdxi(l,xi,eta,mu)
  use common, only : NDIM
  integer l
  double precision dNdxi,xi,eta,mu

  if (NDIM==3) then
     dNdxi=0.125*xic(l)*(1+eta*etac(l))*(1+mu*muc(l))
  else if (NDIM==2) then
     dNdxi=0.250*xic(l)*(1+eta*etac(l))
  else
     dNdxi=0.500*xic(l)
  endif

end function dNdxi

! ######################################################################

function dNdeta(l,xi,eta,mu)  ! 2D or 3D
  use common, only : NDIM
  integer l  
  double precision dNdeta,xi,eta,mu

  if (NDIM==3) then
     dNdeta=0.125*etac(l)*(1+xi*xic(l))*(1+mu*muc(l))
  else
     dNdeta=0.125*etac(l)*(1+xi*xic(l))
  endif

end function dNdeta

! ######################################################################

function dNdmu(l,xi,eta,mu)  ! 3D only
  integer l
  double precision dNdmu,xi,eta,mu
  
  dNdmu=0.125*muc(l)*(1+xi*xic(l))*(1+eta*etac(l))
  
end function dNdmu

! ######################################################################

function Nl(l,xi,eta,mu)
  use common, only : NDIM
  integer l
  double precision Nl,xi,eta,mu

  if (NDIM==3) then
     Nl=0.125*(1+xi*xic(l))*(1+eta*etac(l))*(1+mu*muc(l))
  else if (NDIM==2) then
     Nl=0.250*(1+xi*xic(l))*(1+eta*etac(l))
  else
     Nl=0.500*(1+xi*xic(l))
  endif

end function Nl

! ######################################################################

function N2d(l,xi,eta)
  double precision N2d,xi,eta
  integer l

  N2d=0.25*(1+xi*xic2d(l))*(1+eta*etac2d(l))
end function N2d

! ######################################################################

function dN2d_dxi(l,xi,eta)
  double precision dN2d_dxi,xi,eta
  integer l

  dN2d_dxi=0.25*xic2d(l)*(1+eta*etac2d(l))
end function dN2d_dxi

! ######################################################################

function dN2d_deta(l,xi,eta)
  double precision dN2d_deta,xi,eta
  integer l

  dN2d_deta=0.25*(1+xi*xic2d(l))*etac2d(l)
  
end function dN2d_deta

end module shape
