! Part of Zinc FE package. Author: John Blackburn

module common

! ----------------------------------------------------------------------
! Mesh geometry
! ----------------------------------------------------------------------

  integer, allocatable :: ireg(:,:,:)
  integer, allocatable :: iregup(:,:,:)
  double precision, allocatable :: rnode(:,:,:,:)
  
! ----------------------------------------------------------------------
! FE matrix/vectors Qval, iQ, jQ, RR, kfac and vector solution
! lenQ is length of Q (stored sparse). 
! Left, right and data are for getQ3 and update3 using binary tree
! ----------------------------------------------------------------------
  
  double precision, allocatable :: Qval(:),RR(:),vec(:),vecred(:),Qval0(:),vecp(:)
  integer, allocatable :: iQ(:),jQ(:),irowst(:),irowed(:)
  integer, allocatable :: left(:),right(:),data(:)
  integer, allocatable :: veclookup(:)
  integer lenQ,leniQ,lenjQ

  double precision, allocatable :: kfac(:),uinit(:,:,:,:)
  integer, allocatable :: iunk(:),iunkred(:)

  double precision, allocatable :: Nres(:),dvec(:)

! ----------------------------------------------------------------------
! Various data from .zin and .mtf file
! includes no. nodes in each direction, number of regions and variables
! ----------------------------------------------------------------------
  
  integer imax,jmax,kmax,nvar,nvar_rst,nregel,nregnd,ndof,ndofred,nel,nnod,istep
  integer ng_xi,ng_eta,ng_mu,NDIM
  double precision omega,scale,tstep,tol
  
  integer key_db,key_sim,key_Q,nstride,nstrideNL,nstep,itmax,key_export
  integer idbmax,idbmin,jdbmax,jdbmin,kdbmax,kdbmin         ! debug only

  character(6) zou_format,rst_format,residual,nodecheck,newton,removeFixed

  logical centre_eval,jac_by_element

! ----------------------------------------------------------------------
! The program's name. Numerical constants, variable name strings
! ----------------------------------------------------------------------

  character(len=*), parameter :: name='Zinc 3.6'
  character(20) :: date,time
  character(1000) :: fname,restart
  double precision, parameter :: zero=0d0

  character(30), allocatable :: us(:),dus(:,:),rs(:)
  character(1) :: scoord(3)

  integer, parameter :: EXPRLEN=100   ! Max length of sCC etc strings

! ----------------------------------------------------------------------
!     set gauss quadrature eval points and weights depending on number of
!     points to be evaluated (ng)
! ----------------------------------------------------------------------

  double precision gauss(5,5),wt(5,5)

!!$  data ((gauss(icom,jcom),jcom=1,3),icom=1,3)/ &
!!$       0,0,0, &
!!$       -0.577350269d0,0.577350269d0,0, & ! sqrt(1/3)
!!$       0,-0.774596669d0,0.774596669d0/   ! sqrt(3/5)
!!$
!!$  data ((wt(icom,jcom),jcom=1,3),icom=1,3)/ &
!!$       2,0,0, &
!!$       1,1,0, &
!!$       0.88888888888d0,0.55555555555d0,0.55555555555d0/ ! 8/9,5/9

end module common
