function cfun(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax)

  use fullfuelmod

  implicit none
  character(*) token
  integer nvar,istep,imax,jmax,kmax
  double precision Cfun,x,y,z,ur(nvar),dur(nvar,3)

  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)  ! not used
  double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)            ! not used

  integer lentok,ios,iregion
  character(100) tokenstem
  character(1) cregion

  double precision vEle,p,wH,wO,wH2OCat,vMem,cH3O,cH2O,wH2O,wN
  double precision mass,rho,xH,xO,xH2OCat,xH2O,xN
  double precision mMem,rhoMem,lam, DH2O, DH2O_2, DH2O_3, DH2O_4
  double precision denom,D11Ano,D12Ano,D11Cat,D12Cat,D13Cat,D21Cat,D22Cat,D23Cat
  double precision dxH_dwH,dxO_dwO,dxO_dwH2O,dxH2O_dwO,dxH2O_dwH2O
  double precision c_vEle_vEle,c_p_p,c_wH_p, c_wO_p,c_wH2O_p,c_wH_wH
  double precision c_wO_wO, c_wO_wH2O,c_wH2O_wO,c_wH2O_wH2O
  double precision c_vMem_vMem,c_vMem_cH3O,c_cH2O_cH2O,c_cH3O_cH3O,c_cH3O_vMem

  if (nvar /= 8) then
     print *,"Wrong number of variables!"
     stop
  endif

  call setconst   ! Setup constants if needed

! ----------------------------------------------------------------------
! Split token as name_region, eg C_P_P_1 -> C_P_P and 1
! ----------------------------------------------------------------------

  lentok=len_trim(token)

  cregion=token(lentok:lentok) ! last character "1" "2" etc

  read (cregion,*,iostat=ios) iregion

  if (ios /= 0) then
     print *,'Garbage region specified, got: ',cregion
     stop
  endif

  if (iregion < 1.or.iregion > 5) then
     print *,'Region out of bounds, got:',iregion
     stop
  endif

  tokenstem=token(:lentok-2)

! ----------------------------------------------------------------------
! Unpack variables
! labels = vEle p wH wO wH2OCat vMem cH3O cH2O
! ----------------------------------------------------------------------

  vEle=ur(1)
  p=ur(2)

  wH=ur(3)

  wO=ur(4)
  wH2OCat=ur(5)

  vMem=ur(6)
  cH3O=ur(7)
  cH2O=ur(8)

  wH2O=1-wH       ! Anode side
  wN=1-wO-wH2OCat

  if (iregion == 1.or.iregion == 2) then
     mass=1/(wH/mH+wH2O/mH2O)
  else if (iregion == 4.or.iregion==5) then
     mass=1/(wO/mO+wH2OCat/mH2O+wN/mN)
  else
     mass=0
  endif

  rho=p*mass/(R*T)

  xH=     mass*wH     /mH
  xO=     mass*wO     /mO
  xH2OCat=mass*wH2OCat/mH2O

  xH2O=1-xH            ! on Anode side
  xN=1-xO-xH2OCat

! ----------------------------------------------------------------------
! Nernst Planck Diffusion coefficients
! ----------------------------------------------------------------------

  mMem=1.1       ! kg/mol. Molecular weight of membrane
  rhoMem=2000    ! kg/m^3. Density of membrane

  lam=cH2O*mMem/rhoMem
  DH2O=6e-8*lam*exp(-2436.0/T)

  DH2O_2=0.5*DH2O
  DH2O_3=DH2O
  DH2O_4=0.5*DH2O

! ----------------------------------------------------------------------
! Effective diffusion coefficients and derivative quantities
! ----------------------------------------------------------------------

  denom=xO/(dOH2O*dON)+xH2OCat/(dOH2O*dH2ON)+xN/(dON*dH2ON)

  D11Ano=wH2O**2/(xH*xH2O)*dHH2O
  D12Ano=-wH*(1-wH)/(xH*xH2O)*dHH2O

  D11Cat=((wH2OCat+wN)**2/(xO*dH2ON)+wH2OCat**2/(xH2OCat*dON)+wN**2/(xN*dOH2O))/denom
  D12Cat=(-wO*(wH2OCat+wN)/(xO*dH2ON)-wH2OCat*(wO+wN)/(xH2OCat*dON)+wN**2/(xN*dOH2O))/denom
  D13Cat=(-wO*(wH2OCat+wN)/(xO*dH2ON)+wH2OCat**2/(xH2OCat*dON)-wN*(wO+wH2OCat)/(xN*dOH2O))/denom

  D21Cat=D12Cat
  D22Cat=(wO**2/(xO*dH2ON)+(wO+wN)**2/(xH2OCat*dON)+wN**2/(xN*dOH2O))/denom
  D23Cat=(wO**2/(xO*dH2ON)-wH2OCat*(wO+wN)/(xH2OCat*dON)-wN*(wO+wH2OCat)/(xN*dOH2O))/denom

  dxH_dwH    =-mass**2*(1/mH  -1/mH2O)*wH/mH       +mass/mH
  dxO_dwO    =-mass**2*(1/mO  -1/mN)  *wO/mO       +mass/mO
  dxO_dwH2O  =-mass**2*(1/mH2O-1/mN)  *wO/mO
  dxH2O_dwO  =-mass**2*(1/mO  -1/mN)  *wH2OCat/mH2O
  dxH2O_dwH2O=-mass**2*(1/mH2O-1/mN)  *wH2OCat/mH2O+mass/mH2O

! ----------------------------------------------------------------------
! region-specific expressions
! ----------------------------------------------------------------------

  if (iregion==1) then
     c_vEle_vEle=sig1Ele
  else if (iregion==2) then
     c_vEle_vEle=sig2Ele
  else if (iregion==3) then
     c_vEle_vEle=0
  else if (iregion==4) then
     c_vEle_vEle=sig4Ele
  else if (iregion==5) then
     c_vEle_vEle=sig5Ele
  endif

  if (iregion==1) then
     c_p_p=rho*kap1/eta1
  else if (iregion==2) then
     c_p_p=rho*kap2/eta2
  else if (iregion==3) then
     c_p_p=0
  else if (iregion==4) then
     c_p_p=rho*kap4/eta4
  else if (iregion==5) then
     c_p_p=rho*kap5/eta5
  endif

  if (iregion==1.or.iregion==2) then          ! off-diag anode
     c_wH_p=rho*wH*(xH-wH)*(D11Ano-D12Ano)/p
  else
     c_wH_p=0
  endif

  if (iregion==4.or.iregion==5) then          ! off-diag cathode
     c_wO_p=  rho*wO  *((xO-wO)*(D11Cat-D13Cat)+(xH2OCat-wH2OCat)*(D12Cat-D13Cat))/p
     c_wH2O_p=rho*wH2OCat*((xO-wO)*(D21Cat-D23Cat)+(xH2OCat-wH2OCat)*(D22Cat-D23Cat))/p
  else
     c_wO_p=0
     c_wH2O_p=0
  endif

  if (iregion==1.or.iregion==2) then          ! Anode Maxwell-Stephan
     c_wH_wH=rho*wH*(D11Ano-D12Ano)*dxH_dwH
  else
     c_wH_wH=0
  endif

  if (iregion==4.or.iregion==5) then          ! Cathode M-S
     c_wO_wO  =  rho*wO  *((D11Cat-D13Cat)*dxO_dwO   + (D12Cat-D13Cat)*dxH2O_dwO)
     c_wO_wH2O=  rho*wO  *((D11Cat-D13Cat)*dxO_dwH2O + (D12Cat-D13Cat)*dxH2O_dwH2O)
     c_wH2O_wO=  rho*wH2OCat*((D21Cat-D23Cat)*dxO_dwO   + (D22Cat-D23Cat)*dxH2O_dwO)
     c_wH2O_wH2O=rho*wH2OCat*((D21Cat-D23Cat)*dxO_dwH2O + (D22Cat-D23Cat)*dxH2O_dwH2O)
  else
     c_wO_wO=0
     c_wO_wH2O=0
     c_wH2O_wO=0
     c_wH2O_wH2O=0
  endif

  if (iregion==2) then                   ! Nernst-Planck
     c_vMem_vMem=F*uH3O_2*zH3O*cH3O
     c_vMem_cH3O=DH3O_2
     c_cH2O_cH2O=DH2O_2

  else if (iregion==3) then
     c_vMem_vMem=F*uH3O_3*zH3O*cH3O
     c_vMem_cH3O=DH3O_3
     c_cH2O_cH2O=DH2O_3

  else if (iregion==4) then
     c_vMem_vMem=F*uH3O_4*zH3O*cH3O
     c_vMem_cH3O=DH3O_4
     c_cH2O_cH2O=DH2O_4

  else
     c_vMem_vMem=0
     c_vMem_cH3O=0
     c_cH2O_cH2O=0
  endif

  if (iregion==2.or.iregion==3.or.iregion==4) then  ! NP not region specific
     c_cH3O_cH3O=-D0*zH3O/zM
     c_cH3O_vMem=-F*u0*zH3O*cH3O
  else
     c_cH3O_cH3O=0
     c_cH3O_vMem=0
  endif

! ----------------------------------------------------------------------
! Return matrix components
! ----------------------------------------------------------------------
  
  if (tokenstem=="$C_VELE_VELE") then
     cfun=c_vEle_vEle
  else if (tokenstem=="$C_P_P") then
     cfun=c_p_p

  else if (tokenstem=="$C_WH_P") then
     cfun=c_wH_p
  else if (tokenstem=="$C_WH_WH") then
     cfun=c_wH_wH

  else if (tokenstem=="$C_WO_P") then
     cfun=c_wO_p
  else if (tokenstem=="$C_WO_WO") then
     cfun=c_wO_wO
  else if (tokenstem=="$C_WO_WH2O") then
     cfun=c_wO_wH2O

  else if (tokenstem=="$C_WH2O_P") then
     cfun=c_wH2O_p
  else if (tokenstem=="$C_WH2O_WO") then
     cfun=c_wH2O_wO
  else if (tokenstem=="$C_WH2O_WH2O") then
     cfun=c_wH2O_wH2O

  else if (tokenstem=="$C_VMEM_VMEM") then
     cfun=c_vMem_vMem
  else if (tokenstem=="$C_VMEM_CH3O") then
     cfun=c_vMem_cH3O
  else if (tokenstem=="$C_CH3O_VMEM") then
     cfun=c_cH3O_vMem
  else if (tokenstem=="$C_CH3O_CH3O") then
     cfun=c_cH3O_cH3O

  else if (tokenstem=="$C_CH2O_CH2O") then
     cfun=c_cH2O_cH2O

  else
     print *,'cfun: Unknown token: ',trim(tokenstem)
     stop
  endif

end function cfun

! ######################################################################

function ffun(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax)

  use fullfuelmod

  implicit none
  character(*) token
  integer nvar,istep,imax,jmax,kmax
  double precision ffun,x,y,z,ur(nvar),dur(nvar,3)

  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)  ! not used
  double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)            ! not used

  character(1) cregion
  character(100) tokenstem
  integer lentok,ios,iregion

  double precision vEle,p,wH,wO,wH2OCat,vMem,cH3O,cH2O,wH2O,wN
  double precision mass,rho,xH,xO,xH2OCat,xH2O,xN
  double precision px,py,wHx,wHy,wOx,wOy,wH2OCatx,wH2OCaty,u_Darcy,v_Darcy
  double precision zetAno,vAnoTil,xiAnoOx,xiAnoRed,zH,GAno,cHAgg,iotAno,iAno,Q2
  double precision zetCat,vCatTil,xiCatOx,xiCatRed,zO,GCat,cOAgg,iotCat,iCat,Q4
  double precision kd,alAno,alCat,cEqAno,cEqCat,SHAno,SH2OAno,SOCat,SH2OCat
  double precision RH,RO,RH2O,SH3O_2,SH2O_2,SH3O_4,SH2O_4,FAno,FCat
  double precision f_vEle,f_p,f_wH,f_wO,f_wH2O,f_vMem,f_cH2O
  integer, save :: istepold=-1

  call setconst   ! Setup constants if needed

  if (istep /= istepold) print *,'ffun: istep=',istep
  istepold=istep

! ----------------------------------------------------------------------
! Split token as name_region, eg C_P_P_1 -> C_P_P and 1
! ----------------------------------------------------------------------

  lentok=len_trim(token)

  cregion=token(lentok:lentok) ! last character "1" "2" etc

  read (cregion,*,iostat=ios) iregion

  if (ios /= 0) then
     print *,'Garbage region specified, got: ',cregion
     stop
  endif

  if (iregion < 1.or.iregion > 5) then
     print *,'Region out of bounds, got:',iregion
     stop
  endif

  tokenstem=token(:lentok-2)

! ----------------------------------------------------------------------
! Unpack variables
! labels = vEle p wH wO wH2OCat vMem cH3O cH2O
! ----------------------------------------------------------------------

  vEle=ur(1)
  p=ur(2)

  wH=ur(3)

  wO=ur(4)
  wH2OCat=ur(5)

  vMem=ur(6)
  cH3O=ur(7)
  cH2O=ur(8)

  wH2O=1-wH       ! Anode side
  wN=1-wO-wH2OCat

  if (iregion == 1.or.iregion == 2) then
     mass=1/(wH/mH+wH2O/mH2O)
  else if (iregion == 4.or.iregion==5) then
     mass=1/(wO/mO+wH2OCat/mH2O+wN/mN)
  else
     mass=0
  endif

  rho=p*mass/(R*T)

  xH=     mass*wH     /mH
  xO=     mass*wO     /mO
  xH2OCat=mass*wH2OCat/mH2O

  xH2O=1-xH            ! on Anode side
  xN=1-xO-xH2OCat

! ----------------------------------------------------------------------
! Gradients and Darcy velocity
! ----------------------------------------------------------------------

  px=dur(2,1)
  py=dur(2,2)

  wHx=dur(3,1)
  wHy=dur(3,2)

  wOx=dur(4,1)
  wOy=dur(4,2)

  wH2OCatx=dur(5,1)
  wH2OCaty=dur(5,2)

  if (iregion==1) then
     u_Darcy=-kap1*px/eta1
     v_Darcy=-kap1*py/eta1
  else if (iregion==2) then
     u_Darcy=-kap2*px/eta2
     v_Darcy=-kap2*py/eta2
  else if (iregion==3) then
     u_Darcy=0
     v_Darcy=0
  else if (iregion==4) then
     u_Darcy=-kap4*px/eta4
     v_Darcy=-kap4*py/eta4
  else
     u_Darcy=-kap5*px/eta5
     v_Darcy=-kap5*py/eta5
  endif

!----------------------------------------------------------------------
! Anode side. Calculate Q2. source Q2Ele = -Q2
! ----------------------------------------------------------------------

  zetAno=iAnoRef*rAno**2*sAno/(2*F*dHAgg*cHRef)
  vAnoTil=vEle-vMem-vAnoEq
  xiAnoOx= F*alpAnoOx /(R*T)
  xiAnoRed=F*alpAnoRed/(R*T)

  zH=zetAno*(exp(xiAnoOx*vAnoTil)-exp(-xiAnoRed*vAnoTil))

!  GAno=if(abs(zH)<1e-8,0,if(zH<0,sqrt(-zH)/tan(sqrt(-zH))-1,sqrt(zH)/tanh(sqrt(zH))-1))

  if(abs(zH)<1e-8) then
     GAno=0
  else 
     if (zH<0) then
        GAno=sqrt(-zH)/tan(sqrt(-zH))-1
     else
        GAno=sqrt(zH)/tanh(sqrt(zH))-1
     endif
  endif

  cHAgg=p*xH/kH
  iotAno=6*F*nuAno/rAno**2
  iAno=iotAno*dHAgg*cHAgg*GAno
  Q2=iAno
!  Q2=0

! ----------------------------------------------------------------------
! Cathode side. Calculate Q4. source Q4Ele = +Q4
! ----------------------------------------------------------------------

  zetCat=iCatRef*rCat**2*sCat/(4*F*dOAgg*cORef)
  vCatTil=vEle-vMem-vCatEq
  xiCatOx= F*alpCatOx /(R*T)
  xiCatRed=F*alpCatRed/(R*T)

  zO=zetCat*(-exp(xiCatOx*vCatTil)+exp(-xiCatRed*vCatTil))

! if(abs(zO)<1e-8,0,if(zO<0,sqrt(-zO)/tan(sqrt(-zO))-1,sqrt(zO)/tanh(sqrt(zO))-1))

  if(abs(zO)<1e-8) then
     GCat=0
  else
     if(zO<0) then
        GCat=sqrt(-zO)/tan(sqrt(-zO))-1
     else
        GCat=sqrt(zO)/tanh(sqrt(zO))-1
     endif
  endif

  cOAgg=p*xO/kO
  iotCat=12*F*nuCat/rCat**2
  iCat=iotCat*dOAgg*cOAgg*GCat
  Q4=iCat
!  Q4=0

! ----------------------------------------------------------------------
! alpha factors for use in S and S-bar factors
! ----------------------------------------------------------------------

  kd=4.59e-5/delta*OMH2O*cH2O/(OMMem*cH3O+OMH2O*cH2O)*exp(2416*(1/303.0-1/T))
  alAno=kd/OMMem/Cf
  alCat=alAno

! ----------------------------------------------------------------------
! equilibrium concentrations dependent of molar fractions
! ----------------------------------------------------------------------

  cEqAno=p*xH2O   /kH2OAno       ! xH2O: anode side
  cEqCat=p*xH2OCat/kH2OCat       ! xH2OCat: cathode side

! ----------------------------------------------------------------------
! S factors (based on Q2, Q4) [used for R- and F-factors]
! ----------------------------------------------------------------------

  SHAno=-Q2/(2*F)
  SH2OAno=-lamAno*drag*Q2/F-alAno*(cEqAno-cH2O)
  SOCat=-Q4/(4*F)
  SH2OCat=(2*lamCat*drag+1)*Q4/(2*F)-alCat*(cEqCat-cH2O)

! ----------------------------------------------------------------------
! R factors. Source term of Maxwell Stephan
! ----------------------------------------------------------------------

  RH=(1-wH)*mH*SHAno-wH*mH2O*SH2OAno
  RO=(1-wO)*mO*SOCat-wO*mH2O*SH2OCat
  RH2O=(1-wH2OCat)*mH2O*SH2OCat-wH2OCat*mO*SOCat

! ----------------------------------------------------------------------
! S-bar factors (2=Ano, 4=Cat). Source term of Nernst-Planck
! Note "drag" is "P" in functional spec.
! ----------------------------------------------------------------------

  SH3O_2=Q2/F
  SH2O_2=-(1-lamAno)*drag*Q2/F+alAno*(cEqAno-cH2O)

  SH3O_4=-Q4/F
  SH2O_4= (1-lamCat)*drag*Q4/F+alCat*(cEqCat-cH2O)

! ----------------------------------------------------------------------
! F factors. Source term for Darcy Law
! ----------------------------------------------------------------------

  FAno=mH*SHAno+mH2O*SH2OAno
  FCat=mO*SOCat+mH2O*SH2OCat

! ----------------------------------------------------------------------
! Region dependent
! ----------------------------------------------------------------------

  if (iregion==2) then
     f_vEle = -Q2
     f_p=mH*SHAno+mH2O*SH2OAno
  else if (iregion==4) then
     f_vEle =  Q4
     f_p=mO*SOCat+mH2O*SH2OCat
  else
     f_vEle=0
     f_p=0
  endif

  if (iregion==1) then
     f_wH=  -rho*(u_Darcy*wHx+v_Darcy*wHy)
!    f_wH=  -rho*dot_product(u_Darcy,gradwH)
  else if (iregion==2) then
     f_wH=RH-rho*(u_Darcy*wHx+v_Darcy*wHy)
!    f_wH=RH-rho*dot_product(u_Darcy,gradwH)
  else
     f_wH=0
  endif

  if (iregion==4) then
     f_wO=  RO  -rho*(u_Darcy*wOx     +v_Darcy*wOy)
     f_wH2O=RH2O-rho*(u_Darcy*wH2OCatx+v_Darcy*wH2OCaty)

!     f_wO=  RO  -rho*dot_product(u_Darcy,gradwO)
!     f_wH2O=RH2O-rho*dot_product(u_Darcy,gradwH2O)

  else if (iregion==5) then
     f_wO=      -rho*(u_Darcy*wOx     +v_Darcy*wOy)
     f_wH2O=    -rho*(u_Darcy*wH2OCatx+v_Darcy*wH2OCaty)

!     f_wO=      -rho*dot_product(u_Darcy,gradwO)
!     f_wH2O=    -rho*dot_product(u_Darcy,gradwH2O)

  else
     f_wO=0
     f_wH2O=0
  endif

  if (iregion==2) then
     f_vMem=SH3O_2
     f_cH2O=SH2O_2
  else if (iregion==4) then
     f_vMem=SH3O_4
     f_cH2O=SH2O_4
  else
     f_vMem=0
     f_cH2O=0
  endif

! ----------------------------------------------------------------------
! Return F-vector components
! ----------------------------------------------------------------------

  if (tokenstem=="$F_VELE") then
     ffun=f_vEle
  else if (tokenstem=="$F_P") then
     ffun=f_p
  else if (tokenstem=="$F_WH") then
     ffun=f_wH
  else if (tokenstem=="$F_WO") then
     ffun=f_wO
  else if (tokenstem=="$F_WH2O") then
     ffun=f_wH2O
  else if (tokenstem=="$F_VMEM") then
     ffun=f_vMem
  else if (tokenstem=="$F_CH2O") then
     ffun=f_cH2O
  else
     print *,"ffun: Unknown token: ",trim(tokenstem)
     stop
  endif

  if (isnan(ffun)) then
     print *,trim(token)
     print *,'NAN found'
     print *,'tokenstem,iregion=',trim(tokenstem),iregion
     print *,ur
     print *,dur
     print *,x,y,z
     print *,'Q2=',Q2,iAno,iotAno,dHAgg,cHAgg,GAno
     print *,'Q2 stuff=',zetAno,vAnoTil,xiAnoOx,xiAnoRed,zH
     print *,'vAnoTil parts:',vEle,vMem,vAnoEq
     stop
  endif

!  ffun=(1-exp(-istep/5.0))*ffun
!  ffun=0

end function ffun

! ######################################################################

function BCfun(token,x,y,z,ur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax)

! ----------------------------------------------------------------------
! Boundary concentrations on either side of membrane
! (this is a Dirichlet BC which depends on cH2O)
! ----------------------------------------------------------------------

  use fullfuelmod
  implicit none

  character(*) token
  double precision BCfun,x,y,z,ur(nvar)
  integer nvar,istep,imax,jmax,kmax
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

  double precision K0,H0,KT,cH3OAno,cH3OCat,cH2O

  call setconst
  
  cH2O=ur(8)

  K0=6.2
  H0=-52300
  KT=K0*exp(-H0/R*(1/T-1/298.0))
  cH3OAno=(sqrt((0.5*KT*cH2O/Cf)**2+KT*cH2O/Cf)-0.5*KT*cH2O/Cf)*Cf
  cH3OCat=cH3OAno

  if (token=="$CH3OANO") then
     bcfun=cH3OAno
  else if (token=="$CH3OCAT") then
     bcfun=cH3OCat
  else
     print *,"bcfun: Unknown token: ",trim(token)
     stop
  endif

end function BCfun

! ######################################################################

! labels = vEle p wH wO wH2OCat vMem cH3O cH2O
subroutine dCfun_du(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

  implicit none

  character(*) token
  double precision dCdu(nvar),x,y,z,ur(nvar),dur(nvar,3)
  integer nvar,istep,imax,jmax,kmax,i
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)

  double precision urp(nvar),Cfun
  double precision ur_step(8),Cfun0

  ur_step(1)=0.5            ! central values
  ur_step(2)=111430
  ur_step(3)=0.5
  ur_step(4)=0.5
  ur_step(5)=0.5
  ur_step(6)=0.5
  ur_step(7)=84
  ur_step(8)=27

  ur_step=ur_step*1e-3

  Cfun0=Cfun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)

  do i=1,8
     urp=ur
     urp(i)=urp(i)+ur_step(i)

     dCdu(i)=(Cfun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)-Cfun0)/ur_step(i)
  enddo

end subroutine dCfun_du

! ######################################################################

! labels = vEle p wH wO wH2OCat vMem cH3O cH2O
subroutine dCfun_ddu(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)

  implicit none

  character(*) token
  integer nvar,istep,imax,jmax,kmax
  double precision dCddu(nvar,3),x,y,z,ur(nvar),dur(nvar,3)
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*),Cfun

  double precision durp(nvar,3)
  double precision dur_step(8),Cfun0,dist
  integer i,j

  dist=2e-3                ! characteristic distance

  dur_step(1)=0.1          ! ranges across sim
  dur_step(2)=1000
  dur_step(3)=0.1
  dur_step(4)=0.1
  dur_step(5)=0.1
  dur_step(6)=0.1
  dur_step(7)=5
  dur_step(8)=5

  dur_step=dur_step/dist*1e-3

  Cfun0=Cfun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)

  do i=1,8
     do j=1,3
        durp=dur
        durp(i,j)=durp(i,j)+dur_step(i)
        dCddu(i,j)=(Cfun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)-Cfun0)/dur_step(i)
     enddo
  enddo
     
end subroutine dcfun_ddu

! ######################################################################

! labels = vEle p wH wO wH2OCat vMem cH3O cH2O
subroutine dffun_du(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)

  implicit none

  character(*) token
  double precision dfdu(nvar),x,y,z,ur(nvar),dur(nvar,3)
  integer nvar,istep,imax,jmax,kmax,i
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)

  double precision urp(nvar),ffun
  double precision ur_step(8),ffun0

  ur_step(1)=0.5            ! central values
  ur_step(2)=111430
  ur_step(3)=0.5
  ur_step(4)=0.5
  ur_step(5)=0.5
  ur_step(6)=0.5
  ur_step(7)=84
  ur_step(8)=27

  ur_step=ur_step*1e-3

  ffun0=ffun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)

  do i=1,8
     urp=ur
     urp(i)=urp(i)+ur_step(i)

     dfdu(i)=(ffun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)-ffun0)/ur_step(i)
  enddo

end subroutine dffun_du

! ######################################################################

! labels = vEle p wH wO wH2OCat vMem cH3O cH2O
subroutine dffun_ddu(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)

  implicit none

  character(*) token
  integer nvar,istep,imax,jmax,kmax
  double precision dfddu(nvar,3),x,y,z,ur(nvar),dur(nvar,3)
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)

  double precision durp(nvar,3)
  double precision dur_step(8),ffun0,dist,ffun
  integer i,j

  dist=2e-3                ! characteristic distance

  dur_step(1)=0.1          ! ranges across sim
  dur_step(2)=1000
  dur_step(3)=0.1
  dur_step(4)=0.1
  dur_step(5)=0.1
  dur_step(6)=0.1
  dur_step(7)=5
  dur_step(8)=5

  dur_step=dur_step/dist*1e-3

  ffun0=ffun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)

  do i=1,8
     do j=1,3
        durp=dur
        durp(i,j)=durp(i,j)+dur_step(i)
        dfddu(i,j)=(ffun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax)-ffun0)/dur_step(i)
     enddo
  enddo
     
end subroutine dffun_ddu
