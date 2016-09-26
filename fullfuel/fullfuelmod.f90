module fullfuelmod

  implicit none

  save
  logical :: init=.true.

  double precision pi,kB,F,R,T,TRef,pRef,pRefAno,pRefCat,mH,mN,mH2O,mO
  double precision sig1Ele,sig2Ele,sig4Ele,sig5Ele
  double precision eta1,eta2,eta4,eta5
  double precision kap1,kap2,kap4,kap5
  double precision pAnoH2O,pCatH2O,pH,pAir,pO,pN,rhoAno,rhoCat
  double precision wHIn,wOIn,wH2OIn,xHIn,xOIn,xH2OIn
  double precision CH,kHRef,kH,CO,kORef,kO,kH2OAno,kH2OCat
  double precision cHRef,cORef,pAnoIn,pCatIn,vCellAno,vCellCat,vAnoEq,vCatEq
  double precision drag,delta,rAno,rCat,nuAno,nuCat,sAno,sCat,epsAno,epsCat,epsMic
  double precision Ea,Tstar,IIAno,IICat,iAnoRef,iCatRef,alpAnoOx,alpAnoRed
  double precision alpCatOx,alpCatRed,dHAgg,dOAgg
  double precision dHH2O,dON,dOH2O,dH2ON,D0,u0
  double precision A,B,C,visc,Ri,DH3O_bulk
  double precision ep0,ep,del,muw,Rw,q,tauC
  double precision thI,thF,tauDG,kapG,IG,DH3O_Grotthus,DH3O_2,DH3O_3,DH3O_4
  double precision zH3O, zM, uH3O, uH3O_2, uH3O_3, uH3O_4
  double precision lamAno,lamCat,Cf,OMH2O,OMMem


contains

subroutine setconst

  if (init) then

  init=.false.

! ----------------------------------------------------------------------
! constants
! ----------------------------------------------------------------------

  pi=3.141592654d0
  kB=1.38e-23  
  F=96487
  R=8.314
  T=353
  TRef=298

  pRef=1.013e5
  pRefAno=pRef
  pRefCat=pRef

  mH=2e-3
  mO=32e-3
  mN=28e-3
  mH2O=18e-3

! ----------------------------------------------------------------------
! Conductivities, permeability (kap), viscosity (eta)
! ----------------------------------------------------------------------

  sig1Ele=1000 ! S/m
  sig2Ele=1200 ! S/m
  sig4Ele=1200 ! S/m
  sig5Ele=1000 ! S/m

  kap1=3e-13
  kap2=1e-13
  kap4=1e-13
  kap5=3e-13

  eta1=2.1e-5
  eta2=3.1e-5
  eta4=3.1e-5
  eta5=2.1e-5

! ----------------------------------------------------------------------
! To calculate wHIn etc
! ----------------------------------------------------------------------

  pAnoH2O=47400
  pCatH2O=47400

  pH=  pRefAno-pAnoH2O
  pAir=pRefCat-pCatH2O
  pO=0.21*pAir
  pN=pAir-pO

  rhoAno=(mH2O*pAnoH2O+mH*pH)/(R*T)
  rhoCat=(mH2O*pCatH2O+mO*pO+mN*pN)/(R*T)

  wHIn=mH*pH/(R*T)/rhoAno
  wOIn=mO*pO/(R*T)/rhoCat
  wH2OIn=mH2O*pCatH2O/(R*T)/rhoCat

  print *,'#####',mH2O,pCatH2O,R,T,rhoCat,wH2OIn

! ----------------------------------------------------------------------
! xHIn etc calculated from wHIn etc
! ----------------------------------------------------------------------

  xHIn=wHIn/mH/(wHIn/mH+(1-wHIn)/mH2O)
  xOIn=wOIn/mO/(wOIn/mO+wH2OIn/mH2O+(1-wOIn-wH2OIn)/mN)
  xH2OIn=wH2OIn/mH2O/(wOIn/mO+wH2OIn/mH2O+(1-wOIn-wH2OIn)/mN)

! ----------------------------------------------------------------------
! K factors
! ----------------------------------------------------------------------

  CH=500
  kHRef=129903.7
  kH=kHRef*exp(-CH*(1/T-1/TRef))

  CO=1700
  kORef=77942.2
  kO=kORef*exp(-CO*(1/T-1/TRef))

  kH2OAno=2.0e3
  kH2OCat=2.0e3

! ----------------------------------------------------------------------
! Used for Q2, Q4
! ----------------------------------------------------------------------

  cHRef=xHIn*pRef/kH
  cORef=xOIn*pRef/kO

  pAnoIn=pRef*1.1         ! input output pressures (needs to be in .zin)
  pCatIn=pRef*1.1

  vCellAno=0.0            ! Dirichlet BC
  vCellCat=0.7            ! Dirichlet BC (normally 0.7)

! ----------------------------------------------------------------------
! Equilibrium voltages, used for Q2, Q4
! ----------------------------------------------------------------------

  vAnoEq=0.0
  vCatEq=1.0

! ----------------------------------------------------------------------
! Misc
! ----------------------------------------------------------------------

  drag=1.0
  delta=20e-6

  rAno=1e-7
  rCat=1e-7

  nuAno=0.25
  nuCat=0.25

  sAno=4e7
  sCat=4e7

  epsAno=0.4
  epsCat=0.4
  epsMic=0.2


  Ea=67000         ! J/mol (in comsol input as 67 kJ/mol)
  Tstar=353

  IIAno=1e2        ! called IAno in comsol
  IICat=1e0

  iAnoRef=IIAno*exp(-Ea/R*(1/T-1/Tstar))
  iCatRef=IICat*exp(-Ea/R*(1/T-1/Tstar))

  alpAnoOx=0.5
  alpAnoRed=0.5

  alpCatOx=0.5
  alpCatRed=0.5

  dHAgg=1.2e-10*((1-epsAno)*epsMic)**1.5
  dOAgg=1.2e-10*((1-epsCat)*epsMic)**1.5

! ----------------------------------------------------------------------
! Binary diffusion coefficients
! ----------------------------------------------------------------------

  dHH2O=0.915e-4*(T/307.1)**1.5*epsAno**1.5
  dON  =0.220e-4*(T/293.2)**1.5*epsCat**1.5
  dOH2O=0.282e-4*(T/308.1)**1.5*epsCat**1.5
  dH2ON=0.256e-4*(T/307.5)**1.5*epsCat**1.5

! ----------------------------------------------------------------------
! DM=D0, uM=uO: membrane diffusivity and mobility
! ----------------------------------------------------------------------

  D0=1e-8
  u0=0

! ----------------------------------------------------------------------
! For DH3O_bulk (not used currently)
! ----------------------------------------------------------------------

  A=2.414e-5
  B=570.58
  C=140
  visc=A*exp(B/(T-C))
  Ri=0.141e-9

  DH3O_bulk=kB*T/(6*pi*visc*Ri)

! ----------------------------------------------------------------------
! for tauC (used in Grotthus)
! ----------------------------------------------------------------------

  ep0=8.854e-12
  ep=6
  del=0.143e-9   ! distance between proton...
  muw=2.95*3.336e-30
  Rw=0.141e-9
  q=1.602e-19

  tauC=32*pi**2*visc*ep*ep0*Rw**3*del**2/(muw*1*q)

! ----------------------------------------------------------------------
! For tauDG and hence DH3O_Grotthus
! ----------------------------------------------------------------------

  thI=2*pi/3
  thF=pi/12
  tauDG=tauC*log(tan(thI/2)/tan(thF/2))   ! Check log means log_e in Comsol

  kapG=6
  IG=0.255e-9

  DH3O_Grotthus=IG**2/kapG/tauDG

  DH3O_2=0.5*DH3O_Grotthus
  DH3O_3=    DH3O_Grotthus
  DH3O_4=0.5*DH3O_Grotthus

  lamAno=1.0
  lamCat=1.0

  Cf=1200          ! Combined concentration sulphonate acide + ions
  OMH2O=1.8e-5     ! Molar vol of water
  OMMem=5.5e-4     ! Molar vol of membrane

! ----------------------------------------------------------------------
! Charge factors and mobility
! ----------------------------------------------------------------------

  zH3O=1           ! charge factor of H3O
  zM=-1            ! charge factor of membrane

  uH3O=DH3O_Grotthus/(R*T)

  uH3O_2=DH3O_2/(R*T)
  uH3O_3=DH3O_3/(R*T)
  uH3O_4=DH3O_4/(R*T)

  print *,'This is setconst'
  call dumpconst

  endif

end subroutine setconst

! ######################################################################

subroutine dumpconst

  print *,F
  print *,R
  print *,T
  print *,pRef
  print *,pRefAno
  print *,pAnoH2O
  print *,pCatH2O
  print *,pH
  print *,pAir
  print *,pO
  print *,pN
  print *,rhoAno
  print *,rhoCat

  print *,wHIn
  print *,wOIn
  print *,wH2OIn
  print *,xHIn
  print *,xOIn

  print *,cHRef
  print *,cORef
  print *,pAnoIn
  print *,pCatIn
  print *,vCellAno
  print *,vCellCat
  print *,vAnoEq
  print *,vCatEq
  print *,sig1Ele
  print *,sig2Ele
  print *,sig4Ele
  print *,sig5Ele

  print *,kap1
  print *,kap2
  print *,kap4
  print *,kap5

  print *,eta1
  print *,eta2
  print *,eta4
  print *,eta5

  print *,mH
  print *,mO
  print *,mN
  print *,mH2O
  print *,TRef
  print *,kHRef
  print *,CH
  print *,kH
  print *,kORef
  print *,CO
  print *,kO

  print *,kH2OAno
  print *,kH2OCat
  print *,drag
  print *,delta
  print *,rAno
  print *,rCat
  print *,nuAno
  print *,nuCat
  print *,sAno
  print *,sCat

  print *,epsAno
  print *,epsCat
  print *,epsMic

  print *,Ea
  print *,Tstar
  print *,IIAno
  print *,IICat
  print *,iAnoRef
  print *,iCatRef

  print *,alpAnoOx
  print *,alpAnoRed
  print *,alpCatOx
  print *,alpCatRed

  print *,dHAgg
  print *,dOAgg
  print *,dHH2O
  print *,dON
  print *,dOH2O
  print *,dH2ON

  print *,D0
  print *,u0
  print *,A
  print *,B
  print *,C
  print *,visc
  print *,kB
  print *,Ri
  print *,DH3O_bulk
  print *,'no mMem'
  print *,'no rhoMem'

  print *,kapG
  print *,IG
  print *,ep0
  print *,ep
  print *,del
  print *,muw
  print *,Rw
  print *,q
  print *,tauC
  print *,thI
  print *,thF

  print *,tauDG
  print *,DH3O_Grotthus
  print *,DH3O_2
  print *,DH3O_3
  print *,DH3O_4

  print *,lamAno
  print *,lamCat
  print *,'no K0'
  print *,'no H0'
  print *,Cf
  print *,OMH2O
  print *,OMMem
  print *,zH3O
  print *,zM

  print *,uH3O
  print *,xH2OIn
  print *,uH3O_2
  print *,uH3O_3
  print *,uH3O_4

end subroutine dumpconst

end module fullfuelmod
