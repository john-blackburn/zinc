program ss

! ----------------------------------------------------------------------
! Program to post process results from Zinc simulation:
! planealcu.zin
! key_sim=1: with applied force (use this for comparison)
! key_sim=2: with clamping (not comparable to planalcu.zin)
! ----------------------------------------------------------------------

  use matrix
  implicit none

  integer,parameter :: long=selected_real_kind(12,300)

  real(long), dimension(1000) :: xS,yS,exxS,eyyS,ezzS,exyS,sxxS,syyS,szzS,sxyS,&
       uxS,uyS

  real(long), dimension(6,6) :: sf,cf,sm,cm
  real(long) smdiag,smoff,smtail,sfdiag,sfoff,sftail

  real(long) err,eoo,ezz,srr,soo,szz,urS,uoS,ux,uy,ur,uo
  real(long) sxx,syy,sxy,exx,eyy,exy,dum
  real(long) xmax,xmin,ymax,ymin,xc,yc,x,y,r,theta,s,c,a
  real(long) invkTf,invkTm,kTf,kTm,invlambda,lambda,phi,Af,Am
  real(long) sT,EAf,Etf,nuAf,nuTf,muAf,muTf,EAm,ETm,nuAm,nuTm,muAm,muTm

  real(long) srrS,sooS,errS,eooS,outerad,UU,Vf,Vm

  integer i,j,n,key_sim,Nx,Ny,ix,iy

  character(20) uxfile,uyfile,exxfile,eyyfile,ezzfile,exyfile,sxxfile,syyfile,szzfile,sxyfile

! ----------------------------------------------------------------------

  key_sim=1        ! applied force
  
  a=0.5 ! m
  sT=1 ! N/m2
  
  Nx=30
  Ny=30
  
  cm(1,:)=[1.203e11,5.888e10,5.888e10,0.0,0.0,0.0]           ! Al N/m2
  cm(2,:)=[5.888e10,1.203e11,5.888e10,0.0,0.0,0.0]
  cm(3,:)=[5.888e10,5.888e10,1.203e11,0.0,0.0,0.0]
  cm(4,:)=[0,0,0,1,0,0]*2.8556e10
  cm(5,:)=[0,0,0,0,1,0]*2.8556e10
  cm(6,:)=[0,0,0,0,0,1]*2.8556e10
  
  cf(1,:)=[2.282e11,9.373e10,9.373e10,0.0,0.0,0.0]         ! Cu
  cf(2,:)=[9.373e10,2.282e11,9.373e10,0.0,0.0,0.0]
  cf(3,:)=[9.373e10,9.373e10,2.282e11,0.0,0.0,0.0]
  cf(4,:)=[0,0,0,1,0,0]*4.695e10
  cf(5,:)=[0,0,0,0,1,0]*4.695e10
  cf(6,:)=[0,0,0,0,0,1]*4.695e10

  uxfile='planealcu01.out'
  uyfile='planealcu02.out'
  
  exxfile='planealcu03.out'
  eyyfile='planealcu04.out'
  ezzfile='planealcu05.out'
  exyfile='planealcu06.out'
  
  sxxfile='planealcu07.out'
  syyfile='planealcu08.out'
  szzfile='planealcu09.out'
  sxyfile='planealcu10.out'

! ----------------------------------------------------------------------
! read in computational stress and strain
! ----------------------------------------------------------------------

  open (11,file=uxfile,status='old')
  open (12,file=uyfile,status='old')

  open (13,file=exxfile,status='old')
  open (14,file=eyyfile,status='old')
  open (15,file=ezzfile,status='old')
  open (16,file=exyfile,status='old')

  open (17,file=sxxfile,status='old')
  open (18,file=syyfile,status='old')
  open (19,file=szzfile,status='old')
  open (20,file=sxyfile,status='old')

  do i=11,20
     do j=1,5
        read (i,*)
     enddo
  enddo

  n=0
  do ix=0,Nx
     do iy=0,Ny
        n=n+1

        read (11,*) xS(n),yS(n),dum,uxS(n)
        read (12,*) xS(n),yS(n),dum,uyS(n)

        read (13,*) xS(n),yS(n),dum,exxS(n)
        read (14,*) xS(n),yS(n),dum,eyyS(n)
        read (15,*) xS(n),yS(n),dum,ezzS(n)
        read (16,*) xS(n),yS(n),dum,exyS(n)

        read (17,*) xS(n),yS(n),dum,sxxS(n)
        read (18,*) xS(n),yS(n),dum,syyS(n)
        read (19,*) xS(n),yS(n),dum,szzS(n)
        read (20,*) xS(n),yS(n),dum,sxyS(n)
     enddo
     do i=11,20
        read (i,*)
     enddo
  enddo

  do i=11,20
     close (i)
  enddo

! ----------------------------------------------------------------------
! set material properties and excitation sT
! ----------------------------------------------------------------------

  sf=inverse(cf)

  print *,'sf:'
  do i=1,6
     print '(6f7.3)',(sf(i,j),j=1,6)
  enddo

  sfdiag=(sf(1,1)+sf(2,2)+sf(3,3))/3
  sfoff= (sf(1,2)+sf(1,3)+sf(2,3))/3
  sftail=(sf(4,4)+sf(5,5)+sf(6,6))/3

  EAf=1/sfdiag
  ETf=1/sfdiag

  nuAf=-sfoff/sfdiag
  nuTf=-sfoff/sfdiag

  muAf=1/sftail
  muTf=1/sftail

  sm=inverse(cm)

  print *,'sm:'
  do i=1,6
     print '(6f7.3)',(sm(i,j),j=1,6)
  enddo

  smdiag=(sm(1,1)+sm(2,2)+sm(3,3))/3
  smoff= (sm(1,2)+sm(1,3)+sm(2,3))/3
  smtail=(sm(4,4)+sm(5,5)+sm(6,6))/3

  EAm=1/smdiag
  ETm=1/smdiag

  nuAm=-smoff/smdiag
  nuTm=-smoff/smdiag

  muAm=1/smtail
  muTm=1/smtail

  print *,'Ef=',EAf
  print *,'nuf=',nuAf
  print *,'muf=',muAf

  print *,'Em=',EAm
  print *,'num=',nuAm
  print *,'mum=',muAm

! ----------------------------------------------------------------------
! intermediate quantities, assuming delta T=eps=0
! ----------------------------------------------------------------------

  invkTf=2*(1-nuTf)/ETf-4*nuAf**2/EAf
  invkTm=2*(1-nuTm)/ETm-4*nuAm**2/EAm
  
  kTf=1/invkTf
  kTm=1/invkTm
  
  if (key_sim.eq.1) then
     
     invlambda=0.5*(1/kTf+1/muTm)
     lambda=1/invlambda
     
     phi=lambda*0.5*(1/kTf-1/kTm)*sT
     
     Af=1/(2*kTf)*(sT-phi)
     Am=sT/(2*kTm)
     
  else if (key_sim.eq.2) then

     Vf=(a/outerad)**2
     Vm=1-Vf

     invlambda=0.5*(1/muTm+Vf/kTm+Vm/kTf)
     lambda=1/invlambda

     sT=2*kTm*UU/outerad/(1+0.5*lambda*Vf*(1+kTm/muTm)*(1/kTf-1/kTm))

     phi=a**2*lambda*0.5*(1/kTf-1/kTm)*sT

     Af=1/(2*kTf)*(sT-Vm*phi/a**2)
     Am=1/(2*kTm)*(sT+Vf*phi/a**2)

  endif

  print *,'phi=',phi
  print *,'Af=',Af
  print *,'Am=',Am
  print *,'kTf=',kTf
  print *,'kTm=',kTm
  
! ----------------------------------------------------------------------
! loop over all repatoms. get positiosn and convert into cylindrical
! ----------------------------------------------------------------------

  xmax=maxval(xS)
  xmin=minval(xS)

  ymax=maxval(yS)
  ymin=minval(yS)

  xc=0
  yc=0

  open (1,file='scart.out',status='unknown')
  open (2,file='ecart.out',status='unknown')
  open (3,file='spolar.out',status='unknown')
  open (4,file='epolar.out',status='unknown')
  
  write (1,'(10a13)') 'x','y','sxx','syy','szz','sxy','sxxS','syyS','szzS','sxyS'
  write (2,'(10a13)') 'x','y','exx','eyy','ezz','exy','exxS','eyyS','ezzS','exyS'
  write (3,'(8a13)') 'r','theta','srr','soo','szz','srrS','sooS','szzS'
  write (4,'(8a13)') 'r','theta','err','eoo','ezz','errS','eooS','ezzS'

  do i=1,n
     x=xS(i)
     y=yS(i)

     r=sqrt((x-xc)**2+(y-yc)**2)
     theta=atan2(y-yc,x-xc)

     s=sin(theta)
     c=cos(theta)

! ----------------------------------------------------------------------
! calculate theoretical u, e, s
! Theoretically: uo=0, uz=0 (here with no epsz)
! sro=srz=soz=0 same for e
! ----------------------------------------------------------------------

     if (r.lt.a) then
!        ur=Af*r

        if (key_sim.eq.1) then
           err=Af
           eoo=Af
           ezz=0
           
           srr=sT-phi
           soo=sT-phi
           szz=2*nuAf*(sT-phi)
        else
           err=Af
           eoo=Af
           ezz=0

           srr=sT-Vm*phi/a**2
           soo=sT-Vm*phi/a**2
           szz=2*nuAf*(sT-Vm*phi/a**2)
        endif

     else
!        ur=Am*r+a**2*phi/(2*muTm*r)

        if (key_sim.eq.1) then
           err=Am-a**2*phi/(2*muTm*r**2)
           eoo=Am+a**2*phi/(2*muTm*r**2)
           ezz=0
        
           srr=sT-a**2/r**2*phi
           soo=sT+a**2/r**2*phi
           szz=2*nuAm*sT
        else
           err=Am-phi/(2*muTm*r**2)
           eoo=Am+phi/(2*muTm*r**2)
           ezz=0

           srr=sT+phi*(1/outerad**2-1/r**2)
           soo=sT+phi*(1/outerad**2+1/r**2)
           szz=2*nuAm*(sT+Vf*phi/a**2)
        endif
     endif

     sxx=srr*c**2+soo*s**2
     syy=srr*s**2+soo*c**2
     sxy=srr*c*s-soo*s*c

     exx=err*c**2+eoo*s**2
     eyy=err*s**2+eoo*c**2
     exy=err*c*s-eoo*s*c

     write (1,'(10e13.5)') x,y,sxx,syy,szz,sxy,sxxS(i),syyS(i),szzS(i),sxyS(i)
     write (2,'(10e13.5)') x,y,exx,eyy,ezz,exy,exxS(i),eyyS(i),ezzS(i),exyS(i)

! ----------------------------------------------------------------------
! convert computational result to polars
! sxxS(i) stored computational values. srrC computational srr
! srr, sxx theoretical values
! ----------------------------------------------------------------------

     srrS=sxxS(i)*c**2+syyS(i)*s**2+sxyS(i)*2*s*c
     sooS=sxxS(i)*s**2+syyS(i)*c**2-sxyS(i)*2*s*c

     errS=exxS(i)*c**2+eyyS(i)*s**2+exyS(i)*2*s*c
     eooS=exxS(i)*s**2+eyyS(i)*c**2-exyS(i)*2*s*c

     write (3,'(8e13.5)') r,theta,srr,soo,szz,srrS,sooS,szzS(i)
     write (4,'(8e13.5)') r,theta,err,eoo,ezz,errS,eooS,ezzS(i)

     if (mod(i,Ny+1).eq.0) then
        write (1,*)
        write (2,*)
        write (3,*)
        write (4,*)
     endif

  enddo

  close (1)
  close (2)
  close (3)
  close (4)

! ----------------------------------------------------------------------
! Calculate theoretical displacements
! ----------------------------------------------------------------------

  open (1,file='ucart.out',status='unknown')
  open (2,file='upolar.out',status='unknown')

  write (1,'(6a13)') 'x','y','ux','uy','uxS','uyS'
  write (2,'(6a13)') 'r','theta','ur','uo','urS','uoS'

  do i=1,n
     x=xS(i)
     y=yS(i)

     r=sqrt((x-xc)**2+(y-yc)**2)
     theta=atan2(y-yc,x-xc)

     s=sin(theta)
     c=cos(theta)

     if (r.lt.a) then
        ur=Af*r
        uo=0
     else
        if (key_sim.eq.1) then
           ur=Am*r+a**2*phi/(2*muTm*r)
           uo=0
        else
           ur=Am*r+phi/(2*muTm*r)
           uo=0
        endif
     endif

     ux=ur*c-uo*s
     uy=ur*s+uo*c

     write (1,'(6e13.5)') x,y,ux,uy,uxS(i),uyS(i)

     urS= c*uxS(i)+s*uyS(i)
     uoS=-s*uxS(i)+c*uyS(i)

     write (2,'(6e13.5)') r,theta,ur,uo,urS,uoS

  enddo

  close (1)
  close (2)

end program ss
