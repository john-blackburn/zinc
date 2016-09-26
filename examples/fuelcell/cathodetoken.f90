module nlinc
  integer, save :: init=1
  double precision, save :: R,T,m1,m2,m3,kappa,eta,Pref,epsilon
  double precision, save :: D12tilde,D13tilde,D23tilde,D21tilde,D31tilde,D32tilde
  double precision, save :: x1c,x2c,x3c,pc,RT,m13,m23
end module nlinc

subroutine NLinit

! ----------------------------------------------------------------------
! If the first call, read from cathode.con
! ----------------------------------------------------------------------
  
  use nlinc
  implicit none

  integer i,ios
  character(1000) line,lefts,rights

  if (init.eq.1) then

     print *,'This is the Fortran version. Zinc linked to cathodetoken.f90'
     print *,'NLinit, reading cathodetoken.con'

     init=0

     open (20,file='cathodetoken.con',status='old')
  
     do
        read (20,'(a)',iostat=ios) line
        
        if (ios.ne.0) exit
        
        i=index(line,'=')
        
        if (i.ne.0) then
           lefts=line(:i-1)
           rights=line(i+1:)
           
           if (lefts.eq.'R') then
              read (rights,*) R
              
           else if (lefts.eq.'T') then
              read (rights,*) T
              
           else if (lefts.eq.'m1') then
              read (rights,*) m1
              
           else if (lefts.eq.'m2') then
              read (rights,*) m2
              
           else if (lefts.eq.'m3') then
              read (rights,*) m3
              
           else if (lefts.eq.'kappa') then
              read (rights,*) kappa
              
           else if (lefts.eq.'eta') then
              read (rights,*) eta
              
           else if (lefts.eq.'Pref') then
              read (rights,*) Pref
              
           else if (lefts.eq.'epsilon') then
              read (rights,*) epsilon
              
           else if (lefts.eq.'D12tilde') then
              read (rights,*) D12tilde
           else if (lefts.eq.'D13tilde') then
              read (rights,*) D13tilde
           else if (lefts.eq.'D23tilde') then
              read (rights,*) D23tilde
              
           else if (lefts.eq.'x1c') then
              read (rights,*) x1c
           else if (lefts.eq.'x2c') then
              read (rights,*) x2c
           else if (lefts.eq.'x3c') then
              read (rights,*) x3c

           else if (lefts.eq.'pc') then
              read (rights,*) pc
           endif
        endif
     enddo

     close (20)

     m13=m1-m3
     m23=m2-m3
     RT=R*T

     D21tilde=D12tilde
     D31tilde=D13tilde
     D32tilde=D23tilde

  endif
end subroutine NLinit

! ######################################################################

function Cfun(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

  use nlinc
  implicit none
  character label*(*)
  integer nvar,istep,imax,jmax,kmax
  double precision Cfun,x,y,z,ur(nvar),dur(nvar,3)
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

  double precision D11,D12,D13,D21,D22,D23,D31,D32,D33
  double precision x1,x2,p,x3,w1,w2,w3,m,denom

! ----------------------------------------------------------------------
! Set values appropriately
! ----------------------------------------------------------------------

  call NLinit

  x1=ur(1)
  x2=ur(2)
  p=ur(3)

  x3=1-x1-x2

  m=m1*x1+m2*x2+m3*x3

  w1=m1*x1/m
  w2=m2*x2/m
  w3=m3*x3/m

  denom=x1/(D12tilde*D13tilde)+x2/(D12tilde*D23tilde)+x3/(D13tilde*D23tilde)

  D11=((w2+w3)**2/(x1*D23tilde)+w2**2/(x2*D13tilde)+w3**2/(x3*D12tilde))/denom
  D22=((w1+w3)**2/(x2*D13tilde)+w1**2/(x1*D23tilde)+w3**2/(x3*D12tilde))/denom
  D33=((w1+w2)**2/(x3*D12tilde)+w1**2/(x1*D23tilde)+w2**2/(x2*D13tilde))/denom

  D12=(w1*(w2+w3)/(x1*D23tilde)+w2*(w1+w3)/(x2*D13tilde)-w3**2/(x3*D12tilde))/denom
  D13=(w1*(w2+w3)/(x1*D23tilde)+w3*(w1+w2)/(x3*D12tilde)-w2**2/(x2*D13tilde))/denom
  D23=(w2*(w1+w3)/(x2*D13tilde)+w3*(w1+w2)/(x3*D12tilde)-w1**2/(x1*D23tilde))/denom

  D21=D12
  D31=D13
  D32=D23

  if (label.eq.'$B11') then
     Cfun=p*m1*x1*(D11-D13)/RT

  else if (label.eq.'$B12') then
     Cfun=p*m1*x1*(D12-D13)/RT

  else if (label.eq.'$B21') then
     Cfun=p*m2*x2*(D21-D23)/RT

  else if (label.eq.'$B22') then
     Cfun=p*m2*x2*(D22-D23)/RT

  else if (label.eq.'$Q1') then
     Cfun=m1*x1/RT*(x1*(D11*(1-m1/m)-D13*(1-m3/m)) &
                  + x2*(D12*(1-m2/m)-D13*(1-m3/m)) & 
                                    +D13*(1-m3/m))

  else if (label.eq.'$Q2') then
     Cfun=m2*x2/RT*(x1*(D21*(1-m1/m)-D23*(1-m3/m)) &
                  + x2*(D22*(1-m2/m)-D23*(1-m3/m)) & 
                                    +D23*(1-m3/m))

  else if (label.eq.'$E') then
     Cfun=p*m*kappa/(RT*eta)
  else
     print *,'Error in Cfun: token not recognised:',trim(label)
     stop
  endif
  
end function Cfun

! ######################################################################

function ffun(label,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax)

  use nlinc
  implicit none

  character label*(*)  
  integer nvar,istep,imax,jmax,kmax
  double precision x,y,z,ffun,ur(nvar),dur(nvar,3)
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)


  double precision x1,x2,p,x3,m,x1x,x1y,x1z,x2x,x2y,x2z,px,py,pz
  double precision mx,my,mz

  call NLinit

  x1=ur(1)
  x2=ur(2)
  p=ur(3)

  x3=1-x1-x2

  m=m1*x1+m2*x2+m3*x3

  x1x=dur(1,1)
  x1y=dur(1,2)
  x1z=dur(1,3)

  x2x=dur(2,1)
  x2y=dur(2,2)
  x2z=dur(2,3)
  
  px=dur(3,1)
  py=dur(3,2)
  pz=dur(3,3)

  mx=m1*x1x + m2*x2x + m3*(-x1x-x2x)
  my=m1*x1y + m2*x2y + m3*(-x1y-x2y)
  mz=m1*x1z + m2*x2z + m3*(-x1z-x2z)

  if (label.eq.'$F1') then
     ffun=p*kappa*m1/(eta*RT)*(px*(x1x-x1*mx/m)+py*(x1y-x1*my/m)+pz*(x1z-x1*mz/m))
  else if (label.eq.'$F2') then
     ffun=p*kappa*m2/(eta*RT)*(px*(x2x-x2*mx/m)+py*(x2y-x2*my/m)+pz*(x2z-x2*mz/m))
  else
     print *,'ffun: unknown token encountered'
     stop
  endif

end function ffun

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dCfun_du(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dCdu)

  implicit none

  character(*) token
  double precision dCdu(nvar),x,y,z,ur(nvar),dur(nvar,3)
  integer nvar,istep,imax,jmax,kmax
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)

  double precision urp(nvar),Cfun
  double precision x1cent,x2cent,pcent,dx1,dx2,dp

  x1cent =0.5    ! x1,x2=[0,1]
  x2cent =0.5
  pcent =111430   ! near 1 atm

  dx1=x1cent*1e-3
  dx2=x2cent*1e-3
  dp=  pcent*1e-3

  urp=ur
  urp(1)=urp(1)+dx1
  dCdu(1)=(Cfun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
          -Cfun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/dx1

  urp=ur
  urp(2)=urp(2)+dx2
  dCdu(2)=(Cfun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
          -Cfun(token,x,y,z,ur ,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/dx2

  urp=ur
  urp(3)=urp(3)+dp
  dCdu(3)=(Cfun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
          -Cfun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/dp

end subroutine dCfun_du

! ######################################################################
     

subroutine dCfun_ddu(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dCddu)

  implicit none

  character(*) token
  integer nvar,istep,imax,jmax,kmax
  double precision dCddu(nvar,3),x,y,z,ur(nvar),dur(nvar,3)
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*),Cfun

  double precision durp(nvar,3)
  double precision x1cent,x1range,x2cent,x2range,pcent,prange,dist,ddx1,ddx2,ddp
  integer i,j

  x1cent =0.5    ! x1,x2=[0,1]
  x1range=0.1

  x2cent =0.5
  x2range=0.1

  pcent =111430   ! near 1 atm
  prange= 10000

  dist=2e-3

  ddx1=x1range/dist*0.0001
  ddx2=x2range/dist*0.0001
  ddp=  prange/dist*0.0001

  ! x1
  do j=1,3
     durp=dur
     durp(1,j)=durp(1,j)+ddx1
     dCddu(1,j)=(Cfun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
                -Cfun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/ddx1
  enddo

  ! x2
  do j=1,3
     durp=dur
     durp(2,j)=durp(2,j)+ddx2
     dCddu(2,j)=(Cfun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
                -Cfun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/ddx2
  enddo

  ! p
  do j=1,3
     durp=dur
     durp(3,j)=durp(3,j)+ddp
     dCddu(3,j)=(Cfun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
                -Cfun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/ddp
  enddo

end subroutine dcfun_ddu

! ######################################################################

subroutine dffun_du(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dfdu)
  use nlinc
  implicit none

  character(*) token
  double precision dfdu(nvar),x,y,z,ur(nvar),dur(nvar,3)
  integer nvar,istep,imax,jmax,kmax
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)

  double precision urp(nvar),ffun
  double precision x1cent,x1range,x2cent,x2range,pcent,prange,dist,dx1,dx2,dp

  double precision x1,x2,p,x1x,x1y,x1z,x2x,x2y,x2z,px,py,pz,topx,topy,topz,bot
  logical, parameter :: analytical=.true.

  call NLinit

  if (analytical) then
     x1=ur(1)
     x2=ur(2)
     p=ur(3)

     x1x=dur(1,1)
     x1y=dur(1,2)
     x1z=dur(1,3)

     x2x=dur(2,1)
     x2y=dur(2,2)
     x2z=dur(2,3)
  
     px=dur(3,1)
     py=dur(3,2)
     pz=dur(3,3)

     topx=(m1-m3)*x1x+(m2-m3)*x2x
     topy=(m1-m3)*x1y+(m2-m3)*x2y
     topz=(m1-m3)*x1z+(m2-m3)*x2z
     
     bot =(m1-m3)*x1 +(m2-m3)*x2 +m3

     if (token=='$F1') then

        dfdu(1)=p*kappa*m1/(eta*RT) &
             *(px*(-topx/bot+x1*topx/bot**2*(m1-m3)) &
             + py*(-topy/bot+x1*topy/bot**2*(m1-m3)) &
             + pz*(-topz/bot+x1*topz/bot**2*(m1-m3)))

        dfdu(2)=p*kappa*m1/(eta*RT) &
             *(px*(          x1*topx/bot**2*(m2-m3)) &
             + py*(          x1*topy/bot**2*(m2-m3)) &
             + pz*(          x1*topz/bot**2*(m2-m3)))

        dfdu(3)=kappa*m1/(eta*RT) &
             *(px*(x1x-x1*topx/bot)  &
             + py*(x1y-x1*topy/bot)  &
             + pz*(x1z-x1*topz/bot))

     else if (token=='$F2') then

        dfdu(1)=p*kappa*m2/(eta*RT) &
             *(px*(          x2*topx/bot**2*(m1-m3)) &
             + py*(          x2*topy/bot**2*(m1-m3)) &
             + pz*(          x2*topz/bot**2*(m1-m3)))

        dfdu(2)=p*kappa*m2/(eta*RT) &
             *(px*(-topx/bot+x2*topx/bot**2*(m2-m3)) &
             + py*(-topy/bot+x2*topy/bot**2*(m2-m3)) &
             + pz*(-topz/bot+x2*topz/bot**2*(m2-m3)))

        dfdu(3)=kappa*m2/(eta*RT) &
             *(px*(x2x-x2*topx/bot)  &
             + py*(x2y-x2*topy/bot)  &
             + pz*(x2z-x2*topz/bot))
     else
        print *,'dffun_du: Unknown token'
        stop
     endif
  else

     x1cent =0.5    ! x1,x2=[0,1]
     x2cent =0.5
     pcent =111430   ! near 1 atm
     
     dx1=x1cent*1e-3
     dx2=x2cent*1e-3
     dp=  pcent*1e-3
     
     print *,'dffun_du:'
     print *,ur
     print *,dx1,dx2,dp
     
     urp=ur
     urp(1)=urp(1)+dx1
     dfdu(1)=(ffun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
             -ffun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/dx1
     
     urp=ur
     urp(2)=urp(2)+dx2
     dfdu(2)=(ffun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
             -ffun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/dx2

     urp=ur
     urp(3)=urp(3)+dp
     dfdu(3)=(ffun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
             -ffun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/dp
  endif

end subroutine dffun_du

! ######################################################################

subroutine dffun_ddu(token,x,y,z,ur,dur,nvar,istep, &
     ireg,iregup,rnode,vec,imax,jmax,kmax,dfddu)
  use nlinc
  implicit none

  character(*) token
  integer nvar,istep,imax,jmax,kmax
  double precision dfddu(nvar,3),x,y,z,ur(nvar),dur(nvar,3)
  integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
  double precision rnode(0:imax,0:jmax,0:kmax,*),vec(*)

  double precision durp(nvar,3),ffun
  double precision x1cent,x1range,x2cent,x2range,pcent,prange,dist,ddx1,ddx2,ddp
  integer i,j

  double precision x1,x2,p,x1x,x1y,x1z,x2x,x2y,x2z,px,py,pz,bot
  logical, parameter :: analytical=.true.

  call NLinit

  if (analytical) then
     x1=ur(1)
     x2=ur(2)
     p=ur(3)

     x1x=dur(1,1)
     x1y=dur(1,2)
     x1z=dur(1,3)
     
     x2x=dur(2,1)
     x2y=dur(2,2)
     x2z=dur(2,3)
  
     px=dur(3,1)
     py=dur(3,2)
     pz=dur(3,3)

     bot=(m1-m3)*x1+(m2-m3)*x2+m3

     if (token=='$F1') then
        dfddu(1,:)=p*kappa*m1/(eta*RT)*(1-x1*(m1-m3)/bot)*[px,py,pz]
        dfddu(2,:)=p*kappa*m1/(eta*RT)*( -x1*(m2-m3)/bot)*[px,py,pz]

        dfddu(3,1)=p*kappa*m1/(eta*RT)*(x1x-x1*((m1-m3)*x1x+(m2-m3)*x2x)/bot)
        dfddu(3,2)=p*kappa*m1/(eta*RT)*(x1y-x1*((m1-m3)*x1y+(m2-m3)*x2y)/bot)
        dfddu(3,3)=p*kappa*m1/(eta*RT)*(x1z-x1*((m1-m3)*x1z+(m2-m3)*x2z)/bot)

     else if (token=='$F2') then
        dfddu(1,:)=p*kappa*m2/(eta*RT)*( -x2*(m1-m3)/bot)*[px,py,pz]
        dfddu(2,:)=p*kappa*m2/(eta*RT)*(1-x2*(m2-m3)/bot)*[px,py,pz]

        dfddu(3,1)=p*kappa*m2/(eta*RT)*(x2x-x2*((m1-m3)*x1x+(m2-m3)*x2x)/bot)
        dfddu(3,2)=p*kappa*m2/(eta*RT)*(x2y-x2*((m1-m3)*x1y+(m2-m3)*x2y)/bot)
        dfddu(3,3)=p*kappa*m2/(eta*RT)*(x2z-x2*((m1-m3)*x1z+(m2-m3)*x2z)/bot)
     else
        print *,'dffun_ddu: bad token'
        stop
     endif
  else

     x1cent =0.5    ! x1,x2=[0,1]
     x1range=0.1
     
     x2cent =0.5
     x2range=0.1
     
     pcent =111430   ! near 1 atm
     prange= 10000
     
     dist=2e-3
     
     ddx1=x1range/dist*0.0001
     ddx2=x2range/dist*0.0001
     ddp=  prange/dist*0.0001
     
     ! x1
     do j=1,3
        durp=dur
        durp(1,j)=durp(1,j)+ddx1
        dfddu(1,j)=(ffun(token,x,y,z,ur,durp, nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
                   -ffun(token,x,y,z,ur,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/ddx1
     enddo

     ! x2
     do j=1,3
        durp=dur
        durp(2,j)=durp(2,j)+ddx2
        dfddu(2,j)=(ffun(token,x,y,z,ur,durp, nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
                   -ffun(token,x,y,z,ur,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/ddx2
     enddo

     ! p
     do j=1,3
        durp=dur
        durp(3,j)=durp(3,j)+ddp
        dfddu(3,j)=(ffun(token,x,y,z,ur,durp, nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax) &
                   -ffun(token,x,y,z,ur,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax))/ddp
     enddo
  endif
end subroutine dffun_ddu

