! Part of Zinc FE package. Author: John Blackburn

module tree

! insertQ
! addentry

contains

subroutine insertQ(nentry,iQx,jQx,Qvalx)

! insert nentry matrix elements into tree. All have same row value, iQx
! they have different col values jQx(ientry). Actual value is Qvalx(ientry)

  use common, only : Qval,iQ,jQ,left,right,data,leniQ,lenjQ
  implicit none

  integer nentry,iQx,jQx(8)
  real(8) :: Qvalx(8)

  integer iQp

  if (leniQ.eq.0) then
     leniQ=1
     iQ(1)=iQx
     call addentry(1,nentry,iQx,jQx,Qvalx)
  else
     
     iQp=1
     do
        if (iQ(iQp).eq.iQx) then
           call addentry(iQp,nentry,iQx,jQx,Qvalx)
           exit
        endif
        
        if (iQx.gt.iQ(iQp)) then
           if (right(iQp).eq.0) then
              leniQ=leniQ+1
              iQ(leniQ)=iQx
              right(iQp)=leniQ
              call addentry(leniQ,nentry,iQx,jQx,Qvalx)
              exit
           else
              iQp=right(iQp)
           endif
        else if (iQx.lt.iQ(iQp)) then
           if (left(iQp).eq.0) then
              leniQ=leniQ+1
              iQ(leniQ)=iQx
              left(iQp)=leniQ
              call addentry(leniQ,nentry,iQx,jQx,Qvalx)
              exit
           else
              iQp=left(iQp)
           endif
           
        endif
     enddo
  endif
end subroutine insertQ
     
! ######################################################################

subroutine addentry(iQp,nentry,iQx,jQx,Qvalx)

  use common, only : Qval,iQ,jQ,data,lenjQ
  implicit none

  integer iQp,nentry,iQx,jQx(8)
  integer idone(8),ientry,jQp,ndone
  real(8) Qvalx(8)

  idone=0
  ndone=0

  if (data(iQp).eq.0) then
     
     data(iQp)=lenjQ+1
     do ientry=1,nentry
        lenjQ=lenjQ+1
        jQ(lenjQ)=jQx(ientry)
        Qval(lenjQ)=Qvalx(ientry)
     enddo
     
     lenjQ=lenjQ+1
     jQ(lenjQ)=0
  else
     jQp=data(iQp)
     do
        if (jQ(jQp).lt.0) then
           jQp=abs(jQ(jQp))
        else if (jQ(jQp).eq.0) then
           jQ(jQp)=-(lenjQ+1)
           do ientry=1,nentry
              if (idone(ientry).eq.0) then
                 lenjQ=lenjQ+1
                 jQ(lenjQ)=jQx(ientry)
                 Qval(lenjQ)=Qvalx(ientry)
              endif
           enddo
           lenjQ=lenjQ+1
           jQ(lenjQ)=0
           
           exit
        endif
        
        do ientry=1,nentry
           if (idone(ientry).eq.0) then
              if (jQ(jQp).eq.jQx(ientry)) then
                 Qval(jQp)=Qval(jQp)+Qvalx(ientry)

                 idone(ientry)=1
                 ndone=ndone+1
              endif
           endif
        enddo

        if (ndone.eq.nentry) exit
        
        jQp=jQp+1
     enddo
  endif
end subroutine addentry

end module tree
