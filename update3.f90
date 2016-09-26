! ######################################################################

  subroutine update3

! ----------------------------------------------------------------------
! For when we've stored the Q matrix using tree [not currently used]
! ----------------------------------------------------------------------

    use common
    
    implicit none
    
    integer iQ1,jQ1,iQp,jQp
    double precision tot,Qval1,fac
    
    if (key_sim.eq.1) then
       do iQp=1,leniQ
          iQ1=iQ(iQp)
          
          if (iunk(iQ1).eq.1) then
             
             tot=0
             jQp=data(iQp)
             do
                if (jQ(jQp).eq.0) then
                   exit
                else if (jQ(jQp).lt.0) then
                   jQp=abs(jQ(jQp))
                endif
                
                jQ1=jQ(jQp)
                Qval1=Qval(jQp)
                
                if (jQ1.eq.iQ1) then
                   fac=Qval1
                else
                   tot=tot-Qval1*vec(jQ1)
                endif
                
                jQp=jQp+1
             enddo
             vec2(iQ1)=(tot+RR(iQ1))/fac
          else
             vec2(iQ1)=vec(iQ1)
          endif
       enddo
    else
       do iQp=1,leniQ
          iQ1=iQ(iQp)
          
          if (iunk(iQ1).eq.1) then
             
             tot=0
             jQp=data(iQp)
             do
                if (jQ(jQp).eq.0) then
                   exit
                else if (jQ(jQp).lt.0) then
                   jQp=abs(jQ(jQp))
                endif
              
                jQ1=jQ(jQp)
                Qval1=Qval(jQp)
                
                tot=tot-Qval1*vec(jQ1)
                
                jQp=jQp+1
             enddo
             vec2(iQ1)=vec(iQ1)+tstep*(tot+RR(iQ1))/kfac(iQ1)
          else
             vec2(iQ1)=vec(iQ1)
          endif
       enddo
    endif
  end subroutine update3

