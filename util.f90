! Part of Zinc FE package. Author: John Blackburn

module util

implicit none

contains

! ######################################################################

  function mag(r)
    double precision mag,r(3)
    mag=sqrt(r(1)**2+r(2)**2+r(3)**2)
  end function mag

! ######################################################################
  
  function cross(a,b)
    
    ! a1 a2 a3
    ! b1 b2 b3
    
    double precision cross(3),a(3),b(3)
    
    cross(1)=  a(2)*b(3)-b(2)*a(3)
    cross(2)=-(a(1)*b(3)-b(1)*a(3))
    cross(3)=  a(1)*b(2)-b(1)*a(2)
  end function cross

! ######################################################################
  
  function i2c(i)
    integer i
    character i2c*3
    
    write (i2c,'(i3.3)') i
    
  end function i2c

end module util
