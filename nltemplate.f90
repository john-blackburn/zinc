function cfun(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

   character(*) label
   integer nvar,istep,imax,jmax,kmax
   double precision cfun,x,y,z,ur(nvar),dur(nvar,3)
   integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
   double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

! Your code goes here

end function cfun

! ######################################################################
     
function afun(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

   character(*) label
   integer nvar,istep,imax,jmax,kmax
   double precision afun,x,y,z,ur(nvar),dur(nvar,3)
   integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
   double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

! Your code goes here

end function afun

! ######################################################################

function ffun(label,x,y,z,ur,dur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

   character(*) label  
   integer nvar,istep,imax,jmax,kmax
   double precision x,y,z,ffun,ur(nvar),dur(nvar,3)
   integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
   double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

! Your code goes here

end function ffun

! ######################################################################

function qfun(label,x,y,z,nx,ny,nz,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

   character(*) label
   double precision qfun,x,y,z,nx,ny,nz,ur(nvar)
   integer nvar,istep,imax,jmax,kmax
   integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
   double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

! Your code goes here

end function qfun

! ######################################################################
     
function gfun(label,x,y,z,nx,ny,nz,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

   character(*) label
   double precision gfun,x,y,z,nx,ny,nz,ur(nvar)
   integer nvar,istep,imax,jmax,kmax
   integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
   double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

! Your code goes here

end function gfun

! ######################################################################

function BCfun(label,x,y,z,ur,nvar,istep, &
          ireg,iregup,rnode,vec,imax,jmax,kmax)

   character(*) label
   double precision BCfun,x,y,z,ur(nvar)
   integer nvar,istep,imax,jmax,kmax
   integer ireg(0:imax,0:jmax,0:kmax),iregup(0:imax,0:jmax,0:kmax)
   double precision rnode(0:imax,0:jmax,0:kmax,3),vec(*)

! Your code goes here

end function BCfun

! ######################################################################

function scanfun(label,x,y,z,nx,ny,nz,ur,dur,nvar)

  character(*) label
  integer nvar
  double precision scanfun,x,y,z,nx,ny,nz,ur(nvar),dur(nvar,3)

! Your code goes here

end function scanfun
