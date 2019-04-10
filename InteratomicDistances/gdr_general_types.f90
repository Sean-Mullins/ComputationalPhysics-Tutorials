! Calculates the pair distribution function for a Ag-Au cluster
! Writes the numbers per interval as well as the gdr convoluted with 
! a Gaussian;
! This is done seperately for {all;AgAg,AgAu,AuAu} bonds
!
! the normalization of the Gauss-convoluted output is fine, i.e., 
! integrating over the full range gives the number of all atom pairs 
! (note that here, no cutoff is applied, unlike for the counting of
! the bonds)

! written based on the little code by Christin Mottet, 
! 21.04.2012 Hansi Weissker

program gdr_new
  
  implicit none
  
  integer :: natoms, i, j, typeflag
!  integer :: nij!, nijAgAg, nijAgAu, nijAuAu

  integer file_count

  integer :: l,k, ntypes
!  integer :: nij(10,10)
!  real :: dijsum(10,10)
!  character (len=2) :: Spec(10)
  integer, dimension(:,:), allocatable  :: nij
  real, dimension(:,:), allocatable :: dijsum
  character (len=2), dimension(:), allocatable :: Spec




  real, dimension(:,:), allocatable :: dij,dijtrans,dijAgAg,dijAgAu,dijAuAu
  real :: dijcut
!  real :: dijsum, dijsumAgAg, dijsumAgAu, dijsumAuAu
  real, dimension(:), allocatable :: x, y, z
  character (len=80) :: cdummy,infile
  character (len=2),  dimension(:), allocatable :: type(:)

  integer :: nstep
  real :: d1,d2,dmax,dinterval
  integer, dimension(:), allocatable :: nbox
  real :: sigma, sigma2
  real :: xc, xpc, gauss
  real, dimension(:), allocatable :: gdr
  integer :: ix, ixp


  real, dimension(:), allocatable :: Mass
  real :: Msum, masstemp, r
  real :: Cmass(3),  vsum(3), vtemp(3)
  real, dimension(:,:), allocatable :: rdist
  real, dimension(:), allocatable :: rtrans
  integer, dimension(:), allocatable :: icount
    real, dimension(:), allocatable :: rbonddist

!---------------------------------------------------------------
! changeable input  
  ntypes = 4
  allocate(Spec(ntypes))
  allocate(Mass(ntypes))
  Spec(1) = "Au"
  Mass(1) = 196.97
  Spec(2) = "S"
  Mass(2) = 32.06
  Spec(3) = "C"
  Mass(3) = 12.011
  Spec(4) = "H"
  Mass(4) = 1.0079
! Masses:
! 'Ag'  | 107.8682 
! 'Au'  | 196.97   
!  'H'  | 1.0079   
!  'C'  | 12.011   
!  'S'  | 32.06    

  dijcut =  3.2 !! this determines what a "bond" is; it does
                !! not act in the calculation of the gdr
  dmax = 20 ! max interatomic distance (quasi infinity for small clusters)
  nstep = 20000
  sigma = .001 ! sigma of the (normalized) Gaussian
  infile = 'Au60-144TD-40A-550Ry-ST-SP-Gamma-PDOS.xyz'
!  infile='in.xyz'
!---------------------------------------------------------------
  

  allocate(rdist(ntypes,natoms))
  allocate(rtrans(natoms))
  allocate(rbonddist(nstep))

  allocate(icount(ntypes))
  allocate(gdr(nstep))
  allocate(nbox(nstep))
  open(12,file=infile)

  read(12,*) natoms
  read(12,*) ! haven't understood that -- like this it works even when
  ! there is no entry in this line.
  write(*,*) 'Found',natoms,'atoms:'

   allocate(nij(ntypes,ntypes))
   allocate(dijsum(ntypes,ntypes))


  allocate(x(natoms))
  allocate(y(natoms))
  allocate(z(natoms))
  allocate(dij(natoms,natoms))

  allocate(dijAgAg(natoms,natoms))
  allocate(dijAgAu(natoms,natoms))
  allocate(dijAuAu(natoms,natoms))

  allocate(dijtrans(natoms,natoms))

  allocate(type(natoms))
  x=0.d0; y=0.d0; z=0.d0; dij=0.d0; dijsum=0.d0
!  dijsumAgAg=0.d0;  dijsumAgAu=0.d0;  dijsumAuAu=0.d0;
  nij=0  !; nijAgAg=0; nijAgAu=0; nijAuAu=0

! read coords from input file
  do i=1,natoms    
     read(12,*) type(i),x(i),y(i),z(i)
     write(*,20) i, type(i), x(i), y(i), z(i)
  enddo
20 format(i5,a5,3f12.8)


! calculate center of mass
  Vsum = 0.d0
  Msum = 0.d0
  do i = 1, natoms
     do k = 1,ntypes
        if (type(i).eq.Spec(k))   then 
           masstemp = Mass(k)
           !           write(*,*) ' test' , masstemp
        end if
     end do
     vtemp = (/ x(i),y(i),z(i) /)     
!     write(*,*) ' test', vtemp, masstemp
     
     Vsum = Vsum + masstemp * vtemp
     Msum = Msum + masstemp
  enddo
  Cmass = Vsum / Msum
  write(*,*) 'Cmass =' , Cmass
  
! calc distances from center of mass
  icount = 0
  do i = 1, natoms
     do  k = 1,ntypes
        if (type(i).eq.Spec(k)) then
           r = sqrt((x(i)-Cmass(3))**2 + (y(i)-Cmass(2))**2 + (z(i)-Cmass(1))**2)
           !           write(500+k,*) i, r  ! OK, that works
           icount(k) = icount(k) + 1
           rdist(k,icount(k)) = r
           
        end if
     end do
  end do
  
!------------------------------------------------------------
  do k = 1,ntypes
     do i = 1, icount(k)
!!!        write(*,*) 'check... why do they start with zero? '
        write(600+k,*) i, rdist(k,i)  ! OK, that works
     end do
  end do

!-------------------------------------------------------
  file_count = 1
  do k = 1, ntypes

     write(*,*) 'Species:', Spec(k)!, Spec(l)
     
     rtrans = 0.d0

     cdummy = 'rbonddist_'//TRIM(Spec(k))//'.dat'
     write(*,*) cdummy
     open(700+file_count,file=cdummy)
     write(700+file_count,*) '#', cdummy

!     rtrans(:) = rdist(k,:)
     do i = 1, natoms
        rtrans(i) = rdist(k,i)
        write(*,*) rtrans(i),'this is my little test for type', k, i 
    end do

!- OK now for rtrans
     do i = 1, natoms
        write(*,*) 'check',i, rdist(k,i)  !rtrans
!        write(900+k,*) rdist(k,i)  !rtrans
        write(900+k,*) rtrans(i)  !rtrans
     end do

     call calc_rdist(dmax,nstep,rtrans,type,sigma,natoms,nbox,rbonddist)
     dinterval = dmax / float(nstep)
     
!     do i = 1,nstep
!        write(999,*) i*dinterval, rbonddist(i)
!     end do

     write(*,*) icount(k), file_count
     
     do i = 1, nstep
!        write(*,*) 'i=',i ,700+file_count
        !         write(60+typeflag,*) dinterval*i, nbox(i)
        write(700+file_count,*)  i*dinterval, rbonddist(i)
     enddo
     
     close(700+file_count)
     file_count = file_count + 1
     write(*,*) "file_count=",file_count
     !   end do
  end do
  

!-- calculate all distances dij between all atoms 
  do i=1,natoms-1
     do j=i+1,natoms
        dij(i,j)=sqrt((x(j)-x(i))**2 + (y(j)-y(i))**2 + (z(j)-z(i))**2)
        !- all bonds
!        if (dij(i,j).lt.dijcut) then
!           dijsum=dijsum+dij(i,j)
!           nij=nij+1
!        end if        
     enddo
  enddo

! careful: k,l refer to the type list (Au,S,C,H...), while
! i,j are the types of each individual atom
  write(*,*) 'Cut-off for dij is ', dijcut,'A or whatever unit is in your xyz file'
  do k = 1, ntypes
     do l = k, ntypes-1
        do i=1,natoms-1
           do j=i+1,natoms
              if ((type(i).eq.Spec(k).AND.type(j).eq.Spec(l)).OR.  &
                   &              (type(i).eq.Spec(l).AND.type(j).eq.Spec(k))) then
                 !                 dijAgAu(i,j)= dij(i,j)
!                 write(*,*) Spec(k), Spec(l)
! determine where "bonds" are; can do this later type specific
                 if (dij(i,j).lt.dijcut) then
                    !                    nijAgAu = nijAgAu + 1
                    nij(k,l) = nij(k,l) + 1
                    dijsum(k,l)=dijsum(k,l)+dij(i,j)
                 endif
              endif
           end do
        end do
     end do
  end do



! just for test
  do k = 1, ntypes
     do l = k, ntypes-1
        write(*,*) nij(k,l), dijsum(k,l)/float(nij(k,l)),Spec(k), Spec(l)
     end do
  end do
  
write(*,*) Spec(1),Spec(2),Spec(3),Spec(4)
write(*,*) Spec(1),Spec(2),trim(Spec(3)),Spec(4)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
file_count = 1
do k = 1, ntypes
   do l = k, ntypes

      write(*,*) 'Species:', Spec(k), Spec(l)

      dijtrans = 0

      cdummy = 'dist_'//TRIM(Spec(k))//'_'//TRIM(Spec(l))//'.dat'
      write(*,*) cdummy
      open(300+file_count,file=cdummy)
!      open(300+file_count)
      write(300+file_count,*) '#', cdummy
      do i=1,natoms-1
         do j=i+1,natoms
            if( (type(i).eq.Spec(l).AND.type(j).eq.Spec(k)).OR.    &
                 & (type(j).eq.Spec(l).AND.type(i).eq.Spec(k))  ) then
               dijtrans(i,j) = dij(i,j)
            end if
         end do
      end do

      call calc_gdr(dmax,nstep,dijtrans,type,sigma,natoms,nbox,gdr)
      dinterval = dmax / float(nstep)
  
      do i = 1, nstep
!         write(60+typeflag,*) dinterval*i, nbox(i)
         write(300+file_count,*)  i*dinterval, gdr(i)
      enddo

      close(300+file_count)
      file_count = file_count + 1
      write(*,*) "file_count=",file_count
   end do
end do


end program gdr_new


!!!  open(61,file='out_count_all_bonds.dat')
!!!  open(62,file='out_count_AgAg_bonds.dat')
!!!  open(63,file='out_count_AgAu_bonds.dat')
!!!  open(64,file='out_count_AuAu_bonds.dat')
!!!  open(71,file='out_convoluted_all_bonds.dat')
!!!  open(72,file='out_convoluted_AgAg_bonds.dat')
!!!  open(73,file='out_convoluted_AgAu_bonds.dat')
!!!  open(74,file='out_convoluted_AuAu_bonds.dat')
!!!
!!!  do typeflag = 1,4 ! corresponds to (all, AgAg, AgAu, AuAu)
!!!     select case (typeflag)
!!!     case (1)
!!!        dijtrans = dij
!!!     case (2)
!!!        dijtrans = dijAgAg
!!!     case (3)
!!!        dijtrans = dijAgAu
!!!     case (4)
!!!        dijtrans = dijAuAu
!!!     case DEFAULT
!!!        stop 'something wrong'
!!!     end select
!!!     
!!!     call calc_gdr(dmax,nstep,dijtrans,type,sigma,natoms,nbox,gdr)
!!!     dinterval = dmax / float(nstep)
!!!  
!!!     do i = 1, nstep
!!!        write(60+typeflag,*) dinterval*i, nbox(i)
!!!        write(70+typeflag,*)  i *dinterval, gdr(i)
!!!     enddo
!!!  end do
!!!
!!!
!!!
!!!!!write(*,*)   dijAuAu
!!!
!!!
!!!
!!!end program gdr_new

!-----------------------------------------------------------------------------

subroutine calc_gdr(dmax,nstep,dij,type,sigma,natoms,nbox,gdr)
  
  implicit none
  real :: dmax, sigma
  integer :: nstep, natoms
  real :: dij(natoms,natoms)
  character (len=2)   :: type(natoms)
  
  real :: dinterval,d1,d2,gauss,xc,xpc,sigma2
  integer :: nbox(nstep)
  real :: gdr(nstep)
  integer :: ix, ixp, i,j,k

  REAL, PARAMETER :: Pi = 3.1415927

  dinterval = dmax / float(nstep)
  d1 = dinterval * 0.5d0
  d2 = dinterval * 1.5d0
  nbox=0
  do k = 1, nstep
     do i = 1, natoms - 1
        do j = i, natoms
!       write(*,*) d1, d2
           if(dij(i,j).gt.d1.AND.dij(i,j).le.d2) then
              nbox(k) = nbox(k)+1
           endif
        enddo
     enddo
     d1 = d1 + dinterval
     d2 = d2 + dinterval
! now in main:     write(66,*) dinterval*k, nbox(k)
  enddo

  sigma2 = sigma*sigma

  gdr=0.d0
  
  do ix = 1,nstep
     xc = ix *dinterval !+ dinterval*0.5d0
     !     write(*,*) 'xc=',xc
     do ixp = -nstep, nstep
        xpc = ixp*dinterval !+ dinterval*0.5d0
        !        write(*,*) 'xpc=',xpc
        gauss = 1/(sigma*SQRT(2.d0*Pi))*exp(-0.5d0*((xpc-xc)/sigma)**2)
        if(xpc.gt.0.AND.xpc.lt.dmax) then
           !          write(*,*) 'ixp=',ixp
           gdr(ixp) = gdr(ixp) + nbox(ix)*gauss 
!                      write(*,*) ixp, nbox(ix),gauss, gdr(ixp)
        endif
     enddo
  enddo

! now moved to main:  
!  do ix = 1, nstep
!     !     write(77,*)  ix *dinterval + dinterval*0.5d0, gdr(ix)
!     write(77,*)  ix *dinterval, gdr(ix)
!  enddo
  
  return
  
end subroutine calc_gdr

!-------------------------------------------------------

!-----------------------------------------------------------------------------

subroutine calc_rdist(dmax,nstep,rdist,type,sigma,natoms,nbox,rbonddist)
  
  implicit none
  real :: dmax, sigma
  integer :: nstep, natoms
  !  real :: dij(natoms,natoms)
  real :: rdist(natoms)
  character (len=2)   :: type(natoms)
  
  real :: dinterval,d1,d2,gauss,xc,xpc,sigma2
  integer :: nbox(nstep)
  !  real :: gdr(nstep)
  real :: rbonddist(nstep)
  integer :: ix, ixp, i,j,k
  
  REAL, PARAMETER :: Pi = 3.1415927

  dinterval = dmax / float(nstep)
  d1 = dinterval * 0.5d0
  d2 = dinterval * 1.5d0
  nbox=0
  do k = 1, nstep
     do i = 1, natoms !- 1
!        do j = i, natoms
           if(rdist(i).gt.d1.AND.rdist(i).le.d2) then
              nbox(k) = nbox(k)+1
           endif
        enddo
!     enddo
     d1 = d1 + dinterval
     d2 = d2 + dinterval
! hier stimmt es noch:   write(66,*) dinterval*k, nbox(k)
  enddo

  sigma2 = sigma*sigma

!  gdr=0.d0
  rbonddist = 0.d0
 
  do ix = 1,nstep
     xc = ix *dinterval 
     do ixp = -nstep, nstep
        xpc = ixp*dinterval
        gauss = 1/(sigma*SQRT(2.d0*Pi))*exp(-0.5d0*((xpc-xc)/sigma)**2)
        if(xpc.gt.0.AND.xpc.lt.dmax) then
!           gdr(ixp) = gdr(ixp) + nbox(ix)*gauss 
            rbonddist(ixp) = rbonddist(ixp) + nbox(ix)*gauss 
        endif
     enddo
  enddo

  return
  
end subroutine calc_rdist


