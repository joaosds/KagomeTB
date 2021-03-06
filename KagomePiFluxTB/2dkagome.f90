module variables
 use ISO_FORTRAN_ENV
 integer :: L, rede, nz, coord, sublattice, natoms
 real(real64) :: t
 real(real64), allocatable :: H(:,:), kx(:), ky(:), dif(:), dif2(:),  dif3(:), energykspace(:) ! allowed
 integer(int32), allocatable :: vnx(:),vny(:), nneighbors(:,:), site1n(:,:)
 character*4 :: bc
 ! Variables from external function
 integer :: LDA, INFO
 integer ::   LWORK
 character, external :: JOBZ, UPLO
 real(real64), allocatable :: W(:)
 real(real64), allocatable :: P(:)
 real(real64), allocatable :: k(:)
 real(real64), allocatable :: WORK(:)
 real(real64),parameter :: pi = 4.0d0*DATAN(1.0d0) ! pi in double precision * error at last digit (1, not 2)
end module variables

program kagomePifluxTB
 use variables
 implicit none

 open(unit = 150, file = 'input.dat')

  read(150,*) L
  read(150,*) t, bc

  print*, L
  print*, t, bc

 close(150)

 nz = 4
 rede = 3*L*L
 sublattice = L*L
 LDA = rede
 LWORK = 64*LDA

 coord = 2 ! coordinates; for 2D x and y only

 allocate(energykspace(rede))
 allocate(H(LDA,rede), vnx(rede), vny(rede), nneighbors(rede,nz))
 allocate(W(rede), WORK(LWORK), ky(sublattice), kx(sublattice))
 allocate(dif(rede), dif2(rede), dif3(rede))
 allocate(site1n(nz,coord))

 energykspace = 0.0d0
 vnx = 0
 vny = 0

 H = 0.0d0
 dif = 0.0d0
 dif2 = 0.0d0

 ! W will return the eigenvalues in ascending order

 call twodlattice

 call DSYEV('v','U', rede, H, LDA, W, WORK, LWORK, INFO)

 call energyanalyt

 deallocate(H, vnx, vny, nneighbors, W, WORK, kx, ky, dif, dif2, dif3)

end program kagomePifluxTB

subroutine twodlattice

 ! B.D.: Constructs the Kagome Hamiltonian with appropriate first neighbor
 ! interactions.

 use ISO_FORTRAN_ENV
 use variables
 implicit none
 integer :: j, v, ix , iy, i
 integer :: Bright, Bdown, CupA, CupB, subB, subC

 call neighbors

 do i = 1, rede
  do v = 1, 3, 2 ! One neighbor at a1 and other at a2
    j = nneighbors(i,v)
    H(i,j) = -t
    H(j,i) = H(i,j)
   enddo
enddo

 select case (bc)

 case('MBCY')

   do i = 1, rede
    do v = 1, 3, 2 ! only one neighbor over each direction
      j = nneighbors(i,v)
      ix = abs(vnx(i)-vnx(j))
      iy = abs(vny(i)-vny(j))
      if((iy .gt. 2)) then ! impose ABC over ix or iy only
       H(i,j) = t
      endif
    H(j,i) = H(i,j)
    enddo
  enddo

case('MBCX')

  do i = 1, rede
   do v = 1, 3, 2 ! only one neighbor over each direction
     j = nneighbors(i,v)
     ix = abs(vnx(i)-vnx(j))
     iy = abs(vny(i)-vny(j))
     if((ix .gt. 2)) then ! impose ABC over ix or iy only
      H(i,j) = t
     endif
   H(j,i) = H(i,j)
   enddo
 enddo




 case('ABC')

  do i = 1, rede
   do v = 1, 3, 2 ! only one neighbor over each direction
     j = nneighbors(i,v)
     ix = abs(vnx(i)-vnx(j))
     iy = abs(vny(i)-vny(j))
     if(ix .gt. 2) then ! impose MBC = impose ABC for kagome 0 flux
      H(i,j) = t
     endif
     H(j,i) = H(i,j)
   enddo
 enddo

 end select


 !-------------------------- Pi flux state ------------------------------------!
 do i = 1, sublattice

    if (mod(vnx(i),2) .ne. 0) then ! odd x atom
     subB = i + sublattice
     subC = i + sublattice + sublattice
     Bright = nneighbors(subB,1) ! A right
     Bdown = nneighbors(subB,4)  ! C down/right
     H(subB,Bright) = -H(subB,Bright) ! pi flux
     H(subB,Bdown) = -H(subB,Bdown)   ! pi flux
     H(Bright,subB) = H(subB, Bright) ! impose matrix symmetry
     H(Bdown,subB) = H(subB,Bdown) ! impose matrix symmetry

     CupA = nneighbors(subC,1)  ! A up
     CupB = nneighbors(subC,3)  ! B up
     H(subC,CupA) = -H(subC,CupA) ! pi flux
     H(subC,CupB) = -H(subC,CupB) ! pi flux
     H(CupA,subC) = H(subC, CupA) ! impose matrix symmetry
     H(CupB,subC) = H(subC, CupB) ! impose matrix symmetry
    endif
enddo

 return
end subroutine twodlattice

subroutine neighbors

 ! B.D.: Constructs the Kagome lattice with appropriate neighbors;

 use variables
 implicit none
 integer :: ix, iy, lsite, temp, m, temp1

!     Lattice site spatial localization and nearest neighbor definitions.

 do ix = 1,L !running over A index, or unit cell index
  do iy = 1,L

  ! ------------------- Defining the triangular Unit Cell -------------------- !

   temp = lsite(ix,iy)            ! Unit cell label (or 3 atoms on each site)
   vnx(temp) = ix                  ! Unit cell x-coordinate for A
   vny(temp) = iy                  ! Unit cell y-coordinate for A

   vnx(temp+sublattice) = ix               ! Unit cell x-coordinate for b
   vny(temp+sublattice) = iy               ! Unit cell y-coordinate for b

   vnx(temp+sublattice+sublattice) = ix               ! Unit cell x-coordinate for C
   vny(temp+sublattice+sublattice) = iy               ! Unit cell y-coordinate for C

!  Consider the fact the nearest neighbors are actually different for the sublattices A, B and C

   call neighbors1a(ix,iy,site1n,coord, nz)

   ! --------------------- Neighbors along a1 ------------------------------- !

   do m = 1,nz/2 ! Define the B neighbors for sublattice A
    temp1 = lsite(site1n(m,1),site1n(m,2)) + sublattice
    nneighbors(temp,m) = temp1
   enddo

   ! --------------------- Neighbors along a2 ------------------------------- !

   do m = nz/2+1,nz ! Define the C neighbors for sublattice A
    temp1 = lsite(site1n(m,1),site1n(m,2)) + sublattice + sublattice
    nneighbors(temp,m) = temp1
   enddo

   call neighbors1b(ix,iy,site1n,coord, nz)

   ! --------------------- Neighbors along a1 ------------------------------- !

   do m = 1,nz/2 ! Define the A neighbors for sublattice B
    temp1 = lsite(site1n(m,1),site1n(m,2))
    nneighbors(temp+sublattice,m) = temp1
   enddo

   ! --------------------- Neighbors along a2 ------------------------------- !

   do m = nz/2+1,nz ! Define the C neighbors for sublattice B
    temp1 = lsite(site1n(m,1),site1n(m,2)) + sublattice + sublattice
    nneighbors(temp+sublattice,m) = temp1
   enddo

   call neighbors1c(ix,iy,site1n,coord, nz)

   do m = 1,nz/2 ! Define the A neighbors for sublattice C
    temp1 = lsite(site1n(m,1),site1n(m,2))
    nneighbors(temp+sublattice+sublattice,m) = temp1
   enddo

   do m = nz/2+1,nz ! Define the B neighbors for sublattice C
    temp1 = lsite(site1n(m,1),site1n(m,2)) + sublattice
    nneighbors(temp+sublattice+sublattice,m) = temp1
   enddo
  enddo
 enddo

 return
end subroutine neighbors

subroutine neighbors1a(ix,iy,site,coord1,nzt)
  use variables
  implicit none
  integer :: ix,iy,ix0,iy0,coord1, nzt
  integer :: site(nzt,coord1)

  ! ------------------------- Triangular sublaticce A ------------------------ !

  ! ------------------------------ B neighbors -----------------------------!

  ! B right neighbor
  iy0 = iy
  ix0 = ix
  site(1,1) = ix0
  site(1,2) = iy0

  ! B left neighbor
  ix0 = ix - 1
  if (ix0 .lt. 1) then
   ix0 = ix0 + l
  endif
  iy0 = iy
  site(2,1) = ix0
  site(2,2) = iy0

  ! ------------------------------ C neighbors -----------------------------!

  ! C up neighbor
  iy0 = iy
  ix0 = ix
  site(3,1) = ix0
  site(3,2) = iy0

  ! c down neighbor
  ix0 = ix
  iy0 = iy - 1
  if (iy0 .lt. 1) then
   iy0 = iy0 + L
  endif
  site(4,1) = ix0
  site(4,2) = iy0

  return
end subroutine neighbors1a

subroutine neighbors1b(ix,iy,site,coord1,nzt)
  use variables
  implicit none
  integer :: ix,iy,ix0,iy0,coord1, nzt
  integer :: site(nzt,coord1)

  ! ------------------------- Triangular sublaticce B ------------------------ !

  ! ------------------------------ A neighbors -----------------------------!

  ! A right neighbor
  ix0 = ix + 1
  if (ix0 .gt. L) then
   ix0 = ix0 - L
  endif
  iy0 = iy
  site(1,1) = ix0
  site(1,2) = iy0

  ! A left neighbor
  ix0 = ix
  iy0 = iy
  site(2,1) = ix0
  site(2,2) = iy0

  ! ------------------------------ C neighbors -----------------------------!

  ! C up neighbor
  iy0 = iy
  ix0 = ix
  site(3,1) = ix0
  site(3,2) = iy0

  ! C down neighbor
  ix0 = ix + 1
  if (ix0 .gt. L) then
   ix0 = ix0 - L
  endif
  iy0 = iy - 1
  if (iy0 .lt. 1) then
   iy0 = iy0 + L
  endif
  site(4,1) = ix0
  site(4,2) = iy0

  return
end subroutine neighbors1b

subroutine neighbors1c(ix,iy,site,coord1,nzt)
  use variables
  implicit none
  integer :: ix,iy,ix0,iy0,coord1, nzt
  integer :: site(nzt,coord1)

  ! ------------------------- Triangular sublaticce C ------------------------ !

  ! ------------------------------ A neighbors -----------------------------!

  ! A up neighbor
  iy0 = iy + 1
  ix0 = ix
  if (iy0 .gt. L) then
   iy0 = iy0 - L
  endif
  site(1,1) = ix0
  site(1,2) = iy0

  ! A down neighbor
  ix0 = ix
  iy0 = iy
  site(2,1) = ix0
  site(2,2) = iy0

  ! ------------------------------ B neighbors -----------------------------!

  ! B up neighbor

  ix0 = ix - 1
  if (ix0 .lt. 1) then
   ix0 = ix0 + L
  endif
  iy0 = iy + 1
  if (iy0 .gt. L) then
   iy0 = iy0 - L
  endif
  site(3,1) = ix0
  site(3,2) = iy0

  ! B down neighbor
  iy0 = iy
  ix0 = ix
  site(4,1) = ix0
  site(4,2) = iy0


  return
end subroutine neighbors1c

subroutine energyanalyt

 ! B.D.: Calculates the eigenvalues of the matrices in real space and builds
 ! the energy vector with the complete kagome spectrum for 0-flux;

 use variables
 implicit none
 !real(real64), allocatable :: energy(:), fk, f, f2, n1, n2, n3
 integer :: i

 open(100, file = 'dispanalyt.dat')

 call numspectrum

 do i = 1 ,rede
  !dif(i) = abs(energy(i)-W(i))
  !dif2(i) = abs(energy(i)-energykspace(i))
  dif3(i) = abs(energykspace(i)-W(i))
 enddo


 !write(100,*) 'Max dif between elements (Analytical and Numerical)', maxval(dif)
 !write(100,*)
 !write(100,*) 'Max dif between elements (Real and K-space)', maxval(dif2)
 !write(100,*)
 write(100,*) 'Max dif between elements (Real and K-space - Numerical)', maxval(dif3)
 write(100,*)
 write(100,*) '    ', 'Real Space' , '                ', 'K-Space'
 write(100,*)

 do i = 1, rede
  write(100,*) W(i), energykspace(i)
 end do

 close(100)

 deallocate(energykspace)
 return
end subroutine energyanalyt

subroutine numspectrum

 ! Fazer vetor de momentos pra ser ordenado com energia pra termos banda

 ! B.D.: Calculates the eigenvalues of the matrices in k-space and builds
 ! the energykspace vector with the complete kagome spectrum for 0-flux;

 use variables
 implicit none
 complex*16, allocatable :: Hnm(:,:), expk(:), Workn(:)
 integer :: kys, kxs, temp, jj, i, LRWORK, tp
 real*8 :: kxtemp, kytemp
 real(real64), allocatable :: RWORK(:), Wn(:)
 real(real64) :: f
 complex*16 :: zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0)

 natoms = 6
 LWORK = 64*natoms
 LRWORK = 3*natoms-2

 allocate(Hnm(natoms, natoms), Wn(natoms), Workn(LWORK), expk(natoms/3), RWORK(LRWORK))

 jj = 0

 select case(BC)

  case('PBC')

   f = 2.0d0*pi/real(L,8)

   do kys = 0, L-1
     do kxs = 0, L/2 -1 !doubled unit cell

      Hnm = zero
      expk = zero
      Wn = 0.0d0
      Workn = zero
      Rwork = 0.0d0
      temp = 0

      kxtemp = 2.0d0*(real(kxs,8)*f) ! PBC
      kytemp = real(kys,8)*f ! PBC
      expk(1) = cmplx(dcos(kxtemp), dsin(kxtemp), kind=real64)
      expk(2) = cmplx(dcos(kytemp), dsin(kytemp), kind=real64)  ! PBC

      ! -------------- Hamiltonian Construction  Upper triangle -------------- !
      ! We don't need to find the lower triangle
      do i = 1, natoms-1
       Hnm(i,i+1) = one
      enddo
      Hnm(3,4) = zero

      Hnm(1,3) = one - conjg(expk(2))
      Hnm(2,4) = -one
      Hnm(3,5) = -conjg(expk(1))*expk(2)
      Hnm(4,6) = one + conjg(expk(2))

      Hnm(1,5) = conjg(expk(1))
      Hnm(2,6) = -conjg(expk(2))

      ! -------------- Hamiltonian Diagonalization ------------------- !

      call ZHEEV('N','U', natoms, Hnm, natoms, Wn, Workn, LWORK, RWORK, INFO)

      do tp = 1, natoms
      jj = jj +1
        energykspace(jj) = Wn(tp)
      enddo

     end do
   end do

  case('ABC')

   f = pi/real(L,8)


      do kys = 0, L-1
        do kxs = 0, L/2-1 !doubled unit cell

         Hnm = zero
         expk = zero
         Wn = 0.0d0
         Workn = zero
         Rwork = 0.0d0
         temp = 0

         kxtemp = 2.0d0*(2.0d0*real(kxs, 8)+1.0d0)*f ! ABC
         kytemp = (2.0d0*real(kys, 8)+1.0d0)*f  ! ABC
         expk(1) = cmplx(dcos(kxtemp), dsin(kxtemp), kind=real64)
         expk(2) = cmplx(dcos(kytemp), dsin(kytemp), kind=real64)  ! PBC

         ! -------------- Hamiltonian Construction  Upper triangle -------------- !
         ! We don't need to find the lower triangle
         do i = 1, natoms-1
          Hnm(i,i+1) = one
         enddo
         Hnm(3,4) = zero

         Hnm(1,3) = one - conjg(expk(2))
         Hnm(2,4) = -one
         Hnm(3,5) = -conjg(expk(1))*expk(2)
         Hnm(4,6) = one + conjg(expk(2))

         Hnm(1,5) = conjg(expk(1))
         Hnm(2,6) = -conjg(expk(2))
         Hnm = -1.0d0*Hnm
         ! -------------- Hamiltonian Diagonalization ------------------- !

         call ZHEEV('N','U', natoms, Hnm, natoms, Wn, Workn, LWORK, RWORK, INFO)

         do tp = 1, natoms
           jj = jj + 1
           energykspace(jj) = Wn(tp)
         enddo

        end do
      end do

 case('MBCX')

   f = 2.0d0*pi/real(L,8)

   do kys = 0, L-1
     do kxs = 0, L/2-1 !doubled unit cell

      Hnm = zero
      expk = zero
      Wn = 0.0d0
      Workn = zero
      Rwork = 0.0d0
      temp = 0

      kxtemp = 2.0d0*(2.0d0*real(kxs, 8)+1.0d0)*f/2.0d0 ! ABC
      kytemp = real(kys, 8)*f  ! PBC
      expk(1) = cmplx(dcos(kxtemp), dsin(kxtemp), kind=real64)
      expk(2) = cmplx(dcos(kytemp), dsin(kytemp), kind=real64)  ! PBC

      ! -------------- Hamiltonian Construction  Upper triangle -------------- !
      ! We don't need to find the lower triangle
      do i = 1, natoms-1
       Hnm(i,i+1) = one
      enddo
      Hnm(3,4) = zero

      Hnm(1,3) = one - conjg(expk(2))
      Hnm(2,4) = -one
      Hnm(3,5) = -conjg(expk(1))*expk(2)
      Hnm(4,6) = one + conjg(expk(2))

      Hnm(1,5) = conjg(expk(1))
      Hnm(2,6) = -conjg(expk(2))

      Hnm = -1.0d0*Hnm
      ! -------------- Hamiltonian Diagonalization ------------------- !

      call ZHEEV('N','U', natoms, Hnm, natoms, Wn, Workn, LWORK, RWORK, INFO)

      do tp = 1, natoms
        jj = jj + 1
        energykspace(jj) = Wn(tp)
      enddo

     end do
   end do

 case('MBCY')

   f = 2.0d0*pi/real(L,8)

   do kys = 0, L-1
     do kxs = 0, L/2-1 !doubled unit cell

      Hnm = zero
      expk = zero
      Wn = 0.0d0
      Workn = zero
      Rwork = 0.0d0
      temp = 0

      kxtemp = 2.0d0*(real(kxs, 8))*f    ! PBC
      kytemp = (2.0d0*real(kys, 8)+1.0d0)*f/2.0d0  ! ABC
      expk(1) = cmplx(dcos(kxtemp), dsin(kxtemp), kind=real64)
      expk(2) = cmplx(dcos(kytemp), dsin(kytemp), kind=real64)  ! PBC

      ! -------------- Hamiltonian Construction  Upper triangle -------------- !
      ! We don't need to find the lower triangle
      do i = 1, natoms-1
       Hnm(i,i+1) = one
      enddo
      Hnm(3,4) = zero

      Hnm(1,3) = one - conjg(expk(2))
      Hnm(2,4) = -one
      Hnm(3,5) = -conjg(expk(1))*expk(2)
      Hnm(4,6) = one + conjg(expk(2))

      Hnm(1,5) = conjg(expk(1))
      Hnm(2,6) = -conjg(expk(2))
      Hnm = -1.0d0*Hnm
      ! -------------- Hamiltonian Diagonalization ------------------- !

      call ZHEEV('N','U', natoms, Hnm, natoms, Wn, Workn, LWORK, RWORK, INFO)

      do tp = 1, natoms
        jj = jj + 1
        energykspace(jj) = Wn(tp)
      enddo

     end do
   end do
   end select

  call Order(energykspace, rede)
  deallocate(Wn,Hnm,expk, Workn, RWORK)

 return
end subroutine numspectrum

integer function lsite(x,y)
  !B.D.: index of unit cell
 use variables
 implicit none
 integer :: x,y

 lsite = x+(y-1)*L

end function

subroutine Order(matrix,dims)

 ! B.D.: Order an Array in Ascending order by the bubble sort method.

 use ISO_FORTRAN_ENV
 implicit none
 integer :: INC, dims, i, J
 real(real64) :: vny, matrix(dims)

 INC = 1
 do while(INC .LE. dims)
   INC = 3*INC+1
 enddo
 do while(INC .GT. 1)
  INC = INC/3
  do I = INC+1,dims
   vny = matrix(I)
   J = I
   do while(matrix(J-INC) .GT. vny)
    matrix(J) = matrix(J-INC)
    J = J-INC
    if (J .LE. INC) then
     exit
    endif
   enddo
   matrix(J)=vny
  enddo
 enddo

 return
end subroutine Order

subroutine pbctest(ind)
 use variables
 implicit none
 integer :: ind

 if (ind .gt. L) then
  ind = ind - L
 elseif (ind .lt. 1) then
  ind = ind + L
 endif

 return
end subroutine pbctest
