 !================================================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Huckel's method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!=================================================================================

!-------------------------Author: Carlos Cruz--------------------------------

!********************************************************************************
!* nC = number of Carbons  / nH = number of Hydrogens  / ne= number of electrons *
!* alpha, beta = energy parameters  / bond_distance = threshold for bonds        *
!* huckel_energy = total energy / occ = vector of occupancies                    *
!*   mul = vector containing mulliken analysis / H = hamiltonian                 *
!*  coord, coord C, coordH = matrices of coordinates / adj = adjacency matrix    *
!* Eigenvalues= matrix of eigenvalues / Eigenvectors = matrix of eigenvectors    *
!* BO= matrix of bonding orders / top,topC,topH = vector separating C and H      *
!*                                                                        	 *
!* In case of any doubt contact carloscruzmarin.bdn@gmail.com           	 *
!*                                                                      	 *
!********************************************************************************

program huckel
implicit none
integer::n,i,j,k,nC,nH,charge,ne
real, parameter:: alpha=-11.4 , beta = -0.8
real::bond_distance,huckel_energy
real,dimension(100) :: occ,mull
real,dimension(100,100):: H,coord,coordC,coordH,mat_dist,adj,Eigenvalues,Eigenvectors,BO
character(len=1),dimension(100):: top,topC,topH
character(len=100):: name,input, output
! H is the hamiltonian , nC the number of AO (C atoms) 

! Threshold for considering distance 
write(*,*) "Insert the name of the molecule you want to do the analysis: "
read(*,*) name
input = trim(name) // ".xyz"
output= trim(name) // ".txt"
bond_distance=1.58 !angstroms
open(unit=10, file=adjustl(trim(input)),status='old', action='read') ! Change here the name of the file before running the program
open(unit=20, file=adjustl(trim(output)),status='new', action='write')
open(unit=30, file="coords.txt" ,status="new" , action="write")
open(unit=40, file="coeffs.txt",status="new", action="write")


write(20,*) "Dear user: "
write(20,*) "The program will calculate the Huckel and Mulliken analysis for the desired molecule"
write(20,*) "For doing this, you will need to provide a .xyz file and insert the charge in the prompt"
write(20,*) "--------------------------------------------------------------------------------------------------------------"

write(*,*) "Charge of the molecule:"
read(*,*) charge
read(10,*) n
write(20,*) "Number of Atoms:" ,n
write(20,*) " "
write(20,*) "Charge of the molecule:" ,charge
write(20,*) " "
read(10,*) 

do i=1,n
  read(10,*) top(i),(coord(i,j),j=1,3) 
end do
!We separate the C and the H, because for our calculations we only need the C
nC=0; nH=0
do i=1,n
 if (top(i) == 'C') then
  nC=nC+1
  do j=1,3
   topC(nC)=top(i) 
   coordC(nC,j)=coord(i,j)
  end do
 else
 nH=nH+1
  do j=1,3
   topH(nH)=top(i) 
   coordH(nH,j)=coord(i,j)
  end do
 end if
end do 

do i=1,nC
  write(30,*) (coordC(i,j),j=1,3)
end do

call matr_dist(coordC,nC,mat_dist)
write(20,*) " "
write(20,*) "Distances between atoms :"
write(20,*) " "
do i=1,nC
 write(20,*) (mat_dist(i,j),j=1,nC)
end do

call adjacency_mat(mat_dist,adj,bond_distance,nC)
call hamiltonian (adj,H,nC)
call jacobi(H,Eigenvalues,Eigenvectors,nC)

! Computing the occupancy, where each pz orbital contrinutes with 1 electron and the charge has to be considered. 

ne = nC-charge !Calculation of the number of electrons
huckel_energy = 0
occ = 0
i=1
write(20,*) " "
write(20,*) "Number of pi - electrons: ",ne
do while (ne > 0)
  occ(i)= 2
  ne=ne-2
  i=i+1
end do

if (ne == 1) then
  occ(i+1) = 1
end if

ne = nC-charge
write(20,*) " "
write(20,*) "--------------------------------------------------------------------------------------------------------------"
write(20, '(A,10X,*(I2,:,"               "))') "Functions ψ", (i, i=1,nC)
write(20,*) "Occupancy:    ",(occ(i),i=1,nC)
write(20,*) "Energy:       ",(Eigenvalues(i,i),i=1,nC)
do i=1,nC
  write(20,*) "C",i," ",(Eigenvectors(i,j),j=1,nC)
end do

do i=1,ne
  huckel_energy=huckel_energy + (occ(i)) * Eigenvalues(i,i)
end do
write(20,*) "---------------------------------------------------------------------------------------------------------------"
write(20,*) "Huckel energy of the system is:", huckel_energy
write(20,*) "---------------------------------------------------------------------------------------------------------------"

!Mulliken population
mull = 0
do i=1,nC
 do j=1,nC
  mull(i) = mull(i)+ occ(j)*(Eigenvectors(i,j)**2)
 end do
end do

write(20,*) "---------------------------------------------------------------------------------------------------------------"
write(20,*) "                                        Mulliken analysis"
write(20,*) "---------------------------------------------------------------------------------------------------------------"
write(20,*) " "
write(20, '(A,10X,*(I2,:,"               "))') "Atoms", (i, i=1,nC)
write(20,*) "       " ,(mull(i),i=1,nC)
write(20,*) "---------------------------------------------------------------------------------------------------------------"

  ! Bond order
do i=1,nC
  do j=1,nC
    do k=1,nC
      BO(j,k)=BO(j,k) + occ(i)*Eigenvectors(j,i)*Eigenvectors(k,i)
    end do 
  end do
end do

write(20,*) "---------------------------------------------------------------------------------------------------------------"
write(20,*) "                                             Bond order"
write(20,*) "---------------------------------------------------------------------------------------------------------------"
write(20,*) " "
do i=1,nC
  do j=1,nC
   if (abs(BO(i,j))<10E-3) then
    BO(i,j) = 0
   end if
  end do
end do
do i=1,nC
   write(20,*) (BO(i,j),j=1,nC)
end do

do i=1,nC
  write(40,*) (Eigenvectors(i,j),j=1,nC)
end do

write(20,*) "The program finished grafecully"
close(30) ; close(40)
close(10) ; close(20)
end program huckel


subroutine matr_dist(A,n,dist)
implicit none
real,dimension(100,100):: A,dist
integer:: n,i,j,k
do i=1,n
 do j=1,n
  do k=1,3
   dist(i,j)=dist(i,j)+(A(i,k)-A(j,k))**2
  end do
   dist(i,j)=sqrt(dist(i,j))
 end do
end do
return
end subroutine matr_dist

subroutine adjacency_mat(A,B,limit_bonding,n)
implicit none
integer:: i,j,n
real,dimension(100,100):: A,B
real:: limit_bonding
do i=1,n
 do j=1,n
  if (A(i,j).lt.(limit_bonding)) then
    B(i,j)=1
  else 
    B(i,j)=0
  end if
 end do
end do
end subroutine adjacency_mat


subroutine hamiltonian(A,B,n)
implicit none
integer:: i,j,n
real,dimension(100,100):: A,B
real, parameter:: alpha=-11.4 , beta = -0.8
do i=1,n
 do j=1,n
  if (A(i,j)==1) then
    if (i==j) then
      B(i,j)=alpha
    else 
      B(i,j)=beta
    endif 
   else 
    B(i,j)=0
   end if
 end do
end do
end subroutine hamiltonian

subroutine jacobi(A,B,P_norm,n)
implicit none
integer:: i,j,n,u,v,it
real,dimension(100,100):: A,O,trans_O,B,P,C,P_norm,D,E
real,dimension(100):: norm
real:: maximum,theta
real,parameter:: pi=3.14159265358979323846

P=0
do i=1,n
 P(i,i)=1.0
end do
it=0
maximum = 0.1

!Here the Jacobi algorithm starts

do while (maximum>= 10E-4)
  !We restart the operations matrix
  it=it+1
  do i=1,n
    do j=1,n
      if (i== j) then
        O(i,j) = 1
      else
        O(i,j) = 0
      end if
    end do
  end do

  !We restart the maximum 
  maximum= 0.0
   
  !We find again the new matrix O finding the maximum value of the Hamiltonian
  do i = 1, n-1
    do j = i+1, n
       if (abs(A(i,j)) > maximum) then
        maximum = abs(A(i,j))
        u=i ; v=j
       endif
    enddo
  enddo
  !We create the new matrix of operation O
  
  if (A(u,u) == A(v,v)) then
    theta = pi/4
  else
   theta=abs(0.5*atan((2*A(u,v))/(A(u,u)-A(v,v))))
  end if
   O(u,v)=sin(theta)
   O(v,u)=-sin(theta)
   O(u,u)=cos(theta)
   O(v,v)=O(u,u)

  !We need the transpose of O
  trans_O = transpose(O)
  P=matmul(P,O)
  
  !Then we calculate transpose of O x A (Hamiltonian) and then this product x the Operations matrix (O)
  C=0
  C=matmul(trans_O,A)
  B=matmul(C,O)
   !And we save the new value obtained
  A=B
end do


do i=1,n
  do j=1,n
    if (abs(B(i,j)).lt.10E-4) then
      B(i,j)= 0
    end if
  end do
end do

! normalization just in case
do i=1,n
  do j=1,n
    norm(i)=norm(i)+(P(i,j))**2
    P_norm(i,j)=P(i,j)
  end do
  P_norm(i,:)=P_norm(i,:)/norm(i)
end do

! Now we have to sort the eigenvalues and then its respectives eigenvector

D=0 !D is a matrix of support for sorting eigenvalues
E=0 !E is a matrix of support for sorting eigenvectors depending on eigenvalues' order
do i=1,n
  do j=i,n
   if (abs(B(i,i)) < abs(B(j,j))) then 
    D(j,j) = B(j,j)
    B(j,j) = B(i,i)
    B(i,i) = D(j,j)

    E(:,j) = P_norm(:,j)
    P_norm(:,j) = P_norm(:,i)
    P_norm(:,i) = E(:,j)
   end if
  end do 
end do
return
end subroutine jacobi
