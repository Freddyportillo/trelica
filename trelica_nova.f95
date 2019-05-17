program trelica

    implicit none
    
     ! Este programa calcula as deformacoes de uma trelica calculando a matriz de !rigidez global com loop
     
     real, dimension (3,4) :: prop
     integer, dimension (3,4) :: mat_conect
     integer, dimension (3,2) :: cc
     real, dimension (3,2) :: fa
     real, dimension (3) :: alfa, kprop
     integer, dimension (3) :: gdl_livres, gdl_impostos
     real, dimension (2,2) :: matop, kprop2
     real, dimension (4,4) :: kele
     integer :: ii, nele, nnodes, ngdl
     integer, dimension (6,6) :: ident
     real, dimension (6,6) :: kglobal, kglobal1
     real, dimension (4,6) :: mat_transf
     real, dimension (6,4) :: transmat
     real :: f
     real, dimension (2,4) :: T
     real, dimension (4,2) :: Ttrans
     real*8, dimension (3,3) :: kll, invkll, kli, transkli, kii
     integer :: n
     real, dimension (3) :: u_liv, f_imp, sigma
     REAL, PARAMETER :: Pi = 3.1415927
     
     f= 100000
     prop = 0
     nele = 3
     nnodes = 3
     ngdl = 6
     prop(1,:) = (/32.3/10**4, 6.9*10**10, 2.54, 90.0/)
     prop(2,:) = (/38.7/10**4, 20.7*10**10, 2.54, 180.0/)
     prop(3,:) = (/25.8/10**4, 20.7*10**10, 3.59, 135.0/)
    
     mat_conect(1,:) = (/3,4,5,6/)
     mat_conect(2,:) = (/1,2,3,4/)
     mat_conect(3,:) = (/1,2,5,6/)
     
     cc(1,:) = (/2,0/)
     cc(2,:) = (/5,0/)
     cc(3,:) = (/6,0/)
     
     fa(1,:) = (/1.0,0.0/)
     fa(2,:) = (/3.0,-f*cos(pi/4)/)
     fa(3,:) = (/4.0,-f*sin(pi/4)/)
     
     
     !Construção das matrizes de rigides globais e elementares
    matop(1,:) = (/1,-1/)
    matop(2,:) = (/-1,1/)
    
    ident = 0
    mat_transf = 0
    kele = 0
    
    do ii=1,ngdl
        ident(ii,ii) = 1
    end do
    
    kele = 0
    kglobal = 0    
    kglobal1 = 0
    alfa = prop(:,4)*pi/180  !Alfa em rad
    
     do ii=1,nele 
        kprop = prop(ii,1)*prop(ii,2)/prop(ii,3)
        
        call Tmat(T,alfa(ii))
        kprop2 = kprop(ii)*matop
        Ttrans = transpose(T)
       
        call matmu1(Ttrans, kprop2, T, kele)
         
        mat_transf(1,:) = ident(mat_conect(ii,1),:)
        mat_transf(2,:) = ident(mat_conect(ii,2),:)
        mat_transf(3,:) = ident(mat_conect(ii,3),:)
        mat_transf(4,:) = ident(mat_conect(ii,4),:)
        transmat = transpose(mat_transf)
     
        call matmu2(transmat, kele, mat_transf,kglobal1)
    
        kglobal = kglobal + kglobal1
        end do
    
     print*, kglobal
     
    
    print*, '-----------------'
    
    
    gdl_livres = fa(:,1)
    gdl_impostos = cc(:,1)
    
    kll = kglobal(gdl_livres,gdl_livres)
    kli = kglobal(gdl_livres,gdl_impostos)
    kii = kglobal(gdl_impostos,gdl_impostos)
    n = 3.0
    print*, kll
    call inv(kll, invkll, n)
     u_liv = matmul(invkll,fa(:,2))*1000
            print*, 'os valosres dos deloscamentos no gdl livres sao:'
            print*,u_liv
    print*, 'os valores das reacoes sao:'
            transkli = transpose(kli)
            u_liv = u_liv/1000 !em m
            f_imp = matmul(transkli,u_liv)/1000 !em kN
            print*, f_imp
!        sigma = 0     
!     do i = 1,nele
!         sigma = prop(i,2) * u_liv(i)
     
     end program
     subroutine matmu1 (A,B,C,D)
      !programa para multiplicacao de 3 matrizes
  
  real, dimension (4,2) :: A, mul
  real, dimension (2,2) :: B
  real, dimension (2,4) :: C
  real, dimension (4,4) :: D
  
  mul = matmul(A,B)
  D = matmul (mul,C)
  
    end subroutine matmu1
    
    subroutine matmu2 (A,B,C,D)
      !programa para multiplicacao de 3 matrizes
  
  real, dimension (6,4) :: A, mul
  real, dimension (4,4) :: B
  real, dimension (4,6) :: C
  real, dimension (6,6) :: D
  
  mul = matmul(A,B)
  D = matmul (mul,C)
  
    end subroutine matmu2
    
     subroutine Tmat(T,alfa)
     
     real, dimension (2,4) :: T
     real :: alfa
     T = 0
     T(1,1) = cos(alfa)
     T(1,2) = sin(alfa)
     T(2,3) = cos(alfa)
     T(2,4) = sin(alfa)
     
     end subroutine Tmat
   subroutine inv(a,c,n)
!============================================================
! Inverse matrix
! Method: Based LU decomposition for Ax=b
! Marouf Abderahmane - Strasbourg University 2016 
!-----------------------------------------------------------
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

L=0.0
U=0.0
b=0.0

do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

do i=1,n
  L(i,i) = 1.0
end do
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inv
