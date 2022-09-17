program huckel
  implicit none
  integer::a,b,i,j,nsiti,ne, lda, lwork,nl,n,info,k
  real*8::ac,bc,c,d, energy,q,bo, kxy, check,mux, muy, muz, risene, energygs,domega,intensita,sigma,pi,nomega, omegastart,omegaend,omega
  character*1::jobz, uplo
  integer,allocatable:: bond(:,:)
  real*8,allocatable::hx(:),nz(:),ham(:,:), w(:), work(:), coords(:,:), charge(:), x(:,:), y(:,:), z(:,:)
  ac=0
  bc=-8.760442803
  pi =dacos(-1.d0)

  open(1,file='input.dat')
  open(2,file='output.dat')

  !reading input values and allocating operating matrix
  read(1,*) nsiti, ne, nl
  allocate(hx(nsiti),nz(nsiti),ham(nsiti,nsiti),bond(nl,2), coords(nsiti,3), charge(nsiti),x(nsiti-ne/2,ne/2),y(nsiti-ne/2,ne/2),z(nsiti-ne/2,ne/2))
  do i=1,nsiti
     read(1,*) a, hx(i), nz(i), coords(i,1), coords(i,2), coords(i,3)
     ham(i,i)=ac + hx(i)*bc
  enddo

  do i=1,nl
     read(1,*) a,b, kxy
     ham(a,b)= kxy*bc
     bond(i,1)=a
     bond(i,2)=b
   enddo

  do i=1,nsiti
     do j=1,nsiti
        ham(j,i)=ham(i,j)
     enddo
  enddo
  
  write(2,*) 'HAMILTONIAN'
  do i=1,nsiti
     write(2,'(<nsiti>(2x,f10.5))')(ham(i,j),j=1,nsiti)
  enddo
  write(2,*) 

  jobz='V'
  uplo='U'
  n=nsiti
  lda=nsiti
  lwork=3*N-1
  allocate(w(n), work(lwork))
  call dsyev (jobz,uplo,n,ham,lda,w,work,lwork,info)
  if(info.eq.0) write (*,*) 'successful diagonalization'

  write(2,*) 'EIGENVALUES'
  do i=1,nsiti
     write(2,*) i, w(i)
  enddo
  energy=0.d0
  do i=1,7
     energy=energy+2*w(i)
  enddo
  energygs=energy
  
  write(2,*) 'EIGENVECTORS'
  do i=1,nsiti
     write(2,'(<nsiti>(2x,f10.5))')(ham(i,j),j=1,nsiti)
  enddo
  risene=energy-ne*bc
    
  write(2,*) 'GROUND STATE ENERGY=',energy
  write(2,*) 'RISONANCE ENERGY FOR GS=', dabs(risene)
  energy=0.d0
  do i=1,ne/2-1
     energy=energy+2*w(i)
  enddo
  energy=energy+w(ne/2)+w(ne/2+1)
  write(2,*) 'ENERGY OF FIRST EXCITED STATE=', energy
  energy=-w(ne/2)+w(ne/2+1)
  write(2,*) 'EXITATION ENERGY', energy


  !=========================CARICHE=========================
  write(2,*)
  write(2,*) 'GROUND STATE ATOMIC CHARGES'
  check=0
  do i=1,nsiti
     q=0.d0
     charge=0.d0
     do j=1,7
        q=q+2*(dabs(ham(i,j)))**2
     enddo
     charge(i)=nz(i)-q
     check=check+charge(i)
     write(2,*) i,q,charge(i)
  enddo

  write(2,*) 'TOTAL CHARGE=',check
  write(2,*)
  write(2,*) 'CHARGES OF FIRST EXCITED STATE'
  check=0
  do i=1,nsiti
     q=0.d0
     charge(i)=0.d0
     do j=1,ne/2-1
        q=q+2*ham(i,j)**2
     enddo
     q=q+ham(i,ne/2+1)**2+ham(i,ne/2)**2
     charge(i)=nz(i)-q
     check=check+charge(i)
     write(2,*) i,q, charge(i)
  enddo

  write(2,*) 'TOTAL CHARGE=',check

  write(2,*)
  !=========================ORDINI DI LEGAME=========================

  write(2,*) 'BOND ORDERS'
  do i=1,nl
     BO=0
     do j=1,ne/2
        bo=bo+2*ham(bond(i,1),j)*ham(bond(i,2),j)
     enddo
     write(2,*) bond(i,1), bond(i,2), bo
  enddo
  write(2,*)
  !=========================MOMENTI DI DIPOLO=========================
  write(2,*) 'MOMENTS OF PERMANENT DIPOLE'

  mux=0
  muy=0
  muz=0

  do i=1,nsiti
     mux=mux+charge(i)*coords(i,1)
     muy=muy+charge(i)*coords(i,2)
     muz=muz+charge(i)*coords(i,3)
  enddo
  write(2,*) 'X COMPONENT=',mux
  write(2,*) 'Y COMPONENT=',muy
  write(2,*) 'Z COMPONENT=',muz
  write(2,*)
  !=========================MOMENTI DI DIPOLO DI TRANSIZIONE=========================
  write(2,*) 'MOMENTS OF TRANSITION DIPOLE'
 
  do j=ne/2,1, -1
     do k=ne/2+1, nsiti
         mux=0
         muy=0
         muz=0
        do i=1,nsiti
           mux=mux+ham(i,j)*ham(i,k)*coords(i,1)
           muy=muy+ham(i,j)*ham(i,k)*coords(i,2)
           muz=muz+ham(i,j)*ham(i,k)*coords(i,3)
        enddo
        x(k,j)=mux
        y(k,j)=muy
        z(k,j)=muz
        
        energy=w(k)-w(j)
        
        write(2,*) 'TRANSITION FROM',J, 'TO', K
        write(2,*) 'X COMPONENT=',mux
        write(2,*) 'Y COMPONENT=',muy
        write(2,*) 'Z COMPONENT=',muz
        write(2,*) 'EXITATION ENERGY=', energy
     enddo
  enddo
  !=========================COSTRUZIONE SPETTRO=========================
  open(3,file='spectra.dat')
  sigma=0.1d0
  nomega=100000
  omegastart=0.d0
  omegaend=2
  domega=(omegaend-omegastart)/nomega
  omega=omegastart-domega
  do i=ne/2,ne/2
     do j=1,nomega
        omega = omega + domega
        intensita=0.d0
        do k= ne/2+1, nsiti
           intensita=intensita+(w(k)-w(i)) * (x(k,i)**2+y(k,i)**2) *((2*pi)**(-0.5)*sigma)* dexp(-(omega-(w(k)-w(i)))**2/(2*sigma**2))
        enddo
        write(3,*) omega, intensita
     enddo
  enddo
end program
