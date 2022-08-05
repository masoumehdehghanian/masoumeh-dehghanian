program SplitOperator
	implicit none
	integer :: j,j_Mirror,n,m
	integer,parameter :: Num1 = 70400,Num2 = 2080
	real(8),parameter :: Total1 = 1760,Total2 = 520
	real(8),parameter :: pi = 4*atan(1.0)
	real(8),parameter :: dy = 0.2500,dt = 0.0250
	real(8),parameter :: E0 = 0.119,tou = 110.00,w0 = 2*pi/tou,etta = 0,delta = 385.70/sqrt(2*log(2.0))
	complex(8),parameter :: i_cmplx = cmplx(0,1)
	real(8),dimension(:),allocatable :: t,f
	real(8),dimension(:),allocatable :: y
	real(8),dimension(:),allocatable :: mask,v_ql,Psi_0 
	real(8),dimension(:),allocatable :: E
	real(8),dimension(:),allocatable :: plus
	real(8),dimension(:),allocatable :: a_nl,v_nl,d_nl
	real(8),dimension(:,:),allocatable :: increment
	complex(8),dimension(:),allocatable :: P,Q,denom,Phi 
	complex(8),dimension(:),allocatable :: dft
	complex(8),dimension(:),allocatable :: intermediate 
	complex(8),dimension(:,:),allocatable :: Psi 
	complex(8),dimension(:,:),allocatable :: U0_Right,U0_Left  
	complex(8),dimension(:,:),allocatable :: Uint 
	allocate(t(1:Num1),f(1:Num1))
	allocate(y(1:Num2))
	allocate(mask(1:Num2),v_ql(1:Num2),Psi_0(1:Num2))
	allocate(E(1:Num1))
	allocate(plus(1:Num1))
	allocate(a_nl(1:Num1),v_nl(1:Num1),d_nl(1:Num1))
	allocate(increment(1:Num1,1:Num2))
	allocate(P(1:Num2),Q(1:Num2),denom(1:Num2),Phi(1:Num2))
	allocate(dft(1:Num1))
	allocate(intermediate(1:Num2))
	allocate(Psi(1:Num1,1:Num2))
	allocate(U0_Right(1:Num2,1:Num2),U0_Left(1:Num2,1:Num2))
	allocate(Uint(1:Num2,1:Num2))
	t = 0
	f = 0
	y = 0
	mask = 0
	v_ql = 0
	Psi_0 = 0
	E = 0
	plus = 0
	a_nl = 0
	v_nl = 0
	d_nl = 0
	increment = 0
	P = 0
	Q = 0
	denom = 0
	Phi = 0
	dft = 0
	intermediate = 0
	Psi = 0
	U0_Right = 0
	U0_Left = 0
	Uint = 0	
 	open(unit = 1, file = 'Mydata07_01.dat')
 	open(unit = 2, file = 'Mydata07_02.dat')
 	open(unit = 3, file = 'Mydata07_03.dat')
 	open(unit = 4, file = 'Mydata07_04.dat')
 	open(unit = 5, file = 'Mydata07_05.dat')
 	open(unit = 6, file = 'Mydata07_06.dat')
 	open(unit = 7, file = 'Mydata07_07.dat')
 	do n = 1, Num1
		t(n) = -Total1/2 + n*Total1/Num1
	enddo
	do j = 1, Num2
  		y(j) = -Total2/2 + j*Total2/Num2
 	enddo
 	do m = 1, Num1
 		f(m) = m/Total1
 	enddo
 	do j = 1, Num2
 		v_ql(j) = -1/sqrt(1+y(j)**2)
 	enddo
 	do j = 1, Num2
  		read(1,*) Psi_0(j)
  		Psi(1,j) = Psi_0(j)
 	enddo
 	do n = 1, Num1
		E(n) = E0*exp(-(t(n)/delta)**2)*cos(w0*t(n)+etta)	
 	enddo
 	do j = 1, Num2
 		if (j == 1) then
   			mask(j) = 0
  		elseif (j == Num2) then
   			mask(j) = 0
  		elseif (j < 208) then
   			mask(j) = cos((pi/2)*(y(j)-y(208))/(y(1)-y(208)))**0.125
  		elseif (j > 1872) then
   			mask(j) = cos((pi/2)*(y(j)-y(1872))/(y(Num2)-y(1872)))**0.125
  		else
   			mask(j) = 1.0
  		endif
 	enddo
 	do j = 1, Num2
  		if (j == 1) then 
   			do j_Mirror = 1, Num2
    				if (j_Mirror == j) then  
     					U0_Right(j,j_Mirror) = 1 - i_cmplx*(dt/2)*(1/dy**2)
     					U0_Left(j,j_Mirror) = 1 + i_cmplx*(dt/2)*(1/dy**2)
    				elseif (j_Mirror == j+1) then
     					U0_Right(j,j_Mirror) = -i_cmplx*(dt/2)*(-1/(2*dy**2))
     					U0_Left(j,j_Mirror) = +i_cmplx*(dt/2)*(-1/(2*dy**2))
    				else
     					U0_Right(j,j_Mirror) = 0
     					U0_Left(j,j_Mirror) = 0
    				endif
   			enddo
   		elseif (j == Num2) then
   			do j_Mirror = 1, Num2
    				if (j_Mirror == j) then
     					U0_Right(j,j_Mirror) = 1 - i_cmplx*(dt/2)*(1/dy**2)
     					U0_Left(j,j_Mirror) = 1 + i_cmplx*(dt/2)*(1/dy**2)
    				elseif (j_Mirror == j-1) then
     					U0_Right(j,j_Mirror) = -i_cmplx*(dt/2)*(-1/(2*dy**2))
     					U0_Left(j,j_Mirror) = +i_cmplx*(dt/2)*(-1/(2*dy**2))
    				else
     					U0_Right(j,j_Mirror) = 0
     					U0_Left(j,j_Mirror) = 0
    				endif
   			enddo
   		else
   			do j_Mirror = 1, Num2
    				if (j_Mirror == j) then
     					U0_Right(j,j_Mirror) = 1 - i_cmplx*(dt/2)*(1/dy**2)
     					U0_Left(j,j_Mirror) = 1 + i_cmplx*(dt/2)*(1/dy**2)
    				elseif (j_Mirror == j-1) then
     					U0_Right(j,j_Mirror) = -i_cmplx*(dt/2)*(-1/(2*dy**2))
     					U0_Left(j,j_Mirror) = +i_cmplx*(dt/2)*(-1/(2*dy**2))
    				elseif (j_Mirror == j+1) then
     					U0_Right(j,j_Mirror) = -i_cmplx*(dt/2)*(-1/(2*dy**2))
     					U0_Left(j,j_Mirror) = +i_cmplx*(dt/2)*(-1/(2*dy**2))
    				else
     					U0_Right(j,j_Mirror) = 0
     					U0_Left(j,j_Mirror) = 0
    				endif
   			enddo
  		endif
	enddo
	do n = 1, Num1-1
  		 do j = 1, Num2
   			if (j == 1) then
    				do j_Mirror = 1, Num2
     					if (j_Mirror == j) then  
      						Uint(j,j_Mirror) = exp(-i_cmplx*(dt/2)*(+v_ql(j) -y(j)*E(n)))
     					else
      						Uint(j,j_Mirror) = 0
     					endif
    				enddo
    			elseif (j == Num2) then
    				do j_Mirror = 1, Num2
     					if (j_Mirror == j) then
      						Uint(j,j_Mirror) = exp(-i_cmplx*(dt/2)*(+v_ql(j) -y(j)*E(n)))
     					else
      						Uint(j,j_Mirror) = 0
     					endif
    				enddo
    			else
    				do j_Mirror = 1, Num2
    					if (j_Mirror == j) then
      						Uint(j,j_Mirror) = exp(-i_cmplx*(dt/2)*(+v_ql(j) -y(j)*E(n)))
     					else
      						Uint(j,j_Mirror) = 0
     					endif
    				enddo
   			endif
  		enddo 
  		do j = 1, Num2
   			do j_Mirror = 1, Num2 
    				Phi(j) = Phi(j) + Uint(j,j_Mirror)*Psi(n,j_Mirror)
   			enddo
  		enddo
  		do j = 1, Num2
   			intermediate(j) = mask(j)*Phi(j)
  		enddo
  		deallocate(Phi)
  		allocate(Phi(1:Num2))
  		Phi = 0 
		do j = 1, Num2
   			do j_Mirror = 1, Num2  
    				Phi(j) = Phi(j) + U0_Right(j,j_Mirror)*intermediate(j_Mirror)
   			enddo
  		enddo
  		do j = 1, Num2
   			if (j == 1) then
    				P(j) = -U0_Left(j,j+1)/U0_Left(j,j)
    				Q(j) = Phi(j)/U0_Left(j,j)
   			elseif (j /= Num2) then
    				denom(j) = U0_Left(j,j) + U0_Left(j,j-1)*P(j-1)
    				P(j) = -U0_Left(j,j+1)/denom(j)
    				Q(j) = (Phi(j) - U0_Left(j,j-1)*Q(j-1))/denom(j)
   			else
    				denom(j) = U0_Left(j,j) + U0_Left(j,j-1)*P(j-1)
    				P(j) = 0
    				Q(j) = (Phi(j) - U0_Left(j,j-1)*Q(j-1))/denom(j)
   			endif 
  		enddo
  		deallocate(intermediate)
  		allocate(intermediate(1:Num2))
  		intermediate = 0
  		do j = Num2, 1, -1
   			if (j == Num2) then
    				intermediate(j) = mask(j)*Q(j)
   			else
    				intermediate(j) = mask(j)*(P(j)*intermediate(j+1) + Q(j))
   			endif
  		enddo
  		deallocate(Phi)
  		allocate(Phi(1:Num2))
  		Phi = 0
  		do j = 1, Num2
   			do j_Mirror = 1, Num2 
    				Phi(j) = Phi(j) + Uint(j,j_Mirror)*intermediate(j_Mirror)
   			enddo
  		enddo
  		do j = 1, Num2
   			Psi(n+1,j) = mask(j)*Phi(j)
  		enddo
  		deallocate(Phi)
  		deallocate(Uint)
  		deallocate(intermediate)
  		deallocate(P,Q,denom)
  		allocate(Phi(1:Num2))
  		allocate(Uint(1:Num2,1:Num2))
  		allocate(intermediate(1:Num2))
  		allocate(P(1:Num2),Q(1:Num2),denom(1:Num2))
  		Phi = 0
  		Uint = 0
  		intermediate = 0
  		P = 0
  		Q = 0
  		denom = 0
  	enddo
  	do n = 1, Num1
  		do j = 1, Num2
  			increment(n,j) = (E(n) - y(j)/(sqrt(1+y(j)**2))**3)*abs(Psi(n,j))**2
  			plus(n) = plus(n) + increment(n,j)
  		enddo
  		a_nl(n) = dy*(plus(n)-(increment(n,1)+increment(n,Num2))/2)
  		write(2,*) t(n),a_nl(n)
  	enddo
  	deallocate(increment)
  	deallocate(plus)
  	allocate(increment(1:Num1,1:Num2))
  	allocate(plus(1:Num1))
  	increment = 0
  	plus = 0
  	do m = 1, Num1
		do n = 1, Num1
			dft(m) = dft(m) + a_nl(n)*exp((-2*pi*i_cmplx/Num1)*m*n)
		enddo
	enddo
	do m = 1, Num1
		write(3,*) f(m)/0.01,log10(abs(dft(m))**2)
	enddo
	deallocate(dft)
	allocate(dft(1:Num1))
	dft = 0
  	do n = 1, Num1
  		do j = 1, Num2
  			increment(n,j) = (1/dy)*abs(Psi(n,j))**2
  			plus(n) = plus(n) + increment(n,j)
  		enddo
  		v_nl(n) = dy*(plus(n)-(increment(n,1)+increment(n,Num2))/2)
  		write(4,*) t(n),v_nl(n)
  	enddo
  	deallocate(increment)
  	deallocate(plus)
  	allocate(increment(1:Num1,1:Num2))
  	allocate(plus(1:Num1))
  	increment = 0
  	plus = 0
  	do m = 1, Num1
		do n = 1, Num1
			dft(m) = dft(m) + v_nl(n)*exp((-2*pi*i_cmplx/Num1)*m*n)
		enddo
	enddo
	do m = 1, Num1
		write(5,*) f(m)/0.01,log10(abs(dft(m))**2)
	enddo
	deallocate(dft)
	allocate(dft(1:Num1))
	dft = 0
  	do n = 1, Num1
  		do j = 1, Num2
  			increment(n,j) = (y(j))*abs(Psi(n,j))**2
  			plus(n) = plus(n) + increment(n,j)
  		enddo
  		d_nl(n) = dy*(plus(n)-(increment(n,1)+increment(n,Num2))/2)
  		write(6,*) t(n),d_nl(n)
  	enddo
	do m = 1, Num1
		do n = 1, Num1
			dft(m) = dft(m) + d_nl(n)*exp((-2*pi*i_cmplx/Num1)*m*n)
		enddo
	enddo
	do m = 1, Num1
		write(7,*) f(m)/0.01,log10(abs(dft(m))**2)
	enddo
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(6)
	close(7)
	deallocate(t,f)
	deallocate(y)
	deallocate(mask,v_ql,Psi_0)
	deallocate(E)
	deallocate(plus)
	deallocate(a_nl,v_nl,d_nl)
	deallocate(increment)
	deallocate(P,Q,denom,Phi)
	deallocate(dft)
	deallocate(intermediate)
	deallocate(Psi)
	deallocate(U0_Right,U0_Left)
	deallocate(Uint)
end program SplitOperator
