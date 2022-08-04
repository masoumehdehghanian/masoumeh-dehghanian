program Mordad
	implicit none
	integer,parameter :: Num1 = 70400,Num3 = 704,Num = 5
	real(8),parameter :: Total1 = 1760,Total3 = 176000
	real(8),parameter :: pi = 4*atan(1.0)
	real(8),parameter :: w0 = 600000,landa = 16000
	real(8),parameter :: zR = pi*w0**2/landa
	real(8),parameter :: E0 = 0.119,tou = 110.00,FWHM = 385.70
	real(8),parameter :: w1 = 2*pi/tou,k1 = 2*pi/landa
	real(8),parameter :: c = 137.03599
	real(8),parameter :: ip = 0.67
	real(8),parameter :: n0 = 0.004
	real(8),parameter :: del_t = 0.025,del_z = 250
	complex(8),parameter :: i_cmplx = cmplx(0,1)
	integer :: k,m,n,n_prime,Itr_num
	real(8) Plus
	real(8),dimension(:),allocatable :: t,f
	real(8),dimension(:),allocatable :: ADK_rate,ne,w_Plasma
	real(8),dimension(:),allocatable :: z,w,phiG
	complex(8),dimension(:,:),allocatable :: E,G
	complex(8),dimension(:,:),allocatable :: E_bar,G_bar,Itr_value
	allocate(t(1:Num1),f(1:Num1))
	allocate(ADK_rate(1:Num1),ne(1:Num1),w_Plasma(1:Num1))
	allocate(z(+0:Num3),w(+0:Num3),phiG(+0:Num3))
	allocate(E(1:Num1,+0:Num3),G(1:Num1,+0:Num3))
	allocate(E_bar(1:Num1,+0:Num3),G_bar(1:Num1,+0:Num3),Itr_value(1:Num1,+0:Num3))
	t = 0
	f = 0
	ADK_rate = 0
	Plus = 0
	ne = 0
	w_Plasma = 0
	z = 0
	w = 0
	phiG = 0
	E = 0
	G = 0
	E_bar = 0
	G_bar = 0
	Itr_value = 0
	open(unit = 11, file = 'Mydata01.dat')
	open(unit = 12, file = 'Mydata02.dat')
	open(unit = 13, file = 'Mydata03.dat')
	open(unit = 14, file = 'Mydata04.dat')
	open(unit = 15, file = 'Mydata05.dat')
	open(unit = 16, file = 'Mydata06.dat')
	open(unit = 17, file = 'Mydata07.dat')
	open(unit = 18, file = 'Mydata08.dat')
	open(unit = 19, file = 'Mydata09.dat')
	open(unit = 20, file = 'Mydata10.dat')
	do n = 1, Num1
		t(n) = -Total1/2+n*Total1/Num1
	enddo
	do m = 1, Num1
		f(m) = m/Total1
	enddo
	do k = +0, Num3, +8
		z(k) = -Total3/2+k*Total3/Num3
 		w(k) = w0*sqrt(1+(z(k)/zR)**2)
 		phiG(k) = atan(z(k)/zR)
	enddo
	do n = 1, Num1
		do k = +0, Num3, +8
			E(n,k) = E0*(w0/w(k))*cos(k1*z(k)-w1*t(n)+phiG(k))*exp(-2*log(2.0)*(t(n)/FWHM)**2)
		enddo
	enddo
	do k = +0, Num3, +8
		do m = 1, Num1
			do n = 1, Num1
				E_bar(m,k) = E_bar(m,k) + E(n,k)*exp((-2*pi*i_cmplx/Num1)*m*n)
			enddo
		enddo
	enddo
	do m = 1, Num1
		do k = +0, Num3, +8
			Itr_value(m,k) = E_bar(m,k)
		enddo
	enddo
	do k = +0, Num3, +8
		if ( k == +0 ) then
			do n = 1, Num1
				ADK_rate(n) = 4*sqrt(2*ip)**5*exp(-2*sqrt(2*ip)**3/(3*abs(E(n,k))))/abs(E(n,k))
				Plus = Plus + ADK_rate(n)
				ne(n) = n0*(1-exp(-Plus))
				w_Plasma(n) = sqrt(4*pi*ne(n))
				G(n,k) = E(n,k)*(w_Plasma(n)/c)**2
			enddo
			do m = 1, Num1
				do n = 1, Num1
					G_bar(m,k) = G_bar(m,k) + G(n,k)*exp((-2*pi*i_cmplx/Num1)*m*n)
				enddo
			enddo
		endif
	enddo
	deallocate(ADK_rate,ne,w_Plasma)
	deallocate(E,E_bar)
	allocate(ADK_rate(1:Num1),ne(1:Num1),w_Plasma(1:Num1))
	allocate(E(1:Num1,+0:Num3),E_bar(1:Num1,+0:Num3))
	ADK_rate = 0
	Plus = 0
	ne = 0
	w_Plasma = 0
	E = 0
	E_bar = 0
	do k = +0, Num3, +8
		if ( k == +0 ) then
			do m = 1, Num1
				E_bar(m,k) = Itr_value(m,k)
			enddo
		endif
	enddo
	do k = +0, Num3-8, +8
		do m = 1, Num1
			E_bar(m,k+8) = E_bar(m,k) + (i_cmplx*c*del_z/(2*f(m)))*G_bar(m,k)
		enddo
		do m = 1, Num1
			if ( abs(E_bar(m,k+8) - Itr_value(m,k+8)) > 0.01 ) then
				E_bar(m,k) = E_bar(m,k+8)
			endif
		enddo
		do n = 1, Num1
			E(n,k) = 0
		enddo
		do n = 1, Num1
			do m = 1, Num1
				E(n,k) = E(n,k) + E_bar(m,k)*exp((+2*pi*i_cmplx/Num1)*m*n)/Num1
			enddo
		enddo
		do n = 1, Num1
			ADK_rate(n) = 4*sqrt(2*ip)**5*exp(-2*sqrt(2*ip)**3/(3*abs(E(n,k))))/abs(E(n,k))
			Plus = Plus + ADK_rate(n)
			ne(n) = n0*(1-exp(-Plus))
			w_Plasma(n) = sqrt(4*pi*ne(n))
			G(n,k) = E(n,k)*(w_Plasma(n)/c)**2
		enddo
		do m = 1, Num1
			E_bar(m,k) = 0
		enddo
		do m = 1, Num1
			do n = 1, Num1
				E_bar(m,k) = E_bar(m,k) + E(n,k)*exp((-2*pi*i_cmplx/Num1)*m*n)
			enddo
		enddo
		do m = 1, Num1
			G_bar(m,k) = 0
		enddo
		do m = 1, Num1
			do n = 1, Num1
				G_bar(m,k) = G_bar(m,k) + G(n,k)*exp((-2*pi*i_cmplx/Num1)*m*n)
			enddo
		enddo
		do m = 1, Num1
			E_bar(m,k+8) = E_bar(m,k) + (i_cmplx*c*del_z/(2*f(m)))*G_bar(m,k)
		enddo
		deallocate(ADK_rate,ne,w_Plasma)
		allocate(ADK_rate(1:Num1),ne(1:Num1),w_Plasma(1:Num1))
		ADK_rate = 0
		Plus = 0
		ne = 0
		w_Plasma = 0
		do Itr_num = 1, Num
			do m = 1, Num1
				if ( abs(E_bar(m,k+8) - E_bar(m,k)) > 0.01 ) then
					E_bar(m,k) = E_bar(m,k+8)
				endif
			enddo
			do n = 1, Num1
				E(n,k) = 0
			enddo
			do n = 1, Num1
				do m = 1, Num1
					E(n,k) = E(n,k) + E_bar(m,k)*exp((+2*pi*i_cmplx/Num1)*m*n)/Num1
				enddo
			enddo
			do n = 1, Num1
				ADK_rate(n) = 4*sqrt(2*ip)**5*exp(-2*sqrt(2*ip)**3/(3*abs(E(n,k))))/abs(E(n,k))
				Plus = Plus + ADK_rate(n)
				ne(n) = n0*(1-exp(-Plus))
				w_Plasma(n) = sqrt(4*pi*ne(n))
				G(n,k) = E(n,k)*(w_Plasma(n)/c)**2
			enddo
			do m = 1, Num1
				E_bar(m,k) = 0
			enddo
			do m = 1, Num1
				do n = 1, Num1
					E_bar(m,k) = E_bar(m,k) + E(n,k)*exp((-2*pi*i_cmplx/Num1)*m*n)
				enddo
			enddo
			do m = 1, Num1
				G_bar(m,k) = 0
			enddo
			do m = 1, Num1
				do n = 1, Num1
					G_bar(m,k) = G_bar(m,k) + G(n,k)*exp((-2*pi*i_cmplx/Num1)*m*n)
				enddo
			enddo
			do m = 1, Num1
				E_bar(m,k+8) = E_bar(m,k) + (i_cmplx*c*del_z/(2*f(m)))*G_bar(m,k)
			enddo
			deallocate(ADK_rate,ne,w_Plasma)
			allocate(ADK_rate(1:Num1),ne(1:Num1),w_Plasma(1:Num1))
			ADK_rate = 0
			Plus = 0
			ne = 0
			w_Plasma = 0
		enddo
	enddo
	do k = +0, Num3, +8
		if ( k == Num3 ) then
			do n = 1, Num1
				do m = 1, Num1
					E(n,k) = E(n,k) + E_bar(m,k)*exp((+2*pi*i_cmplx/Num1)*m*n)/Num1
				enddo
			enddo
		endif
	enddo
	do n = 1, Num1
		write(11,*) t(n),real(E(n,+0))
		write(12,*) t(n),real(E(n,+8))
		write(13,*) t(n),real(E(n,+16))
		write(14,*) t(n),real(E(n,+64))
		write(15,*) t(n),real(E(n,+80))
		write(16,*) t(n),real(E(n,+88))
		write(17,*) t(n),real(E(n,+176))
		write(18,*) t(n),real(E(n,+320))
		write(19,*) t(n),real(E(n,+352))
		write(20,*) t(n),real(E(n,+704))
	enddo
	close(unit = 11)
	close(unit = 12)
	close(unit = 13)
	close(unit = 14)
	close(unit = 15)
	close(unit = 16)
	close(unit = 17)
	close(unit = 18)
	close(unit = 19)
	close(unit = 20)
	deallocate(t,f)
	deallocate(ADK_rate,ne,w_Plasma)
	deallocate(z,w,phiG)
	deallocate(E,G)
	deallocate(E_bar,G_bar,Itr_value)
end program Mordad
