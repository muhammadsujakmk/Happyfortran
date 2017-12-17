	program PD2_Runge_Kutta
	
	implicit none
	integer :: i,N
	real*8 :: u1,x0,y0,y_0,u0,t,y,tmin,tmax,h,f0,f1,t0,u
	real*8, external :: f
	
	y0=-0.99/2! Initialize
	y_0=0.d0
	
	tmin=0.d0  ! x range.
	tmax=100d0 !
	N=100000
	h=(tmax-tmin)/real(N,8)
	
	
	! Menghitung y'=u(x,y) dan u'=f(x,y,u)
	u0=y_0
	DO i=0,N
	
	t0=tmin+h*i
	
	f0=f(t0,y0,u0)
	u1=u0+h*f0
	f1=f(t0+h, y0+(h*u0), u1)
	
	y=y0+(0.5d0*h*(u0+u1))
	u=u0+(0.5d0*h*(f0+f1))
	
	open(unit=10,file="PD2_RK.dat",status="unknown")
	write(10,*) t0,"                      ",y
	t0=t0
	u0=u
	y0=y

	END DO
	close(10)

	
	
	write(*,*) "Execution Success!." 
	write(*,*) "Take your file in PD2_RK.dat and Plot ."
	
	end program
	
	function f(z,x,y)
	implicit none
	real*8 :: x,y,f,z,w,T,g,pi
	w=195 ! kg
	T=2 ! periode
	pi=3.14
	g = 9.8 ! kg/m^2
	f=(-62.5*pi*g*x/w)-((y/(10*T))) !Equation for Damping Harmonic Motion
	return 
	end
 