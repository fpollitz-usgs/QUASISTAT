	function lapl(f,t)
c	Calculate inverse Laplace transform of function f at time t.
c	The array of complex s-values where f(s) was computed
c	is already stored in common.
c	Assume that only the imaginary part of s-values vary,
c	not the real part.
c	Further assume that the imaginary part of s varies quadratically.
	real*8 lapl
	real*8 t,v
	real*8 s0
	complex*16 ui,s,f(15),ds,fn,gn,val
c	real*8 wt(15)
	real*8 sarr(7)
	real*8 bigb(8,8),c(8),a(8)
	real*8 emax,eig(8),vm(8,8)
	complex*16 bcal,ccal,acal
	complex*16 sval
	common/flapl/sval(15)
	ui=dcmplx(0.,1.)
	s0=5.6d0*dble(sval(1)) 
c		write(6,*)'LAPL: s0=',s0
c	Assign weights so that the static displacement response in f(15) is fit exactly.
c	do n=1,15
c	if(n.lt.15) then
c	wt(n)=1.d0
c	else 
c	wt(n)=1.d-4
c	endif
c	enddo

	sarr(1)=s0*2.d0
	sarr(2)=s0
	sarr(3)=s0/2.d0
	sarr(4)=s0/10.d0
	sarr(5)=s0/50.d0
	sarr(6)=s0/100.d0
	sarr(7)=s0/500.d0
c		write(6,*)'LAPL: sarr=',sarr
		do l=1,8
	ccal=0.d0
	do k=1,8
	bcal=0.
	do n=1,15
	
	if(l.le.7) fn=1.d0/(sval(n)*(sval(n)+sarr(l)))
	if(l.eq.8) fn=1.d0/sval(n) 
	if(k.le.7) gn=1.d0/(sval(n)*(sval(n)+sarr(k)))
	if(k.eq.8) gn=1.d0/sval(n) 
	bcal=bcal+fn*gn
	enddo
	bigb(k,l)=real(bcal)
	enddo
	do n=1,15
	if(l.le.7) fn=1./(sval(n)*(sval(n)+sarr(l)))
	if(l.eq.8) fn=1./sval(n) 
	ccal=ccal+f(n)*fn
	enddo
	c(l)=real(ccal)
		enddo
c			bigb(8,8)=bigb(8,8)+1.e+6
cNOTE
c		write(6,*)'entering ainver: bigb=',bigb
c		write(6,*)'c=',c
c	call ainver(a,bigb,c,8)
c		write(6,*)'a=',a
c		write(6,*)'bigb=',bigb
	m=8
	mf1=m
c		write(6,*) 'entering svdcmp: m=',m
	call svdcmp(bigb,m,mf1,8,8,eig,vm)
		write(6,*)'LAPL: initial eigenvalues are'
		write(6,*) (eig(j), j=1,m)
	emax=0.d0
	do i=1,m
	if(eig(i).gt.emax) emax=eig(i)
	enddo
	do i=1,m
	eig(i)=eig(i)+(2.d-6)*emax
	enddo
		write(6,*)'LAPL: revised eigenvalues are'
		write(6,*) (eig(j), j=1,m)
c		pause
	call svdksb(bigb,eig,vm,m,mf1,8,8,c,a)
		write(6,*)'a=',a
c		pause
	lapl=0.
	do k=1,7
	v=-sarr(k)*t
	if(v.lt.-20.d0) v=-20.d0
	lapl=lapl+a(k)*(1.d0-dexp(v))/sarr(k)
	enddo
	lapl=lapl+a(8)
c--------------
cNOTE
c	Evaluate misfit
	do n=1,15
	val=0.
	do l=1,8
	if(l.le.7) fn=1./(sval(n)*(sval(n)+sarr(l)))
	if(l.eq.8) fn=1./sval(n) 
	val=val+fn*a(l)
c		write(6,*)'LAPL: fn=',fn,'sval(',n,')=',sval(n)
	enddo
	write(6,*)'LAPL: f(',n,')=',f(n),'calc=',val
	enddo
c-------
c		pause
	return
	end	
