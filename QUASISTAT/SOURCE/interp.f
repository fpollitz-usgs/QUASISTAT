	subroutine interp(mus,lams,denss)
c	Give back mus ,lams, denss at radius r0 interpolated between layers.
	double complex svalue
	double complex mus,lams,kaps
	real*8 denss,r0,f
        real*8 tau1,tau2
c	real*8 pi,twopi
	real*4 kappa,mu,mu2
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),mu2(200),eta(200),
     &	eta1(200),svalue,r0
	common/matl/l
c	pi=3.14159265358979d0
	j=0
5	j=j+1
c		if(j.eq.34) write(6,*)'INTERP34: j=',j,'r0=',r0,'rt(',j,')=',rt(j)
c		if(j.eq.35) write(6,*)'INTERP35: j=',j,'r0=',r0,'rt(',j,')=',rt(j)
	if(rt(j).lt.real(r0)) go to 5
c		write(6,*)'INTERP: j=',j,'r0=',r0,'rt(',j,')=',rt(j)
cOLD	f=(r0-dble(rb(j)))/dble(rt(j)-rb(j))
c		write(6,*)'INTERP: qbet(',j,')=',qbet(j)
c	The following is adapted from visco2pt5dsubs.f
c	tau1=mu0/et1val
c          tau2=mu2val/etaval
c          mu=mu0*s*(s+tau2)/((s+tau2)*(s+tau1)+mu0*s/etaval)
c
	tau1=dble(mu(j)/eta1(j))
	tau2=dble(mu2(j)/eta(j))
c		write(6,*)'INTERP: tau1,tau2=',tau1,tau2
c		write(6,*)'INTERP: mu(',j,')=',mu(j),'svalue=',svalue
c		write(6,*)'INTERP: mu2(',j,')=',mu2(j)
c		write(6,*)'INTERP: eta1(',j,')=',eta1(j),'eta=',eta(j)
	mus=dble(mu(j))*svalue*(svalue+tau2)/((svalue+tau2)*(svalue+tau1)+dble(mu(j))*svalue/eta(j)) 
c		write(6,*)'INTERP: mus=',mus
	kaps=dble(kappa(j))
	lams=kaps-(2.d0/3.d0)*mus
	denss=dble(dens(j))
	return
	end
