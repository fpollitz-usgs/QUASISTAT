	subroutine lbound(zyst,j,csys)
c
* Calculates STARTING VALUES for integration of the spheroidal
* and toroidal motion upwards to the source level.
* The values are derived from a homogeneous isotropic Earth model
* assumed beneath the starting radius r1. 
c
	character*1 csys
	real*4 kappa,mu,mu2 
	double complex zyst(j)
	double complex svalue
	double complex mus,lams,biga
	real*8 r,r2,flp,flp3
	real*8 g0
        real*8 tau1,tau2
        double complex a1p,a2p,a3p,a4p
        double complex b1p,b2p,b3p,b4p
	real*8 rlp
	real*8 denss,rmus
      	real*8 r0
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),mu2(200),eta(200),
     &	eta1(200),svalue,r0
	common/gval/g0
	common/matl/l
	double complex aj(4,2)
c	Assume that we're at bottom of layer 1. r0=rb(1) is already in the calling program.
	tau1=dble(mu(1)/eta1(1))
        tau2=dble(mu2(1)/eta(1))
        r=dble(rb(1))
c	Obtain complex-valued elastic parameters at radius r0
	call interp(mus,lams,denss)
		write(6,*)'LBOUND: r=',r,'r0=',r0
		write(6,*)'LBOUND: mus=',mus

      if(csys.eq.'s') then            !   Spheroidal motion
      if (l.eq.0) then
c	Use eqns 8.167-8.168 of Dahlen and Tromp for to specify the (U,R,V,S) displacement-stress vector
c	for the l=0 case, noting that j_0(gamma r) tends to 1 and j_1(gamma r) tends to (1/3)*gamma*r
c	as gamma tends to zero.
c	Divide both equations by j_1(gamma r) to get the elements below.
        zyst(1)=1.d0
        zyst(2)=(3.d0*lams+2.d0*mus)/r0
	zyst(3)=0.d0
	zyst(4)=0.d0
	zyst(5)=0.d0
        return
      else
	biga=lams+2.d0*mus
c       Determine matrix elements (Takeuchi and Saito, 1972, p. 243-244).
        flp=dble(l)
        flp3=2.d0*(2.d0*flp+3.d0)
        r2=r*r
	aj(1,1)=flp*(flp+1.d0)
        aj(2,1)=2.d0*mus*(flp*(flp*flp-1.d0))
        aj(3,1)=(flp+1.d0)
        aj(4,1)=2.d0*mus*(flp*flp-1.d0)
        a1p=-(flp+2.d0)*r2/flp3
        a2p=-(lams+2.d0*mus*(flp+2.d0)*(flp+1.d0)/flp3)*r2
        a3p=-r2/flp3
        a4p=-2.d0*mus*(flp+1.d0)*r2/flp3
        b1p=flp*(flp+1.d0)*r2/flp3
        b2p=2.d0*mus*flp*(flp+1.d0)*(flp+1.d0)*r2/flp3
        b3p=(flp+3.d0)*r2/flp3
        b4p=mus*2.d0*flp*(flp+2.d0)*r2/flp3
        aj(1,2)=mus*(flp+1.d0)*a1p+biga*b1p
        aj(2,2)=mus*(flp+1.d0)*a2p+biga*b2p
        aj(3,2)=mus*(flp+1.d0)*a3p+biga*b3p
        aj(4,2)=mus*(flp+1.d0)*a4p+biga*b4p
c       Normalize matrix elements.
        rlp=1.d0
        do m=1,2
        aj(1,m)=aj(1,m)*rlp/r
        aj(2,m)=aj(2,m)*rlp/r2
        aj(3,m)=aj(3,m)*rlp/r
        aj(4,m)=aj(4,m)*rlp/r2
	enddo

c	Construct second order minors.
	zyst(1)=aj(1,1)*aj(2,2)-aj(2,1)*aj(1,2)
        zyst(2)=aj(1,1)*aj(3,2)-aj(3,1)*aj(1,2)
        zyst(3)=aj(1,1)*aj(4,2)-aj(4,1)*aj(1,2)
        zyst(4)=aj(2,1)*aj(3,2)-aj(3,1)*aj(2,2)
        zyst(5)=aj(2,1)*aj(4,2)-aj(4,1)*aj(2,2)
cTE
c	When gravitation is included, impose a zero-displacement boundary condition
c	at the bottom of the model in order to stabilize the solution at small l.
		if(g0.gt.0.d0) then
		zyst(1)=0.d0
		zyst(2)=0.d0
		zyst(3)=0.d0
		zyst(4)=0.d0
		zyst(5)=1.d0
		endif
      endif
	
      else if(csys.eq.'t') then   !    Toroidal motion
	zyst(1)=1.d0
        zyst(2)=mus*dble(l-1)*(1./r)

      else
       stop 'Error in <stavani>: No motion specified!'
      endif

      return
      end 


	subroutine derivs2(x,yn,yout)
	real*8 r0,r2,facl,fac
	real*8 x,rintp
	double complex yn(2),yout(2),aj(2,2)
	double complex svalue
	double complex mus,lams
	real*8 denss
	real*4 kappa,mu,mu2
	common/radlev/rintp,bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),mu2(200),eta(200),
     &	eta1(200),svalue,r0
	common/matl/l
	r0=rintp
c		write(6,*)'entering DERIVS2: r0=',r0
	call interp(mus,lams,denss)
c		write(6,*)'out of DERIVS2: r0=',r0,'mus=',mus,'l=',l
	r0=x
	r2=r0*r0
	facl=dble((l-1)*(l+2))
c	facl=(l-1.d0)*(l+2.d0)
	  aj(1,1)=1.d0/r0
	  aj(1,2)=1.d0/mus
  	  aj(2,1)=facl*mus/r2
	  aj(2,2)=-3.d0/r0
	call prod(aj,yn,2,yout) 
	return
	end

	subroutine derivs4(x,yn,yout)
	real*8 r0
	real*8 x,rintp
	double complex yn(4),yout(4),aj(4,4)
	double complex svalue
	real*4 kappa,mu,mu2
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),mu2(200),eta(200),
     &	eta1(200),svalue,r0
	common/matl/l
	common/radlev/rintp,bigr
c		write(6,*)'derivs4: x=',x
c		write(6,*)'derivs4: yn=',yn
	r0=x
c		write(6,*)'entering amat4: r0=',r0
	call amat4(aj,1)
c		write(6,*)'derivs4: out of amat4: aj=',aj
	call prod(aj,yn,4,yout) 
c		write(6,*)'derivs4: out of prod: dydr=',yout
	return
	end

	subroutine derivs5(x,yn,yout)
	real*8 r0
	real*8 x,rintp
	double complex yn(5),yout(5),aj(5,5)
	double complex svalue
	real*4 kappa,mu,mu2
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),mu2(200),eta(200),
     &	eta1(200),svalue,r0
	common/matl/l
	common/radlev/rintp,bigr
c		write(6,*)'derivs5: x[input]=',x
	r0=x
c		write(6,*)'derivs5: yn=',yn
c		write(6,*)'derivs5: r0=',r0,'rintp=',rintp
c		write(6,*)'derivs5: qbet(1),qbet(2)=',qbet(1),qbet(2)
	call amat5(aj,1)
c		write(6,*)'out of amat5: aj=',aj
	call prod(aj,yn,5,yout) 
c		write(6,*)'derivs5: yout=',yout
c		pause
	return
	end

    	subroutine amat4(aj,i)
c 	Use to determine A(r) (returned in aj) in layer at radius r=r0.
c	dY(r)/dr = A(r)Y.
	real*4 kappa,mu,mu2 
	double complex aj(4,4)
	double complex svalue
	double complex mus,lams,biga
	real*8 rintp
	real*8 fl21,r,r0,denss,r2
	real*8 g0,g
	common/radlev/rintp,bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),mu2(200),eta(200),
     &	eta1(200),svalue,r0
	common/gval/g0
	common/matl/l
c	mus=dble(mu(i))
c	lams=dble(kappa(i)-2.*mu(i)/3.)
c	denss=dble(dens(i)) 
c		write(6,*)'amat4: r0=',r0,'rintp=',rintp
	r=r0
	r0=rintp
	call interp(mus,lams,denss)
	r0=r
	biga=lams+2.d0*mus
	fl21=dble(l*(l+1))
	r=r0
	r2=r*r
c	Employ the Cowling approximation
c	We also assume that gravitational acceleration varies as 1/r,
c	which is a good approximation for the upper mantle.
	g=(g0*1.d+2)*dble(bigr*1.e+5)*1.d-11
	aj(1,1)=-2.d0*lams/(r*biga)
	aj(1,2)=1.d0/biga
	aj(1,3)=fl21*lams/(r*biga)
	aj(1,4)=0.d0
	aj(2,1)=(4.d0*(lams+mus)/r)/r+(2.d0/r)*lams*aj(1,1)-4.d0*denss*g/r2
	aj(2,2)=-2.d0/r+(2.d0/r)*lams*aj(1,2)
	aj(2,3)=-fl21*(2.d0*(lams+mus)/r)/r+(2.d0/r)*lams*aj(1,3)
     &	+denss*g*fl21/r2
	aj(2,4)=fl21/r
	aj(3,1)=-1.d0/r
	aj(3,2)=0.d0
	aj(3,3)=1.d0/r
	aj(3,4)=1.d0/mus
	aj(4,1)=-(lams/r)*aj(1,1)+2.d0*(mus-biga)/r2+denss*g/r2 
	aj(4,2)=-(lams/r)*aj(1,2)
	aj(4,3)=-(lams/r)*aj(1,3)+(fl21*biga-2.d0*mus)/r2
	aj(4,4)=-3.d0/r
c	aj contains A(r).  
	return
	end 

    	subroutine amat5(aj,i0)
c 	Use to determine A(r) (returned in aj) in layer at radius r=r0.
c	dY(r)/dr = A(r)Y.
	real*4 kappa,mu,mu2 
	dimension id(4,4)
	double complex aj(5,5),a(4,4)
	double complex svalue
	double complex mus,lams,biga
	real*8 rintp
	real*8 fac
	real*8 fl21,r,r0,denss,r2
	real*8 g0,g
	common/radlev/rintp,bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),mu2(200),eta(200),
     &	eta1(200),svalue,r0
	common/gval/g0
	common/matl/l
c	mus=dble(mu(i0))
c	lams=dble(kappa(i0)-2.*mu(i0)/3.)
c	denss=dble(dens(i0)) 
	r=r0
	r0=rintp
	call interp(mus,lams,denss)
c		write(6,*)'AMAT5: mus=',mus
c		write(6,*)'amat5: out of interp: mus,lams=',mus,lams
c		if(r0.ge.0.999572) write(6,*)'out of interp: r0=',r0,'denss=',denss
	r0=r
c		write(6,*)'amat5: mus,lams,denss=',mus,lams,denss
c		pause
	biga=lams+2.d0*mus
c	fl21=l*(l+1.d0)
	fl21=dble(l*(l+1))
	r=r0
	r2=r*r
c	Below commented-out lines are for the non-gravitational case.
c	aj(1,1)=-2.d0/r
c	aj(1,2)=-2.d0*fl21*(biga-lams*lams/biga-mus)/r2
c	aj(1,3)=fl21/r
c	aj(1,4)=-fl21*(lams/biga)/r
c	aj(1,5)=0.d0
c	aj(2,1)=0.d0
c	aj(2,2)=(1.d0-2.d0*lams/biga)/r
c	aj(2,3)=1.d0/mus
c	aj(2,4)=1.d0/biga
c	aj(2,5)=0.d0
c	aj(3,1)=-2.d0*(lams/biga)/r
c	aj(3,2)=fl21*(biga-lams*lams/biga)/r2-2.d0*mus/r2
c	aj(3,3)=-(3.d0+2.d0*lams/biga)/r
c	aj(3,4)=0.d0
c	aj(3,5)=1.d0/biga
c	aj(4,1)=2.d0/r
c	aj(4,2)=4.d0*(biga-lams*lams/biga-mus)/r2
c	aj(4,3)=0.d0
c	aj(4,4)=-(1.d0-2.d0*lams/biga)/r
c	aj(4,5)=1.d0/mus
c	aj(5,1)=4.*(biga-lams*lams/biga-mus)/r2
c	aj(5,2)=0.d0
c	aj(5,3)=4.*(biga-lams*lams/biga-mus)/r2
c	aj(5,4)=fl21*(biga-lams*lams/biga)/r2-2.d0*mus/r2
c	aj(5,5)=(2.d0*lams/biga-5.d0)/r
c	aj contains A(r).  
c	  write(6,*)'B: aj=',aj
c	  pause
c	Employ the Cowling approximation
	g=(g0*1.d+2)*dble(bigr*1.e+5)*1.d-11
	do i=1,5
        do j=1,5
       	aj(i,j)=0.d0
	enddo
	enddo
c       specify id-array, which identifies the index of the second order minor.
        id(1,1)=0
        id(1,2)=1
        id(1,3)=2
        id(1,4)=3
        id(2,1)=id(1,2)
        id(2,2)=0
        id(2,3)=4
        id(2,4)=5
        id(3,1)=id(1,3)
        id(3,2)=id(2,3)
        id(3,3)=0
        id(3,4)=6
        id(4,1)=id(1,4)
        id(4,2)=id(2,4)
        id(4,3)=id(3,4)
        id(4,4)=0
c
	a(1,1)=-2.d0*lams/(r*biga)
	a(1,2)=1.d0/biga
	a(1,3)=fl21*lams/(r*biga)
	a(1,4)=0.d0
	a(2,1)=(4.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,1)-4.d0*denss*g/r2
	a(2,2)=-2.d0/r+(2.d0/r)*lams*a(1,2)
	a(2,3)=-fl21*(2.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,3)
     &	+denss*g*fl21/r2
	a(2,4)=fl21/r
	a(3,1)=-1.d0/r
	a(3,2)=0.d0
	a(3,3)=1.d0/r
	a(3,4)=1.d0/mus
	a(4,1)=-(lams/r)*a(1,1)+2.d0*(mus-biga)/r2+denss*g/r2 
	a(4,2)=-(lams/r)*a(1,2)
	a(4,3)=-(lams/r)*a(1,3)+(fl21*biga-2.d0*mus)/r2
	a(4,4)=-3.d0/r
        do 5 j=1,4
        do 10 k=j,4
        do 15 n=1,4
        m0=id(j,k)
        if(m0.eq.0.or.m0.eq.6) go to 15
        m=id(n,k)
        if(m.eq.0) go to 16
        fac=1.d0
        if(m.eq.6) fac=-1.d0/fl21
        if(m.eq.6) m=1
        if(n.gt.k) fac=-fac
        aj(m0,m)=aj(m0,m)+a(j,n)*fac
c         write(6,*)'m0,m,=',m0,m,'j,n=',j,n,'fac=',fac
16      m=id(j,n)
        if(m.eq.0) go to 15
        fac=1.d0
        if(m.eq.6) fac=-1.d0/fl21
        if(m.eq.6) m=1
        if(j.gt.n) fac=-fac
        aj(m0,m)=aj(m0,m)+a(k,n)*fac
c         write(6,*)'m0,m,=',m0,m,'k,n=',k,n,'fac=',fac
15      continue
10      continue
5       continue
	return
	end 

	subroutine prod(a,c,n,b)
c	Form matrix product b=(A)c
	double complex a,b,c 
	dimension a(n,n),c(n),b(n)
	do 5 i=1,n
	b(i)=0.d0
	do 4 j=1,n
4	b(i)=b(i)+a(i,j)*c(j)
5	continue
	return
	end 

