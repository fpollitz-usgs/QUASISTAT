c	Program QSTAT1A
c***
c	Determine static and postseismic strains using the
c	response functions computed by QSTAT0A and stored in 'qstat0.out'.
c	The source lower and upper depths and observation-pt depth are previously specified
c	in the input to QSTAT0A.  QSTAT1A essentially supplies the remaining source
c	information and writes out the static displacements and strains and
c	time series of postseismic displacements and strains at a number of observation pts.
c***	INPUT
c	This program assumes a source uniformly distributed in a specified depth range.  
c	It handles dipping faults with any shear dislocation, or tensile opening or closing.
c	The source location (more precisely, the endpoint of the lower edge
c	closest to the strike direction) is [flat],[flon].
c	The seismic moment tensor is computed for a finite fault plane with
c	input parameters 
c	fault length ([fbigl]), dip ([dip]), rake ([frake]), strike [fstr],
c	min and max depth range (dmin,dmax), and  slip([fwt]).
c	A number of such planes [iseg], each with different strike, etc.,
c	but the same depth range (dmin,dmax) may be read in.
c	All fault segments have the same dip ([dip]).
c	The dip, strike, and rake follow
c	the convention of Ben Menahem and Singh...
c       If the rake is >181 (<-181) deg., a tensile (compressional) dislocation
c       is assumed.
c	There are [ipts] observation points with locations [olat(i)],[olon(i)],
c	for i=1,...,ipts.
c***	OUTPUT
c	Static and postseismic displacements and strains written out to 
c	'qstat1.out'
c	Output displacement units are cm
c	Output strains in units of 10**(-6) .
c***
	character*80 aread
	real*4 mrr,mtt,mpp,mrt,mrp,mtp  
	parameter (maxsrc=30)
	real*4 mrr0(maxsrc),mtt0(maxsrc),mpp0(maxsrc),mrt0(maxsrc),mrp0(maxsrc),mtp0(maxsrc)  
	real*4 kappa,mu,lam
cOLD	common/source/mrr,mtt,mpp,mrt,mrp,mtp,depth,slat,slon
	real*8 cl,sl,xl0(25001),dxl0(25001),xl1(25001),dxl1(25001)
	real*8 xl2(25001),dxl2(25001)
	parameter (maxlen=15)
	double complex y1(7*maxlen,25000)
	double complex y2(7*maxlen,25000)
	double complex y3(5*maxlen,25000)
	double complex y4(5*maxlen,25000)
c	double complex y1(7*maxlen,11500)
c	double complex y2(7*maxlen,11500)
c	double complex y3(5*maxlen,11500)
c	double complex y4(5*maxlen,11500)
	double complex bfacl,bfacm 
	double complex dfac,dfac1
	real*8 fm(6)
	real*8 x1s(25000),x2s(25000),x3s(25000)
	real*8 dx1s(25000),dx2s(25000),dx3s(25000)
	real*8 ddx1s(25000),ddx2s(25000),ddx3s(25000)
	real*8 x1,x2,x3,dx1,dx2,dx3,ddx1,ddx2,ddx3,faclk,cotd
	parameter (maxpts=10000)
	dimension olat(maxpts),olon(maxpts)
	dimension flat(maxsrc),flon(maxsrc)
     	dimension fbigl(maxsrc),fstr(maxsrc),fdip(maxsrc)
	dimension fwt(maxsrc),frake(maxsrc)
cN
	dimension deltr(40)
	complex*16 rd1S(6,15,40),r1hS(6,15,40),r1tS(6,15,40)
	complex*16 rttS(6,15,40),rppS(6,15,40),rtrS(6,15,40)
	complex*16 rtpS(6,15,40),rprS(6,15,40),rrrS(6,15,40)
	complex*16 romT(6,15,40)
	complex*16 r1hT(6,15,40),r1tT(6,15,40)
	complex*16 rttT(6,15,40),rppT(6,15,40),rtrT(6,15,40)
	complex*16 rtpT(6,15,40),rprT(6,15,40)
	complex*16 vd1S(6),v1hS(6),v1tS(6),vttS(6),vppS(6)
	complex*16 vtrS(6),vtpS(6),vprS(6),vrrS(6)
	complex*16 vomT(6)
	complex*16 v1hT(6),v1tT(6),vttT(6),vppT(6)
	complex*16 vtrT(6),vtpT(6),vprT(6)
	complex*16 yd1S(40),y1hS(40),y1tS(40),yttS(40),yppS(40)
	complex*16 ytrS(40),ytpS(40),yprS(40),yrrS(40)
	complex*16 yomT(40)
	complex*16 y1hT(40),y1tT(40),yttT(40),yppT(40)
	complex*16 ytrT(40),ytpT(40),yprT(40)
	complex*16 sarr(40)
	complex*16 sd1S(6,15,40),s1hS(6,15,40),s1tS(6,15,40),sttS(6,15,40),
     &	sppS(6,15,40)
	complex*16 strS(6,15,40),stpS(6,15,40),sprS(6,15,40),srrS(6,15,40)
	complex*16 somT(6,15,40)
	complex*16 s1hT(6,15,40),s1tT(6,15,40),sttT(6,15,40),sppT(6,15,40)
	complex*16 strT(6,15,40),stpT(6,15,40),sprT(6,15,40)
c--
	real*8 dep,odep
	real*8 ypoint
	complex*16 evaleq
	double complex wkv1(7),wkv2(7),wkv3(5),wkv4(5)
	double complex z1,z2,z3,z4,z5,dz1,dz2,dz3,dz4,dz5
	double complex z1h,z2h,z3h,z4h,z5h,dz1h,dz2h,dz3h,dz4h,dz5h
	double complex disp1,disp1h,disp1t
	double complex ett,epp,err,etp,etr,epr
	double complex dispx1,dispy1,dispz1,exx1,eyy1,ezz1,exy1,exz1,eyz1
	double complex dispx(maxlen),dispy(maxlen),dispz(maxlen)
	double complex exx(maxlen),eyy(maxlen),ezz(maxlen)
	double complex exy(maxlen),exz(maxlen),eyz(maxlen)
cOLD	real*8 ymu,ylam
	real*8 deltp
	real*4 len,kmin
	common/maxl/lmax
	complex*16 oome,ui
	parameter (dstmax=11700.)
c       dstmax has the maximum fault-observation pt. distance.
c
	complex*16 sval
	real*8 lapl
	real*8 tm0,tm1,tm2,tm1o,tm2o,tm1p,tm2p
	real*8 corrf(10)
	common/tmvals/tm1,tm2
	real*8 dspx0(maxpts,11),dspy0(maxpts,11),dspz0(maxpts,11)
	real*8 xtxx0(maxpts,11),xtxy0(maxpts,11),xtxz0(maxpts,11)
	real*8 xtyy0(maxpts,11),xtyz0(maxpts,11),xtzz0(maxpts,11)
	real*8 dspx1(maxpts,11),dspy1(maxpts,11),dspz1(maxpts,11)
	real*8 xtxx1(maxpts,11),xtxy1(maxpts,11),xtxz1(maxpts,11)
	real*8 xtyy1(maxpts,11),xtyz1(maxpts,11),xtzz1(maxpts,11)
c---------------
	parameter (ndep=16)
c---------------
	common/flapl/sval(15)
	real*8 efac
	parameter (efac=dlog(100.d0)/9.d0)
c 
c	ethr=earth's radius in km.
c	bigr=earth radius in 10**6 cm.
	pi=3.1415926535  
	twopi=2.*pi
	ui=cmplx(0.d0,1.0d0)
	rad=180./3.1415926
c	Initialize coseismic and postseismic deformation output arrays
	do i=1,maxpts
	do np=1,11
	dspx0(i,np)=0.d0
        dspy0(i,np)=0.d0
        dspz0(i,np)=0.d0
	xtxx0(i,np)=0.d0
	xtxy0(i,np)=0.d0
        xtxz0(i,np)=0.d0
        xtyy0(i,np)=0.d0
        xtyz0(i,np)=0.d0
        xtzz0(i,np)=0.d0
	dspx1(i,np)=0.d0
        dspy1(i,np)=0.d0
        dspz1(i,np)=0.d0
	xtxx1(i,np)=0.d0
	xtxy1(i,np)=0.d0
        xtxz1(i,np)=0.d0
        xtyy1(i,np)=0.d0
        xtyz1(i,np)=0.d0
        xtzz1(i,np)=0.d0
	enddo
	enddo
c
	write(6,*)'# year of earthquake, year obs.#1, year obs.#2 (yrs), viscosity multiplier'
	read(5,5) aread
	read(5,*) tm0,tm1,tm2,vmult
	write(6,*) tm0,tm1,tm2,vmult
	tm1=(tm1-tm0)/vmult
	tm2=(tm2-tm0)/vmult
	write(6,*)'  '
	write(6,*)'  '
	write(6,*)'NOTE: The source depth runs over ',ndep,' discrete values'
	write(6,*)'processed by QSTAT0A that spans the depth range of each fault plane'
c 
	write(6,*)'finite fault with [iseg] # segments.  iseg=?'
	read(5,5) aread
	read(5,*) iseg
c	iobs=0
	read(5,5) aread
	do 39 i=1,iseg
	write(6,*)'segment #',i,'Lower edge lat,lon(deg.),length(km)'
	write(6,*) 'strike(deg.),dip(deg.),rake(deg.),slip (cm)?'
	read(5,*) flat(i),flon(i),fbigl(i),fstr(i),fdip(i),frake(i),fwt(i)
	flat(i)=(pi/2.-flat(i)/rad)
	flon(i)=flon(i)/rad
	fstr(i)=fstr(i)/rad 
	frake(i)=frake(i)/rad 
	if(fwt(i).lt.0.0) then
	write(6,*)'Read in moment tensor components'
	write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp (10^20 N m)?'
	write(6,*)'where r=r hat=Up, t=theta hat=South, p=phi hat=East'
	read(5,*) mrr0(i),mtt0(i),mpp0(i),mrt0(i),mrp0(i),mtp0(i)
	endif
39	continue
c 
	write(6,*)'number of observation points [ipts]?'
	read(5,5) aread
	read(5,*) ipts
	write(6,*)'ipts=',ipts
	do 40 i=1,ipts 
c	write(6,*)'point #',i,'lat,lon(deg.)?'
	read(5,*) olat(i),olon(i)
	olat(i)=(pi/2.-olat(i)/rad)
	olon(i)=olon(i)/rad
40	continue
c
	open(2,file='qstat0.out',form='unformatted')
	rewind(2)
	read(2) dmax,dmin
c
	do 45 j=1,ndep
c*-*-*-*-*-*-*-*-*-*-*
	read(2) ethr,sdep
	  write(6,*)'ethr=',ethr
	rmin=1.-dmax/ethr
        rmax=1.-dmin/ethr
	  write(6,*)'dmax,dmin=',dmax,dmin,'rmax,rmin=',rmax,rmin
	bigr=ethr/10.
	read(2) lmin,lmax
	  write(6,*)'lmin,lmax=',0,lmax
	kmin=twopi/(real(lmax)+0.5)
	read(2) odep
	read(2) kappa,mu
	  write(6,*)'kappa,mu=',kappa,mu
	lam=kappa-2.*mu/3.
c	ymu=dble(mu)
c	ylam=dble(lam)
c	  write(6,*)'ymu,ylam=',ymu,ylam
	read(2) sval
		write(6,*)'sval=',sval
c--
	write(6,*)'Reading in spheroidal and toroidal motion depth functions'
	do 55 iom=1,15
c		write(6,*)'reading iom=',iom,'out of',15
	do 50 lr=lmin,lmax+1
	l=lr
c	Do l=0 case last.
	if(lr.eq.(lmax+1)) l=0
c		write(6,*)'l,lr=',l,lr
	do 69 idec=1,5
	k=15*idec-15+iom
	read(2) ldum,y1(k,lr),y2(k,lr),y3(k,lr),y4(k,lr)
	read(2) bfacl,bfacm
69	continue
	if(lr.eq.(lmax+1)) go to 50
	do 58 idec=6,7
	k=15*idec-15+iom
	read(2) ldum,y1(k,lr),y2(k,lr)
c--
58	continue
c	  pause
c	Note: y1 - y4(1 thru 80,l) are Spheroidal motion coeff.
c	      y1 - y2(81 thru 112,l) are Toroidal motion coeff.
50	continue
55	continue
	write(6,*)'Done reading in spheroidal and toroidal motion depth functions'
cTEST
c		lmax=2000
	write(6,*)'finished reading displacement coefficients'
c	Divide moment tensor by nmesh (# elements in fault subdivision).
	nmesh2=41
c	nmesh2-1=# horizontal length points 
5	format(a80)
c*****
c	There is coupling between static and postseismic (exponentially decaying) terms
c	when LAPL is called to evaluate the inverse Laplace transform.  Evaluate the
c	effect of the static term on the postseismic terms.
	do iom=1,15
	dispx(iom)=1.d0/sval(iom)
	enddo
c       Do postseismic displacements for 10 time intervals
	do np=1,10
	tm1p=tm1*(0.129155d0*dexp(efac*dble(np-1)))/3.16881d0
	tm2p=tm2*(0.129155d0*dexp(efac*dble(np-1)))/3.16881d0
	corrf(np)=lapl(dispx,tm2p)-lapl(dispx,tm1p)
	write(6,*)'corrf(',np,')=',corrf(np)
	enddo
c*****	
c	Next, calculate response Greens functions for m=0, 1, and 2 at 
c	specific values of DELTA.
	write(6,*)'calculate response Greens functions for m=0, 1, and 2 at' 
	write(6,*)'specific values of DELTA'
	do i=1,40
	deltr(i)=dstmax**(0.025*real(i))/ethr
	enddo
	do 105 i=1,40
	deltpf=deltr(i)
	write(6,*)'i=',i,' DELTA=',deltpf
	cotd=cos(deltpf)/sin(deltpf)
	slf=sin(deltpf)
	call lgndr0(deltpf,xl0,dxl0)
	call lgndr1(deltpf,xl1,dxl1)
	call lgndr2(deltpf,xl2,dxl2) 
c	Store Legendre functions for later use.
	do 364 lk=lmin,lmax+1
	faclk=sqrt(real(lk*(lk+1)))
	x1=real(xl0(lk))
	x2=real(xl1(lk))
	x3=real(xl2(lk))
	dx1=real(dxl0(lk))
	dx2=real(dxl1(lk))
	dx3=real(dxl2(lk))
	ddx1=(-faclk*faclk)*x1-cotd*dx1
	ddx2=(1./(slf*slf)-faclk*faclk)*x2-cotd*dx2 
	ddx3=(4./(slf*slf)-faclk*faclk)*x3-cotd*dx3
	x1s(lk)=x1
	x2s(lk)=x2
	x3s(lk)=x3
	dx1s(lk)=dx1
	dx2s(lk)=dx2
	dx3s(lk)=dx3
	ddx1s(lk)=ddx1
	ddx2s(lk)=ddx2
	ddx3s(lk)=ddx3
364	continue 
	do 107 im=1,6
	fm(1)=0.
	fm(2)=0.
	fm(3)=0.
	fm(4)=0.
	fm(5)=0.
	fm(6)=0.
	fm(im)=1.	
	do 106 iom=1,15
	rd1S(im,iom,i)=0.
	r1hS(im,iom,i)=0.
	r1tS(im,iom,i)=0.
	rttS(im,iom,i)=0.
	rppS(im,iom,i)=0.
	rtrS(im,iom,i)=0.
	rtpS(im,iom,i)=0.
	rprS(im,iom,i)=0.
	rrrS(im,iom,i)=0.
	romT(im,iom,i)=0.
	r1hT(im,iom,i)=0.
	r1tT(im,iom,i)=0.
	rttT(im,iom,i)=0.
	rppT(im,iom,i)=0.
	rtrT(im,iom,i)=0.
	rtpT(im,iom,i)=0.
	rprT(im,iom,i)=0.
	cp=1.
	sp=1.
	cp2=1.
	sp2=1.
	do 371 lr=lmin,lmax+1
	l=lr
c	Do l=0 case last.
	if(lr.eq.(lmax+1)) l=0
c	With facf apply tapering over last fifth of frequency range.
	fl=real(l)
	flmax=real(lmax)
	iw0m=lmax/5
	iw0x=lmax-iw0m+1
	fiw0m=real(iw0m)
	facf=1.0
	if(l.ge.iw0x) facf=1.0-cos((flmax-fl)/(fiw0m-1.) * real(pi)/2.)
c
	facl=sqrt(real(l*(l+1)))
c
	dfac=facf*(1.e+6)/ethr**2
	dfac1=dfac 
c	Now account for the fact that the observation radius is generally
c	different from the surface radius.  This affects displ calculations
c	(for example, y2 and y4-arrays have functions
c	y1(r)/r and y3(r)/r, where r is dimensionless radius,
c	which equals 1 only at the surface), but not strain calculations.
	dfac=dfac*real(odep)
c	SPHEROIDAL MODES
	do 375 idec=1,5
	k=15*idec-15+iom
	wkv1(idec)=y1(k,lr)
	wkv2(idec)=y2(k,lr)
	wkv3(idec)=y3(k,lr)
	wkv4(idec)=y4(k,lr)
375	continue
c		write(6,*)'wkv4=',wkv4
c	Note extra factor of 2 in z2 - z5 components in order to account
c	for m=-1 and m=-2 contributions.
c * *	Moment tensor excitation
	z1=fm(1)*(wkv2(1))+fm(2)*(wkv2(2))
	z2=2.*fm(3)*(wkv2(3))
	z3=-2.*fm(4)*(wkv2(3))
	z4=2.*fm(5)*(wkv2(4))
	z5=-2.*fm(6)*(wkv2(4))
c		write(6,*)'fm(5),wkv2(4),z4=',fm(5),wkv2(4),z4
c		write(6,*)'fm(6),wkv2(4),z5=',fm(6),wkv2(4),z5
	dz1=fm(1)*(wkv1(1))+fm(2)*(wkv1(2))
	dz2=2.*fm(3)*(wkv1(3))
	dz3=-2.*fm(4)*(wkv1(3))
	dz4=2.*fm(5)*(wkv1(4))
	dz5=-2.*fm(6)*(wkv1(4))
	dz1h=fm(1)*(wkv3(1))+fm(2)*(wkv3(2))
	dz2h=2.*fm(3)*(wkv3(3))
	dz3h=-2.*fm(4)*(wkv3(3))
	dz4h=2.*fm(5)*(wkv3(4))
	dz5h=-2.*fm(6)*(wkv3(4))
	z1h=fm(1)*(wkv4(1))+fm(2)*(wkv4(2))
	z2h=2.*fm(3)*(wkv4(3))
	z3h=-2.*fm(4)*(wkv4(3))
	z4h=2.*fm(5)*(wkv4(4))
	z5h=-2.*fm(6)*(wkv4(4))
c * *
c	  write(6,*)'fm(6),wkv4(4),z5h=',fm(6),wkv4(4),z5h
	if(l.eq.0) then
	x1=sqrt(0.25/pi)
	x2=0.
	x3=0.
	dx1=0.
	dx2=0.
	dx3=0.
	ddx1=0.
	ddx2=0.
	ddx3=0.
	else
	x1=x1s(lr)
	x2=x2s(lr)
	x3=x3s(lr)
	dx1=dx1s(lr)
	dx2=dx2s(lr)
	dx3=dx3s(lr)
	ddx1=ddx1s(lr)
	ddx2=ddx2s(lr)
	ddx3=ddx3s(lr)
	endif
c 
	rd1S(im,iom,i)=rd1S(im,iom,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac
	r1hS(im,iom,i)=r1hS(im,iom,i)+(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac
c		write(6,*)'r1hS(',im,iom,i,')=',r1hS(im,iom,i)
	r1tS(im,iom,i)=r1tS(im,iom,i)+(-z2h*x2*sp
     &	+z3h*x2*cp -2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac/slf 
	rttS(im,iom,i)=rttS(im,iom,i)+(z1h*ddx1 +z2h*ddx2*cp
     &	+z3h*ddx2*sp +z4h*ddx3*cp2 +z5h*ddx3*sp2)*dfac1/bigr 
	rttS(im,iom,i)=rttS(im,iom,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	rppS(im,iom,i)=rppS(im,iom,i)+(-z2h*x2*cp -z3h*x2*sp 
     &	-4.*z4h*x3*cp2 -4.*z5h*x3*sp2)*(dfac1/slf)/(slf*bigr) 
	rppS(im,iom,i)=rppS(im,iom,i)+(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1*cotd/bigr 
	rppS(im,iom,i)=rppS(im,iom,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	rtrS(im,iom,i)=rtrS(im,iom,i)+(dz1h*dx1 +dz2h*dx2*cp
     &	+dz3h*dx2*sp +dz4h*dx3*cp2 +dz5h*dx3*sp2)*dfac1/(2.*bigr) 
	rtrS(im,iom,i)=rtrS(im,iom,i)-(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1/(2.*bigr) 
	rtrS(im,iom,i)=rtrS(im,iom,i)+(z1*dx1 +z2*dx2*cp
     &	+z3*dx2*sp +z4*dx3*cp2 +z5*dx3*sp2)*dfac1/(2.*bigr) 
	rtpS(im,iom,i)=rtpS(im,iom,i)+(-z2h*dx2*sp
     &	+z3h*dx2*cp -2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*bigr*slf) 
	rtpS(im,iom,i)=rtpS(im,iom,i)-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rtpS(im,iom,i)=rtpS(im,iom,i)+(-z2h*dx2*sp +z3h*dx2*cp 
     &	-2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*slf*bigr)
	rprS(im,iom,i)=rprS(im,iom,i)+(-z2*x2*sp
     &	+z3*x2*cp -2.*z4*x3*sp2 +2.*z5*x3*cp2)*dfac1/(2.*slf*bigr) 
	rprS(im,iom,i)=rprS(im,iom,i)+(-dz2h*x2*sp +dz3h*x2*cp 
     &	-2.*dz4h*x3*sp2 +2.*dz5h*x3*cp2)*dfac1/(2.*bigr*slf) 
	rprS(im,iom,i)=rprS(im,iom,i)-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1/(2.*bigr*slf) 
c		if(l.gt.7490) write(6,*)'rprS(',im,iom,i,')=',rprS(im,iom,i)
	rrrS(im,iom,i)=rrrS(im,iom,i)+(dz1*x1 +dz2*x2*cp
     &	+dz3*x2*sp +dz4*x3*cp2 +dz5*x3*sp2)*dfac1/bigr 
c	Done with Spheroidal motion component of degree l.
c	TOROIDAL MODES
	do 475 idec=6,7
	k=15*idec-15+iom
	wkv1(idec)=y1(k,lr)
	wkv2(idec)=y2(k,lr)
475	continue
	z1=-2.*fm(4)*(wkv2(6))
	z2=-2.*fm(3)*(wkv2(6))
	z3=-2.*fm(6)*(wkv2(7))
	z4=-2.*fm(5)*(wkv2(7))
	dz1=-2.*fm(4)*(wkv1(6)+wkv2(6))
	dz2=-2.*fm(3)*(wkv1(6)+wkv2(6))
	dz3=-2.*fm(6)*(wkv1(7)+wkv2(7))
	dz4=-2.*fm(5)*(wkv1(7)+wkv2(7))
	romT(im,iom,i)=romT(im,iom,i)-facl*facl*((-z1*cp-z2*sp)*x2
     &	-(z3*cp2+z4*sp2)*x3)*dfac/(2.*bigr)
	r1tT(im,iom,i)=r1tT(im,iom,i)+(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac
	r1hT(im,iom,i)=r1hT(im,iom,i)-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac/slf 
	rttT(im,iom,i)=rttT(im,iom,i)-(-z1*dx2*sp +z2*dx2*cp
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf) 
	rttT(im,iom,i)=rttT(im,iom,i)+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rppT(im,iom,i)=rppT(im,iom,i)+(-z1*dx2*sp +z2*dx2*cp 
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf)
	rppT(im,iom,i)=rppT(im,iom,i)-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rtpT(im,iom,i)=rtpT(im,iom,i)+(z1*ddx2*cp +z2*ddx2*sp
     &	+z3*ddx3*cp2 +z4*ddx3*sp2)*dfac1/(2.*bigr)
	rtpT(im,iom,i)=rtpT(im,iom,i)-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1*cotd/(2.*bigr)
	rtpT(im,iom,i)=rtpT(im,iom,i)-(-z1*x2*cp -z2*x2*sp
     &	-4.*z3*x3*cp2 -4.*z4*x3*sp2)*(dfac1/slf)/(2.*bigr*slf)
	rtrT(im,iom,i)=rtrT(im,iom,i)-(-dz1*x2*sp +dz2*x2*cp
     &	-2.*dz3*x3*sp2 +2.*dz4*x3*cp2)*dfac1/(2.*bigr*slf) 
	rtrT(im,iom,i)=rtrT(im,iom,i)+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1/(2.*bigr*slf) 
	rprT(im,iom,i)=rprT(im,iom,i)+(dz1*dx2*cp +dz2*dx2*sp
     &	+dz3*dx3*cp2 +dz4*dx3*sp2)*dfac1/(2.*bigr)
	rprT(im,iom,i)=rprT(im,iom,i)-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1/(2.*bigr)
c	Done with Toroidal motion component of degree l.
371	continue
106	continue
107	continue
105	continue
	write(6,*)'Finished determining response Greens functions'
	write(6,*)'rd1S(1,1,35)=',rd1S(1,1,35)
	write(6,*)'rd1S(1,8,35)=',rd1S(1,8,35)
	write(6,*)'rd1S(1,15,35)=',rd1S(1,15,35)
c	pause
c*****
c	Set up spline interpolation arrays for the response functions.
	do 111 im=1,6
	write(6,*)'setting up spline interpolation arrays for response functions'
	write(6,*)'im=',im
	do 112 iom=1,15
	write(6,*)'iom=',iom
	do 113 i=1,40
	yd1S(i)=dcmplx(rd1S(im,iom,i))
	y1hS(i)=dcmplx(r1hS(im,iom,i))
	y1tS(i)=dcmplx(r1tS(im,iom,i))
	yttS(i)=dcmplx(rttS(im,iom,i))
	yppS(i)=dcmplx(rppS(im,iom,i))
	ytrS(i)=dcmplx(rtrS(im,iom,i))
	ytpS(i)=dcmplx(rtpS(im,iom,i))
	yprS(i)=dcmplx(rprS(im,iom,i))
	yrrS(i)=dcmplx(rrrS(im,iom,i))
	yomT(i)=dcmplx(romT(im,iom,i))
	y1hT(i)=dcmplx(r1hT(im,iom,i))
	y1tT(i)=dcmplx(r1tT(im,iom,i))
	yttT(i)=dcmplx(rttT(im,iom,i))
	yppT(i)=dcmplx(rppT(im,iom,i))
	ytrT(i)=dcmplx(rtrT(im,iom,i))
	ytpT(i)=dcmplx(rtpT(im,iom,i))
	yprT(i)=dcmplx(rprT(im,iom,i))
113	continue
	call splneq(40,yd1S,sarr)
	do i=1,40
	sd1S(im,iom,i)=sarr(i)
	enddo
	call splneq(40,y1hS,sarr)
	do i=1,40
	s1hS(im,iom,i)=sarr(i)
c		write(6,*)'s1hS(',im,iom,i,')=',s1hS(im,iom,i)
	enddo
	call splneq(40,y1tS,sarr)
	do i=1,40
	s1tS(im,iom,i)=sarr(i)
	enddo
	call splneq(40,yttS,sarr)
	do i=1,40
	sttS(im,iom,i)=sarr(i)
	enddo
	call splneq(40,yppS,sarr)
	do i=1,40
	sppS(im,iom,i)=sarr(i)
	enddo
	call splneq(40,ytrS,sarr)
	do i=1,40
	strS(im,iom,i)=sarr(i)
	enddo
	call splneq(40,ytpS,sarr)
	do i=1,40
	stpS(im,iom,i)=sarr(i)
	enddo
	call splneq(40,yprS,sarr)
	do i=1,40
	sprS(im,iom,i)=sarr(i)
	enddo
	call splneq(40,yrrS,sarr)
	do i=1,40
	srrS(im,iom,i)=sarr(i)
	enddo
c
	call splneq(40,yomT,sarr)
	do i=1,40
	somT(im,iom,i)=sarr(i)
	enddo
	call splneq(40,y1hT,sarr)
	do i=1,40
	s1hT(im,iom,i)=sarr(i)
	enddo
	call splneq(40,y1tT,sarr)
	do i=1,40
	s1tT(im,iom,i)=sarr(i)
	enddo
	call splneq(40,yttT,sarr)
	do i=1,40
	sttT(im,iom,i)=sarr(i)
	enddo
	call splneq(40,yppT,sarr)
	do i=1,40
	sppT(im,iom,i)=sarr(i)
	enddo
	call splneq(40,ytrT,sarr)
	do i=1,40
	strT(im,iom,i)=sarr(i)
	enddo
	call splneq(40,ytpT,sarr)
	do i=1,40
	stpT(im,iom,i)=sarr(i)
	enddo
	call splneq(40,yprT,sarr)
	do i=1,40
	sprT(im,iom,i)=sarr(i)
	enddo
c
112	continue
111	continue
	write(6,*)'Finished determining spline interpolations of response functions'
c	pause
c*****	
c	Begin sweep over observation points.
	do 60 ip=1,ipts 
		write(6,*)'Doing obs pt #',i,' out of ',ipts
	do 84 iom=1,15
		write(6,*)'Doing Laplace transform parameter # ',iom,' out of ',15
	dispx(iom)=0.
	dispy(iom)=0.
	dispz(iom)=0.
	exx(iom)=0.
	eyy(iom)=0.
	exy(iom)=0.
	exz(iom)=0.
	eyz(iom)=0.
	ezz(iom)=0.
c	Compute displacements and strains
c	at points specified by
c	arrays olat and olon.  Begin summation over fault segments.
	jf=0
65	jf=jf+1
	if(jf.gt.iseg) go to 84
c		write(6,*)'jf=',jf
	write(6,*) flat(jf),flon(jf),fbigl(jf),fstr(jf)
	sstr=sin(fstr(jf))
	cstr=cos(fstr(jf))
	s2str=2.*sstr*cstr
	c2str=cstr*cstr-sstr*sstr 
	cdip=cos(fdip(jf)/rad)
	sdip=sin(fdip(jf)/rad)
	c2dip=cdip*cdip-sdip*sdip
		 dh=(rmax-rmin)*(cdip/sdip)/real(ndep)
		 bigh=-dh/2.+dh*real(j)
	srak=sin(frake(jf))
	crak=cos(frake(jf))
	dlen=(fbigl(jf)/ethr)/real(nmesh2-1)
	dlon=olon(ip)-flon(jf)
c	Find angular distance and azimuth to observation point from
c	the initial point on the fault segment.
c	Note: delta is the angular distance from the earthquake.
c	Note: phi is the azimuth meazsured positive counterclockwise
c	from south.
	cdelt=cos(flat(jf))*cos(olat(ip))+sin(flat(jf))*
     &	sin(olat(ip))*cos(dlon)
	if (cdelt.le.0.9999) delta=acos(cdelt)
	if(cdelt.gt.0.9999) delta=sqrt((flat(jf)-olat(ip))**2+
     &	(dlon*sin(flat(jf)))**2)
	spsi=sin(dlon)*sin(olat(ip))/sin(delta)
	cpsi=(cos(olat(ip))-cos(flat(jf))*cdelt)/(sin(flat(jf))*sin(delta))
	phi=pi-atan2(spsi,cpsi) 
	dwrite=rad*delta
	pwrite=rad*phi
	  write(6,*)'olat(',i,')=',olat(ip),'olon=',olon(ip)
	  write(6,*)'flat(',jf,')=',flat(jf),'flon=',flon(jf)
	  write(6,*)'delta,phi=',dwrite,'deg.',pwrite 
	cphi=cos(phi)
	sphi=sin(phi)
cNEW
	spsi=sin(dlon)*sin(flat(jf))/sin(delta)
	cpsi=(cos(flat(jf))-cos(olat(ip))*cdelt)/(sin(olat(ip))*sin(delta))
	peps=atan2(spsi,cpsi)
	cphio=cos(peps)
	sphio=sin(peps)
c--
	if(jf.gt.1) go to 59
	cphi1=cphi
	sphi1=sphi
	delt1=delta 
59	delta=delta*(6371./ethr)
	  write(6,*)'ethr=',ethr
	  write(6,*)'delta (after mult)=',delta
c	Angular distances are a factor of (6371./ethr) larger on the
c	earth with radius ethr km.
	  write(6,*)'after 59, delta=',delta
c	Integrate displacements over fault elements.
c ***	Determine moment tensor elements for fault element (ilen).
c	Units of moment tensor are 10**20 N-m.
	  write(6,*)'mu,lam=',mu,lam
	u1=fwt(jf)
c       u1=magnitude of slip.
	amesh=real(nmesh2-1)*real(ndep)
c	Use input moment tensor elements if scalar moment < 0
	if(fwt(jf).lt.0.0) then
	mrr=mrr0(jf)/amesh
	mtt=mtt0(jf)/amesh
	mpp=mpp0(jf)/amesh
	mrt=mrt0(jf)/amesh
	mrp=mrp0(jf)/amesh
	mtp=mtp0(jf)/amesh
	go to 51
	endif
c	Note that shear modulus [mu] is already input in units of 10**10 Pa.
c	Moment tensor from Ben Menahem and Singh, eqn. 4.115b for
c	shear dislocation, and derived from eqn. (4.101), (4.110), and
c	(4.113) for a tensile dislocation.
	if(abs(frake(jf)).gt.3.159046) go to 49
c	Next line is shear moment.
	shrm1=u1*fbigl(jf)*((dmax-dmin)/sdip)*mu*(1.e-6)/amesh
c--
	p1=srak*sdip*cdip*s2str+crak*sdip*c2str
	p2=crak*sdip*s2str-srak*sdip*cdip*c2str
	p3=-crak*cdip*sstr+srak*c2dip*cstr
	p4=srak*c2dip*sstr+crak*cdip*cstr
	p5=srak*sdip*cdip
	mrr=shrm1*2.*p5
	mtt=-shrm1*(p2+p5)
	mpp=shrm1*(p2-p5)
	mrt=-shrm1*p4
	mrp=-shrm1*p3
	mtp=-shrm1*p1
	go to 51
49      srak=1.
        if(frake(jf).lt.-3.159046) srak=-1.
          if(srak.eq.1.0) write(6,*)'tensile dislocation'
          if(srak.eq.-1.0) write(6,*)'compressional dislocation'
        shrm1=2.*u1*fbigl(jf)*((dmax-dmin)/sdip)*srak*mu*(1.e-6)/amesh
        mrr=shrm1*cdip**2
        mtt=shrm1*(sdip*sstr)**2
        mpp=shrm1*(sdip*cstr)**2
        mrt=shrm1*sdip*cdip*sstr
        mrp=shrm1*sdip*cdip*cstr
        mtp=shrm1*sdip**2*sstr*cstr
        shrm2=0.5*shrm1*(lam/mu)
        mrr=mrr+shrm2
        mtt=mtt+shrm2
        mpp=mpp+shrm2
51	continue
	  write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp='
	  write(6,*) mrr,mtt,mpp,mrt,mrp,mtp
c ***
c	Determine multiplying factors for the cases m=0,1,2 
	fm(1)=dble((1./(2.*mu*mu+3.*lam*mu))*((lam+mu)*mrr-(lam/2.)*(mtt+mpp)))
	fm(2)=dble((1./(2.*mu*mu+3.*lam*mu))*(-lam*mrr+(lam/2.+mu)*(mtt+mpp)))
	fm(3)=dble((1./(2.*mu))*mrt)
	fm(4)=dble(-(1./(2.*mu))*mrp)
	fm(5)=dble((1./(2.*mu))*(mtt-mpp))
	fm(6)=dble(-(1./mu)*mtp)
	len=-dlen/2.
	ilen=0
83	ilen=ilen+1
	if(ilen.eq.nmesh2) go to 65
c	  write(6,*)'ilen=',ilen
	len=len+dlen
	bigx=delta*sphi+len*sstr+bigh*cstr
	bigy=-delta*cphi+len*cstr-bigh*sstr
	deltp=dble(sqrt(bigx*bigx+bigy*bigy))
	phip=atan2(bigy,bigx)+pi/2.	
	  deltpw=real(deltp)*rad
	  phipw=real(phip)*rad 
c	  write(6,*)'bigx,bigy,deltp,phip=',bigx,bigy,deltp,phipw 
	cp=cos(phip)
	sp=sin(phip)
	cp2=cos(2.*phip)
	sp2=sin(2.*phip)
	cl=dcos(deltp)
	sl=dsin(deltp)
	cotd=cl/sl
	slf=real(sl)
	deltpf=real(deltp)
c*****
c	Now interpolate response functions at this deltpf.
	ypoint=dble(40.*log(deltpf*ethr)/log(dstmax))
	do 114 im=1,6
	do i=1,40
	sarr(i)=sd1S(im,iom,i)
	yd1S(i)=rd1S(im,iom,i)
	enddo
	vd1S(im)=evaleq(ypoint,40,yd1S,sarr)
	do i=1,40
	sarr(i)=s1hS(im,iom,i)
	y1hS(i)=r1hS(im,iom,i)
	enddo
	v1hS(im)=evaleq(ypoint,40,y1hS,sarr)
	do i=1,40
	sarr(i)=s1tS(im,iom,i)
	y1tS(i)=r1tS(im,iom,i)
	enddo
	v1tS(im)=evaleq(ypoint,40,y1tS,sarr)
	do i=1,40
	sarr(i)=sttS(im,iom,i)
	yttS(i)=rttS(im,iom,i)
	enddo
	vttS(im)=evaleq(ypoint,40,yttS,sarr)
	do i=1,40
	sarr(i)=sppS(im,iom,i)
	yppS(i)=rppS(im,iom,i)
	enddo
	vppS(im)=evaleq(ypoint,40,yppS,sarr)
	do i=1,40
	sarr(i)=strS(im,iom,i)
	ytrS(i)=rtrS(im,iom,i)
	enddo
	vtrS(im)=evaleq(ypoint,40,ytrS,sarr)
	do i=1,40
	sarr(i)=stpS(im,iom,i)
	ytpS(i)=rtpS(im,iom,i)
	enddo
	vtpS(im)=evaleq(ypoint,40,ytpS,sarr)
	do i=1,40
	sarr(i)=sprS(im,iom,i)
	yprS(i)=rprS(im,iom,i)
	enddo
	vprS(im)=evaleq(ypoint,40,yprS,sarr)
	do i=1,40
	sarr(i)=srrS(im,iom,i)
	yrrS(i)=rrrS(im,iom,i)
	enddo
	vrrS(im)=evaleq(ypoint,40,yrrS,sarr)
c
	do i=1,40
	sarr(i)=somT(im,iom,i)
	yomT(i)=romT(im,iom,i)
	enddo
	vomT(im)=evaleq(ypoint,40,yomT,sarr)
	do i=1,40
	sarr(i)=s1hT(im,iom,i)
	y1hT(i)=r1hT(im,iom,i)
	enddo
	v1hT(im)=evaleq(ypoint,40,y1hT,sarr)
	do i=1,40
	sarr(i)=s1tT(im,iom,i)
	y1tT(i)=r1tT(im,iom,i)
	enddo
	v1tT(im)=evaleq(ypoint,40,y1tT,sarr)
	do i=1,40
	sarr(i)=sttT(im,iom,i)
	yttT(i)=rttT(im,iom,i)
	enddo
	vttT(im)=evaleq(ypoint,40,yttT,sarr)
	do i=1,40
	sarr(i)=sppT(im,iom,i)
	yppT(i)=rppT(im,iom,i)
	enddo
	vppT(im)=evaleq(ypoint,40,yppT,sarr)
	do i=1,40
	sarr(i)=strT(im,iom,i)
	ytrT(i)=rtrT(im,iom,i)
	enddo
	vtrT(im)=evaleq(ypoint,40,ytrT,sarr)
	do i=1,40
	sarr(i)=stpT(im,iom,i)
	ytpT(i)=rtpT(im,iom,i)
	enddo
	vtpT(im)=evaleq(ypoint,40,ytpT,sarr)
	do i=1,40
	sarr(i)=sprT(im,iom,i)
	yprT(i)=rprT(im,iom,i)
	enddo
	vprT(im)=evaleq(ypoint,40,yprT,sarr)
114	continue
c*****
cNEW
	bigx=delta*sphio+len*sstr+bigh*cstr
	bigy=-delta*cphio+len*cstr-bigh*sstr
	deltp=dble(sqrt(bigx*bigx+bigy*bigy))
	phip=atan2(bigy,bigx)+pi/2.	
	  gamma1=sin(phip)
	  gamma2=-cos(phip)
c--
	disp1=vd1S(1)*fm(1)+vd1S(2)*fm(2)+vd1S(3)*fm(3)*cp
     &	+vd1S(4)*fm(4)*sp+vd1S(5)*fm(5)*cp2+vd1S(6)*fm(6)*sp2
	disp1h=v1hS(1)*fm(1)+v1hS(2)*fm(2)+(v1hS(3)+v1hT(3))*fm(3)*cp
     &	+(v1hS(4)+v1hT(4))*fm(4)*sp+(v1hS(5)+v1hT(5))*fm(5)*cp2
     &	+(v1hS(6)+v1hT(6))*fm(6)*sp2
		write(6,*)'disp1h=',disp1h
	disp1t=(v1tS(3)+v1tT(3))*fm(3)*sp
     &	+(v1tS(4)+v1tT(4))*fm(4)*cp+(v1tS(5)+v1tT(5))*fm(5)*sp2
     &	+(v1tS(6)+v1tT(6))*fm(6)*cp2
	ett=vttS(1)*fm(1)+vttS(2)*fm(2)+(vttS(3)+vttT(3))*fm(3)*cp
     &	+(vttS(4)+vttT(4))*fm(4)*sp+(vttS(5)+vttT(5))*fm(5)*cp2
     &	+(vttS(6)+vttT(6))*fm(6)*sp2
	epp=vppS(1)*fm(1)+vppS(2)*fm(2)+(vppS(3)+vppT(3))*fm(3)*cp
     &	+(vppS(4)+vppT(4))*fm(4)*sp+(vppS(5)+vppT(5))*fm(5)*cp2
     &	+(vppS(6)+vppT(6))*fm(6)*sp2
	etr=vtrS(1)*fm(1)+vtrS(2)*fm(2)+(vtrS(3)+vtrT(3))*fm(3)*cp
     &	+(vtrS(4)+vtrT(4))*fm(4)*sp+(vtrS(5)+vtrT(5))*fm(5)*cp2
     &	+(vtrS(6)+vtrT(6))*fm(6)*sp2
	etp=(vtpS(3)+vtpT(3))*fm(3)*sp
     &	+(vtpS(4)+vtpT(4))*fm(4)*cp+(vtpS(5)+vtpT(5))*fm(5)*sp2
     &	+(vtpS(6)+vtpT(6))*fm(6)*cp2
	epr=(vprS(3)+vprT(3))*fm(3)*sp
     &	+(vprS(4)+vprT(4))*fm(4)*cp+(vprS(5)+vprT(5))*fm(5)*sp2
     &	+(vprS(6)+vprT(6))*fm(6)*cp2
	err=vrrS(1)*fm(1)+vrrS(2)*fm(2)+vrrS(3)*fm(3)*cp
     &	+vrrS(4)*fm(4)*sp+vrrS(5)*fm(5)*cp2+vrrS(6)*fm(6)*sp2
c	Rotate strain tensor into x-y-z Cartesian coordinates.
c	Also rotate displacements into x-y coordinates.
	dispx1=gamma1*disp1h-gamma2*disp1t
	dispy1=gamma2*disp1h+gamma1*disp1t
	dispz1=disp1
	exx1=gamma1**2*ett+gamma2**2*epp-2.*gamma1*gamma2*etp
	exy1=gamma1*gamma2*(ett-epp)+(gamma1**2-gamma2**2)*etp
	eyy1=gamma2**2*ett+gamma1**2*epp+2.*gamma1*gamma2*etp
	exz1=etr*gamma1-epr*gamma2
	eyz1=etr*gamma2+epr*gamma1
	ezz1=err
c
70	dispx(iom)=dispx(iom)+dispx1 
	dispy(iom)=dispy(iom)+dispy1 
	dispz(iom)=dispz(iom)+dispz1 
	exx(iom)=exx(iom)+exx1
	eyy(iom)=eyy(iom)+eyy1
	exy(iom)=exy(iom)+exy1
	exz(iom)=exz(iom)+exz1
	eyz(iom)=eyz(iom)+eyz1
	ezz(iom)=ezz(iom)+ezz1
	go to 83
84	continue
		write(6,*)'After 84'
c--
	do iom=1,15
c	Multiply everything by 1/s to account for step-function
c	character of source.
	oome=1.d0/sval(iom)
	dispx(iom)=dispx(iom)*oome
	dispy(iom)=dispy(iom)*oome
	dispz(iom)=dispz(iom)*oome
	exx(iom)=exx(iom)*oome
	eyy(iom)=eyy(iom)*oome
	ezz(iom)=ezz(iom)*oome
	exz(iom)=exz(iom)*oome
	eyz(iom)=eyz(iom)*oome
	exy(iom)=exy(iom)*oome
	enddo
c	First update static displacements and strains
	tm0=0.d0
	dspx0(ip,1)=dspx0(ip,1)+lapl(dispx,tm0)
        dspy0(ip,1)=dspy0(ip,1)+lapl(dispy,tm0)
c		write(6,*)'dspy0(',i,1,')=',dspy0(ip,1)
        dspz0(ip,1)=dspz0(ip,1)+lapl(dispz,tm0)
	xtxx0(ip,1)=xtxx0(ip,1)+lapl(exx,tm0)
	xtxy0(ip,1)=xtxy0(ip,1)+lapl(exy,tm0)
        xtxz0(ip,1)=xtxz0(ip,1)+lapl(exz,tm0)
        xtyy0(ip,1)=xtyy0(ip,1)+lapl(eyy,tm0)
        xtyz0(ip,1)=xtyz0(ip,1)+lapl(eyz,tm0)
        xtzz0(ip,1)=xtzz0(ip,1)+lapl(ezz,tm0)
c       Update postseismic displacements for 10 time intervals
	do np=1,10
	tm1p=tm1*(0.129155d0*dexp(efac*dble(np-1)))/3.16881d0
	tm2p=tm2*(0.129155d0*dexp(efac*dble(np-1)))/3.16881d0
c		write(6,*)'tm1p,tm2p=',tm1p,tm2p
	dspx1(ip,np)=dspx1(ip,np)+lapl(dispx,tm2p)-lapl(dispx,tm1p) - corrf(np)*lapl(dispx,tm0)
	dspy1(ip,np)=dspy1(ip,np)+lapl(dispy,tm2p)-lapl(dispy,tm1p) - corrf(np)*lapl(dispx,tm0)
	dspz1(ip,np)=dspz1(ip,np)+lapl(dispz,tm2p)-lapl(dispz,tm1p) - corrf(np)*lapl(dispz,tm0)
	xtxx1(ip,np)=xtxx1(ip,np)+lapl(exx,tm2p)-lapl(exx,tm1p) - corrf(np)*lapl(exx,tm0)
	xtxy1(ip,np)=xtxy1(ip,np)+lapl(exy,tm2p)-lapl(exy,tm1p) - corrf(np)*lapl(exy,tm0)
	xtxz1(ip,np)=xtxz1(ip,np)+lapl(exz,tm2p)-lapl(exz,tm1p) - corrf(np)*lapl(exz,tm0)
	xtyy1(ip,np)=xtyy1(ip,np)+lapl(eyy,tm2p)-lapl(eyy,tm1p) - corrf(np)*lapl(eyy,tm0)
	xtyz1(ip,np)=xtyz1(ip,np)+lapl(eyz,tm2p)-lapl(eyz,tm1p) - corrf(np)*lapl(eyz,tm0)
	xtzz1(ip,np)=xtzz1(ip,np)+lapl(ezz,tm2p)-lapl(ezz,tm1p) - corrf(np)*lapl(ezz,tm0)
	enddo
60	continue
c*-*-*-*-*-*-*-*-*-*-*
45	continue
	close(2)

	open(2,file='qstat1.out')
	do i=1,ipts 
c	First write out static displacements and strains
	tm0=0.d0
	write(2,76) tm0,tm0,dspx0(i,1),dspy0(i,1),dspz0(i,1),xtxx0(i,1),xtxy0(i,1),xtxz0(i,1),
     &	xtyy0(i,1),xtyz0(i,1),xtzz0(i,1)
c       Do postseismic displacements for 10 time intervals
	do np=1,10
	tm1p=tm1*(0.129155d0*dexp(efac*dble(np-1)))/3.16881d0
	tm2p=tm2*(0.129155d0*dexp(efac*dble(np-1)))/3.16881d0
	write(2,76) tm1p*3.16881,tm2p*3.16881,dspx1(i,np),dspy1(i,np),dspz1(i,np),xtxx1(i,np),xtxy1(i,np),xtxz1(i,np),
     &	xtyy1(i,np),xtyz1(i,np),xtzz1(i,np)
	enddo
	enddo
	close(2)

76      format(2f10.3,9e13.5e2)
c
	write(6,*)'end of qstat1A'
	end  
	 
