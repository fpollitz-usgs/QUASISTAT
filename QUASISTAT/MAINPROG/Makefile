# Makefile for QSTAT0, QSTAT1
#

sdir = ../SOURCE

#FC=ifort
#FFLAGS=-O2 -ftz -ip -ipo -axP -align all -extend_source

FC=gfortran
FFLAGS=-O2 -ffixed-form -ffixed-line-length-132 

all : qstat0 qstat1

qstat0  : ${sdir}/qstat0.f ${sdir}/msubs.f ${sdir}/bsstep.f ${sdir}/rzextr.f ${sdir}/mmid.f ${sdir}/interp.f
	${FC} -o $@ ${FFLAGS} ${sdir}/qstat0.f ${sdir}/msubs.f ${sdir}/bsstep.f ${sdir}/rzextr.f ${sdir}/mmid.f ${sdir}/interp.f 

qstat1  : ${sdir}/qstat1.f ${sdir}/lgndrl.f ${sdir}/lapl.f ${sdir}/svdcmp8.f 
	${FC} -o $@ ${FFLAGS} ${sdir}/qstat1.f ${sdir}/lgndrl.f ${sdir}/lapl.f ${sdir}/svdcmp8.f

