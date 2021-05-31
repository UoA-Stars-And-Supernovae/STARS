# Warrick's attempt at a new Makefile.
# Drawn from
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

# FC = f77
FC = gfortran

#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15 -fbounds-check
#FFLAGS = -ffixed-li	ne-length-none -finit-local-zero -fno-automatic -O15 -m32
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O2
#FFLAGS = -extend_source -fast # -fpe3
#FFLAGS = -extend_source -g
FFLAGS = -g -ffixed-line-length-none -fno-automatic

#ifort options
#FFLAGS = -e -fast
#FFLAGS = -extend_source -ftz -static -C -fast

ODIR=obj
SDIR=src

_OBJ = main.o compos.o difrns.o divide.o elimn8.o equns1.o equns2.o \
funcs1.o funcs2.o nucrat.o nucrat2.o pressi.o printa.o printb.o printc.o \
remesh.o neutron.o xopac.o massloss.o diffusion.o diffusion2.o \
solver.o statef.o statel.o fdirac.o consts.o opacty.o opspln.o spline.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

VPATH = $(SDIR)

$(ODIR)/%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS)

bs: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS)

.PHONY: clean

clean:
	touch fort.zorkzork
	rm -f fort.*
	rm -f $(ODIR)/*.o *~ bs
