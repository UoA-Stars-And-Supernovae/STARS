# Drawn from
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

FC = gfortran

FFLAGS = -g -ffixed-line-length-none -fno-automatic -march=native -fimplicit-none


ODIR=obj
SDIR=src

_OBJ = main.o compos.o difrns.o divide.o elimn8.o equns1.o equns2.o \
funcs1.o funcs2.o nucrat.o nucrat2.o pressi.o printa.o printb.o printc.o \
remesh.o neutron.o xopac.o massloss.o diffusion.o diffusion2.o \
solver.o statef.o statel.o fdirac.o consts.o opacty.o opspln.o spline.o \
overflow.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

VPATH = $(SDIR)

$(ODIR)/%.o: %.f
	@mkdir -p $(ODIR)
	$(FC) -c -o $@ $< $(FFLAGS)

bs: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS)

clean:
	touch fort.tmp
	rm -f fort.*
	rm -f $(ODIR)/*.o *~ bs
