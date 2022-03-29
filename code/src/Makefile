FC = gfortran
OBJ = detcovid.o stochcovid.o mcmc.o tanhrt.o ignbin.o toms343.o utils.o
FFLAGS = -c -O3 
all: $(OBJ)

#R CMD SHLIB $(OBJ)

detcovid.o:  detcovid.f
	$(FC) $(FFLAGS) $(PKG_LIBS) detcovid.f

stochcovid.o:  stochcovid.f
	$(FC) $(FFLAGS) $(PKG_LIBS) stochcovid.f

mcmc.o:  mcmc.f
	$(FC) $(FFLAGS) $(PKG_LIBS) mcmc.f

tanhrt.o:  tanhrt.f
	$(FC) $(FFLAGS) $(PKG_LIBS) tanhrt.f

ignbin.o: ignbin.f
	$(FC) $(FFLAGS) $(PKG_LIBS) ignbin.f

toms343.o: toms343.f
	$(FC) $(FFLAGS) $(PKG_LIBS) toms343.f

utils.o: utils.f
	$(FC) $(FFLAGS) $(PKG_LIBS) utils.f
