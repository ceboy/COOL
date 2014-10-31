SVPATH = .
include $(SVPATH)/Mdefs.mk

EXEC = sv.exe

OBJ1 = m_zeromachine.f90 m_physics.f90 m_cell.f90 m_data.f90 svini.f90 sv.f90
OBJ2 = conv.py plot.py Mdefs.mk Makefile

all: $(EXEC)

# sv.exe: svini.f90 sv.f90
sv.exe: $(OBJ1)
	$(FC) -finit-local-zero $(OBJ1) -o $@ 
# sv.exe: sv.f90
# 	$(FC) sv.f90 -o $@ 

clean:
	$(RM) *.res *.png 

#dt.res h.res Q.res u.res P.res phi.res

cleanall: clean
	$(RM) *.mod $(EXEC) *~ 

cleanallfiles: cleanall
	$(RM) $(OBJ1) $(OBJ2)

