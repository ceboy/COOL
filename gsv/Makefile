SVPATH = .
include $(SVPATH)/Mdefs.mk

EXEC = sv.exe

OBJ1 = m_cell.f90 m_data.f90 svini.f90 sv.f90

all: $(EXEC)

# sv.exe: svini.f90 sv.f90
sv.exe: $(OBJ1)
	$(FC) $(OBJ1) -o $@ 
# sv.exe: sv.f90
# 	$(FC) sv.f90 -o $@ 

clean:
	$(RM) h.res Q.res u.res P.res

cleanall: clean $(EXEC)