CXX=g++  -O2 

all:	$(SRC)
	make $(OBJ)
	make theta_naive
	make psi
	make theta

SRC=$(wildcard *.cc)
OBJ=$(patsubst %.cc, %.o, $(SRC))


%.o: %.cc
	$(CXX) -c $^

theta_naive: theta_naive.cpp
	$(CXX) -o theta_naive theta_naive.cpp $(OBJ) $(LFLAGS) -lgmpxx -lgmp -lmpfr

psi: tst_psi.cpp
	$(CXX) -o psi tst_psi.cpp $(OBJ) $(LFLAGS)  -lgmpxx -lgmp -lmpfr

theta: tst_theta.cpp
	$(CXX) -o theta tst_theta.cpp $(OBJ) $(LFLAGS) -lgmpxx -lgmp -lmpfr

clean:
	rm -f *.o psi theta_naive theta



