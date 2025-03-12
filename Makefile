CPP = g++
FLAGS = -O3 -m64 -std=c++11 -Wall -Wextra -DNDEBUG -g
GRBPATH = /home/ling/opt/gurobi912/linux64


SG: SingleTarget.cpp
	$(CPP) $(FLAGS) SingleTarget.cpp -o SingleTarget -I$(GRBPATH)/include -L$(GRBPATH)/lib -lgurobi_c++ -lgurobi91 -lm

MT: MultiTarget.cpp
	$(CPP) $(FLAGS) MultiTarget.cpp -o MultiTarget -I$(GRBPATH)/include -L$(GRBPATH)/lib -lgurobi_c++ -lgurobi91 -lm


clean: 
	rm -f code *.o *~ls

