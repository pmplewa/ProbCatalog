CC = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
LIBS = -ldnest4 -lpthread

OPTIONS = -DNoiseModel_UniformNoise -DPSFModel_GaussianPSF

default:
	$(CC) -I$(DNEST4_PATH) $(CXXFLAGS) $(OPTIONS) -c *.cpp
	$(CC) -pthread -L$(DNEST4_PATH)/DNest4/code -o main *.o $(LIBS)
	rm *.o

clean:
	rm -f main levels.txt sample.txt sample_info.txt posterior_sample.txt weights.txt
