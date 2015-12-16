All:
#	g++-5 -o bin/tunnelShock tunnelShock.C -std=c++11 -O3
#	g++ -o bin/tunnelShock tunnelShock.C -std=c++11 -O3
	mpiCC tunnelShock.C -o bin/tunnelShock -std=c++11 -O3
