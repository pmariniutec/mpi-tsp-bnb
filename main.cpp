#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <mpi.h>

int n = 8, lrank, numProcess;

struct Location {
  int visited;
  float x, y;
};
std::vector<Location> locations;

int main(int argc, char** argv) {
  float totalTime;

  std::ifstream fp{ "in.txt" };
  float x, y;
  while (fp >> x >> y) {
	locations.emplace_back(0, x, y);
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &lrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcess);

  double start = MPI_Wtime();
  tsp();
  double end = MPI_Wtime();
  totalTime = end - start;

  MPI_Finalize();
  return 0;
}
