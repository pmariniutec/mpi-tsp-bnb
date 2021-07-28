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

void swap(int* p1, int* p2) {
  int temp;
  temp = *p1;
  *p1 = *p2;
  *p2 = temp;
}

float distance(Location a, Location b) {
  return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

float calculateCost(int* a) {
  int i;
  float cost = 0.0f;
  for (i = 0; i < n - 1; i++)
	cost += distance(locations[a[i] - 1], locations[a[i + 1] - 1]);
  return cost;
}

int getNearest(int curr, int start, int end) {
  int i, index = -1;
  float min = std::numeric_limits<float>::max();
  float dist;
  for (i = start; i <= end; ++i) {
	dist = distance(locations[curr], locations[i]);
	if (dist < min && i != curr && locations[i].visited == 0) {
	  min = dist;
	  index = i;
	}
  }
  return index;
}

void tsp() {
  int index, sloc, eloc, next;
  int locpn = n / numProcess;
  int* inm;
  int fpath[n];
  double start, end, dt;
  float min = std::numeric_limits<float>::max();
  float dist;
  float cost = 0.0f;

  inm = new int[numProcess];

  sloc = locpn * lrank;
  eloc = sloc + locpn - 1;

  if (lrank == numProcess - 1)
	eloc += n % numProcess;

  next = 0;
  fpath[0] = 0;

  for (int i = 0; i < n - 1; i++) {
	MPI_Bcast(&next, 1, MPI_INT, 0, MPI_COMM_WORLD);
	locations[next].visited = 1;
	auto index = getNearest(next, sloc, eloc);
	MPI_Gather(&index, 1, MPI_INT, inm, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (lrank == 0) {
	  index = inm[0];
	  min = std::numeric_limits<float>::max();
	  for (int j = 0; j < numProcess; ++j) {
		if (inm[j] < 0)
		  continue;
		dist = distance(locations[next], locations[inm[j]]);
		if (dist < min) {
		  min = dist;
		  index = inm[j];
		}
	  }
	  next = index;
	  fpath[i + 1] = index;
	}

	MPI_Barrier(MPI_COMM_WORLD);
  }

  delete inm;
}

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

  auto start = MPI_Wtime();
  tsp();
  auto end = MPI_Wtime();
  totalTime = end - start;

  MPI_Finalize();
  return 0;
}
