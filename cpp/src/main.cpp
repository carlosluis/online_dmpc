#include <iostream>
#include <Eigen/Dense>
#include "simulator.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;



int main() {
	cout << "Hello world!" << endl;

    std::ifstream my_config_file("../config/config.json");
    assert(my_config_file && "Couldn't find the config file");

    Simulator sim(my_config_file);

    int T = 20; // simulation duration
    sim.run(T);

    // Save data to file
    char const *file = "/home/carlos/repos/bezier_dmpc/cpp/results/trajectories.txt";
    sim.saveDataToFile(file);

	return 0;
}


