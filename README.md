# Online Multi-Robot Motion Planning 
This code accompanies the RA-L/ICRA 2020 paper

>  C. E. Luis, M. Vukosavljev, and A. P. Schoellig, “Online trajectory generation with distributed model predictive control for multi-robot motion
planning,” IEEE Robot. Autom. Lett., vol. 5, no. 2, pp. 604–611, Jan. 2020.

## Citation
If you use this library for your own work, consider citing:
```
@article{luis2020online,
  title={Online trajectory generation with distributed model predictive control for multi-robot motion planning},
  author={Luis, Carlos E and Vukosavljev, Marijan and Schoellig, Angela P},
  journal={IEEE Robotics and Automation Letters},
  volume={5},
  number={2},
  pages={604--611},
  year={2020},
  publisher={IEEE}
}
```

## What's included
- Standalone C++ library implementing the algorithm.
- MATLAB code for running the benchmark and visualize data.

## Usage
Below you will find instructions on how to use the two main pieces of software included in this repo.

## C++ Library

Dependencies:
- C++14
- Cmake >= 3.0
- Eigen >= 3.0

### Installation

1. Initialize qpOASES submodule
```
cd <path-to-repo>
git submodule init && git submodule update
```

2. Build the project
```
mkdir build
cd build
cmake ..
make
```
3. Test installation by running example scenario 
```
cd ../bin
./run
```

4. You should see a stream of data in the console. The expected last lines of the console output are:
```
No collisions found!
All the vehicles reached their goals!
Writing solution to text file...
```

5. The generated simulation data is in `cpp/results/trajectories.txt`. You can run the MATLAB script `plot_results.m` for a 3D visualization of the generated trajectories.

### Running your own scenarios
The entry point of the code is `src/main.cpp`. If you want to run your own transition scenarios, or play around with the (many) hyperparameters of the algorithm, the main configuration file is in `cpp/config/config.json`. You can find an explanation of each hyperparameter in `cpp/config/help.txt`. 

## MATLAB
The code was written and executed with **MATLAB2018a**. There's no guarantees it will work in other versions.

### Running benchmark
To run the benchmark presented in the paper against the Buffered Voronoi Cells method, execute `matlab/tests/comp_dmpc_bvc.m`. At the top of the file there's several parameters to change the test characteristics.

### Plotting
- `plot_2agent_video.m`: dynamic 3D visualization of experiment with 2 drones exchanging positions
- `plot_comp_allCA.m`: summarize results from benchmark against BVC method
- `plot_disturbance_exp.m`: plots in paper for continuous vs event-based replanning
- `plot_hoop_test_paper.m`: static 3D visualization of 10 drones passing through a hula-hoop
- `plot_hoop_test_video.m`: dynamic 3D visualization of 10 drones passing through a hula-hoop
- `plot_obstaclefree_video.m`: dynamic 3D visualization of 20 drones randomly transitioning between goal points