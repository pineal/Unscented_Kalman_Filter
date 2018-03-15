# Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` Previous versions use i/o from text files.  The current state uses i/o
from the simulator.

## Other Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Accuracy
The result meet the requirement that RMSE <= [.09, .10, .40, .30] for Dataset 1 after it becomes stable.

## Follows the Correct Algorithm
Sensor Fusion algorithm follows the general processing flow as taught in the preceding lessons. The Entry is in function processMeasurement() from line 82 to line 123 in ukf.cpp, it processed lidar and radar measurement package, it also token care of first measurement initilization, and the UKF did the predict and update steps as correct algorithm.

## Code Efficiency
No significant time complexity improvement has been spotted following UKF algorithm in this project. 


