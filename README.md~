# Extended Kalman Filter
Self-Driving Car Engineer Nanodegree Program

This project implements an Extended Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. It will use RMES to measure the results.

This project involves the Term 2 Simulator to which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository uses [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) to communicate with the simulator.

## Installation

Once the install for uWebSocketIO is complete, the main program can be built and run by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ExtendedKF

## INPUT / OUTPUT

INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

---

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

## SRC structure

`main.cpp` Initialize the objects necessary for the calculation, and creates the server that reads the input (the lidar and radar measurements), process the data and sends the estimations to the simulator.

`FusionEKF` Class that contains the EKF object, and acts as the interface to process the measurement.

`KalmanFilter` Implementation of the Kalman Filter. It has an update function which acts on the vector data of the measurement, and a predict function.

`Tools` Contains thee functions to calculate the Jacobian and the RMSE.

## Results

![Image](./Images/Image1.png?raw=true)
