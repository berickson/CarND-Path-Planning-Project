# CarND-Path-Planning-Project - Brian Erickson
Self-Driving Car Engineer Nanodegree Program

### Abstract
This is my implementation of the Udacity Self Driving Car Nanodegree, Term 3 Path Planning Project.

The goal of this project is to create a path planner that can navigate a car through safely and efficiently through traffic along a busy highway.

The project is based on a fork of the Udacity project at https://github.com/udacity/CarND-Path-Planning-Project.

There rubric can be found at https://review.udacity.com/#!/rubrics/1020/view.

### Overview
##### Path planning
I used basic metrics for determining what lane was best based on occupancy data from the sensor fusion.  Lane changes take into account the velocity and position of cars in the current and neighboring lanes.  Preference is given to getting in lanes with faster cars ahead.

##### Track Smoothing
One of the big issues I had to start was that the track was given in large piecewise segments.  Boundary conditions near track segements like this can cause exessive accelerations, jerk, and speed around the segement ends.  One solution in the walkthrough was to use a "look ahead" strategy.  I found this solution to be too jerky and instead decided to represent the track using smooth parameterized splines.  I implemented this in the SmoothTrack class.  This class handles providing locations in both Frenet and Cartesian coordinates.  To go from Cartesian to Frenet, I used Newton's method of solving and minimization.  See the "newton_solve" and "newton_minimize" functions in main.cpp.

Note: I do  occasionaly see "out of lane" conditions reported in the simulator when the car is in the outside lane and is clearly travelling smoothly within the lane markings.  I believe this may be caused by the simulator reporting the condition based on distance from straight line track segments instead of the better, smooth trajectory.

#### Simulator.
The code is run designed for use with the Udacity Term3 Simulator from https://github.com/udacity/self-driving-car-sim/releases.

### Goals

#### Compilation

##### The code is to compile under cmake using CMakeLists.txt
I didn't change the structure of CMakeLists.txt, and the code compiles cleanly with cmake.

#### Valid Trajectories

##### The car is able to drive at least 4.32 miles without incident
I typically see it going much greather than this distance without incident.

##### The car drives according to the speed limit.

This was one of the more difficult problems. The issue was I wanted to keep as close to the speed limit without going over.  When I set speeds of 49.5 mph, it would go over the speed limit around corners and during lane changes.

Simply using the Frenet s coordinate for speeds wasn't working since speeds would be too great when going around corners or when doing lane changes.  The supplied S coordinates are for lane center, not for the different lanes.  This made the speed higher around left turns.  With my SmoothTrack class, I solved this by computing distances in Cartesian coordinates and mapping back to Frenet.

For lane changes, I had to ensure that lateral speed was accounted for, not just along S.  I solved this with simple use of Pythagorean theorem to find the maximum s velocity given a d velocity.

##### Max Acceleration and Jerk are not Exceeded.

I solved this goal by using smooth routes, and smooth lane changes.  The smooth lane changes use a jerk minimizing function as suggested in the course.  See the "jerk_minimizing_trajectory" function in main.cp.


##### Car does not have collisions.

##### The car stays in its lane, except for the time between changing lanes.

##### The car is able to change lanes
