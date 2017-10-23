# CarND-Path-Planning-Project - Brian Erickson
Self-Driving Car Engineer Nanodegree Program



<a href="http://www.youtube.com/watch?feature=player_embedded&v=Sl-rqcI8tUs
" target="_blank"><img src="http://img.youtube.com/vi/Sl-rqcI8tUs/0.jpg" 
alt="Project running in simulator" width="240" height="180" border="10" /></a>

[![project in simulator](https://img.youtube.com/vi/Sl-rqcI8tUs/0.jpg)](https://www.youtube.com/watch?v=Sl-rqcI8tUs)

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

I solved this goal by using smooth routes, and smooth lane changes.  The smooth lane changes use a jerk minimizing function as suggested in the course.  See the "jerk_minimizing_trajectory" function in main.cpp.

##### Car does not have collisions.

The car avoid collisions by only changing lanes when it is safe and by tracking speed to the car ahead.

Lane change safety is maintained by calculating an expected clearance gap at the begin and end of the lane change.  Both the current lane and the candidate lane are considered. Speeds of cars in front of the car and behind the car are assumed to keep constant velocity during the lane change.

Tracking speed of the car ahead is fairly simple and doesn't really need a PID.  Speed is reduced if the car is too close, and if the car is a good distance, the algorithm simply sets a target speed to matcht the speed of the car ahead.

##### The car stays in its lane, except for the time between changing lanes.

Smooth tracks are computed for the center of lanes using Frenet coordinates and the bezier based SmoothTrack class.  Lane centers are kept perfectly except when a lane change maneuver is executed.  Lane change manuevers always start and end in lane centers.

##### The car is able to change lanes

The car uses a cost function to know when to change lanes.  It considers the speed and distance of the car ahead for each lane. Faster speeds and greater distances are rewarded with lower costs.  Additionally, the center lane is preferred and given a lower cost.  The tanh function is used to keep scores in a manageble range.


### Reflection

Generating paths has two distict modes.  One is for lane keeping, and one is for lane following.

For lane keeping, a small number of points are generated.  Goal velocities are determined based on speed limit, speed and distance of the car ahead.  The car accelerates or decelerates to reach goal velocity.  I chose a rather aggressive 5 m/s^3 for the acceleration.  This keeps a good amount of acceleration available for the lateral acceleration needed to keep lanes during curves.

For lane changing, a complete lane change set of points is created and committed to.  That is why it is so important to analyze the safety of the lane change in advance. It this model once it is started, it is always completed.

Here are some places that could be improved for a more reliable system:
- Consider more cases for safety, for example, when another car enters your lane space.
- Look more steps ahead in path planning
- Have the ability to abort lane change manouvers.
