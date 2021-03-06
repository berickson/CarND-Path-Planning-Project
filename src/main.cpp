#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "math.h"

#include "Eigen-3.3/Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// for convenience
using json = nlohmann::json;

const double mph_to_m_s = 0.44704;

double clamp(double a, double low, double high) {
  if (a < low) return low;
  if (a > high) return high;
  return a;
}

double netwon_solve(std::function<double(double)> f, double guess) {
    const double dx = 1E-6;
    const double epsilon = 1E-9;
    int max_iters = 20;
    for(int i = 0; i < max_iters; i++) {
        double x = f(guess);
        if(x<epsilon) break;
        double x_left = f(guess-dx);
        double x_right = f(guess+dx);
        double x_dot = (x_right-x_left)/(2*dx);
        guess -= x/x_dot;
    }
    return guess;
}

double newton_minimize(std::function<double(double)> f, double guess) {
    const double dx = 1E-6;
    const double epsilon = 1E-9;
    int max_iters = 20;
    for(int i = 0; i < max_iters; i++) {
        double x = f(guess);
        double x_left = f(guess-dx);
        double x_right = f(guess+dx);
        double x_dot = (x_right-x_left)/(2*dx);
        if(fabs(x_dot)<epsilon) break;
        double x_double_dot = ((x_right-x)/dx - (x-x_left)/ dx)/dx;
        guess -= x_dot/x_double_dot;
    }
    return guess;
}

// structs, because they're easier to debug than vectors
struct Point {
  Point(double x=NAN, double y=NAN):x(x),y(y){}
  double x = NAN;
  double y = NAN;
  double distance_to(const Point & p2) const {
    const Point & p1 = *this;
    double dx = p2.x-p1.x;
    double dy = p2.y-p1.y;
    return sqrt(dx*dx+dy*dy);
  }
};

struct Frenet {
  Frenet(double s, double d): s(s), d(d) {}
  double s = NAN;
  double d = NAN;
};


// For converting back and forth between radians and degrees.

constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
Frenet getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

  return Frenet(frenet_s,frenet_d);

}

// Transform from Frenet s,d coordinates to Cartesian x,y
Point getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

  return Point(x,y);

}


class SmoothTrack {
  tk::spline s_to_x;
  tk::spline s_to_y;
  tk::spline s_to_dx;
  tk::spline s_to_dy;
  double max_s;
public:
  void init(
      vector<double> s,
      vector<double> x,
      vector<double> y,
      vector<double> dx,
      vector<double> dy,
      double _max_s)
  {
    vector<double> s_,x_,y_,dx_,dy_;

    // add elements to make a smooth loop.
    s_=s;
    x_=x;
    y_=y;
    dx_=dx;
    dy_=dy;
    s_.push_back(_max_s + s[0]);
    s_.push_back(_max_s + s[1]);
    x_.push_back(x[0]);
    x_.push_back(x[1]);
    y_.push_back(y[0]);
    y_.push_back(y[1]);
    dx_.push_back(dx[0]);
    dx_.push_back(dx[1]);
    dy_.push_back(dy[0]);
    dy_.push_back(dy[1]);
    s_to_x.set_points(s_,x_);
    s_to_y.set_points(s_,y_);
    s_to_dx.set_points(s_,dx_);
    s_to_dy.set_points(s_,dy_);
    max_s = _max_s;
  }

  Point get_point(double s, double d=0) {
    Point p;
    s = fmod(s,max_s);
    p.x = s_to_x(s) + d * s_to_dx(s);
    p.y = s_to_y(s) + d * s_to_dy(s);
    return p;
  }

  double get_s_ahead(double s, double d, double desired_distance) {
    assert(desired_distance>=0);
    Point p1 = get_point(s,d);

    // returns point ds in front of (s,d)
    std::function<double(double)> error = [this,d,s,p1,desired_distance](double ds){
      Point p2 = get_point(s+fabs(ds), d);
      double delta = p1.distance_to(p2)-desired_distance;
      return delta*delta;
    };
    double ds = newton_minimize(error, desired_distance);
    return s + ds;
  }

  Frenet get_frenet(Point p, double approx_s) {
      std::function<double(double)> d2 = [this,p](double s){
          Point p2 = get_point(s);
          return p2.x * p.x + p2.y * p.y;
      };
      double s = newton_minimize(d2, approx_s);
      double d = sqrt(d2(s));
      return Frenet(s, d);
      // todo: find sign of d, currently we only need positive
  }
};

SmoothTrack smooth_track;
double last_s_sent = NAN;

enum PlanningMode {
  keep_lane,
  lane_change
};

struct CarState {
  double x,y,s,d,m_s,acceleration_s, acceleration_d,sequence;
  PlanningMode mode;

};

struct FusionCar {
  double x,y,vx,vy,s,d,speed_m_s;
  int id,lane;
   FusionCar(json::reference j) {
     id=j[0];
     x=j[1];
     y=j[2];
     vx=j[3];
     vy=j[4];
     s=j[5];
     d=j[6];
     speed_m_s = sqrt(vx*vx+vy*vy);
     lane=int(d)/4;
   }
};

vector<CarState> path;

struct LaneStatus{
  int lane_number = 0;
  double closest_d_ahead = 999;
  double v_ahead = NAN;
  double v_behind = NAN;
  double closest_d_behind = -999;
  bool has_car_ahead = false;
  bool has_car_behind = false;
};

class Polynomial {
public:
  vector<double> coefs;
  Polynomial(vector<double>coefs_) {
    coefs = coefs_;
  }

  double eval(double t) {
    double rv = 0;
    for(int i = 0; i < coefs.size(); ++i) {
      rv += coefs[i]*pow(t,i);
    }
    return rv;
  }

};

vector<double> jerk_minimizing_trajectory(vector< double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT
    an array of length 6, each value corresponding to a coefficent in the polynomial
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */

    MatrixXd A = MatrixXd(3, 3);
  A << T*T*T, T*T*T*T, T*T*T*T*T,
          3*T*T, 4*T*T*T,5*T*T*T*T,
          6*T, 12*T*T, 20*T*T*T;

  MatrixXd B = MatrixXd(3,1);
  B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
          end[1]-(start[1]+start[2]*T),
          end[2]-start[2];

  MatrixXd Ai = A.inverse();

  MatrixXd C = Ai*B;

  vector <double> result = {start[0], start[1], .5*start[2]};
  for(int i = 0; i < C.size(); i++)
  {
      result.push_back(C.data()[i]);
  }

    return result;
}

bool is_lane_change_safe(LaneStatus & lane_status, double seconds_required, double current_speed) {
  double min_front_gap = 10;
  double min_rear_gap = 20;

  // note: closest_d values are negative when behind self

  // not safe if the cars are currently too close
  if(lane_status.has_car_ahead) {
    if(lane_status.closest_d_ahead < min_front_gap) {
      return false;
    }
  }

  if(lane_status.has_car_behind) {
    if(lane_status.closest_d_behind > -min_rear_gap) {
      return false;
    }
  }

  // not safe if the cars will be too close, assuming constant velocities
  if(lane_status.has_car_ahead) {
    double end_d_ahead = lane_status.closest_d_ahead + seconds_required * (lane_status.v_ahead - current_speed);
    if(end_d_ahead < min_front_gap) {
      return false;
    }
  }
  if(lane_status.has_car_behind) {
    double end_d_behind = lane_status.closest_d_behind + seconds_required * (lane_status.v_behind - current_speed);
    if(end_d_behind > -min_rear_gap) {
      return false;
    }
  }

  return true;
}



int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  smooth_track.init(map_waypoints_s, map_waypoints_x, map_waypoints_y, map_waypoints_dx, map_waypoints_dy, max_s);

  double ref_vel_m_s = 0.0;
  int lane = 1;
  h.onMessage([&lane, &ref_vel_m_s, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
            double car_speed_mph = j[1]["speed"];
            double car_speed_m_s = car_speed_mph * mph_to_m_s;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];


            LaneStatus lane_status[3];

            for(int i = 0; i < sensor_fusion.size(); i++) {
              FusionCar check_car(sensor_fusion[i]);
              LaneStatus & check_lane = lane_status [check_car.lane];
              double d_ahead = check_car.s - car_s;
              if(d_ahead > 0 && d_ahead < check_lane.closest_d_ahead) {
                check_lane.closest_d_ahead = d_ahead;
                check_lane.v_ahead = check_car.speed_m_s;
                check_lane.has_car_ahead = true;
              }
              if(d_ahead < 0 && d_ahead > check_lane.closest_d_behind) {
                check_lane.closest_d_behind = d_ahead;
                check_lane.v_behind = check_car.speed_m_s;
                check_lane.has_car_behind = true;
              }
            }

            for(int i=0;i<3;i++) {
              LaneStatus & s = lane_status[i];
              s.lane_number = 0;
              cout << "L" << i << " (" << std::fixed << std::setprecision(1) << std::setw(5)<< s.closest_d_ahead << ", " << s.closest_d_behind << ") ";
            }
            cout << endl;

            double speed_limit = 49.5;
            double speed_limit_m_s = speed_limit * mph_to_m_s;
            cout << "car_s: " << car_s << endl;

            double min_front_gap = 20;
            double min_back_gap = 10;
            double min_follow_gap = 15;
            double max_follow_gap = 18;

            LaneStatus & my_lane = lane_status[lane];

            double target_speed_m_s = speed_limit_m_s;

            bool too_close = false;
            // try to follow car ahead
            if(my_lane.has_car_ahead) {
              too_close = my_lane.closest_d_ahead < min_follow_gap;
              bool too_far = my_lane.closest_d_ahead > max_follow_gap;
              if(too_close) {
                target_speed_m_s = my_lane.v_ahead - 2;
              } else if (too_far ){
                target_speed_m_s = speed_limit_m_s;
              } else {
                target_speed_m_s = my_lane.v_ahead;
              }

            }

            target_speed_m_s = clamp(target_speed_m_s, 0, speed_limit_m_s);
            // set path to car position if path is empty, make first car state
            if(path.size() == 0) {
              CarState car;
              car.s = car_s;
              car.d = 2+4.*lane;//car_d;
              Point p = smooth_track.get_point(car.s, car.d);
              car.x = p.x;
              car.y = p.y;
              car.m_s = ref_vel_m_s;
              car.acceleration_s = 0;
              car.mode = keep_lane;
              path.push_back(car);
            }
            CarState car_state = path[path.size()-1];
            const double min_lane_change_m_s = 5;     // speed to be going before considering a lane change
            double seconds_to_execute_lane_change = 3;

            LaneStatus & r = lane_status[lane+1];
            LaneStatus & l = lane_status[lane-1];

            // see which lane is best
            bool right_lane_change_safe = false;
            if(lane < 2) {
              right_lane_change_safe =
                  is_lane_change_safe(r, seconds_to_execute_lane_change, car_state.m_s)
                  && is_lane_change_safe(my_lane, seconds_to_execute_lane_change, car_state.m_s);
            }


            bool left_lane_change_safe = false;
            if(lane > 0) {
              left_lane_change_safe =
                  is_lane_change_safe(l, seconds_to_execute_lane_change, car_state.m_s)
                  && is_lane_change_safe(my_lane, seconds_to_execute_lane_change, car_state.m_s);
            }

            // get scores for each line
            vector<double> lane_costs={0,0,0};

            for(int i= 0; i< 3; i++) {
              LaneStatus l = lane_status[i];
              if(l.has_car_ahead) {
                // prefer cars to be far away
                lane_costs[i] += tanh((50.-l.closest_d_ahead)/100.);
                // prefer cars in front to be going fast
                lane_costs[i] += tanh((50.-l.v_ahead)/30.);
              } else {
                lane_costs[i] -= 2;
              }
              // keep options open by preferring middle lane
              if(i != 1) {
                lane_costs[i] += 0.2;
              }
            }

            // safety overrides above
            double keep_cost, left_lane_change_cost, right_lane_change_cost;

            if(right_lane_change_safe) {
              right_lane_change_cost = lane_costs[lane+1];
            } else {
              right_lane_change_cost = 999;
            }
            if(left_lane_change_safe) {
              left_lane_change_cost = lane_costs[lane-1];
            } else {
              left_lane_change_cost = 999;
            }
            keep_cost = lane_costs[lane];


            double lane_delta = 0;
            if(left_lane_change_cost < keep_cost && left_lane_change_cost < right_lane_change_cost) {
              lane_delta = -1;
            } else if (right_lane_change_cost < keep_cost && right_lane_change_cost < left_lane_change_cost) {
              lane_delta = 1;
            }

            // remove path points already visited
            while(path.size() > previous_path_x.size()) {
              path.erase(path.begin());
            }

            vector<double> next_x_vals, next_y_vals;

            double dt = 0.02;
            const int number_of_points_to_send = 20;

            if(path.size() < number_of_points_to_send) {
              if(lane_delta == 0 || car_state.m_s < min_lane_change_m_s) {
                // generate trajector for keep lane
                while(path.size() < number_of_points_to_send) {
                  // try to reach target speed
                  if(car_state.m_s < target_speed_m_s) {
                    car_state.m_s += std::min(0.1, target_speed_m_s - car_state.m_s);
                  } else if (car_state.m_s > target_speed_m_s) {
                    car_state.m_s -= std::min(0.1, car_state.m_s - target_speed_m_s);
                  }
                  car_state.m_s = clamp(car_state.m_s, 0, speed_limit_m_s);

                  car_state.d = 2+4.*lane;
                  double dx = dt * car_state.m_s;
                  car_state.s = smooth_track.get_s_ahead(car_state.s, car_state.d, dx);
                  Point p = smooth_track.get_point(car_state.s, car_state.d);
                  car_state.x = p.x;
                  car_state.y = p.y;
                  car_state.mode = keep_lane;
                  path.push_back(car_state);
                }
              } else {
                // generate trajectory for changing lanes
                Polynomial trajectory(jerk_minimizing_trajectory({car_state.d,0,0},{car_state.d+4*lane_delta,0,0},seconds_to_execute_lane_change));
                for(double t = dt; t<=seconds_to_execute_lane_change; t+= dt) {
                  double d_normal = trajectory.eval(t)-car_state.d;
                  double dx = dt * car_state.m_s;
                  double d_tangent = sqrt(dx*dx-d_normal*d_normal);
                  car_state.s = smooth_track.get_s_ahead(car_state.s, car_state.d, d_tangent);
                  car_state.d = trajectory.eval(t);
                  Point p = smooth_track.get_point(car_state.s, car_state.d);
                  car_state.x = p.x;
                  car_state.y = p.y;
                  path.push_back(car_state);
                  car_state.mode = lane_change;
                }
                lane += lane_delta;
              }
            }

            for(CarState & car : path) {
              next_x_vals.push_back(car.x);
              next_y_vals.push_back(car.y);
            }

            json msgJson;
            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}


/*

todo:
- add scores per lane ( calc scores per state )
*/
