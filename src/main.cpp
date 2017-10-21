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

using namespace std;

// for convenience
using json = nlohmann::json;

const double mph_to_m_s = 0.44704;


// structs, because they're easier to debug than vectors
struct Point {
  Point(double x=NAN, double y=NAN):x(x),y(y){}
  double x = NAN;
  double y = NAN;
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
    s_to_x.set_points(s,x);
    s_to_y.set_points(s,y);
    s_to_dx.set_points(s,dx);
    s_to_dy.set_points(s,dy);
    max_s = _max_s;
  }

  Point get_point(double s, double d=0) {
    Point p;
    s = fmod(s,max_s);
    p.x = s_to_x(s) + d * s_to_dx(s);
    p.y = s_to_y(s) + d * s_to_dy(s);
    return p;
  }
};

SmoothTrack smooth_track;
double last_s_sent = NAN;


struct CarState {
  double x,y,s,d,speed,acceleration_s, acceleration_d,sequence;
};

vector<CarState> path;

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

  double ref_vel = 0.0;
  int lane = 1;
  h.onMessage([&lane, &ref_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

            bool too_close = false;
            double car_ahead_speed;
            for(int i = 0; i < sensor_fusion.size(); i++) {
              // car is in my lane
              float d = sensor_fusion[i][6]; // why 6?
              if((d < 2. + 4. * lane + 2) && d > (2. + 4. * lane-2)) {
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double check_speed = sqrt(vx*vx+vy*vy);
                double check_car_s = sensor_fusion[i][5];
                // predict where car will be
                //check_car_s += ((double)prev_size*0.2*check_speed);
                double d_ahead = check_car_s - car_s;
                //cout << "car ahead " << d_ahead << " meters." << endl;
                if( d_ahead > 0 && d_ahead <30) {
                  too_close = true;
                  car_ahead_speed = check_speed;
                }
              }
            }

            double speed_limit = 190.0;

            cout << "car_s: " << car_s << endl;

            if(too_close) {
              if(car_ahead_speed < car_speed_mph * mph_to_m_s) {
                cout << "too close, slowing down" << endl;
                ref_vel -= 0.45;
              }
            } else {
              ref_vel += 0.45;
            }

            if(ref_vel > speed_limit) {
              ref_vel = speed_limit;
            }

            /*
            if(too_close ) {
              if(lane>0) {
                lane--;
              } else {
                lane++;
              }
            }
            */

            // remove path points already visited
            while(path.size() > previous_path_x.size()) {
              path.erase(path.begin());
            }

            // set path to car position if path is empty
            if(path.size() == 0) {
              CarState car;
              car.s = car_s;
              car.d = 2+4.*lane;//car_d;
              Point p = smooth_track.get_point(car.s, car.d);
              car.x = p.x;
              car.y = p.y;
              car.speed = car_speed_mph;
              car.acceleration_s = 0;
              path.push_back(car);
            }

            vector<double> next_x_vals, next_y_vals;

            double dt = 0.02;
            const int number_of_points_to_send = 20;
            CarState car_state = path[path.size()-1];
            while(path.size() < number_of_points_to_send) {
              car_state.s += 0.02 * ref_vel * mph_to_m_s;
              car_state.d = 2+4.*lane;
              Point p = smooth_track.get_point(car_state.s, car_state.d);
              car_state.x = p.x;
              car_state.y = p.y;
              path.push_back(car_state);
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
