#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>
#include <thread>
#define Pi M_PI
#define T 3.0

inline double Sin(double x){return sin(x);}
inline double Cos(double x){return cos(x);}
inline double Power(double x, double y){return pow(x,y);}

inline double fx(double t){return 4.5*Sin(Pi/2. + (2*Pi*t)/T);}
inline double fx_d(double t) {return (-28.274333882308138*Sin((2*Pi*t)/T))/T;}
inline double fx_dd(double t){return (-177.65287921960845*Cos((2*Pi*t)/T))/Power(T,2);}

inline double fy(double t){return 1.5*Sin((6*Pi*t)/T);}
inline double fy_d(double t){return (28.274333882308138*Cos((6*Pi*t)/T))/T;}
inline double fy_dd(double t){return (-532.9586376588253*Sin((6*Pi*t)/T))/Power(T,2);}

double vFunction(double t){
    double xd2 = pow(fx_d(t),2.0);
    double yd2 = pow(fy_d(t),2.0);
    double result = sqrt(xd2 + yd2);
    return result;
}

double wFunction(double t){
    double up = (fx_d(t) * fy_dd(t) - fx_dd(t) * fy_d(t));
    double down = pow(vFunction(t),2.0);
    return up / down;
}

int main() {
    double t = 0;
    using namespace std::this_thread;
    using namespace std::chrono_literals;
    using std::chrono::system_clock;

    std::cout << "Lineal\t\tAngular" << std::endl;
    while(t < 6) {
        std::cout << vFunction(t) << "\t\t" << wFunction(t) << std::endl;
        t += 0.1;
        sleep_for(100ms);
    }
    return 0;
}