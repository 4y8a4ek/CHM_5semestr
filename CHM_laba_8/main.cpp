#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;
const double a = 0.7;       
const double gamma = 0.7;
const double x_eq = 1624; 
const double y_eq = 540;
const double beta = a / y_eq; 
const double delta = gamma / x_eq;
double dx(double x, double y)
{
    return (a - beta * y) * x;
}
double dy(double x, double y)
{
    return (-gamma + delta * x) * y;
}
void runge_kutta(double t0, double xzero, double yzero, double h, int steps)
{
    ofstream x_file("x_values_runge_kutta.txt");
    ofstream y_file("y_values_runge_kutta.txt");
    double t = t0, x = xzero, y = yzero;
    for (int i = 0; i < steps; ++i)
    {
        double k1x = h * dx(x, y);
        double k1y = h * dy(x, y);
        double k2x = h * dx(x + k1x / 2, y + k1y / 2);
        double k2y = h * dy(x + k1x / 2, y + k1y / 2);
        double k3x = h * dx(x + k2x / 2, y + k2y / 2);
        double k3y = h * dy(x + k2x / 2, y + k2y / 2);
        double k4x = h * dx(x + k3x, y + k3y);
        double k4y = h * dy(x + k3x, y + k3y);
        x += (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        t += h;
            x_file << fixed << setprecision(4) << x << "\n";
            y_file << fixed << setprecision(4) << y << "\n";
   }
    x_file.close();
    y_file.close();
}

void adams_bashforth(double t0, double xzero, double yzero, double h, int steps)
{
    ofstream x_file("x_values_adams_bashforth.txt");
    ofstream y_file("y_values_adams_bashforth.txt");
    vector<double> t(steps + 1), x(steps + 1), y(steps + 1);
    t[0] = t0;
    x[0] = xzero;
    y[0] = yzero;
    for (int i = 0; i < 4; ++i)
    {
        double k1x = h * dx(x[i], y[i]);
        double k1y = h * dy(x[i], y[i]);
        double k2x = h * dx(x[i] + k1x / 2, y[i] + k1y / 2);
        double k2y = h * dy(x[i] + k1x / 2, y[i] + k1y / 2);
        double k3x = h * dx(x[i] + k2x / 2, y[i] + k2y / 2);
        double k3y = h * dy(x[i] + k2x / 2, y[i] + k2y / 2);
        double k4x = h * dx(x[i] + k3x, y[i] + k3y);
        double k4y = h * dy(x[i] + k3x, y[i] + k3y);
        x[i + 1] = x[i] + (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y[i + 1] = y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        t[i + 1] = t[i] + h;
    }
    for (int i = 3; i < steps; ++i)
    {
        x[i + 1] = x[i] + h / 24 * (55 * dx(x[i], y[i]) - 59 * dx(x[i - 1], y[i - 1]) + 37 * dx(x[i - 2], y[i - 2]) - 9 * dx(x[i - 3], y[i - 3]));
        y[i + 1] = y[i] + h / 24 * (55 * dy(x[i], y[i]) - 59 * dy(x[i - 1], y[i - 1]) + 37 * dy(x[i - 2], y[i - 2]) - 9 * dy(x[i - 3], y[i - 3]));
        t[i + 1] = t[i] + h;
    }
    for (int i = 0; i <= steps; ++i)
    {
            x_file << fixed << setprecision(4) << x[i] << "\n";
            y_file << fixed << setprecision(4) << y[i] << "\n";
    }
    x_file.close();
    y_file.close();
}

int main()
{
    double t0 = 0.0;  
    double xzero = 2189; 
    double yzero = 667; 
    double h = 0.1;   
    int steps = 3650;
    cout << "Вычисленные параметры:\n";
    cout << "a = " << a << ", beta = " << beta << ", gamma = " << gamma << ", delta = " << delta << "\n";
    runge_kutta(t0, xzero, yzero, h, steps); 
    adams_bashforth(t0, xzero, yzero, h, steps);
    return 0;
}
