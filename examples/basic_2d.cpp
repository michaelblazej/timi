#include <timi/timi.hpp>
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "=== timi: 2D Interpolation Example ===\n\n";

    // Define a 3x3 grid for temperature data
    // Axes: x = [0, 5, 10] meters, y = [0, 5, 10] meters
    std::vector<std::vector<double>> axes = {
        {0.0, 5.0, 10.0},  // x-axis
        {0.0, 5.0, 10.0}   // y-axis
    };

    // Temperature values at each grid point (row-major order)
    // Row-major: x varies fastest
    //           x=0    x=5   x=10
    // y=0  :    20.0,  22.0, 25.0
    // y=5  :    21.0,  24.0, 28.0
    // y=10 :    23.0,  27.0, 32.0
    std::vector<double> temperatures = {
        20.0, 22.0, 25.0,  // y=0
        21.0, 24.0, 28.0,  // y=5
        23.0, 27.0, 32.0   // y=10
    };

    // Create a stateful interpolator
    timi::InterpolatorND<double, double> temp_interp(axes, temperatures);

    std::cout << "Temperature field (interpolated at 0.5m intervals):\n\n";
    std::cout << std::fixed << std::setprecision(2);

    // Print header
    std::cout << "      ";
    for (double x = 0.0; x <= 10.0; x += 2.5) {
        std::cout << std::setw(6) << x;
    }
    std::cout << "  (x)\n";

    // Print interpolated grid
    for (double y = 10.0; y >= 0.0; y -= 2.5) {
        std::cout << "y=" << std::setw(4) << y << " ";
        for (double x = 0.0; x <= 10.0; x += 2.5) {
            double t = temp_interp({x, y});
            std::cout << std::setw(6) << t;
        }
        std::cout << "\n";
    }

    std::cout << "\n--- 1D Interpolation Example ---\n\n";

    // 1D example: pressure vs altitude
    std::vector<double> altitude = {0.0, 1000.0, 5000.0, 10000.0};  // meters
    std::vector<double> pressure = {101.3, 89.9, 54.0, 26.5};       // kPa

    timi::Interpolator1D<double, double> pressure_interp(altitude, pressure);

    std::cout << "Atmospheric pressure at various altitudes:\n";
    for (double h = 0.0; h <= 10000.0; h += 2000.0) {
        std::cout << "  " << std::setw(6) << h << " m: "
                  << std::setw(6) << pressure_interp(h) << " kPa\n";
    }

    std::cout << "\n--- Custom Type Example ---\n\n";

    // Define a simple 2D vector type
    struct Point {
        double x, y;
    };

    // Required operators for interpolation
    auto add_points = [](const Point& a, const Point& b) -> Point {
        return {a.x + b.x, a.y + b.y};
    };

    // We can't use lambdas directly with the library, so let's use functional API
    // with a struct that has the operators defined

    std::cout << "See test.cpp for custom type (Vec2) interpolation example.\n";

    return 0;
}
