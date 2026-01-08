# TIMI

**T**emplated **I**nterpolation for **M**ulti-dimensional **I**nputs. A modern C++17 header-only library for n-dimensional linear interpolation.

## Overview

TIMI provides efficient multilinear interpolation on rectilinear grids. It supports both 1D interpolation and arbitrary n-dimensional grids, with a clean API that works with custom data types.

**Use cases:**
- Scientific computing and numerical analysis
- Graphics and color interpolation
- Lookup tables and data approximation
- Signal processing and sensor data

## Features

- **Header-only** - Single `#include`, no linking required
- **No dependencies** - Uses only the C++ Standard Library
- **1D and N-dimensional** - From simple curves to multidimensional fields
- **Dual API** - Functional (one-shot) and stateful (reusable) interfaces
- **Custom types** - Works with any type supporting `+` and `*` operators
- **Thread-safe** - Const sampling operations safe for concurrent use
- **Boundary handling** - Automatic clamping for out-of-range queries

## Requirements

- C++17 compatible compiler (GCC, Clang, MSVC)
- CMake 3.14+ (for building tests and examples)

## Installation

### Header-only (recommended)

Copy the `include/timi` directory to your project:

```bash
cp -r include/timi /path/to/your/project/include/
```

Then include in your code:

```cpp
#include <timi/timi.hpp>
```

### CMake Integration

Add as a subdirectory:

```cmake
add_subdirectory(timi)
target_link_libraries(your_target PRIVATE timi)
```

Or use FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(
    timi
    GIT_REPOSITORY https://github.com/your-username/timi.git
    GIT_TAG main
)
FetchContent_MakeAvailable(timi)
target_link_libraries(your_target PRIVATE timi)
```

## Quick Start

### 1D Interpolation

```cpp
#include <timi/timi.hpp>

// Sample data: pressure vs altitude
std::vector<double> altitude = {0.0, 1000.0, 5000.0, 10000.0};
std::vector<double> pressure = {101.3, 89.9, 54.0, 26.5};

// Create a reusable interpolator
timi::Interpolator1D<double, double> interp(altitude, pressure);

// Query at any point
double p = interp(2500.0);  // Interpolated pressure at 2500m
```

### 2D Grid Interpolation

```cpp
// Define a 3x3 temperature grid
std::vector<std::vector<double>> axes = {
    {0.0, 5.0, 10.0},  // x-axis (meters)
    {0.0, 5.0, 10.0}   // y-axis (meters)
};

// Values in row-major order (x varies fastest)
std::vector<double> temperatures = {
    20.0, 22.0, 25.0,  // y=0
    21.0, 24.0, 28.0,  // y=5
    23.0, 27.0, 32.0   // y=10
};

timi::InterpolatorND<double, double> temp_field(axes, temperatures);

// Query at any point
double t = temp_field({2.5, 7.5});  // Temperature at (2.5, 7.5)
```

### Custom Types

Any type with `operator+` and `operator*` works:

```cpp
struct Color {
    double r, g, b;
};

Color operator+(const Color& a, const Color& b) {
    return {a.r + b.r, a.g + b.g, a.b + b.b};
}

Color operator*(const Color& c, double s) {
    return {c.r * s, c.g * s, c.b * s};
}

// Create a color gradient
std::vector<double> positions = {0.0, 0.5, 1.0};
std::vector<Color> colors = {
    {1.0, 0.0, 0.0},  // Red
    {0.0, 1.0, 0.0},  // Green
    {0.0, 0.0, 1.0}   // Blue
};

timi::Interpolator1D<double, Color> gradient(positions, colors);
Color c = gradient(0.25);  // Interpolated color
```

## API Reference

### Functional API

One-shot interpolation without storing data:

```cpp
// 1D interpolation
template<typename Scalar, typename T>
T timi::interpolate(
    const std::vector<Scalar>& x_vals,  // Sorted ascending coordinates
    const std::vector<T>& y_vals,       // Corresponding values
    Scalar x                            // Query point
);

// N-D interpolation
template<typename Scalar, typename T>
T timi::interpolate(
    const std::vector<std::vector<Scalar>>& axes,  // Coordinate vectors per dimension
    const std::vector<T>& values,                  // Flat array in row-major order
    const std::vector<Scalar>& point               // Query point
);
```

### Stateful API

Reusable interpolators that store data:

```cpp
// 1D interpolator
template<typename Scalar, typename T>
class timi::Interpolator1D {
    Interpolator1D(std::vector<Scalar> x_vals, std::vector<T> y_vals);
    T operator()(Scalar x) const;

    const std::vector<Scalar>& x_values() const;
    const std::vector<T>& y_values() const;
    std::size_t size() const;
};

// N-D interpolator
template<typename Scalar, typename T>
class timi::InterpolatorND {
    InterpolatorND(std::vector<std::vector<Scalar>> axes, std::vector<T> values);
    T operator()(const std::vector<Scalar>& point) const;

    std::size_t num_dimensions() const;
    const std::vector<std::vector<Scalar>>& axes() const;
    const std::vector<T>& values() const;
};
```

### Data Layout

For N-dimensional grids, values are stored in row-major order where the **first axis varies fastest**:

```
2D example with axes x=[0,1], y=[0,1]:
values = {f(0,0), f(1,0), f(0,1), f(1,1)}
         └─x=0─┘  └─x=1─┘  └─x=0─┘  └─x=1─┘
          y=0       y=0      y=1      y=1
```

### Constraints

- Coordinate vectors must be sorted in ascending order
- Each axis requires at least 2 points
- Value count must equal the product of axis sizes
- Maximum 16 dimensions supported

## Building

```bash
# Configure and build
cmake -B build
cmake --build build

# Run tests
cd build && ctest --output-on-failure

# Or run directly
./build/timi_tests

# Run examples
./build/basic_2d
./build/custom_class
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) fore the details.