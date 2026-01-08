# Design Prompt: `timi` — C++17 N-Dimensional Linear Interpolation Library

This document defines the design and implementation prompt for **`timi`**, a modern C++ library for **n-dimensional linear interpolation** over templated input and output types. The library supports both **functional** and **stateful (object-oriented)** usage patterns.

---

## 1. Goals & Constraints

- Language: **C++17**
- Compiler: **GCC**
- No external dependencies beyond the C++ standard library
- Header-only preferred (but not required)
- Focus on:
  - Minimal API surface
  - Clear type requirements
  - Practical usage for scientific / numeric code

---

## 2. Core Concept

`timi` operates on **known sampled coordinate/value pairs** and interpolates values at **new coordinates**.

### 1D Concept

Given known sample points:

```
(x0, y0), (x1, y1), ..., (xn, yn)
```

Interpolating at a new coordinate:

```cpp
y = timi::interpolate(x_vals, y_vals, x_new);
```

Returns the linearly interpolated value.

---

### N-D Concept

Generalize this to **N dimensions**:

- Samples lie on a **rectilinear grid**
- Each axis has its own coordinate vector
- Values exist at all grid intersections

Example (2D):

```cpp
axes = {
  {x0, x1},
  {y0, y1}
};

values = {
  f(x0,y0), f(x1,y0),
  f(x0,y1), f(x1,y1)
};
```

Interpolating:

```cpp
value = timi::interpolate(axes, values, {x, y});
```

---

## 3. Functional API

### 1D Interpolation

```cpp
template<typename Scalar, typename T>
T interpolate(
    const std::vector<Scalar>& x_vals,
    const std::vector<T>& y_vals,
    Scalar x
);
```

Requirements:
- `x_vals.size() == y_vals.size()`
- `x_vals` sorted ascending
- Linear interpolation between nearest neighbors
- Out-of-range behavior must be documented (assert or clamp)

---

### N-D Interpolation

```cpp
template<typename Scalar, typename T>
T interpolate(
    const std::vector<std::vector<Scalar>>& axes,
    const std::vector<T>& values,
    const std::vector<Scalar>& point
);
```

Where:
- `axes.size() == point.size() == N`
- `values.size() == Π axes[i].size()`
- Values are stored in **row-major order**
- Interpolation is multilinear

---

## 4. Stateful Interpolator Objects

`timi` must support **stateful interpolator objects** that encapsulate sample data and can be reused for repeated queries.

### 1D Interpolator

```cpp
template<typename Scalar, typename T>
class Interpolator1D {
public:
    Interpolator1D(
        std::vector<Scalar> x_vals,
        std::vector<T> y_vals
    );

    T operator()(Scalar x) const;
};
```

Behavior:
- Stores copies of `x_vals` and `y_vals`
- Validates sizes and ordering on construction
- Sampling is `const`

---

### N-D Interpolator

```cpp
template<typename Scalar, typename T>
class InterpolatorND {
public:
    InterpolatorND(
        std::vector<std::vector<Scalar>> axes,
        std::vector<T> values
    );

    T operator()(const std::vector<Scalar>& point) const;
};
```

Behavior:
- Dimensionality inferred from `axes.size()`
- Precomputes strides and axis sizes
- No allocations during sampling

---

## 5. Type Requirements for `T`

The **minimum required operations** for a value type `T`:

```cpp
T operator+(const T&, const T&);
T operator*(const T&, Scalar);
```

No other requirements:
- No inheritance
- No traits
- No virtual functions

---

## 6. Usage Examples

### 1D Functional Usage

```cpp
std::vector<double> x = {0.0, 1.0, 2.0};
std::vector<double> y = {0.0, 10.0, 20.0};

double v = timi::interpolate(x, y, 0.5);
```

---

### 1D Stateful Usage

```cpp
timi::Interpolator1D<double, double> interp(x, y);

double a = interp(0.25);
double b = interp(1.75);
```

---

### N-D Stateful Usage

```cpp
timi::InterpolatorND<double, double> interp({
    {0.0, 1.0},
    {0.0, 2.0}
}, {
    0.0, 10.0,
    20.0, 30.0
});

double v = interp({0.5, 1.0});
```

---

### Custom Type Example

```cpp
struct Vec2 {
    double x, y;
};

Vec2 operator+(const Vec2& a, const Vec2& b) {
    return {a.x + b.x, a.y + b.y};
}

Vec2 operator*(const Vec2& v, double s) {
    return {v.x * s, v.y * s};
}

Vec2 result = interp(0.5);
```

---

## 7. Tests (Required)

- No third-party frameworks
- Use `assert`
- Include:
  - 1D, 2D, 3D interpolation
  - Boundary cases (`0.0`, `1.0`)
  - Custom type interpolation
  - Repeated sampling of the same interpolator
  - Consistency between functional and stateful APIs

Tests must compile and run with:

```bash
g++ -std=gnu++17 test.cpp && ./a.out
```

---

## 8. Documentation Requirements

Include documentation explaining:
- Axis ordering and value layout
- Functional vs stateful API usage
- Construction cost vs sampling cost
- Thread safety (sampling must be thread-safe if immutable)

---

## 9. Non-Goals / Optional Future Work

These should be mentioned but **not implemented**:

- Extrapolation policies
- Compile-time fixed dimensions
- `std::span`-based views
- SIMD acceleration

---

## 10. Output Expectations

Produce:

1. `timi.hpp`
2. `test.cpp`
3. Example usage (inline or separate file)

All code must compile cleanly under **GCC with C++17**.

