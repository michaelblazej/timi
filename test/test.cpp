#include <timi/timi.hpp>
#include <gtest/gtest.h>
#include <cmath>

// =============================================================================
// Custom Type Definitions
// =============================================================================

struct Vec2 {
    double x, y;
};

Vec2 operator+(const Vec2& a, const Vec2& b) {
    return {a.x + b.x, a.y + b.y};
}

Vec2 operator*(const Vec2& v, double s) {
    return {v.x * s, v.y * s};
}

struct Vec3 {
    double x, y, z;
};

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator*(const Vec3& v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}

struct Color {
    double r, g, b;
};

Color operator+(const Color& a, const Color& b) {
    return {a.r + b.r, a.g + b.g, a.b + b.b};
}

Color operator*(const Color& c, double s) {
    return {c.r * s, c.g * s, c.b * s};
}

struct Complex {
    double re, im;
};

Complex operator+(const Complex& a, const Complex& b) {
    return {a.re + b.re, a.im + b.im};
}

Complex operator*(const Complex& c, double s) {
    return {c.re * s, c.im * s};
}

struct FloatVec2 {
    float x, y;
};

FloatVec2 operator+(const FloatVec2& a, const FloatVec2& b) {
    return {a.x + b.x, a.y + b.y};
}

FloatVec2 operator*(const FloatVec2& v, float s) {
    return {v.x * s, v.y * s};
}

// =============================================================================
// Custom Assertion Helpers for Custom Types
// =============================================================================

void ExpectVec2Near(const Vec2& actual, const Vec2& expected, double eps = 1e-9) {
    EXPECT_NEAR(actual.x, expected.x, eps) << "Vec2.x mismatch";
    EXPECT_NEAR(actual.y, expected.y, eps) << "Vec2.y mismatch";
}

void ExpectVec3Near(const Vec3& actual, const Vec3& expected, double eps = 1e-9) {
    EXPECT_NEAR(actual.x, expected.x, eps) << "Vec3.x mismatch";
    EXPECT_NEAR(actual.y, expected.y, eps) << "Vec3.y mismatch";
    EXPECT_NEAR(actual.z, expected.z, eps) << "Vec3.z mismatch";
}

void ExpectColorNear(const Color& actual, const Color& expected, double eps = 1e-9) {
    EXPECT_NEAR(actual.r, expected.r, eps) << "Color.r mismatch";
    EXPECT_NEAR(actual.g, expected.g, eps) << "Color.g mismatch";
    EXPECT_NEAR(actual.b, expected.b, eps) << "Color.b mismatch";
}

void ExpectComplexNear(const Complex& actual, const Complex& expected, double eps = 1e-9) {
    EXPECT_NEAR(actual.re, expected.re, eps) << "Complex.re mismatch";
    EXPECT_NEAR(actual.im, expected.im, eps) << "Complex.im mismatch";
}

void ExpectFloatVec2Near(const FloatVec2& actual, const FloatVec2& expected, float eps = 1e-5f) {
    EXPECT_NEAR(actual.x, expected.x, eps) << "FloatVec2.x mismatch";
    EXPECT_NEAR(actual.y, expected.y, eps) << "FloatVec2.y mismatch";
}

// =============================================================================
// 1D Interpolation Tests
// =============================================================================

TEST(Interp1D, Basic) {
    // Linear function y = 2x
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 2.0, 4.0};

    // Interior points
    EXPECT_NEAR(timi::interpolate(x, y, 0.5), 1.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 1.5), 3.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 0.25), 0.5, 1e-9);

    // Exact grid points
    EXPECT_NEAR(timi::interpolate(x, y, 0.0), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 1.0), 2.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 2.0), 4.0, 1e-9);
}

TEST(Interp1D, OutOfRange) {
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 10.0, 20.0};

    // Out of range (clamped)
    EXPECT_NEAR(timi::interpolate(x, y, -1.0), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, -100.0), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 3.0), 20.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 100.0), 20.0, 1e-9);
}

TEST(Interp1D, Stateful) {
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 10.0, 20.0};

    timi::Interpolator1D<double, double> interp(x, y);

    // Multiple samples
    EXPECT_NEAR(interp(0.25), 2.5, 1e-9);
    EXPECT_NEAR(interp(0.5), 5.0, 1e-9);
    EXPECT_NEAR(interp(1.5), 15.0, 1e-9);
    EXPECT_NEAR(interp(0.0), 0.0, 1e-9);
    EXPECT_NEAR(interp(2.0), 20.0, 1e-9);
}

TEST(Interp1D, StatefulRepeated) {
    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {0.0, 100.0};

    timi::Interpolator1D<double, double> interp(x, y);

    // Repeated sampling at same point
    for (int i = 0; i < 1000; ++i) {
        EXPECT_NEAR(interp(0.5), 50.0, 1e-9);
    }
}

TEST(Interp1D, NonuniformGrid) {
    // Non-uniform spacing
    std::vector<double> x = {0.0, 0.1, 0.5, 2.0};
    std::vector<double> y = {0.0, 1.0, 5.0, 20.0};

    // Interpolate in each interval
    EXPECT_NEAR(timi::interpolate(x, y, 0.05), 0.5, 1e-9);   // halfway in [0, 0.1]
    EXPECT_NEAR(timi::interpolate(x, y, 0.3), 3.0, 1e-9);    // halfway in [0.1, 0.5]
    EXPECT_NEAR(timi::interpolate(x, y, 1.25), 12.5, 1e-9);  // halfway in [0.5, 2.0]
}

// =============================================================================
// 2D Interpolation Tests
// =============================================================================

TEST(Interp2D, Basic) {
    // 2x2 grid
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0}
    };
    // Values in row-major: (0,0), (1,0), (0,1), (1,1)
    // f(x,y) = x + 2*y: f(0,0)=0, f(1,0)=1, f(0,1)=2, f(1,1)=3
    std::vector<double> values = {0.0, 1.0, 2.0, 3.0};

    // Center point: (0.5, 0.5) -> f = 0.5 + 2*0.5 = 1.5
    double v = timi::interpolate(axes, values, {0.5, 0.5});
    EXPECT_NEAR(v, 1.5, 1e-9);

    // Corner points
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 0.0}), 1.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 1.0}), 2.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0}), 3.0, 1e-9);

    // Edge midpoints
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, 0.0}), 0.5, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, 1.0}), 2.5, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.5}), 1.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 0.5}), 2.0, 1e-9);
}

TEST(Interp2D, Stateful) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 2.0}
    };
    std::vector<double> values = {0.0, 10.0, 20.0, 30.0};

    timi::InterpolatorND<double, double> interp(axes, values);

    // Multiple samples
    EXPECT_NEAR(interp({0.5, 1.0}), 15.0, 1e-9);
    EXPECT_NEAR(interp({0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(interp({1.0, 2.0}), 30.0, 1e-9);
}

TEST(Interp2D, Grid3x3) {
    // 3x3 grid: f(x,y) = x + y
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0, 2.0},
        {0.0, 1.0, 2.0}
    };
    // Row-major: (0,0), (1,0), (2,0), (0,1), (1,1), (2,1), (0,2), (1,2), (2,2)
    std::vector<double> values = {
        0.0, 1.0, 2.0,  // y=0
        1.0, 2.0, 3.0,  // y=1
        2.0, 3.0, 4.0   // y=2
    };

    // Test various points
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, 0.5}), 1.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.5, 1.5}), 3.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 2.0}), 2.0, 1e-9);
}

// =============================================================================
// 3D Interpolation Tests
// =============================================================================

TEST(Interp3D, Basic) {
    // 2x2x2 grid
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    // 8 values in row-major order
    // Index = x + 2*y + 4*z
    std::vector<double> values = {
        0.0, 1.0, 2.0, 3.0,  // z=0 plane: (0,0,0), (1,0,0), (0,1,0), (1,1,0)
        4.0, 5.0, 6.0, 7.0   // z=1 plane: (0,0,1), (1,0,1), (0,1,1), (1,1,1)
    };

    // Center point: average of all 8 corners = (0+1+2+3+4+5+6+7)/8 = 3.5
    double v = timi::interpolate(axes, values, {0.5, 0.5, 0.5});
    EXPECT_NEAR(v, 3.5, 1e-9);

    // Corner points
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0, 1.0}), 7.0, 1e-9);
}

TEST(Interp3D, Stateful) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<double> values = {
        0.0, 1.0, 2.0, 3.0,
        4.0, 5.0, 6.0, 7.0
    };

    timi::InterpolatorND<double, double> interp(axes, values);

    EXPECT_NEAR(interp({0.5, 0.5, 0.5}), 3.5, 1e-9);
    EXPECT_NEAR(interp({0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(interp({1.0, 1.0, 1.0}), 7.0, 1e-9);
}

// =============================================================================
// Custom Type Tests - Vec2
// =============================================================================

TEST(CustomTypes, Vec2_1D) {
    std::vector<double> x = {0.0, 1.0};
    std::vector<Vec2> y = {{0.0, 0.0}, {2.0, 4.0}};

    Vec2 result = timi::interpolate(x, y, 0.5);
    EXPECT_NEAR(result.x, 1.0, 1e-9);
    EXPECT_NEAR(result.y, 2.0, 1e-9);

    // Using stateful API
    timi::Interpolator1D<double, Vec2> interp(x, y);
    Vec2 result2 = interp(0.25);
    EXPECT_NEAR(result2.x, 0.5, 1e-9);
    EXPECT_NEAR(result2.y, 1.0, 1e-9);
}

TEST(CustomTypes, Vec2_2D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<Vec2> values = {
        {0.0, 0.0}, {1.0, 0.0},
        {0.0, 1.0}, {1.0, 1.0}
    };

    Vec2 result = timi::interpolate(axes, values, {0.5, 0.5});
    EXPECT_NEAR(result.x, 0.5, 1e-9);
    EXPECT_NEAR(result.y, 0.5, 1e-9);
}

// =============================================================================
// Custom Type Tests - Vec3
// =============================================================================

TEST(CustomTypes, Vec3_1D) {
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<Vec3> y = {
        {0.0, 0.0, 0.0},
        {1.0, 2.0, 3.0},
        {2.0, 4.0, 6.0}
    };

    Vec3 result = timi::interpolate(x, y, 0.5);
    ExpectVec3Near(result, Vec3{0.5, 1.0, 1.5});

    result = timi::interpolate(x, y, 1.5);
    ExpectVec3Near(result, Vec3{1.5, 3.0, 4.5});

    // Stateful API
    timi::Interpolator1D<double, Vec3> interp(x, y);
    result = interp(0.25);
    ExpectVec3Near(result, Vec3{0.25, 0.5, 0.75});
}

TEST(CustomTypes, Vec3_2D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<Vec3> values = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
    };

    Vec3 result = timi::interpolate(axes, values, {0.5, 0.5});
    ExpectVec3Near(result, Vec3{0.5, 0.5, 0.0});

    // Corner points
    ExpectVec3Near(timi::interpolate(axes, values, {0.0, 0.0}), Vec3{0.0, 0.0, 0.0});
    ExpectVec3Near(timi::interpolate(axes, values, {1.0, 1.0}), Vec3{1.0, 1.0, 0.0});
}

TEST(CustomTypes, Vec3_3D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    // 8 Vec3 values for 2x2x2 grid
    std::vector<Vec3> values = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}
    };

    Vec3 result = timi::interpolate(axes, values, {0.5, 0.5, 0.5});
    ExpectVec3Near(result, Vec3{0.5, 0.5, 0.5});

    // Stateful API
    timi::InterpolatorND<double, Vec3> interp(axes, values);
    ExpectVec3Near(interp({0.0, 0.0, 0.0}), Vec3{0.0, 0.0, 0.0});
    ExpectVec3Near(interp({1.0, 1.0, 1.0}), Vec3{1.0, 1.0, 1.0});
}

// =============================================================================
// Custom Type Tests - Color
// =============================================================================

TEST(CustomTypes, Color_1D) {
    // Color gradient from black to red to white
    std::vector<double> x = {0.0, 0.5, 1.0};
    std::vector<Color> y = {
        {0.0, 0.0, 0.0},  // black
        {1.0, 0.0, 0.0},  // red
        {1.0, 1.0, 1.0}   // white
    };

    Color result = timi::interpolate(x, y, 0.25);
    ExpectColorNear(result, Color{0.5, 0.0, 0.0});

    result = timi::interpolate(x, y, 0.75);
    ExpectColorNear(result, Color{1.0, 0.5, 0.5});

    // Stateful API
    timi::Interpolator1D<double, Color> interp(x, y);
    ExpectColorNear(interp(0.0), Color{0.0, 0.0, 0.0});
    ExpectColorNear(interp(1.0), Color{1.0, 1.0, 1.0});
}

TEST(CustomTypes, Color_2D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0}
    };
    // 2x2 color grid
    std::vector<Color> values = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},  // black, red
        {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}   // green, yellow
    };

    Color result = timi::interpolate(axes, values, {0.5, 0.5});
    ExpectColorNear(result, Color{0.5, 0.5, 0.0});

    // Edge midpoints
    ExpectColorNear(timi::interpolate(axes, values, {0.5, 0.0}), Color{0.5, 0.0, 0.0});
    ExpectColorNear(timi::interpolate(axes, values, {0.0, 0.5}), Color{0.0, 0.5, 0.0});
}

TEST(CustomTypes, Color_3D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    // RGB color cube corners
    std::vector<Color> values = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}
    };

    // Center of color cube should be gray
    Color result = timi::interpolate(axes, values, {0.5, 0.5, 0.5});
    ExpectColorNear(result, Color{0.5, 0.5, 0.5});

    // Stateful API
    timi::InterpolatorND<double, Color> interp(axes, values);
    ExpectColorNear(interp({0.0, 0.0, 0.0}), Color{0.0, 0.0, 0.0});
    ExpectColorNear(interp({1.0, 1.0, 1.0}), Color{1.0, 1.0, 1.0});
}

// =============================================================================
// Custom Type Tests - Complex
// =============================================================================

TEST(CustomTypes, Complex_1D) {
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<Complex> y = {
        {0.0, 0.0},
        {1.0, 1.0},
        {2.0, 0.0}
    };

    Complex result = timi::interpolate(x, y, 0.5);
    ExpectComplexNear(result, Complex{0.5, 0.5});

    result = timi::interpolate(x, y, 1.5);
    ExpectComplexNear(result, Complex{1.5, 0.5});

    // Stateful API
    timi::Interpolator1D<double, Complex> interp(x, y);
    ExpectComplexNear(interp(1.0), Complex{1.0, 1.0});
}

TEST(CustomTypes, Complex_2D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<Complex> values = {
        {0.0, 0.0}, {1.0, 0.0},
        {0.0, 1.0}, {1.0, 1.0}
    };

    Complex result = timi::interpolate(axes, values, {0.5, 0.5});
    ExpectComplexNear(result, Complex{0.5, 0.5});

    // Stateful API
    timi::InterpolatorND<double, Complex> interp(axes, values);
    ExpectComplexNear(interp({0.0, 0.0}), Complex{0.0, 0.0});
    ExpectComplexNear(interp({1.0, 1.0}), Complex{1.0, 1.0});
}

// =============================================================================
// Custom Type Tests - FloatVec2
// =============================================================================

TEST(CustomTypes, FloatVec2_1D) {
    std::vector<float> x = {0.0f, 1.0f};
    std::vector<FloatVec2> y = {{0.0f, 0.0f}, {2.0f, 4.0f}};

    FloatVec2 result = timi::interpolate(x, y, 0.5f);
    ExpectFloatVec2Near(result, FloatVec2{1.0f, 2.0f});

    // Stateful API
    timi::Interpolator1D<float, FloatVec2> interp(x, y);
    result = interp(0.25f);
    ExpectFloatVec2Near(result, FloatVec2{0.5f, 1.0f});
}

TEST(CustomTypes, FloatVec2_2D) {
    std::vector<std::vector<float>> axes = {
        {0.0f, 1.0f},
        {0.0f, 1.0f}
    };
    std::vector<FloatVec2> values = {
        {0.0f, 0.0f}, {1.0f, 0.0f},
        {0.0f, 1.0f}, {1.0f, 1.0f}
    };

    FloatVec2 result = timi::interpolate(axes, values, {0.5f, 0.5f});
    ExpectFloatVec2Near(result, FloatVec2{0.5f, 0.5f});

    // Stateful API
    timi::InterpolatorND<float, FloatVec2> interp(axes, values);
    ExpectFloatVec2Near(interp({0.0f, 0.0f}), FloatVec2{0.0f, 0.0f});
    ExpectFloatVec2Near(interp({1.0f, 1.0f}), FloatVec2{1.0f, 1.0f});
}

// =============================================================================
// Float Type Tests
// =============================================================================

TEST(FloatType, Basic1D) {
    std::vector<float> x = {0.0f, 1.0f, 2.0f};
    std::vector<float> y = {0.0f, 10.0f, 20.0f};

    float result = timi::interpolate(x, y, 0.5f);
    EXPECT_NEAR(result, 5.0f, 1e-5f);

    timi::Interpolator1D<float, float> interp(x, y);
    EXPECT_NEAR(interp(1.5f), 15.0f, 1e-5f);
}

TEST(FloatType, Basic2D) {
    std::vector<std::vector<float>> axes = {
        {0.0f, 1.0f},
        {0.0f, 1.0f}
    };
    std::vector<float> values = {0.0f, 1.0f, 2.0f, 3.0f};

    float v = timi::interpolate(axes, values, {0.5f, 0.5f});
    EXPECT_NEAR(v, 1.5f, 1e-5f);

    // Corner points
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0f, 0.0f}), 0.0f, 1e-5f);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0f, 1.0f}), 3.0f, 1e-5f);

    // Stateful API
    timi::InterpolatorND<float, float> interp(axes, values);
    EXPECT_NEAR(interp({0.5f, 0.5f}), 1.5f, 1e-5f);
}

TEST(FloatType, Basic3D) {
    std::vector<std::vector<float>> axes = {
        {0.0f, 1.0f},
        {0.0f, 1.0f},
        {0.0f, 1.0f}
    };
    std::vector<float> values = {
        0.0f, 1.0f, 2.0f, 3.0f,
        4.0f, 5.0f, 6.0f, 7.0f
    };

    float v = timi::interpolate(axes, values, {0.5f, 0.5f, 0.5f});
    EXPECT_NEAR(v, 3.5f, 1e-5f);

    EXPECT_NEAR(timi::interpolate(axes, values, {0.0f, 0.0f, 0.0f}), 0.0f, 1e-5f);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0f, 1.0f, 1.0f}), 7.0f, 1e-5f);
}

// =============================================================================
// Long Double Type Tests
// =============================================================================

TEST(LongDoubleType, Basic1D) {
    std::vector<long double> x = {0.0L, 1.0L, 2.0L};
    std::vector<long double> y = {0.0L, 10.0L, 20.0L};

    long double result = timi::interpolate(x, y, 0.5L);
    EXPECT_NEAR(static_cast<double>(result), 5.0, 1e-12);

    result = timi::interpolate(x, y, 1.5L);
    EXPECT_NEAR(static_cast<double>(result), 15.0, 1e-12);

    // Stateful API
    timi::Interpolator1D<long double, long double> interp(x, y);
    EXPECT_NEAR(static_cast<double>(interp(0.25L)), 2.5, 1e-12);
}

TEST(LongDoubleType, Basic2D) {
    std::vector<std::vector<long double>> axes = {
        {0.0L, 1.0L},
        {0.0L, 1.0L}
    };
    std::vector<long double> values = {0.0L, 1.0L, 2.0L, 3.0L};

    long double v = timi::interpolate(axes, values, {0.5L, 0.5L});
    EXPECT_NEAR(static_cast<double>(v), 1.5, 1e-12);

    // Stateful API
    timi::InterpolatorND<long double, long double> interp(axes, values);
    EXPECT_NEAR(static_cast<double>(interp({0.0L, 0.0L})), 0.0, 1e-12);
    EXPECT_NEAR(static_cast<double>(interp({1.0L, 1.0L})), 3.0, 1e-12);
}

// =============================================================================
// Higher Dimension Tests
// =============================================================================

TEST(HigherDims, Basic4D) {
    // 2x2x2x2 grid = 16 values
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };

    // Values are just indices 0-15
    std::vector<double> values(16);
    for (int i = 0; i < 16; ++i) {
        values[i] = static_cast<double>(i);
    }

    // Center point: average of all 16 corners = (0+1+...+15)/16 = 7.5
    double v = timi::interpolate(axes, values, {0.5, 0.5, 0.5, 0.5});
    EXPECT_NEAR(v, 7.5, 1e-9);

    // Corner tests
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0, 1.0, 1.0}), 15.0, 1e-9);
}

TEST(HigherDims, Stateful4D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };

    std::vector<double> values(16);
    for (int i = 0; i < 16; ++i) {
        values[i] = static_cast<double>(i);
    }

    timi::InterpolatorND<double, double> interp(axes, values);

    EXPECT_NEAR(interp({0.5, 0.5, 0.5, 0.5}), 7.5, 1e-9);
    EXPECT_NEAR(interp({0.0, 0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(interp({1.0, 1.0, 1.0, 1.0}), 15.0, 1e-9);

    // Test some edge midpoints
    EXPECT_NEAR(interp({0.5, 0.0, 0.0, 0.0}), 0.5, 1e-9);
    EXPECT_NEAR(interp({0.0, 0.5, 0.0, 0.0}), 1.0, 1e-9);
}

TEST(HigherDims, Vec2_4D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };

    // 16 Vec2 values
    std::vector<Vec2> values(16);
    for (int i = 0; i < 16; ++i) {
        values[i] = {static_cast<double>(i), static_cast<double>(i * 2)};
    }

    Vec2 result = timi::interpolate(axes, values, {0.5, 0.5, 0.5, 0.5});
    ExpectVec2Near(result, Vec2{7.5, 15.0});

    timi::InterpolatorND<double, Vec2> interp(axes, values);
    ExpectVec2Near(interp({0.0, 0.0, 0.0, 0.0}), Vec2{0.0, 0.0});
    ExpectVec2Near(interp({1.0, 1.0, 1.0, 1.0}), Vec2{15.0, 30.0});
}

TEST(HigherDims, Basic5D) {
    // 2x2x2x2x2 grid = 32 values
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };

    std::vector<double> values(32);
    for (int i = 0; i < 32; ++i) {
        values[i] = static_cast<double>(i);
    }

    // Center: average of 0-31 = 15.5
    double v = timi::interpolate(axes, values, {0.5, 0.5, 0.5, 0.5, 0.5});
    EXPECT_NEAR(v, 15.5, 1e-9);

    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0, 0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0, 1.0, 1.0, 1.0}), 31.0, 1e-9);
}

TEST(HigherDims, Stateful5D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };

    std::vector<double> values(32);
    for (int i = 0; i < 32; ++i) {
        values[i] = static_cast<double>(i);
    }

    timi::InterpolatorND<double, double> interp(axes, values);

    EXPECT_NEAR(interp({0.5, 0.5, 0.5, 0.5, 0.5}), 15.5, 1e-9);
    EXPECT_NEAR(interp({0.0, 0.0, 0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(interp({1.0, 1.0, 1.0, 1.0, 1.0}), 31.0, 1e-9);
}

TEST(HigherDims, Basic6D) {
    // 2x2x2x2x2x2 grid = 64 values
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };

    std::vector<double> values(64);
    for (int i = 0; i < 64; ++i) {
        values[i] = static_cast<double>(i);
    }

    // Center: average of 0-63 = 31.5
    double v = timi::interpolate(axes, values, {0.5, 0.5, 0.5, 0.5, 0.5, 0.5});
    EXPECT_NEAR(v, 31.5, 1e-9);

    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}), 63.0, 1e-9);
}

TEST(HigherDims, MaxDims16D) {
    // Test with maximum 16 dimensions (2^16 = 65536 values)
    // Use smaller grid (all dims size 2) to keep test fast
    std::vector<std::vector<double>> axes(16, {0.0, 1.0});

    std::vector<double> values(65536);
    for (size_t i = 0; i < 65536; ++i) {
        values[i] = static_cast<double>(i);
    }

    // Center point
    std::vector<double> center(16, 0.5);
    double v = timi::interpolate(axes, values, center);
    // Average of 0 to 65535 = 32767.5
    EXPECT_NEAR(v, 32767.5, 1e-9);

    // Corner points
    std::vector<double> origin(16, 0.0);
    std::vector<double> far_corner(16, 1.0);
    EXPECT_NEAR(timi::interpolate(axes, values, origin), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, far_corner), 65535.0, 1e-9);
}

// =============================================================================
// Boundary Clamping Tests
// =============================================================================

TEST(BoundaryClamping, Clamping2D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<double> values = {0.0, 1.0, 2.0, 3.0};

    // Out of range should clamp
    EXPECT_NEAR(timi::interpolate(axes, values, {-1.0, 0.5}), 1.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {2.0, 0.5}), 2.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, -1.0}), 0.5, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, 2.0}), 2.5, 1e-9);

    // Both out of range
    EXPECT_NEAR(timi::interpolate(axes, values, {-1.0, -1.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {2.0, 2.0}), 3.0, 1e-9);
}

TEST(BoundaryClamping, Clamping3D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<double> values = {
        0.0, 1.0, 2.0, 3.0,
        4.0, 5.0, 6.0, 7.0
    };

    // Out of range in various dimensions
    EXPECT_NEAR(timi::interpolate(axes, values, {-1.0, 0.5, 0.5}), 3.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, -1.0, 0.5}), 2.5, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, 0.5, -1.0}), 1.5, 1e-9);

    // All out of range
    EXPECT_NEAR(timi::interpolate(axes, values, {-1.0, -1.0, -1.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {2.0, 2.0, 2.0}), 7.0, 1e-9);
}

// =============================================================================
// Edge Case Tests
// =============================================================================

TEST(EdgeCases, TwoPointGrid) {
    // Minimum valid grid: 2 points
    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {10.0, 20.0};

    EXPECT_NEAR(timi::interpolate(x, y, 0.0), 10.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 0.5), 15.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 1.0), 20.0, 1e-9);
}

TEST(EdgeCases, ExactGridPoints) {
    // Test that exact grid points return exact values
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0, 2.0},
        {0.0, 1.0}
    };
    std::vector<double> values = {
        0.0, 1.0, 2.0,
        3.0, 4.0, 5.0
    };

    // All exact grid points
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 0.0}), 1.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {2.0, 0.0}), 2.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 1.0}), 3.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0}), 4.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {2.0, 1.0}), 5.0, 1e-9);
}

TEST(EdgeCases, NegativeCoordinates) {
    // Grid with negative coordinates
    std::vector<double> x = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<double> y = {-20.0, -10.0, 0.0, 10.0, 20.0};

    EXPECT_NEAR(timi::interpolate(x, y, -1.5), -15.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, 0.5), 5.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(x, y, -0.5), -5.0, 1e-9);

    // 2D with negative coordinates
    std::vector<std::vector<double>> axes = {
        {-1.0, 1.0},
        {-1.0, 1.0}
    };
    std::vector<double> values = {0.0, 1.0, 2.0, 3.0};

    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0}), 1.5, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {-1.0, -1.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0}), 3.0, 1e-9);
}

TEST(EdgeCases, NonuniformSpacing2D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 0.1, 1.0},   // non-uniform x
        {0.0, 0.9, 1.0}    // non-uniform y
    };
    // Values: f(x,y) = x + y at grid points
    std::vector<double> values = {
        0.0, 0.1, 1.0,    // y=0
        0.9, 1.0, 1.9,    // y=0.9
        1.0, 1.1, 2.0     // y=1
    };

    // Test midpoints within non-uniform cells
    double v = timi::interpolate(axes, values, {0.05, 0.45});
    // Midpoint between (0, 0) and (0.1, 0.9) cell
    EXPECT_NEAR(v, 0.5, 1e-6);

    // Exact grid points
    EXPECT_NEAR(timi::interpolate(axes, values, {0.1, 0.9}), 1.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0}), 2.0, 1e-9);
}

TEST(EdgeCases, NonuniformSpacing3D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 0.25, 1.0},
        {0.0, 0.5, 1.0},
        {0.0, 1.0}
    };

    // 3x3x2 = 18 values
    std::vector<double> values(18);
    for (int i = 0; i < 18; ++i) {
        values[i] = static_cast<double>(i);
    }

    // Test corners
    EXPECT_NEAR(timi::interpolate(axes, values, {0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(timi::interpolate(axes, values, {1.0, 1.0, 1.0}), 17.0, 1e-9);

    // Stateful API
    timi::InterpolatorND<double, double> interp(axes, values);
    EXPECT_NEAR(interp({0.0, 0.0, 0.0}), 0.0, 1e-9);
    EXPECT_NEAR(interp({1.0, 1.0, 1.0}), 17.0, 1e-9);
}

TEST(EdgeCases, SingleCellGrid) {
    // Minimum 2x2x2 grid (single cell)
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<double> values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

    // Center
    EXPECT_NEAR(timi::interpolate(axes, values, {0.5, 0.5, 0.5}), 3.5, 1e-9);

    // Test that it works with 4D single cell too
    std::vector<std::vector<double>> axes4d = {
        {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}
    };
    std::vector<double> values4d(16);
    for (int i = 0; i < 16; ++i) values4d[i] = static_cast<double>(i);

    EXPECT_NEAR(timi::interpolate(axes4d, values4d, {0.5, 0.5, 0.5, 0.5}), 7.5, 1e-9);
}

// =============================================================================
// API Consistency Tests
// =============================================================================

TEST(APIConsistency, Interp1D) {
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 10.0, 15.0, 25.0, 40.0};

    timi::Interpolator1D<double, double> interp(x, y);

    // Functional and stateful should give same results
    for (double t = -0.5; t <= 4.5; t += 0.1) {
        double func_result = timi::interpolate(x, y, t);
        double stat_result = interp(t);
        EXPECT_NEAR(func_result, stat_result, 1e-12);
    }
}

TEST(APIConsistency, InterpND) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0, 2.0},
        {0.0, 1.0}
    };
    std::vector<double> values = {
        0.0, 1.0, 2.0,
        3.0, 4.0, 5.0
    };

    timi::InterpolatorND<double, double> interp(axes, values);

    // Test multiple points
    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0}, {1.0, 0.5}, {0.5, 0.5}, {2.0, 1.0}, {1.5, 0.75}
    };

    for (const auto& pt : test_points) {
        double func_result = timi::interpolate(axes, values, pt);
        double stat_result = interp(pt);
        EXPECT_NEAR(func_result, stat_result, 1e-12);
    }
}

TEST(APIConsistency, Interp3D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0, 2.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<double> values(12);
    for (int i = 0; i < 12; ++i) {
        values[i] = static_cast<double>(i);
    }

    timi::InterpolatorND<double, double> interp(axes, values);

    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0, 0.0}, {1.0, 0.5, 0.5}, {0.5, 0.5, 0.5},
        {2.0, 1.0, 1.0}, {1.5, 0.75, 0.25}
    };

    for (const auto& pt : test_points) {
        double func_result = timi::interpolate(axes, values, pt);
        double stat_result = interp(pt);
        EXPECT_NEAR(func_result, stat_result, 1e-12);
    }
}

TEST(APIConsistency, Interp4D) {
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<double> values(16);
    for (int i = 0; i < 16; ++i) {
        values[i] = static_cast<double>(i);
    }

    timi::InterpolatorND<double, double> interp(axes, values);

    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0, 0.0, 0.0}, {0.5, 0.5, 0.5, 0.5},
        {1.0, 1.0, 1.0, 1.0}, {0.25, 0.75, 0.5, 0.5}
    };

    for (const auto& pt : test_points) {
        double func_result = timi::interpolate(axes, values, pt);
        double stat_result = interp(pt);
        EXPECT_NEAR(func_result, stat_result, 1e-12);
    }
}

TEST(APIConsistency, AccessorMethods) {
    // Test 1D accessor methods
    std::vector<double> x1d = {0.0, 1.0, 2.0};
    std::vector<double> y1d = {10.0, 20.0, 30.0};
    timi::Interpolator1D<double, double> interp1d(x1d, y1d);

    EXPECT_EQ(interp1d.size(), 3u);
    EXPECT_EQ(interp1d.x_values().size(), 3u);
    EXPECT_EQ(interp1d.y_values().size(), 3u);
    EXPECT_NEAR(interp1d.x_values()[0], 0.0, 1e-9);
    EXPECT_NEAR(interp1d.y_values()[2], 30.0, 1e-9);

    // Test ND accessor methods
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0, 2.0}
    };
    std::vector<double> values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    timi::InterpolatorND<double, double> interpnd(axes, values);

    EXPECT_EQ(interpnd.num_dimensions(), 2u);
    EXPECT_EQ(interpnd.axes().size(), 2u);
    EXPECT_EQ(interpnd.axes()[0].size(), 2u);
    EXPECT_EQ(interpnd.axes()[1].size(), 3u);
    EXPECT_EQ(interpnd.values().size(), 6u);
}
