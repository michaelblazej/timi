#include <timi/timi.hpp>
#include <cassert>
#include <cmath>
#include <iostream>

// Helper for approximate equality
template<typename T>
bool approx_equal(T a, T b, T eps = T(1e-9)) {
    return std::abs(a - b) < eps;
}

// =============================================================================
// 1D Tests
// =============================================================================

void test_1d_basic() {
    std::cout << "test_1d_basic... ";

    // Linear function y = 2x
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 2.0, 4.0};

    // Interior points
    assert(approx_equal(timi::interpolate(x, y, 0.5), 1.0));
    assert(approx_equal(timi::interpolate(x, y, 1.5), 3.0));
    assert(approx_equal(timi::interpolate(x, y, 0.25), 0.5));

    // Exact grid points
    assert(approx_equal(timi::interpolate(x, y, 0.0), 0.0));
    assert(approx_equal(timi::interpolate(x, y, 1.0), 2.0));
    assert(approx_equal(timi::interpolate(x, y, 2.0), 4.0));

    std::cout << "PASSED\n";
}

void test_1d_out_of_range() {
    std::cout << "test_1d_out_of_range... ";

    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 10.0, 20.0};

    // Out of range (clamped)
    assert(approx_equal(timi::interpolate(x, y, -1.0), 0.0));
    assert(approx_equal(timi::interpolate(x, y, -100.0), 0.0));
    assert(approx_equal(timi::interpolate(x, y, 3.0), 20.0));
    assert(approx_equal(timi::interpolate(x, y, 100.0), 20.0));

    std::cout << "PASSED\n";
}

void test_1d_stateful() {
    std::cout << "test_1d_stateful... ";

    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 10.0, 20.0};

    timi::Interpolator1D<double, double> interp(x, y);

    // Multiple samples
    assert(approx_equal(interp(0.25), 2.5));
    assert(approx_equal(interp(0.5), 5.0));
    assert(approx_equal(interp(1.5), 15.0));
    assert(approx_equal(interp(0.0), 0.0));
    assert(approx_equal(interp(2.0), 20.0));

    std::cout << "PASSED\n";
}

void test_1d_stateful_repeated() {
    std::cout << "test_1d_stateful_repeated... ";

    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {0.0, 100.0};

    timi::Interpolator1D<double, double> interp(x, y);

    // Repeated sampling at same point
    for (int i = 0; i < 1000; ++i) {
        assert(approx_equal(interp(0.5), 50.0));
    }

    std::cout << "PASSED\n";
}

void test_1d_nonuniform_grid() {
    std::cout << "test_1d_nonuniform_grid... ";

    // Non-uniform spacing
    std::vector<double> x = {0.0, 0.1, 0.5, 2.0};
    std::vector<double> y = {0.0, 1.0, 5.0, 20.0};

    // Interpolate in each interval
    assert(approx_equal(timi::interpolate(x, y, 0.05), 0.5));  // halfway in [0, 0.1]
    assert(approx_equal(timi::interpolate(x, y, 0.3), 3.0));   // halfway in [0.1, 0.5]
    assert(approx_equal(timi::interpolate(x, y, 1.25), 12.5)); // halfway in [0.5, 2.0]

    std::cout << "PASSED\n";
}

// =============================================================================
// 2D Tests
// =============================================================================

void test_2d_basic() {
    std::cout << "test_2d_basic... ";

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
    assert(approx_equal(v, 1.5));

    // Corner points
    assert(approx_equal(timi::interpolate(axes, values, {0.0, 0.0}), 0.0));
    assert(approx_equal(timi::interpolate(axes, values, {1.0, 0.0}), 1.0));
    assert(approx_equal(timi::interpolate(axes, values, {0.0, 1.0}), 2.0));
    assert(approx_equal(timi::interpolate(axes, values, {1.0, 1.0}), 3.0));

    // Edge midpoints
    assert(approx_equal(timi::interpolate(axes, values, {0.5, 0.0}), 0.5));
    assert(approx_equal(timi::interpolate(axes, values, {0.5, 1.0}), 2.5));
    assert(approx_equal(timi::interpolate(axes, values, {0.0, 0.5}), 1.0));
    assert(approx_equal(timi::interpolate(axes, values, {1.0, 0.5}), 2.0));

    std::cout << "PASSED\n";
}

void test_2d_stateful() {
    std::cout << "test_2d_stateful... ";

    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 2.0}
    };
    std::vector<double> values = {0.0, 10.0, 20.0, 30.0};

    timi::InterpolatorND<double, double> interp(axes, values);

    // Multiple samples
    assert(approx_equal(interp({0.5, 1.0}), 15.0));
    assert(approx_equal(interp({0.0, 0.0}), 0.0));
    assert(approx_equal(interp({1.0, 2.0}), 30.0));

    std::cout << "PASSED\n";
}

void test_2d_3x3_grid() {
    std::cout << "test_2d_3x3_grid... ";

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
    assert(approx_equal(timi::interpolate(axes, values, {0.5, 0.5}), 1.0));
    assert(approx_equal(timi::interpolate(axes, values, {1.5, 1.5}), 3.0));
    assert(approx_equal(timi::interpolate(axes, values, {0.0, 2.0}), 2.0));

    std::cout << "PASSED\n";
}

// =============================================================================
// 3D Tests
// =============================================================================

void test_3d_basic() {
    std::cout << "test_3d_basic... ";

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
    assert(approx_equal(v, 3.5));

    // Corner points
    assert(approx_equal(timi::interpolate(axes, values, {0.0, 0.0, 0.0}), 0.0));
    assert(approx_equal(timi::interpolate(axes, values, {1.0, 1.0, 1.0}), 7.0));

    std::cout << "PASSED\n";
}

void test_3d_stateful() {
    std::cout << "test_3d_stateful... ";

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

    assert(approx_equal(interp({0.5, 0.5, 0.5}), 3.5));
    assert(approx_equal(interp({0.0, 0.0, 0.0}), 0.0));
    assert(approx_equal(interp({1.0, 1.0, 1.0}), 7.0));

    std::cout << "PASSED\n";
}

// =============================================================================
// Custom Type Tests
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

void test_custom_type_1d() {
    std::cout << "test_custom_type_1d... ";

    std::vector<double> x = {0.0, 1.0};
    std::vector<Vec2> y = {{0.0, 0.0}, {2.0, 4.0}};

    Vec2 result = timi::interpolate(x, y, 0.5);
    assert(approx_equal(result.x, 1.0));
    assert(approx_equal(result.y, 2.0));

    // Using stateful API
    timi::Interpolator1D<double, Vec2> interp(x, y);
    Vec2 result2 = interp(0.25);
    assert(approx_equal(result2.x, 0.5));
    assert(approx_equal(result2.y, 1.0));

    std::cout << "PASSED\n";
}

void test_custom_type_2d() {
    std::cout << "test_custom_type_2d... ";

    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},
        {0.0, 1.0}
    };
    std::vector<Vec2> values = {
        {0.0, 0.0}, {1.0, 0.0},
        {0.0, 1.0}, {1.0, 1.0}
    };

    Vec2 result = timi::interpolate(axes, values, {0.5, 0.5});
    assert(approx_equal(result.x, 0.5));
    assert(approx_equal(result.y, 0.5));

    std::cout << "PASSED\n";
}

// =============================================================================
// API Consistency Tests
// =============================================================================

void test_consistency_1d() {
    std::cout << "test_consistency_1d... ";

    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 10.0, 15.0, 25.0, 40.0};

    timi::Interpolator1D<double, double> interp(x, y);

    // Functional and stateful should give same results
    for (double t = -0.5; t <= 4.5; t += 0.1) {
        double func_result = timi::interpolate(x, y, t);
        double stat_result = interp(t);
        assert(approx_equal(func_result, stat_result, 1e-12));
    }

    std::cout << "PASSED\n";
}

void test_consistency_nd() {
    std::cout << "test_consistency_nd... ";

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
        assert(approx_equal(func_result, stat_result, 1e-12));
    }

    std::cout << "PASSED\n";
}

// =============================================================================
// Edge Case Tests
// =============================================================================

void test_two_point_grid() {
    std::cout << "test_two_point_grid... ";

    // Minimum valid grid: 2 points
    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {10.0, 20.0};

    assert(approx_equal(timi::interpolate(x, y, 0.0), 10.0));
    assert(approx_equal(timi::interpolate(x, y, 0.5), 15.0));
    assert(approx_equal(timi::interpolate(x, y, 1.0), 20.0));

    std::cout << "PASSED\n";
}

void test_float_type() {
    std::cout << "test_float_type... ";

    std::vector<float> x = {0.0f, 1.0f, 2.0f};
    std::vector<float> y = {0.0f, 10.0f, 20.0f};

    float result = timi::interpolate(x, y, 0.5f);
    assert(approx_equal(result, 5.0f, 1e-5f));

    timi::Interpolator1D<float, float> interp(x, y);
    assert(approx_equal(interp(1.5f), 15.0f, 1e-5f));

    std::cout << "PASSED\n";
}

// =============================================================================
// Main
// =============================================================================

int main() {
    std::cout << "Running timi tests...\n\n";

    // 1D tests
    test_1d_basic();
    test_1d_out_of_range();
    test_1d_stateful();
    test_1d_stateful_repeated();
    test_1d_nonuniform_grid();

    // 2D tests
    test_2d_basic();
    test_2d_stateful();
    test_2d_3x3_grid();

    // 3D tests
    test_3d_basic();
    test_3d_stateful();

    // Custom type tests
    test_custom_type_1d();
    test_custom_type_2d();

    // Consistency tests
    test_consistency_1d();
    test_consistency_nd();

    // Edge case tests
    test_two_point_grid();
    test_float_type();

    std::cout << "\nAll tests PASSED!\n";
    return 0;
}
