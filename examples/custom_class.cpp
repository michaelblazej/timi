#include <timi/timi.hpp>
#include <iostream>
#include <iomanip>

// =============================================================================
// Custom Class: RGB Color
// =============================================================================

/// RGB color with components in [0, 1] range
struct Color {
    double r, g, b;

    // Helper to print as RGB values (0-255)
    void print() const {
        std::cout << "rgb("
                  << std::setw(3) << static_cast<int>(r * 255) << ", "
                  << std::setw(3) << static_cast<int>(g * 255) << ", "
                  << std::setw(3) << static_cast<int>(b * 255) << ")";
    }
};

// Required operator for timi: addition
Color operator+(const Color& a, const Color& b) {
    return {a.r + b.r, a.g + b.g, a.b + b.b};
}

// Required operator for timi: scalar multiplication
Color operator*(const Color& c, double s) {
    return {c.r * s, c.g * s, c.b * s};
}

// =============================================================================
// Main
// =============================================================================

int main() {
    std::cout << "=== timi: Custom Class (RGB Color) Example ===\n\n";

    // -------------------------------------------------------------------------
    // 1D Color Gradient: Temperature scale (cold to hot)
    // -------------------------------------------------------------------------
    std::cout << "--- 1D Color Gradient: Temperature Scale ---\n\n";

    // Sample points: blue (cold) -> cyan -> green -> yellow -> red (hot)
    std::vector<double> temp = {0.0, 25.0, 50.0, 75.0, 100.0};
    std::vector<Color> colors = {
        {0.0, 0.0, 1.0},  // Blue   (0%)
        {0.0, 1.0, 1.0},  // Cyan   (25%)
        {0.0, 1.0, 0.0},  // Green  (50%)
        {1.0, 1.0, 0.0},  // Yellow (75%)
        {1.0, 0.0, 0.0}   // Red    (100%)
    };

    // Create stateful interpolator
    timi::Interpolator1D<double, Color> temp_gradient(temp, colors);

    std::cout << "Temperature gradient (0-100):\n";
    for (double t = 0.0; t <= 100.0; t += 10.0) {
        Color c = temp_gradient(t);
        std::cout << "  " << std::setw(3) << static_cast<int>(t) << "%: ";
        c.print();
        std::cout << "\n";
    }

    // -------------------------------------------------------------------------
    // 2D Color Field: Interpolate across a surface
    // -------------------------------------------------------------------------
    std::cout << "\n--- 2D Color Field ---\n\n";

    // 2x2 grid with corners: Red, Green, Blue, White
    std::vector<std::vector<double>> axes = {
        {0.0, 1.0},  // x-axis
        {0.0, 1.0}   // y-axis
    };

    // Corner colors in row-major order:
    // (0,0)=Red, (1,0)=Green, (0,1)=Blue, (1,1)=White
    std::vector<Color> field = {
        {1.0, 0.0, 0.0},  // (0,0) Red
        {0.0, 1.0, 0.0},  // (1,0) Green
        {0.0, 0.0, 1.0},  // (0,1) Blue
        {1.0, 1.0, 1.0}   // (1,1) White
    };

    timi::InterpolatorND<double, Color> color_field(axes, field);

    std::cout << "Color field (x=0..1, y=0..1):\n\n";
    std::cout << "         x=0.0       x=0.5       x=1.0\n";

    for (double y = 1.0; y >= 0.0; y -= 0.5) {
        std::cout << "y=" << std::fixed << std::setprecision(1) << y << "  ";
        for (double x = 0.0; x <= 1.0; x += 0.5) {
            Color c = color_field({x, y});
            c.print();
            std::cout << "  ";
        }
        std::cout << "\n";
    }

    // -------------------------------------------------------------------------
    // Functional API example
    // -------------------------------------------------------------------------
    std::cout << "\n--- Functional API Example ---\n\n";

    // One-shot interpolation without creating an interpolator object
    std::vector<double> x_vals = {0.0, 1.0};
    std::vector<Color> y_vals = {
        {0.0, 0.0, 0.0},  // Black
        {1.0, 1.0, 1.0}   // White
    };

    std::cout << "Grayscale gradient (black to white):\n";
    for (double x = 0.0; x <= 1.0; x += 0.25) {
        Color c = timi::interpolate(x_vals, y_vals, x);
        std::cout << "  " << std::setw(4) << static_cast<int>(x * 100) << "%: ";
        c.print();
        std::cout << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
