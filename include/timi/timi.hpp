#ifndef TIMI_HPP
#define TIMI_HPP

#include <vector>
#include <cstddef>
#include <cassert>
#include <algorithm>

namespace timi {

namespace detail {

// Linear interpolation: result = a * (1 - t) + b * t
// Uses only T + T and T * Scalar operations
template<typename Scalar, typename T>
T lerp(const T& a, const T& b, Scalar t) {
    return a * (Scalar(1) - t) + b * t;
}

// Find the interval [i, i+1] where x falls using binary search
// Returns index i such that x_vals[i] <= x < x_vals[i+1]
// Clamps to valid range: returns 0 if x <= x_vals[0], n-2 if x >= x_vals[n-1]
template<typename Scalar>
std::size_t find_interval(const std::vector<Scalar>& x_vals, Scalar x) {
    assert(x_vals.size() >= 2 && "Need at least 2 points for interpolation");

    // Binary search for first element > x
    auto it = std::upper_bound(x_vals.begin(), x_vals.end(), x);

    if (it == x_vals.begin()) {
        // x <= x_vals[0], clamp to first interval
        return 0;
    }
    if (it == x_vals.end()) {
        // x >= x_vals[n-1], clamp to last interval
        return x_vals.size() - 2;
    }

    // upper_bound returns first element > x
    // We want interval [i, i+1] where x_vals[i] <= x < x_vals[i+1]
    return static_cast<std::size_t>(it - x_vals.begin() - 1);
}

// 1D interpolation implementation
template<typename Scalar, typename T>
T interpolate_1d(
    const std::vector<Scalar>& x_vals,
    const std::vector<T>& y_vals,
    Scalar x
) {
    std::size_t i = find_interval(x_vals, x);

    Scalar x0 = x_vals[i];
    Scalar x1 = x_vals[i + 1];

    // Compute interpolation parameter t in [0, 1]
    Scalar t = (x1 != x0) ? (x - x0) / (x1 - x0) : Scalar(0);

    // Clamp t to [0, 1] for robustness
    if (t < Scalar(0)) t = Scalar(0);
    if (t > Scalar(1)) t = Scalar(1);

    return lerp<Scalar, T>(y_vals[i], y_vals[i + 1], t);
}

// Compute strides from axis sizes where first axis varies fastest
// stride[0] = 1, stride[i] = stride[i-1] * sizes[i-1]
// This matches the spec's row-major layout: values = {f(x0,y0), f(x1,y0), f(x0,y1), f(x1,y1)}
inline std::vector<std::size_t> compute_strides(const std::vector<std::size_t>& sizes) {
    std::size_t n = sizes.size();
    std::vector<std::size_t> strides(n);

    if (n == 0) return strides;

    strides[0] = 1;
    for (std::size_t i = 1; i < n; ++i) {
        strides[i] = strides[i - 1] * sizes[i - 1];
    }

    return strides;
}

// Convert N-D indices to flat row-major index
inline std::size_t flat_index(
    const std::size_t* indices,
    const std::size_t* strides,
    std::size_t n_dims
) {
    std::size_t idx = 0;
    for (std::size_t d = 0; d < n_dims; ++d) {
        idx += indices[d] * strides[d];
    }
    return idx;
}

// Maximum supported dimensions for stack allocation
constexpr std::size_t MAX_DIMS = 16;

// N-D multilinear interpolation implementation using iterative blending
template<typename Scalar, typename T>
T interpolate_nd_impl(
    const std::vector<std::vector<Scalar>>& axes,
    const std::vector<T>& values,
    const std::vector<Scalar>& point,
    const std::vector<std::size_t>& strides
) {
    const std::size_t n_dims = axes.size();
    assert(n_dims <= MAX_DIMS && "Too many dimensions");

    // Per-dimension data: interval indices and interpolation parameters
    std::size_t lower_indices[MAX_DIMS];
    Scalar t_values[MAX_DIMS];

    for (std::size_t d = 0; d < n_dims; ++d) {
        std::size_t i = find_interval(axes[d], point[d]);
        lower_indices[d] = i;

        Scalar x0 = axes[d][i];
        Scalar x1 = axes[d][i + 1];
        Scalar t = (x1 != x0) ? (point[d] - x0) / (x1 - x0) : Scalar(0);

        // Clamp to [0, 1]
        if (t < Scalar(0)) t = Scalar(0);
        if (t > Scalar(1)) t = Scalar(1);
        t_values[d] = t;
    }

    // Number of corners = 2^n_dims
    const std::size_t n_corners = std::size_t(1) << n_dims;

    // Gather corner values into a contiguous buffer
    // We use a vector here since T may not be default constructible
    // and we need exactly n_corners elements
    std::vector<T> corners;
    corners.reserve(n_corners);

    std::size_t corner_indices[MAX_DIMS];
    for (std::size_t c = 0; c < n_corners; ++c) {
        // Decode corner: bit d of c indicates lower (0) or upper (1) index
        for (std::size_t d = 0; d < n_dims; ++d) {
            corner_indices[d] = lower_indices[d] + ((c >> d) & 1);
        }
        std::size_t flat_idx = flat_index(corner_indices, strides.data(), n_dims);
        corners.push_back(values[flat_idx]);
    }

    // Iteratively blend along each dimension
    // After blending dimension d, we have 2^(n_dims - d - 1) values
    std::size_t count = n_corners;
    for (std::size_t d = 0; d < n_dims; ++d) {
        Scalar t = t_values[d];
        std::size_t new_count = count / 2;

        for (std::size_t i = 0; i < new_count; ++i) {
            // Blend pair (2*i, 2*i + 1) along dimension d
            corners[i] = lerp<Scalar, T>(corners[2 * i], corners[2 * i + 1], t);
        }
        count = new_count;
    }

    return corners[0];
}

} // namespace detail

// =============================================================================
// Stateful Interpolators
// =============================================================================

/// 1D interpolator that stores sample data for repeated queries
/// Thread-safe: operator() is const and performs only read operations
template<typename Scalar, typename T>
class Interpolator1D {
public:
    Interpolator1D(std::vector<Scalar> x_vals, std::vector<T> y_vals)
        : x_vals_(std::move(x_vals))
        , y_vals_(std::move(y_vals))
    {
        assert(x_vals_.size() >= 2 && "Need at least 2 sample points");
        assert(x_vals_.size() == y_vals_.size() && "x and y sizes must match");
        assert(std::is_sorted(x_vals_.begin(), x_vals_.end()) &&
               "x_vals must be sorted ascending");
    }

    T operator()(Scalar x) const {
        return detail::interpolate_1d(x_vals_, y_vals_, x);
    }

    const std::vector<Scalar>& x_values() const { return x_vals_; }
    const std::vector<T>& y_values() const { return y_vals_; }
    std::size_t size() const { return x_vals_.size(); }

private:
    std::vector<Scalar> x_vals_;
    std::vector<T> y_vals_;
};

/// N-dimensional interpolator that stores grid data for repeated queries
/// Thread-safe: operator() is const and performs only read operations
template<typename Scalar, typename T>
class InterpolatorND {
public:
    InterpolatorND(
        std::vector<std::vector<Scalar>> axes,
        std::vector<T> values
    )
        : axes_(std::move(axes))
        , values_(std::move(values))
    {
        assert(!axes_.empty() && "Need at least one axis");

        std::size_t expected_size = 1;
        for (const auto& axis : axes_) {
            assert(axis.size() >= 2 && "Each axis needs at least 2 points");
            assert(std::is_sorted(axis.begin(), axis.end()) &&
                   "Each axis must be sorted ascending");
            sizes_.push_back(axis.size());
            expected_size *= axis.size();
        }

        assert(values_.size() == expected_size &&
               "Value count must equal product of axis sizes");

        strides_ = detail::compute_strides(sizes_);
    }

    T operator()(const std::vector<Scalar>& point) const {
        assert(point.size() == axes_.size() &&
               "Point dimension must match grid dimension");
        return detail::interpolate_nd_impl(axes_, values_, point, strides_);
    }

    std::size_t num_dimensions() const { return axes_.size(); }
    const std::vector<std::vector<Scalar>>& axes() const { return axes_; }
    const std::vector<T>& values() const { return values_; }

private:
    std::vector<std::vector<Scalar>> axes_;
    std::vector<T> values_;
    std::vector<std::size_t> sizes_;
    std::vector<std::size_t> strides_;
};

// =============================================================================
// Functional API
// =============================================================================

/// 1D linear interpolation
/// @param x_vals Sorted ascending coordinate values (size >= 2)
/// @param y_vals Corresponding sample values (same size as x_vals)
/// @param x Coordinate to interpolate at
/// @return Interpolated value (clamped if x is out of range)
template<typename Scalar, typename T>
T interpolate(
    const std::vector<Scalar>& x_vals,
    const std::vector<T>& y_vals,
    Scalar x
) {
    assert(x_vals.size() >= 2 && "Need at least 2 sample points");
    assert(x_vals.size() == y_vals.size() && "Size mismatch");
    return detail::interpolate_1d(x_vals, y_vals, x);
}

/// N-dimensional multilinear interpolation
/// @param axes Vector of coordinate vectors, one per dimension (each sorted ascending)
/// @param values Flat array of values at grid intersections in row-major order
/// @param point Coordinates to interpolate at (size must match axes.size())
/// @return Interpolated value (clamped if point is out of range on any axis)
template<typename Scalar, typename T>
T interpolate(
    const std::vector<std::vector<Scalar>>& axes,
    const std::vector<T>& values,
    const std::vector<Scalar>& point
) {
    assert(!axes.empty() && "Need at least one axis");
    assert(point.size() == axes.size() && "Point dimension mismatch");

    // Compute sizes and strides
    std::vector<std::size_t> sizes;
    std::size_t expected_size = 1;
    for (const auto& axis : axes) {
        assert(axis.size() >= 2 && "Each axis needs at least 2 points");
        sizes.push_back(axis.size());
        expected_size *= axis.size();
    }
    assert(values.size() == expected_size && "Value count mismatch");

    std::vector<std::size_t> strides = detail::compute_strides(sizes);
    return detail::interpolate_nd_impl(axes, values, point, strides);
}

} // namespace timi

#endif // TIMI_HPP
