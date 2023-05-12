#include "TextureIntegral.h"

#include "helpers.h"

#include <wmtk/quadrature/ClippedQuadrature.h>
#include <wmtk/quadrature/TriangleQuadrature.h>

#include <lagrange/utils/assert.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace wmtk {

namespace {

template <typename T, typename Func>
void sequential_for(T start, T end, Func func)
{
    for (T i = start; i < end; ++i) {
        func(i);
    }
}

struct QuadratureCache
{
    wmtk::Quadrature quad;
    wmtk::Quadrature tmp;
};

// Reference implementation copied over from Displacement.h
template <typename DisplacementFunc>
double get_error_per_triangle_exact(
    const wmtk::Image& m_image,
    const Eigen::Matrix<double, 3, 2, Eigen::RowMajor>& triangle,
    tbb::enumerable_thread_specific<QuadratureCache>& m_cache,
    DisplacementFunc get)
{
    using T = double;
    constexpr int Degree = 4;
    const int order = 2 * (Degree - 1);

    Eigen::Matrix<T, 3, 1> p1 = get(triangle(0, 0), triangle(0, 1));
    Eigen::Matrix<T, 3, 1> p2 = get(triangle(1, 0), triangle(1, 1));
    Eigen::Matrix<T, 3, 1> p3 = get(triangle(2, 0), triangle(2, 1));
    Eigen::AlignedBox2d bbox;
    Eigen::Matrix<double, 3, 2, 1> triangle_double;
    for (auto i = 0; i < 3; i++) {
        Eigen::Vector2d p_double;
        p_double << get_value(triangle(i, 0)), get_value(triangle(i, 1));
        triangle_double.row(i) << p_double(0), p_double(1);
        bbox.extend(p_double);
    }
    auto squared_norm_T = [&](const Eigen::Matrix<T, 3, 1>& row_v) -> T {
        T ret = T(0.);
        for (auto i = 0; i < row_v.rows(); i++) {
            auto debug_rowv = row_v(i, 0);
            ret += pow(row_v(i, 0), 2);
        }
        return ret;
    };
    auto get_coordinate = [&](const double& x, const double& y) -> std::pair<int, int> {
        auto [xx, yy] = m_image.get_pixel_index(get_value(x), get_value(y));
        return {
            m_image.get_coordinate(xx, m_image.get_wrapping_mode_x()),
            m_image.get_coordinate(yy, m_image.get_wrapping_mode_y())};
    };
    auto bbox_min = bbox.min();
    auto bbox_max = bbox.max();
    auto [xx1, yy1] = get_coordinate(bbox_min(0), bbox_min(1));
    auto [xx2, yy2] = get_coordinate(bbox_max(0), bbox_max(1));
    auto num_pixels = std::max(abs(xx2 - xx1), abs(yy2 - yy1)) + 1;
    assert(num_pixels > 0);
    const double pixel_size = bbox.diagonal().maxCoeff() / num_pixels;

    // calculate the barycentric coordinate of the a point using u, v cooridnates
    // returns the 3d coordinate on the current mesh
    auto get_p_tri = [&](const T& u, const T& v) -> Eigen::Matrix<T, 3, 1> {
        /*
        λ1 = ((v2 - v3)(u - u3) + (u3 - u2)(v - v3)) / ((v2 - v3)(u1 - u3) + (u3 - u2)(v1 - v3))
        λ2 = ((v3 - v1)(u - u3) + (u1 - u3)(v - v3)) / ((v2 - v3)(u1 - u3) + (u3 - u2)(v1 - v3))
        λ3 = 1 - λ1 - λ2
        z = λ1 * z1 + λ2 * z2 + λ3 * z3
        */
        auto u1 = triangle(0, 0);
        auto v1 = triangle(0, 1);
        auto u2 = triangle(1, 0);
        auto v2 = triangle(1, 1);
        auto u3 = triangle(2, 0);
        auto v3 = triangle(2, 1);

        auto lambda1 = ((v2 - v3) * (u - u3) + (u3 - u2) * (v - v3)) /
                       ((v2 - v3) * (u1 - u3) + (u3 - u2) * (v1 - v3));
        auto lambda2 = ((v3 - v1) * (u - u3) + (u1 - u3) * (v - v3)) /
                       ((v2 - v3) * (u1 - u3) + (u3 - u2) * (v1 - v3));
        auto lambda3 = 1 - lambda1 - lambda2;
        Eigen::Matrix<T, 3, 1> p_tri = (lambda1 * p1 + lambda2 * p2 + lambda3 * p3);
        return p_tri;
    };

    auto check_degenerate = [&]() -> bool {
        auto u1 = triangle(0, 0);
        auto v1 = triangle(0, 1);
        auto u2 = triangle(1, 0);
        auto v2 = triangle(1, 1);
        auto u3 = triangle(2, 0);
        auto v3 = triangle(2, 1);
        if ((v2 - v3) * (u1 - u3) + (u3 - u2) * (v1 - v3) == 0.) return false;
        if ((v2 - v3) * (u1 - u3) + (u3 - u2) * (v1 - v3) == 0.) return false;
        return true;
    };

    T value = T(0.);
    for (auto x = 0; x < num_pixels; ++x) {
        for (auto y = 0; y < num_pixels; ++y) {
            Eigen::AlignedBox2d box;
            box.extend(bbox.min() + Eigen::Vector2d(x * pixel_size, y * pixel_size));
            box.extend(bbox.min() + Eigen::Vector2d((x + 1) * pixel_size, (y + 1) * pixel_size));
            wmtk::ClippedQuadrature rules;
            rules.clipped_triangle_box_quadrature(
                order,
                triangle_double,
                box,
                m_cache.local().quad,
                &m_cache.local().tmp);
            for (auto i = 0; i < m_cache.local().quad.size(); ++i) {
                auto tmpu = T(m_cache.local().quad.points()(i, 0));
                auto tmpv = T(m_cache.local().quad.points()(i, 1));
                if (!check_degenerate())
                    continue;
                else {
                    Eigen::Matrix<T, 3, 1> tmpp_displaced = get(tmpu, tmpv);
                    Eigen::Matrix<T, 3, 1> tmpp_tri = get_p_tri(tmpu, tmpv);
                    Eigen::Matrix<T, 3, 1> diffp = tmpp_displaced - tmpp_tri;
                    value += squared_norm_T(diffp) * T(m_cache.local().quad.weights()[i]);
                }
            }
        }
    }
    return value;
}

// bool point_in_triangle(
//     const Eigen::Matrix<double, 3, 2, Eigen::RowMajor>& triangle,
//     const Eigen::Vector2d& point)
// {
//     Eigen::Vector2d v0 = triangle.row(0);
//     Eigen::Vector2d v1 = triangle.row(1);
//     Eigen::Vector2d v2 = triangle.row(2);
//     Eigen::Vector2d v0v1 = v1 - v0;
//     Eigen::Vector2d v0v2 = v2 - v0;
//     Eigen::Vector2d v0p = point - v0;

//     double dot00 = v0v2.dot(v0v2);
//     double dot01 = v0v2.dot(v0v1);
//     double dot02 = v0v2.dot(v0p);
//     double dot11 = v0v1.dot(v0v1);
//     double dot12 = v0v1.dot(v0p);

//     double inv_denom = 1 / (dot00 * dot11 - dot01 * dot01);
//     double u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
//     double v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

//     return (u >= 0) && (v >= 0) && (u + v < 1);
// }

bool point_in_triangle(
    const Eigen::Matrix<double, 3, 2, Eigen::RowMajor>& triangle,
    const Eigen::Vector2d& point)
{
    const auto& a = triangle.row(0);
    const auto& b = triangle.row(1);
    const auto& c = triangle.row(2);

    const double as_x = point.x() - a.x();
    const double as_y = point.y() - a.y();

    bool s_ab = (b.x() - a.x()) * as_y - (b.y() - a.y()) * as_x > 0;

    if ((c.x() - a.x()) * as_y - (c.y() - a.y()) * as_x > 0 == s_ab) {
        return false;
    }
    if ((c.x() - b.x()) * (point.y() - b.y()) - (c.y() - b.y()) * (point.x() - b.x()) > 0 != s_ab) {
        return false;
    }
    return true;
}

// Optimized implementation that switches between nearest and bilinear interpolation
template <typename DisplacementFunc>
double get_error_per_triangle_adaptive(
    const std::array<wmtk::Image, 3>& images,
    const Eigen::Matrix<double, 3, 2, Eigen::RowMajor>& triangle_uv,
    tbb::enumerable_thread_specific<QuadratureCache>& m_cache,
    DisplacementFunc get)
{
    using T = double;
    const int order = 1;
    // constexpr int Degree = 2;
    // const int order = 2 * (Degree - 1);

    Eigen::Matrix3d triangle_3d;
    triangle_3d.row(0) = get(triangle_uv(0, 0), triangle_uv(0, 1));
    triangle_3d.row(1) = get(triangle_uv(1, 0), triangle_uv(1, 1));
    triangle_3d.row(2) = get(triangle_uv(2, 0), triangle_uv(2, 1));
    Eigen::AlignedBox2d bbox_uv;
    for (const auto& p : triangle_uv.rowwise()) {
        bbox_uv.extend(p.transpose());
    }

    const int w = images[0].width();
    const int h = images[0].height();
    auto get_coordinate = [&](double u, double v) -> Eigen::Vector2i {
        auto x = u * static_cast<float>(w);
        auto y = v * static_cast<float>(h);
        const auto sx = std::clamp(static_cast<int>(x), 0, w - 1);
        const auto sy = std::clamp(static_cast<int>(y), 0, h - 1);
        return {sx, sy};
    };

    const auto min_pixel = get_coordinate(bbox_uv.min()(0), bbox_uv.min()(1));
    const auto max_pixel = get_coordinate(bbox_uv.max()(0), bbox_uv.max()(1));
    const Eigen::Vector2d pixel_size(1.0 / w, 1.0 / h);

    const double u1 = triangle_uv(0, 0);
    const double v1 = triangle_uv(0, 1);
    const double u2 = triangle_uv(1, 0);
    const double v2 = triangle_uv(1, 1);
    const double u3 = triangle_uv(2, 0);
    const double v3 = triangle_uv(2, 1);
    const double denom = ((v2 - v3) * (u1 - u3) + (u3 - u2) * (v1 - v3));
    if (denom < std::numeric_limits<double>::denorm_min()) {
        // Degenerate triangle
        return 0.;
    }

    // calculate the barycentric coordinate of the a point using u, v coordinates
    // returns the 3d coordinate on the current mesh
    auto get_p_interpolated = [&](double u, double v) -> Eigen::Matrix<T, 3, 1> {
        /*
        λ1 = ((v2 - v3)(u - u3) + (u3 - u2)(v - v3)) / ((v2 - v3)(u1 - u3) + (u3 - u2)(v1 - v3))
        λ2 = ((v3 - v1)(u - u3) + (u1 - u3)(v - v3)) / ((v2 - v3)(u1 - u3) + (u3 - u2)(v1 - v3))
        λ3 = 1 - λ1 - λ2
        z = λ1 * z1 + λ2 * z2 + λ3 * z3
        */
        auto lambda1 = ((v2 - v3) * (u - u3) + (u3 - u2) * (v - v3)) / denom;
        auto lambda2 = ((v3 - v1) * (u - u3) + (u1 - u3) * (v - v3)) / denom;
        auto lambda3 = 1 - lambda1 - lambda2;
        return lambda1 * triangle_3d.row(0) + lambda2 * triangle_3d.row(1) +
               lambda3 * triangle_3d.row(2);
    };

    T value = T(0.);
    size_t num_inside = 0;
    size_t num_total = 0;
    size_t num_boundary = 0;
    for (int x = min_pixel.x(); x <= max_pixel.x(); ++x) {
        for (int y = min_pixel.y(); y <= max_pixel.y(); ++y) {
            ++num_total;
            Eigen::Vector2i pixel_coord(x, y);
            Eigen::AlignedBox2d box;
            box.extend(pixel_coord.cast<double>().cwiseProduct(pixel_size));
            box.extend(
                (pixel_coord + Eigen::Vector2i::Ones()).cast<double>().cwiseProduct(pixel_size));
            bool all_inside = true;
            for (int k = 0; k < 4; ++k) {
                bool inside = point_in_triangle(
                    triangle_uv,
                    box.corner(static_cast<Eigen::AlignedBox2d::CornerType>(k)));
                if (!inside) {
                    all_inside = false;
                    break;
                }
            }
            if (all_inside) {
                ++num_inside;
                Eigen::Matrix<T, 3, 1> p_tri =
                    get_p_interpolated(box.center().x(), box.center().y());
                Eigen::Matrix<T, 3, 1> p_displaced;
                for (size_t k = 0; k < 3; ++k) {
                    p_displaced[k] = images[k].get_raw_image()(x, y);
                }
                value += (p_displaced - p_tri).squaredNorm() * pixel_size(0) * pixel_size(1);
                if (0) {
                    wmtk::ClippedQuadrature rules;
                    auto& quadr = m_cache.local().quad;
                    rules.clipped_triangle_box_quadrature(
                        1,
                        triangle_uv,
                        box,
                        quadr,
                        &m_cache.local().tmp);
                    la_runtime_assert(quadr.size() == 2);
                    for (size_t i = 0; i < quadr.size(); ++i) {
                        double u = quadr.points()(i, 0);
                        double v = quadr.points()(i, 1);
                        Eigen::Matrix<T, 3, 1> p_displaced2 =
                            internal::sample_nearest<3>(images, u, v).cast<T>();
                        Eigen::Matrix<T, 3, 1> p_tri2 = get_p_interpolated(u, v);
                        la_runtime_assert(p_displaced2 == p_displaced);
                    }
                    double w = quadr.weights().sum();
                    double diff = std::abs(quadr.weights().sum() - pixel_size(0) * pixel_size(1));
                    la_runtime_assert(diff < 1e-10);
                }
            } else {
                wmtk::ClippedQuadrature rules;
                auto& quadr = m_cache.local().quad;
                rules.clipped_triangle_box_quadrature(
                    order,
                    triangle_uv,
                    box,
                    quadr,
                    &m_cache.local().tmp);
                if (quadr.size() > 0) {
                    ++num_boundary;
                }
                for (size_t i = 0; i < quadr.size(); ++i) {
                    double u = quadr.points()(i, 0);
                    double v = quadr.points()(i, 1);
                    Eigen::Matrix<T, 3, 1> p_displaced = get(u, v);
                    Eigen::Matrix<T, 3, 1> p_tri = get_p_interpolated(u, v);
                    value += (p_displaced - p_tri).squaredNorm() * quadr.weights()[i];
                }
            }
        }
    }
    // logger().info("Num inside: {}, boundary: {}, total: {}", num_inside, num_boundary,
    // num_total);
    return value;
}

} // namespace

struct TextureIntegral::Cache
{
    // Data for exact error computation
    tbb::enumerable_thread_specific<QuadratureCache> quadrature_cache;
};

TextureIntegral::TextureIntegral(std::array<wmtk::Image, 3> data)
    : m_data(std::move(data))
    , m_cache(lagrange::make_value_ptr<Cache>())
{}

TextureIntegral::~TextureIntegral() = default;

template <TextureIntegral::SamplingMethod Sampling, TextureIntegral::IntegrationMethod Integration>
void TextureIntegral::get_error_per_triangle_internal(
    lagrange::span<const std::array<float, 6>> input_triangles,
    lagrange::span<float> output_errors)
{
    assert(input_triangles.size() == output_errors.size());
    tbb::parallel_for(size_t(0), input_triangles.size(), [&](size_t i) {
        Eigen::Matrix<double, 3, 2, Eigen::RowMajor> triangle;
        triangle.row(0) << input_triangles[i][0], input_triangles[i][1];
        triangle.row(1) << input_triangles[i][2], input_triangles[i][3];
        triangle.row(2) << input_triangles[i][4], input_triangles[i][5];
        auto sampling_func = [&](double u, double v) -> Eigen::Matrix<double, 3, 1> {
            if constexpr (Sampling == SamplingMethod::Bicubic) {
                return internal::sample_bicubic(m_data, u, v).cast<double>();
            } else if constexpr (Sampling == SamplingMethod::Nearest) {
                return internal::sample_nearest(m_data, u, v).cast<double>();
            } else if constexpr (Sampling == SamplingMethod::Bilinear) {
                return internal::sample_bilinear(m_data, u, v).cast<double>();
            }
        };
        if constexpr (Integration == IntegrationMethod::Exact)
            output_errors[i] = get_error_per_triangle_exact(
                m_data[0],
                triangle,
                m_cache->quadrature_cache,
                sampling_func);
        else if constexpr (Integration == IntegrationMethod::Adaptive) {
            output_errors[i] = get_error_per_triangle_adaptive(
                m_data,
                triangle,
                m_cache->quadrature_cache,
                sampling_func);
        }
    });
}

void TextureIntegral::get_error_per_triangle(
    lagrange::span<const std::array<float, 6>> input_triangles,
    lagrange::span<float> output_errors)
{
    assert(input_triangles.size() == output_errors.size());
    switch (m_sampling_method) {
    case SamplingMethod::Bicubic:
        if (m_integration_method == IntegrationMethod::Exact) {
            get_error_per_triangle_internal<SamplingMethod::Bicubic, IntegrationMethod::Exact>(
                input_triangles,
                output_errors);

        } else {
            get_error_per_triangle_internal<SamplingMethod::Bicubic, IntegrationMethod::Adaptive>(
                input_triangles,
                output_errors);
        }
        break;
    case SamplingMethod::Nearest:
        if (m_integration_method == IntegrationMethod::Exact) {
            get_error_per_triangle_internal<SamplingMethod::Nearest, IntegrationMethod::Exact>(
                input_triangles,
                output_errors);

        } else {
            get_error_per_triangle_internal<SamplingMethod::Nearest, IntegrationMethod::Adaptive>(
                input_triangles,
                output_errors);
        }
        break;
    case SamplingMethod::Bilinear:
        if (m_integration_method == IntegrationMethod::Exact) {
            get_error_per_triangle_internal<SamplingMethod::Bilinear, IntegrationMethod::Exact>(
                input_triangles,
                output_errors);

        } else {
            get_error_per_triangle_internal<SamplingMethod::Bilinear, IntegrationMethod::Adaptive>(
                input_triangles,
                output_errors);
        }
        break;
    }
}

void TextureIntegral::get_integral_per_triangle(
    lagrange::span<const std::array<float, 6>> input_triangles,
    lagrange::span<std::array<float, 3>> output_integrals)
{}

} // namespace wmtk
