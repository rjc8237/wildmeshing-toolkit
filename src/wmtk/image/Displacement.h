#pragma once

#include <lagrange/utils/fpe.h>
#include <tbb/enumerable_thread_specific.h>
#include <wmtk/image/Sampling.h>
#include <wmtk/quadrature/ClippedQuadrature.h>
#include <wmtk/quadrature/TriangleQuadrature.h>
#include <wmtk/utils/PolygonClipping.h>
#include <wmtk/utils/autodiff.h>
#include <iostream>
#include <type_traits>
#include <wmtk/utils/LineQuadrature.hpp>
#include "Image.h"
namespace wmtk {

enum class DISPLACEMENT_MODE { MESH_3D = 0, PLANE = 1, VECTOR = 2 };

class Displacement
{
public:
    using DScalar = DScalar2<double, Eigen::Vector2d, Eigen::Matrix2d>;
    virtual ~Displacement() = default;

public:
    virtual Eigen::Matrix<double, 3, 1> get(double x, double y) const = 0;
    virtual Eigen::Matrix<DScalar, 3, 1> get(const DScalar& x, const DScalar& y) const = 0;
    virtual double get_error_per_edge(
        const Eigen::Matrix<double, 2, 1>& p1,
        const Eigen::Matrix<double, 2, 1>& p2) const = 0;
    virtual DScalar get_error_per_edge(
        const Eigen::Matrix<DScalar, 2, 1>& p1,
        const Eigen::Matrix<DScalar, 2, 1>& p2) const = 0;
    virtual double get_error_per_triangle(
        const Eigen::Matrix<double, 3, 2, Eigen::RowMajor>& triangle) = 0;
    virtual DScalar get_error_per_triangle(
        const Eigen::Matrix<DScalar, 3, 2, Eigen::RowMajor>& triangle) = 0;
};

template <typename Derived>
class DisplacementImage : public Displacement
{
protected:
    struct QuadrCache
    {
        wmtk::Quadrature quad;
        wmtk::Quadrature tmp;
    };
    tbb::enumerable_thread_specific<QuadrCache> m_cache;

public:
    DisplacementImage() = default;
    virtual ~DisplacementImage() = default;

    // prevent DisplacementImage being accidentally moved or copied
    DisplacementImage(DisplacementImage&&) = delete;
    DisplacementImage& operator=(DisplacementImage&&) = delete;
    DisplacementImage(const DisplacementImage&) = delete;
    DisplacementImage& operator=(const DisplacementImage&) = delete;

public:
    Eigen::Matrix<double, 3, 1> get(double u, double v) const override
    {
        return static_cast<const Derived*>(this)->get(u, v);
    }

    Eigen::Matrix<DScalar, 3, 1> get(const DScalar& u, const DScalar& v) const override
    {
        return static_cast<const Derived*>(this)->get(u, v);
    }

protected:
    virtual std::pair<int, int> get_coordinate(double x, double y) const = 0;

public:
    template <class T>
    inline T triangle_3d_area(
        const Eigen::Matrix<T, 3, 1>& a,
        const Eigen::Matrix<T, 3, 1>& b,
        const Eigen::Matrix<T, 3, 1>& c)
    {
        const T n0 = (a(1) - b(1)) * (a(2) - c(2)) - (a(1) - c(1)) * (a(2) - b(2));
        const T n1 = -(a(0) - b(0)) * (a(2) - c(2)) + (a(0) - c(0)) * (a(2) - b(2));
        const T n2 = (a(0) - b(0)) * (a(1) - c(1)) - (a(0) - c(0)) * (a(1) - b(1));

        return sqrt(n0 * n0 + n1 * n1 + n2 * n2) * static_cast<T>(0.5);
    };

    // template get 3d tri area
    template <class T>
    inline T triangle_2d_area(
        const Eigen::Matrix<T, 2, 1>& A,
        const Eigen::Matrix<T, 2, 1>& B,
        const Eigen::Matrix<T, 2, 1>& C)
    {
        auto B_A = B - A;
        auto C_A = C - A;
        T area = static_cast<T>(0.5) * abs(B_A.x() * C_A.y() - B_A.y() * C_A.x());
        return area;
    };
    inline double get_error_per_edge(
        const Eigen::Matrix<double, 2, 1>& p1,
        const Eigen::Matrix<double, 2, 1>& p2) const override
    {
        return get_error_per_edge_T<double, 5>(p1, p2);
    }

    inline DScalar get_error_per_edge(
        const Eigen::Matrix<DScalar, 2, 1>& p1,
        const Eigen::Matrix<DScalar, 2, 1>& p2) const override
    {
        return get_error_per_edge_T<DScalar, 5>(p1, p2);
    }

    template <class T, int order>
    inline T get_error_per_edge_T(
        const Eigen::Matrix<T, 2, 1>& uv1,
        const Eigen::Matrix<T, 2, 1>& uv2) const
    {
        throw std::runtime_error("get_error_per_edge should not be called");
        auto p1_displaced = get(uv1(0), uv1(1));
        auto p2_displaced = get(uv2(0), uv2(1));
        // get the pixel index of p1 and p2
        auto [xx1, yy1] = get_coordinate(get_value(uv1(0)), get_value(uv1(1)));
        auto [xx2, yy2] = get_coordinate(get_value(uv2(0)), get_value(uv2(1)));
        // get all the pixels in between p1 and p2
        auto pixel_num = std::max(abs(xx2 - xx1), abs(yy2 - yy1));
        if (pixel_num <= 0) return T(0.);
        assert(pixel_num > 0);
        T error = T(0.0);
        auto norm_T = [&](const Eigen::Matrix<T, Eigen::Dynamic, 1>& row_v) -> T {
            T ret = T(0.);
            for (auto i = 0; i < row_v.rows(); i++) {
                ret += pow(row_v(i, 0), 2);
            }
            return sqrt(ret);
        };
        for (int i = 0; i < pixel_num; i++) {
            const auto r0 = T(i) / pixel_num;
            const auto r1 = T(i + 1) / pixel_num;
            Eigen::Matrix<T, 2, 1> pixel_uv1 = uv1 * (1. - r0) + uv2 * r0;
            Eigen::Matrix<T, 2, 1> pixel_uv2 = uv1 * (1. - r1) + uv2 * r1;
            Eigen::Matrix<T, 3, 1> pixel_p1 = p1_displaced * (1. - r0) + p2_displaced * r0;
            Eigen::Matrix<T, 3, 1> pixel_p2 = p1_displaced * (1. - r1) + p2_displaced * r1;
            wmtk::LineQuadrature quad;
            quad.get_quadrature(order);

            T unweighted_error = T(0.0);
            for (int j = 0; j < quad.points.rows(); j++) {
                Eigen::Matrix<T, 2, 1> tmpuv =
                    T(1 - quad.points(j, 0)) * pixel_uv1 + T(quad.points(j, 0)) * pixel_uv2;
                Eigen::Matrix<T, 3, 1> tmpp_displaced = get(tmpuv(0), tmpuv(1));

                Eigen::Matrix<T, 3, 1> tmpp_tri =
                    T(1 - quad.points(j, 0)) * pixel_p1 + T(quad.points(j, 0)) * pixel_p2;
                Eigen::Matrix<T, 3, 1> diffp = tmpp_tri - tmpp_displaced;
                unweighted_error += T(quad.weights(j)) * norm_T(diffp);
            }
            auto diffuv = pixel_uv1 - pixel_uv2;
            error += unweighted_error * norm_T(diffuv);
        }
        assert(error >= 0);
        return error;
    }

    inline double get_error_per_triangle(
        const Eigen::Matrix<double, 3, 2, Eigen::RowMajor>& triangle) override
    {
        return get_error_per_triangle_T(triangle);
    }

    inline DScalar get_error_per_triangle(
        const Eigen::Matrix<DScalar, 3, 2, Eigen::RowMajor>& triangle) override
    {
        return get_error_per_triangle_T(triangle);
    }

    /**
     * @brief Get the error per triangle T object
     *
     * @tparam T
     * @param triangle 3 x 2 matrix of triangle 3 vertices 2d uv coordinates
     * @return T
     */
    template <class T>
    inline T get_error_per_triangle_T(const Eigen::Matrix<T, 3, 2, Eigen::RowMajor>& triangle)
    {
        constexpr int Degree = 4;
        const int order = 2 * (Degree - 1);

        lagrange::enable_fpe();

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
            for (int i = 0; i < row_v.rows(); i++) {
                ret += pow(row_v(i, 0), 2);
            }
            return ret;
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
            if (triangle_2d_area<T>(
                    triangle.row(0).transpose(),
                    triangle.row(1).transpose(),
                    triangle.row(2).transpose()) < 1e-10)
                return false;
            return true;
        };

        T value = T(0.);
        for (auto y = 0; y < num_pixels; ++y) {
            for (auto x = 0; x < num_pixels; ++x) {
                Eigen::AlignedBox2d box;
                box.extend(bbox.min() + Eigen::Vector2d(x * pixel_size, y * pixel_size));
                box.extend(
                    bbox.min() + Eigen::Vector2d((x + 1) * pixel_size, (y + 1) * pixel_size));
                wmtk::ClippedQuadrature rules;
                rules.clipped_triangle_box_quadrature(
                    order,
                    triangle_double,
                    box,
                    m_cache.local().quad,
                    &m_cache.local().tmp);
                for (size_t i = 0; i < m_cache.local().quad.size(); ++i) {
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
        value = value * triangle_3d_area<T>(p1, p2, p3);
        value = value / triangle_2d_area<T>(
                            triangle.row(0).transpose(),
                            triangle.row(1).transpose(),
                            triangle.row(2).transpose());

        return value;
    }
};

template <typename Derived>
class DisplacementHeight : public DisplacementImage<Derived>
{
protected:
    wmtk::Image m_image; // default, used for height map displaced value
    std::unique_ptr<wmtk::Sampling> m_sampler;

public:
    DisplacementHeight(wmtk::Image img, const SAMPLING_MODE sampling_mode)
        : m_image(std::move(img))
    {
        assert(m_image.width() == m_image.height());
        set_sampling_mode(sampling_mode);
    }

    ~DisplacementHeight() = default;

    // prevent DisplacementHeight being accidentally moved or copied
    DisplacementHeight(DisplacementHeight&&) = delete;
    DisplacementHeight& operator=(DisplacementHeight&&) = delete;
    DisplacementHeight(const DisplacementHeight&) = delete;
    DisplacementHeight& operator=(const DisplacementHeight&) = delete;

protected:
    virtual void set_sampling_mode(const SAMPLING_MODE sampling_mode)
    {
        m_sampler = create_sampler(m_image, sampling_mode);
    }

    virtual std::pair<int, int> get_coordinate(double x, double y) const override
    {
        auto [xx, yy] = m_image.get_pixel_index(get_value(x), get_value(y));
        return {m_image.get_coordinate(xx, m_image.get_wrapping_mode_x()),
                m_image.get_coordinate(yy, m_image.get_wrapping_mode_y())};
    }
};

class DisplacementMesh : public DisplacementHeight<DisplacementMesh>
{
public:
    using DScalar = DScalar2<double, Eigen::Vector2d, Eigen::Matrix2d>;
    using Super = DisplacementHeight<DisplacementMesh>;

public:
    DisplacementMesh(
        wmtk::Image img,
        std::array<wmtk::Image, 6> position_normal_images,
        const SAMPLING_MODE sampling_mode,
        double scale,
        Eigen::Matrix<double, 3, 1> offset)
        : DisplacementHeight(std::move(img), sampling_mode)
        , m_normalization_scale(scale)
        , m_normalization_offset(offset)
    {
        set_position_normal_images(position_normal_images);
        set_sampling_mode(sampling_mode);
        assert(m_image.width() != 0);
        assert(m_position_image[0].width() != 0);
        assert(m_normal_image[0].width() != 0);
    };

    ~DisplacementMesh() = default;

protected:
    std::array<wmtk::Image, 3> m_position_image;
    std::array<std::unique_ptr<wmtk::Sampling>, 3> m_position_sampler;
    std::array<wmtk::Image, 3> m_normal_image;
    std::array<std::unique_ptr<wmtk::Sampling>, 3> m_normal_sampler;
    double m_normalization_scale = 1.0;
    Eigen::Vector3d m_normalization_offset = Eigen::Vector3d::Zero();

protected:
    virtual void set_sampling_mode(const wmtk::SAMPLING_MODE sampling_mode) override
    {
        Super::set_sampling_mode(sampling_mode);
        for (auto i = 0; i < 3; i++) {
            m_position_sampler[i] = create_sampler(m_position_image[i], sampling_mode);
            m_normal_sampler[i] = create_sampler(m_normal_image[i], sampling_mode);
        }
    }

public:
    void set_position_normal_images(std::array<wmtk::Image, 6> images)
    {
        // first 3 of the input are position maps
        // last 3 of the input are normal maps
        for (auto i = 0; i < 3; i++) {
            m_position_image[i] = images[i];
            m_normal_image[i] = images[3 + i];
        }
    }
    Eigen::Matrix<double, 3, 1> get(double u, double v) const override
    {
        double z = m_sampler->sample(u, v);
        Eigen::Matrix<double, 3, 1> displace_3d;
        for (auto i = 0; i < 3; i++) {
            double p = m_position_sampler[i]->sample(u, v);

            double d = 2 * m_normal_sampler[i]->sample(u, v) - 1;

            displace_3d(i, 0) = p * m_normalization_scale - m_normalization_offset(i, 0) + z * d;
        }
        return displace_3d;
    }

    Eigen::Matrix<DScalar, 3, 1> get(const DScalar& u, const DScalar& v) const override
    {
        DScalar z = m_sampler->sample(u, v);
        Eigen::Matrix<DScalar, 3, 1> displace_3d;
        for (auto i = 0; i < 3; i++) {
            DScalar p = m_position_sampler[i]->sample(u, v);

            DScalar d = 2 * m_normal_sampler[i]->sample(u, v) - 1;

            displace_3d(i, 0) = p * m_normalization_scale - m_normalization_offset(i, 0) + z * d;
        }
        return displace_3d;
    }
};

class DisplacementVector : public DisplacementImage<DisplacementMesh>
{
public:
    using DScalar = DScalar2<double, Eigen::Vector2d, Eigen::Matrix2d>;
    using Super = DisplacementImage<DisplacementVector>;

public:
    DisplacementVector(
        std::array<wmtk::Image, 3> displaced_positions,
        const SAMPLING_MODE sampling_mode,
        double scale,
        Eigen::Matrix<double, 3, 1> offset)
        : m_displaced_positions_image(displaced_positions)
        , m_normalization_scale(scale)
        , m_normalization_offset(offset)
    {
        set_sampling_mode(sampling_mode);
        assert(m_displaced_positions_image[0].width() != 0);
    };

    ~DisplacementVector() = default;

protected:
    std::array<wmtk::Image, 3> m_displaced_positions_image;
    std::array<std::unique_ptr<wmtk::Sampling>, 3> m_displaced_positions_sampler;
    double m_normalization_scale = 1.0;
    Eigen::Vector3d m_normalization_offset = Eigen::Vector3d::Zero();

protected:
    virtual std::pair<int, int> get_coordinate(double x, double y) const override
    {
        const auto& image = m_displaced_positions_image[0];
        auto [xx, yy] = image.get_pixel_index(get_value(x), get_value(y));
        return {image.get_coordinate(xx, image.get_wrapping_mode_x()),
                image.get_coordinate(yy, image.get_wrapping_mode_y())};
    }

    void set_sampling_mode(const wmtk::SAMPLING_MODE sampling_mode)
    {
        for (auto i = 0; i < 3; i++) {
            m_displaced_positions_sampler[i] =
                create_sampler(m_displaced_positions_image[i], sampling_mode);
        }
    }

public:
    Eigen::Matrix<double, 3, 1> get(double u, double v) const override
    {
        Eigen::Matrix<double, 3, 1> displace_3d;
        for (auto i = 0; i < 3; i++) {
            double p = m_displaced_positions_sampler[i]->sample(u, v);
            displace_3d[i] = p * m_normalization_scale - m_normalization_offset[i];
        }
        return displace_3d;
    }

    Eigen::Matrix<DScalar, 3, 1> get(const DScalar& u, const DScalar& v) const override
    {
        Eigen::Matrix<DScalar, 3, 1> displace_3d;
        for (auto i = 0; i < 3; i++) {
            DScalar p = m_displaced_positions_sampler[i]->sample(u, v);
            displace_3d[i] = p * m_normalization_scale - m_normalization_offset[i];
        }
        return displace_3d;
    }
};

class DisplacementPlane : public DisplacementHeight<DisplacementPlane>
{
    using DScalar = DScalar2<double, Eigen::Vector2d, Eigen::Matrix2d>;

public:
    DisplacementPlane(const wmtk::Image img, const SAMPLING_MODE sampling_mode)
        : DisplacementHeight(img, sampling_mode){};

public:
    Eigen::Matrix<double, 3, 1> get(double u, double v) const override
    {
        auto z = m_sampler->sample(u, v);
        Eigen::Matrix<double, 3, 1> coord_3d;
        coord_3d << u, v, z;
        return coord_3d;
    }

    Eigen::Matrix<DScalar, 3, 1> get(const DScalar& u, const DScalar& v) const override
    {
        auto z = m_sampler->sample(u, v);
        Eigen::Matrix<DScalar, 3, 1> coord_3d;
        coord_3d << u, v, z;
        return coord_3d;
    }
};
} // namespace wmtk
