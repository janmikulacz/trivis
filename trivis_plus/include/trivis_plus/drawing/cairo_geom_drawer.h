/**
 * File:    cairo_geom_drawer.h
 *
 * Date:    18.08.2020
 * Author:  Jan Mikula
 * E-mail:  jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_DRAWING_CAIRO_GEOM_DRAWER_H_
#define TRIVIS_PLUS_DRAWING_CAIRO_GEOM_DRAWER_H_

#include <memory>
#include <iostream>

#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>

#include "trivis_plus/drawing/colors.h"

#include "trivis/geom/generic_geom_types.h"

namespace trivis_plus::drawing {

class CairoGeomDrawer {

public:

    static constexpr auto kDefaultSizeX = 1920;
    static constexpr auto kDefaultSizeY = 1080;
    static constexpr auto kDefaultRes = 1.0;

    explicit CairoGeomDrawer(
        double scene_size_x = kDefaultSizeX,
        double scene_size_y = kDefaultSizeY,
        double res = kDefaultRes,
        double scene_offset_x = 0.0,
        double scene_offset_y = 0.0
    )
        : _res(res),
          _scene_size_x(scene_size_x),
          _scene_size_y(scene_size_y),
          _scene_off_x(scene_offset_x),
          _scene_off_y(scene_offset_y),
          _pic_size_x(_scene_size_x * _res),
          _pic_size_y(_scene_size_y * _res),
          _pic_size_x_int(static_cast<int>(std::floor(_pic_size_x))),
          _pic_size_y_int(static_cast<int>(std::floor(_pic_size_y))) {
    }

    [[nodiscard]] static RGB SaturateColor(const RGB &color) {
        auto r = color.r, g = color.g, b = color.b;
        auto rgb_min = std::min(std::min(r, g), b);
        r -= rgb_min, g -= rgb_min, b -= rgb_min;
        auto rgb_max = static_cast<double>(std::max(std::max(r, g), b));
        r = static_cast<int>(std::floor(static_cast<double>(r) * 255.0 / rgb_max));
        g = static_cast<int>(std::floor(static_cast<double>(g) * 255.0 / rgb_max));
        b = static_cast<int>(std::floor(static_cast<double>(b) * 255.0 / rgb_max));
        return RGB{r, g, b};
    }

    [[nodiscard]] bool is_open() const { return _is_open; }

    void OpenPDF(const std::string &pdf_file) {
        _is_open = true;
        _surface_ptr = cairo_pdf_surface_create(pdf_file.c_str(), _pic_size_x, _pic_size_y);
    }

    void OpenImage() {
        _is_open = true;
        _surface_ptr = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, _pic_size_x_int, _pic_size_y_int);
    }

    void SaveToPng(const std::string &png_file) const {
        cairo_surface_write_to_png(_surface_ptr, png_file.c_str());
    }

    void Close() {
        _is_open = false;
        cairo_surface_flush(_surface_ptr);
        cairo_surface_destroy(_surface_ptr);
    }

    void DrawPlane(const RGB &fill_color = kColorWhite) const {
        auto fill_color_01 = RGB01(fill_color);
        auto *cr = cairo_create(_surface_ptr);
        cairo_move_to(cr, 0, 0);
        cairo_line_to(cr, _pic_size_x, 0);
        cairo_line_to(cr, _pic_size_x, _pic_size_y);
        cairo_line_to(cr, 0, _pic_size_y);
        cairo_close_path(cr);
        cairo_set_source_rgb(cr, fill_color_01.r, fill_color_01.g, fill_color_01.b);
        cairo_fill(cr);
        cairo_destroy(cr);
    }

    template<typename T>
    void DrawPoint(
        const trivis::geom::Point<T> &p,
        const T &radius = T(1),
        const RGB &fill_color = kColorBlack,
        double opacity = 1.0
    ) const {
        auto fill_color_01 = RGB01(fill_color);
        auto *cr = cairo_create(_surface_ptr);
        cairo_arc(
            cr,
            (p.x + _scene_off_x) * _res,
            _pic_size_y - (p.y + _scene_off_y) * _res,
            radius * _res,
            0.0,
            2.0 * M_PI
        );
        cairo_set_source_rgba(cr, fill_color_01.r, fill_color_01.g, fill_color_01.b, opacity);
        cairo_fill(cr);
        cairo_stroke(cr);
        cairo_destroy(cr);
    }

    template<typename T>
    void DrawText(
        const std::string &text,
        const trivis::geom::Point<T> &p,
        double font_size = 1.0,
        const RGB &fill_color = kColorBlack,
        double opacity = 1.0
    ) const {
        auto fill_color_01 = RGB01(fill_color);
        auto *cr = cairo_create(_surface_ptr);
        cairo_set_source_rgba(cr, fill_color_01.r, fill_color_01.g, fill_color_01.b, opacity);
        cairo_set_font_size(cr, font_size * _res);
        cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_text_extents_t te;
        cairo_text_extents(cr, text.c_str(), &te);
        cairo_move_to(
            cr,
            (p.x + _scene_off_x) * _res - te.x_bearing - te.width / 2,
            _pic_size_y - (p.y + _scene_off_y) * _res + te.height / 2
        );
        cairo_show_text(cr, text.c_str());
        cairo_destroy(cr);
    }

    template<typename T>
    void DrawPoints(
        const trivis::geom::Points<T> &points,
        const T &radius = T(1),
        const RGB &fill_color = kColorBlack,
        double opacity = 1.0
    ) const {
        for (auto &p : points) DrawPoint(p, radius, fill_color, opacity);
    }

    template<typename T>
    void DrawLine(
        const trivis::geom::Point<T> &p1,
        const trivis::geom::Point<T> &p2,
        const T &line_width = T(1),
        const RGB &line_color = kColorBlack,
        double opacity = 1.0
    ) const {
        auto line_color_01 = RGB01(line_color);
        auto *cr = cairo_create(_surface_ptr);
        cairo_move_to(cr, (p1.x + _scene_off_x) * _res, _pic_size_y - (p1.y + _scene_off_y) * _res);
        cairo_line_to(cr, (p2.x + _scene_off_x) * _res, _pic_size_y - (p2.y + _scene_off_y) * _res);
        cairo_set_line_width(cr, line_width * _res);
        cairo_set_source_rgba(cr, line_color_01.r, line_color_01.g, line_color_01.b, opacity);
        cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
        cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
        cairo_stroke(cr);
        cairo_destroy(cr);
    }

    template<typename T>
    void DrawPath(
        const trivis::geom::Points<T> &path,
        const T &line_width = T(1),
        const RGB &line_color = kColorBlack,
        double opacity = 1.0
    ) const {
        for (int i = 0; i < path.size() - 1; ++i) {
            DrawLine(path[i], path[i + 1], line_width, line_color, opacity);
        }
    }

    template<typename T>
    void DrawArc(
        const trivis::geom::Point<T> &p,
        const T &radius = T(1),
        double angle1 = 0.0,
        double angle2 = M_PI,
        const T &line_width = T(1),
        const RGB &line_color = kColorBlack,
        double opacity = 1.0
    ) const {
        auto line_color_01 = RGB01(line_color);
        auto *cr = cairo_create(_surface_ptr);
        cairo_arc(
            cr,
            (p.x + _scene_off_x) * _res,
            _pic_size_y - (p.y + _scene_off_y) * _res,
            radius * _res,
            angle1,
            angle2
        );
        cairo_set_line_width(cr, line_width * _res);
        cairo_set_source_rgba(cr, line_color_01.r, line_color_01.g, line_color_01.b, opacity);
        cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
        cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
        cairo_stroke(cr);
        cairo_destroy(cr);
    }

    template<typename T>
    void DrawPolygon(
        const trivis::geom::Polygon<T> &poly,
        const RGB &fill_color = kColorBlack,
        double opacity = 1.0
    ) const {
        DrawPolygon(poly, true, fill_color, false, T(0), kColorBlack, opacity);
    }

    template<typename T>
    void DrawPolygon(
        const trivis::geom::Polygon<T> &poly,
        const T &line_width,
        const RGB &line_color = kColorBlack,
        double opacity = 1.0
    ) const {
        DrawPolygon(poly, false, kColorBlack, true, line_width, line_color, opacity);
    }

    template<typename T>
    void DrawPolygon(
        const trivis::geom::Polygon<T> &poly,
        const T &line_width,
        const RGB &line_color,
        const RGB &fill_color,
        double opacity = 1.0
    ) const {
        DrawPolygon(poly, true, fill_color, true, line_width, line_color, opacity);
    }

    template<typename T>
    void DrawPolygons(
        const trivis::geom::Polygons<T> &poly,
        const RGB &fill_color = kColorBlack,
        double opacity = 1.0
    ) const {
        DrawPolygons(poly, true, fill_color, false, T(0), kColorBlack, opacity);
    }

    template<typename T>
    void DrawPolygons(
        const trivis::geom::Polygons<T> &poly,
        const T &line_width,
        const RGB &line_color = kColorBlack,
        double opacity = 1.0
    ) const {
        DrawPolygons(poly, false, kColorBlack, true, line_width, line_color, opacity);
    }

    template<typename T>
    void DrawPolygons(
        const trivis::geom::Polygons<T> &poly,
        const T &line_width,
        const RGB &line_color,
        const RGB &fill_color,
        double opacity = 1.0
    ) const {
        DrawPolygons(poly, true, fill_color, true, line_width, line_color, opacity);
    }

protected:

    struct RGB01 {
        double r = 1.0;
        double g = 1.0;
        double b = 1.0;
        explicit RGB01(const RGB &rgb)
            : r(static_cast<double>(rgb.r) / 255.0),
              g(static_cast<double>(rgb.g) / 255.0),
              b(static_cast<double>(rgb.b) / 255.0) {}
    };

private:
    bool _is_open = false;
    double _res;
    double _scene_size_x;
    double _scene_size_y;
    double _scene_off_x;
    double _scene_off_y;
    double _pic_size_x;
    double _pic_size_y;
    int _pic_size_x_int;
    int _pic_size_y_int;
    cairo_surface_t *_surface_ptr = nullptr;

    template<typename T>
    void DrawPolygon(
        const trivis::geom::Polygon<T> &poly,
        bool fill,
        const RGB &fill_color,
        bool stroke,
        const T &line_width,
        const RGB &line_color,
        double opacity = 1.0
    ) const {
        if (!poly.empty() && (fill || stroke)) {
            auto line_color_01 = RGB01(line_color);
            auto fill_color_01 = RGB01(fill_color);
            auto *cr = cairo_create(_surface_ptr);
            cairo_move_to(cr, (poly[0].x + _scene_off_x) * _res, _pic_size_y - (poly[0].y + _scene_off_y) * _res);
            for (int i = 1; i < poly.size(); ++i) {
                cairo_line_to(cr, (poly[i].x + _scene_off_x) * _res, _pic_size_y - (poly[i].y + _scene_off_y) * _res);
            }
            cairo_close_path(cr);
            if (fill) {
                cairo_set_source_rgba(cr, fill_color_01.r, fill_color_01.g, fill_color_01.b, opacity);
                cairo_fill_preserve(cr);
            }
            if (stroke) {
                cairo_set_line_width(cr, line_width * _res);
                cairo_set_source_rgba(cr, line_color_01.r, line_color_01.g, line_color_01.b, opacity);
                cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
                cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
                cairo_stroke(cr);
            }
            cairo_destroy(cr);
        }
    }

    template<typename T>
    void DrawPolygons(
        const trivis::geom::Polygons<T> &polys,
        bool fill,
        const RGB &fill_color,
        bool stroke,
        const T &line_width,
        const RGB &line_color,
        double opacity = 1.0
    ) const {
        for (auto &poly : polys) {
            DrawPolygon(poly, fill, fill_color, stroke, line_width, line_color, opacity);
        }
    }

};

}

#endif //TRIVIS_PLUS_DRAWING_CAIRO_GEOM_DRAWER_H_
