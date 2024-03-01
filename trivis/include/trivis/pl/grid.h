/**
 * File:   grid.h
 *
 * Date:   24.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PL_GRID_H_
#define TRIVIS_PL_GRID_H_

#include <vector>

namespace trivis::pl {

template<typename Cell>
class Grid {

public:

    Grid() = default;

    Grid(int n_row, int n_col) : _n_row(n_row), _n_col(n_col), _cells(n_row * n_col) {}

    void Resize(int n_row, int n_col) {
        _n_row = n_row;
        _n_col = n_col;
        _cells.resize(n_row * n_col);
    }

    void Clear() {
        _n_row = 0;
        _n_col = 0;
        _cells.clear();
    }

    [[nodiscard]] int n_row() const { return _n_row; }

    [[nodiscard]] int n_col() const { return _n_col; }

    [[nodiscard]] const auto &data() const { return _cells; }

    [[nodiscard]] const Cell &operator()(int row, int col) const noexcept(true) { return _cells[col * _n_col + row]; }

    [[nodiscard]] Cell &operator()(int row, int col) noexcept(true) { return _cells[col * _n_col + row]; }

private:

    int _n_row = 0;
    int _n_col = 0;
    std::vector<Cell> _cells;

};

}

#endif //TRIVIS_PL_GRID_H_
