/**
 * File:   grid.cc
 *
 * Date:   21.07.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/bucketing/grid.h"

using namespace trivis;
using namespace trivis::bucketing;

void Grid::Clear() {
    for (auto &cell: _data) {
        cell.clear();
    }
}
