/**
 * File:   status.h
 *
 * Date:   21.11.2024
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PATHFINDER_STATUS_H_
#define TRIVIS_PATHFINDER_STATUS_H_

namespace trivis_pathfinder::utils {

enum class Status : int {
    kOk = 0,
    kErrorSourceOutside = 1,
    kErrorTargetOutside = 2,
    kErrorNoReflexVisibilityGraph = 3,
    kErrorNoCitiesVisibilityGraph = 4,
};

template<typename T>
struct StatusWithResult {
    Status status;
    T result = T(0.0);
};

}

#endif //TRIVIS_PATHFINDER_STATUS_H_
