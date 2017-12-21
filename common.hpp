#pragma once

using Float = double;

enum Type {
        GHOST,
        FLUID,
        WALL,
        DUMMY
};

namespace rw {
	Float r_in, r_out;
	Float theta, w;
}
