#pragma once

using Float = double;

enum Type {
        GHOST,
        FLUID,
        WALL,
        DUMMY,
	NUM_TYPE
};

namespace rw {
	Float r_in, r_out;
	Float theta, w;
}

const int NUM_PER_DST = 30;
