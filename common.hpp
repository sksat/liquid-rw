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

namespace pump {
	Float capacity;
	Float v;
	int num, num_once, times;
}

const int NUM_PER_DST = 20;
