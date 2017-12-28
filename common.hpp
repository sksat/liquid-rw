#pragma once
#include <vector>
#include <string>
#include <sksat/math/vector.hpp>

using Float = double;

#define PRINT(v) {std::cout<< #v <<": "<<v<<std::endl;}

namespace simulation {

	size_t time_step	= 0;
	size_t dim		= 3;
	Float time		= .0;

enum Type {
        GHOST,
        FLUID,
        WALL,
        DUMMY,
	NUM_TYPE
};

namespace particle {
	size_t num;
	Float v;
}

std::vector< sksat::math::vector<Float> > acc, vel, pos;
std::vector<Float> press, pav;
std::vector<int> type;

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

void load_file(const std::string &fname, size_t add_num=0, bool print_flg=false){
	std::ifstream f;

	std::cout<<"loading file \""<<fname<<"\"...";
        f.open(fname);
	if(!f){
		std::cout<<"\tfailed."<<std::endl;
		throw std::runtime_error("cannot open file");
	}

	f >> time_step
		>> time
		>> dim
		>> rw::r_in >> rw::r_out
		>> rw::theta >> rw::w
		>> pump::capacity
		>> pump::v
		>> particle::num;

	particle::num += add_num;

	if(print_flg){
		PRINT(particle::num);
		PRINT(time_step);
		PRINT(time);
		PRINT(dim);
		PRINT(rw::r_in);
		PRINT(rw::r_out);
		PRINT(rw::theta);
		PRINT(rw::w);
		PRINT(pump::capacity);
		PRINT(particle::num);
	}

        pos.reserve(particle::num);
        vel.reserve(particle::num);
        acc.reserve(particle::num);
        press.reserve(particle::num);
        type.reserve(particle::num);
        pav.reserve(particle::num);

        for(auto i=0;i<particle::num;i++){
		f >> type[i]
			>> pos[i].x
			>> pos[i].y
			>> pos[i].z
			>> vel[i].x
			>> vel[i].y
			>> vel[i].z
			>> press[i]
			>> pav[i];
		if(f.fail()){
			for(;i<particle::num;i++)
				type[i] = GHOST;
		}
	}
}

}
