#include <iostream>
#include <vector>
#include <cmath>
#include <sksat/math/vector.hpp>
#include "common.hpp"

namespace simulation {
	size_t time_step = 0;
	Float time=.0, dt=0.001, finish_time=2.0;

	size_t particle_number=0;
	std::vector< sksat::math::vector<Float> > acc, vel, pos;

	size_t output_interval = 20;

	void main_loop();
}

int main(int argc, char **argv){
	std::cout<<" *** START SIMULATION *** "<<std::endl;

	simulation::main_loop();

	std::cout<<" ***  END SIMULATION  *** "<<std::endl;
	return 0;
}

void simulation::main_loop(){
	while(true){
		time_step++;
		time += dt;
		if( (time_step % output_interval) == 0 ){
			std::cout
				<<"time step: "<<time_step<<"  "
				<<"time: "<<time<<"  "
				<<"particle number: "<<particle_number
				<<std::endl;
		}
		if(time >= finish_time) break;
	}
}
