#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <sksat/math/vector.hpp>
#include "common.hpp"

namespace simulation {
	size_t time_step = 0;
	Float time=.0, dt=0.001, finish_time=2.0;

	size_t particle_number=0;
	std::vector< sksat::math::vector<Float> > acc, vel, pos;

	size_t output_interval = 20;
	size_t file_number = 0;

	void main_loop();
	void write_file(const size_t &step, const Float &time);
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
				<< "time step: " << time_step << "  "
				<< "time: " << std::setw(5) << time << "  "
				<< "particle number: " << particle_number
				<< std::endl;
			write_file(time_step, time);
		}
		if(time >= finish_time) break;
	}
}

void simulation::write_file(const size_t &step, const Float &time){
	using std::endl;
	std::stringstream fname;
	fname << "out/"
		<< "output_"
		<< std::setfill('0')
		<< std::setw(5)
		<< file_number
		<< ".prof";
	std::ofstream f(fname.str());
	if(f.fail()) throw std::runtime_error("cannot open file.\nyou should make \'out\' dir.");

	f << time_step << endl
		<< time << endl
		<< particle_number << endl
		<< rw::r_in << "," << rw::r_out << rw::w << endl;

	for(auto i=0;i<particle_number;i++){
		f << pos[i].x
			<< pos[i].y
			<< vel[i].x
			<< vel[i].y
			<< endl;
	}

	file_number++;
}
