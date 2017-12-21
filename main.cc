#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <iomanip>
#include <cmath>
#include <sksat/math/vector.hpp>
#include "common.hpp"

namespace simulation {
	size_t time_step = 0;
	Float time=.0, dt=0.001, finish_time=2.0;
	size_t dim=2;

	const Float pcl_dst = 0.02; // 平均粒子間距離(今は決め打ち)
	const Float crt_num = 0.1; // クーラン条件数

	size_t particle_number=0;
	std::vector< sksat::math::vector<Float> > acc, vel, pos;
	std::vector<int> type;
	sksat::math::vector<Float> min, max;
	Float r, r2; // 影響半径
	namespace bucket {
		Float width, width2, width_inv;
		size_t nx, ny;
		size_t nxy;
		std::unique_ptr<int[]> bfst, blst, nxt;
	}

	size_t output_interval = 20;
	size_t file_number = 0;

	void main_loop();
	void make_bucket();
	void move_particle();
	void move_body();

	void load_file(const std::string &fname);
	void alloc_bucket();
	void write_file(const size_t &step, const Float &time);
}

int usage(const char *s){
	std::cout << "> "<< s << " input_file" << std::endl;
	return -1;
}

int main(int argc, char **argv){
	if(argc != 2)
		return usage(argv[0]);

	simulation::load_file(argv[1]);
	simulation::alloc_bucket();

	std::cout<<" *** START SIMULATION *** "<<std::endl;
	simulation::main_loop();
	std::cout<<" ***  END SIMULATION  *** "<<std::endl;
	return 0;
}

#define PRINT(v) {std::cout<< #v <<": "<<v<<std::endl;}

void simulation::load_file(const std::string &fname){
	std::cout<<"loading file \""<<fname<<"\"..."<<std::endl;
	std::ifstream f(fname);
	f >> time_step
		>> time
		>> dim
		>> rw::r_in >> rw::r_out
		>> rw::theta >> rw::w
		>> particle_number;

	PRINT(time_step);
	PRINT(time);
	PRINT(dim);
	PRINT(rw::r_in);
	PRINT(rw::r_out);
	PRINT(rw::theta);
	PRINT(rw::w);
	PRINT(particle_number);

	acc.reserve(particle_number);
	vel.reserve(particle_number);
	pos.reserve(particle_number);
	type.reserve(particle_number);

	for(auto i=0;i<particle_number;i++){
		f >> type[i]
			>> pos[i].x
			>> pos[i].y
			>> vel[i].x
			>> vel[i].y;
		if(f.eof()){
			std::cerr<<"stop: "<<i<<std::endl;
			throw std::runtime_error("");
		}
	}

	auto tmp = rw::r_out - rw::r_in;
	min.x = min.y = (rw::r_out + tmp) * -1.0;
	max.x = max.y = (rw::r_out + tmp);

	PRINT(min.x);
	PRINT(min.y);
	PRINT(max.x);
	PRINT(max.y);
}

void simulation::alloc_bucket(){
	using namespace bucket;

	r = pcl_dst * 2.1;
	r2 = r*r;

	width = r * (1.0 + crt_num); // バケット一辺の長さ
	width2 = width * width;
	width_inv = 1.0 / width;

	nx = (int)((max.x - min.x) * width_inv) + 1;
	ny = (int)((max.y - min.y) * width_inv) + 1;
	nxy = bucket::nx * bucket::ny;

	PRINT(bucket::nx);
	PRINT(bucket::ny);
	PRINT(bucket::nxy);

	bfst = std::make_unique<int[]>(nxy);
	blst = std::make_unique<int[]>(nxy);
	nxt  = std::make_unique<int[]>(particle_number);
}

void simulation::main_loop(){
	if(time_step == 0)
		write_file(0, time);
	while(true){
		make_bucket();
		move_particle();
		move_body();

		time_step++;
		time += dt;
		if( (time_step % output_interval) == 0 ){
			std::cout
				<< "time step: " << std::setw(5) << time_step << "  "
				<< "time: " << std::setw(5) << time << "  "
				<< "particle number: " << particle_number
				<< std::endl;
			write_file(time_step, time);
		}
		if(time >= finish_time) break;
	}
}

void simulation::make_bucket(){
	using namespace bucket;
	for(auto i=0;i<nxy;i++){ bfst[i] = -1; blst[i] = -1; }
	for(auto i=0;i<particle_number;i++){ nxt[i] = -1; }
	for(auto i=0;i<particle_number;i++){
		if(type[i] == GHOST) continue;
		int ix = static_cast<int>((pos[i].x - min.x) * width_inv) + 1;
		int iy = static_cast<int>((pos[i].y - min.y) * width_inv) + 1;
		int ib = iy * nx + ix;
		int j = blst[ib];
		blst[ib] = i;
		if(j == -1) bfst[ib] = i;
		else nxt[j] = i;
	}
}

void simulation::move_particle(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] != FLUID) continue;
		vel[i].x += acc[i].x;
		vel[i].y += acc[i].y;
		pos[i].x += vel[i].x;
		pos[i].y += vel[i].y;
	}
}

void simulation::move_body(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] != WALL) continue;
		auto dtheta = rw::w *dt;
		auto s = sin(dtheta);
		auto c = cos(dtheta);
		auto x = (pos[i].x * c) - (pos[i].y * s);
		auto y = (pos[i].x * s) + (pos[i].y * c);
		rw::theta += dtheta;
		pos[i].x = x;
		pos[i].y = y;
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
		<< rw::r_in << " " << rw::r_out << endl
		<< rw::theta << " " << rw::w << endl;

	for(auto i=0;i<particle_number;i++){
		f << type[i] << " "
			<< pos[i].x << " "
			<< pos[i].y << " "
			<< vel[i].x << " "
			<< vel[i].y << endl;
	}

	file_number++;
}
