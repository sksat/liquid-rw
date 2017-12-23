#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <sksat/math/vector.hpp>
#include "common.hpp"

#define CHECK_DT
#define CHECK_SAME_POS

namespace simulation {
	size_t time_step = 0;
	Float time=.0, dt=0.000000001, finish_time=0.5;
	size_t dim=2;

	// 定数たち
	Float pcl_dst	= 0.02;		// 平均粒子間距離(今は決め打ち)
	const Float crt_num	= 0.1;		// クーラン条件数
	const Float knm_vsc	= 0.000001;	// 動粘性係数
	const Float dens_fluid	= 1000;		// 液体粒子の密度
	const Float dens_wall	= 1000;		// 壁粒子の密度
	const Float dst_lmt_rat	= 0.9;		// これ以上の粒子間の接近を許さない距離の係数
	const Float col_rat	= 0.2;		// 接近した粒子の反発率
	const Float col		= 1.0 + col_rat;
	const Float sound_vel	= 22.0;		// 流体の音速(?)

	size_t particle_number=0;
	std::vector< sksat::math::vector<Float> > acc, vel, pos;
	std::vector<Float> press;
	std::vector<int> type;
	sksat::math::vector<Float> min, max;

	Float r, r2;			// 影響半径
	Float A1;			// 粘性項の計算に用いる係数
	Float A2;			// 圧力の計算に用いる係数
	Float A3;			// 圧力勾配項の計算に用いる係数
	Float n0;			// 初期粒子数密度
	Float lmd;			// ラプラシアンモデルの係数λ
	Float dens[NUM_TYPE];		// 粒子種類ごとの密度
	Float dens_inv[NUM_TYPE];
	Float rlim, rlim2;		// これ以上の粒子間の接近を許さない距離

	namespace bucket {
		Float width, width2, width_inv;
		size_t nx, ny;
		size_t nxy;
		std::unique_ptr<int[]> first, last, next;
	}

	size_t output_interval = 100;
	size_t file_number = 0;

	template<typename T>
	inline T weight(T dist, T re){ return ((re/dist) - 1.0); } // 重み関数

	void main_loop();
	void make_bucket();
	void calc_tmpacc();
	bool check_overflow(size_t i);
	void move_particle_tmp();
	void check_collision();
	void make_press();
	void calc_press_grad(); // 圧力勾配項
	void move_particle();
	void move_body();

	void load_file(const std::string &fname);
	void alloc_bucket();
	void set_param();
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
	simulation::set_param();

	std::cout<<" *** START SIMULATION *** "<<std::endl;

	auto start = std::chrono::system_clock::now();
	simulation::main_loop();
	auto end = std::chrono::system_clock::now();

	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count(); ;
	std::cout<<"Total : "<<elapsed<<"sec"<<std::endl;

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

	pcl_dst = ((rw::r_out-rw::r_in)+rw::r_out)/static_cast<Float>(NUM_PER_DST);
	PRINT(pcl_dst);

	acc.reserve(particle_number);
	vel.reserve(particle_number);
	pos.reserve(particle_number);
	press.reserve(particle_number);
	type.reserve(particle_number);

	for(auto i=0;i<particle_number;i++){
		f >> type[i]
			>> pos[i].x
			>> pos[i].y
			>> vel[i].x
			>> vel[i].y
			>> press[i];
		if(f.eof()){
			std::cerr<<"stop: "<<i<<std::endl;
			throw std::runtime_error("");
		}
	}

	auto tmp = rw::r_out - rw::r_in;
	min.x = min.y = (rw::r_out + tmp) * -1.0 - pcl_dst;
	max.x = max.y = (rw::r_out + tmp) + pcl_dst;

	PRINT(min.x);
	PRINT(min.y);
	PRINT(max.x);
	PRINT(max.y);

	for(auto i=0;i<particle_number;i++){
		if(pos[i].x <= min.x || max.x <= pos[i].x){
			std::cout<<"error: ";
			PRINT(pos[i].x);
		}
		if(pos[i].y <= min.y || max.y <= pos[i].y){
			std::cout<<"error: ";
			PRINT(pos[i].y);
		}
	}
}

void simulation::alloc_bucket(){
	using namespace bucket;

	r = pcl_dst * 2.1;
	r2 = r*r;

	width = r * (1.0 + crt_num); // バケット一辺の長さ
	width2 = width * width;
	width_inv = 1.0 / width;

	nx = static_cast<int>((max.x - min.x) * width_inv) + 3;
	ny = static_cast<int>((max.y - min.y) * width_inv) + 3;
	nxy = nx * ny;

	PRINT(bucket::nx);
	PRINT(bucket::ny);
	PRINT(bucket::nxy);

	bucket::first = std::make_unique<int[]>(nxy*4);
	bucket::last = std::make_unique<int[]>(nxy*4);
	bucket::next  = std::make_unique<int[]>(particle_number);
}

void simulation::set_param(){
	Float tn0 = 0.0;
	Float tlmd= 0.0;

	for(int ix=-4;ix<5;ix++){
		for(int iy=-4;iy<5;iy++){
			Float x = pcl_dst * static_cast<double>(ix);
			Float y = pcl_dst * static_cast<double>(iy);
			Float dist2 = x*x + y*y;
			if(dist2 <= r2){
				if(dist2 == 0.0) continue;
				Float dist = sqrt(dist2);
				Float tmp = weight(dist, r);
				tn0 += tmp;
				tlmd+= dist2 * tmp;
			}
		}
	}
	n0	= tn0;
	lmd	= tlmd / tn0;
	A1 = 2.0 * knm_vsc * static_cast<Float>(dim) / (n0 * lmd);
	A2 = sound_vel * sound_vel / n0;
	A3 = -1.0 * static_cast<Float>(dim) / n0;

	PRINT(n0);
	PRINT(lmd);
	PRINT(A1);
	PRINT(A2);
	PRINT(A3);

	dens[FLUID]	= dens_fluid;
	dens[WALL]	= dens_wall;
	dens_inv[FLUID]	= 1.0/dens_fluid;
	dens_inv[WALL]	= 1.0/dens_wall;

	rlim	= pcl_dst * dst_lmt_rat;
	rlim2	= rlim * rlim;

	PRINT(dens[FLUID]);
	PRINT(dens[WALL]);
	PRINT(rlim);

	auto D = 0.9;
	if(
		dt < (D * (pcl_dst * pcl_dst) / knm_vsc) &&
		dt < (pcl_dst / sound_vel)
	){
		std::cout<<"dt is ok!"<<std::endl;
	}else{
		throw std::runtime_error("dt is not ok.");
	}
}

void simulation::main_loop(){
//	for(auto i=0;i<particle_number;i++){
//		if(type[i] == WALL)
//			type[i] = GHOST;
//	}
	if(time_step == 0)
		write_file(0, time);
	while(true){
		make_bucket();

		calc_tmpacc();
		Float max_vel = 0.0;
#ifdef CHECK_DT
		for(auto i=0;i<particle_number;i++){
			if(-rw::r_out/100 < pos[i].y && pos[i].y < rw::r_out/100){
//				if(pos[i].x > 0.0)
//					vel[i].y = 0.001;
			}
			if(type[i] != FLUID) continue;
			auto v = vel[i].x*vel[i].x + vel[i].y*vel[i].y;
			if(max_vel < v) max_vel = v;
#ifdef CHECK_SAME_POS
			for(auto j=0;j<particle_number;j++){
				if(i == j) continue;
				if(type[i] != FLUID) continue;
				if(pos[i].x == pos[j].x && pos[i].y == pos[j].y){
					PRINT(i);
					PRINT(j);
					PRINT(pos[i].x);
					PRINT(pos[j].x);
					PRINT(pos[i].y);
					PRINT(pos[j].y);
//					throw std::runtime_error("fuck NVIDIA!");
				}
			}
#endif
		}
#endif
		if(dt < (0.2 * pcl_dst / max_vel)){}
		else
			throw std::runtime_error("dt is not ok");
		move_particle_tmp(); // temporary
		check_collision();
		make_press();
		calc_press_grad();
		move_particle();
		make_press();

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
	for(int i=0;i<nxy*4;i++){ first[i] = -1; last[i] = -1; }
	for(int i=0;i<particle_number;i++){ next[i] = -1; }
	for(int i=0;i<particle_number;i++){
		if(type[i] == GHOST) continue;
		int ix = static_cast<int>((pos[i].x - min.x) * width_inv) + 1;
		int iy = static_cast<int>((pos[i].y - min.y) * width_inv) + 1;
		int ib = iy * nx + ix;
		int j = last[ib];
		last[ib] = i;
		if(j == -1){ first[ib] = i; }
		else{ next[j] = i; }
	}
}

void simulation::calc_tmpacc(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] != FLUID) continue;
		sksat::math::vector<Float> Acc;
		Acc.x = Acc.y = 0.0;
		int ix = static_cast<int>((pos[i].x - min.x) * bucket::width_inv) + 1;
		int iy = static_cast<int>((pos[i].y - min.y) * bucket::width_inv) + 1;
		for(int jy=iy-1;jy<=iy+1;jy++){
			for(int jx=ix-1;jx<=ix+1;jx++){
				int jb = jy*bucket::nx + jx;
				if(jb >= bucket::nxy*4){
					PRINT(i);
					PRINT(type[i]);
					PRINT(pos[i].x);
					PRINT(pos[i].y);
					PRINT(ix);
					PRINT(iy);
					PRINT(jx);
					PRINT(jy);
					PRINT(jb);
					getchar();
				}
				int j = bucket::first[jb];
				if(j == -1) continue;
				for(;;){
					Float v0 = pos[j].x - pos[i].x;
					Float v1 = pos[j].y - pos[i].y;
					Float dist2 = v0*v0 + v1*v1;
					if(dist2 < r2){
						if(j!=i && type[j]!=GHOST){
							Float dist = sqrt(dist2);
							Float w = weight(dist, r);
							Acc += (vel[j] - vel[i]) * w;
						}
					}
					j = bucket::next[j];
					if(j == -1) break;
				}
			}
		}
		acc[i] = Acc * A1;
//		acc[i].y += -9.8;
	}
}

bool simulation::check_overflow(size_t i){
	if(
		pos[i].x > max.x || pos[i].x < min.x ||
		pos[i].y > max.y || pos[i].y < min.y
	  ){
		type[i] = GHOST;
		press[i] = vel[i].x = vel[i].y = 0.0;
		return true;
	}
	return false;
}

void simulation::move_particle_tmp(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] != FLUID) continue;
		vel[i] += acc[i] * dt;
		pos[i] += vel[i] * dt;
		acc[i].x = acc[i].y = acc[i].z = 0.0;
		check_overflow(i);
//		if(check_overflow(i)) std::cout<<"move_particle_tmp: "<<i<<std::endl;
	}
}

void simulation::check_collision(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] != FLUID) continue;
		Float mi = dens[type[i]];
		sksat::math::vector<Float> vec_i = vel[i], vec_i2 = vel[i];
		int ix = static_cast<int>((pos[i].x - min.x) * bucket::width_inv) + 1;
		int iy = static_cast<int>((pos[i].y - min.y) * bucket::width_inv) + 1;
		for(int jy=iy-1; jy<=iy+1; jy++){
			for(int jx=ix-1; jx<=ix+1; jx++){
				int jb = jy * bucket::nx + ix;
				int j = bucket::first[jb];
				if(j == -1) continue;
				for(;;){
					Float v0 = pos[j].x - pos[i].x;
					Float v1 = pos[j].y - pos[i].y;
					Float dist2 = v0*v0 + v1*v1;
					if(dist2 < rlim2){
						if(j!=i && type[i]!= GHOST){
							Float fDT = (vec_i.x-vel[j].x)*v0+(vec_i.y-vel[j].y)*v1;
							if(fDT > 0.0){
								Float mj = dens[type[j]];
								fDT *= col * mj / (mi+mj) / dist2;
								vec_i2.x -= v0 * fDT;
								vec_i2.y -= v1 * fDT;
							}
						}
					}
					j = bucket::next[j]; // バケット内の次の粒子へ
					if(j == -1) break;
				}
			}
		}
		acc[i] = vec_i2;
	}
	for(auto i=0;i<particle_number;i++){
		vel[i] = acc[i];
	}
}

void simulation::make_press(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] == GHOST) continue;
//		if(-0.1 < pos[i].y && pos[i].y < 0.1){
//			press[i] = 100.0;
//			std::cout<<"a:"; getchar();
//			continue;
//		}
		Float ni = 0.0;
		int ix = static_cast<int>((pos[i].x - min.x)*bucket::width_inv) + 1;
		int iy = static_cast<int>((pos[i].y - min.y)*bucket::width_inv) + 1;
		for(int jy=iy-1;jy<=iy+1;jy++){
			for(int jx=ix-1;jx<=ix+1;jx++){
				int jb = jy * bucket::nx + ix;
				int j = bucket::first[jb];
				if(j == -1) continue;
				for(;;){
					Float v0 = pos[j].x - pos[i].x;
					Float v1 = pos[j].y - pos[i].y;
					Float dist2 = v0*v0 + v1*v1;
					if(dist2 < r2){
						if(j!=i && type[j]!=GHOST){
							Float dist = sqrt(dist2);
							Float w = weight(dist, r);
							ni += w;
						}
					}
					j = bucket::next[j];
					if(j == -1) break;
				}
			}
		}
		Float mi = dens[type[i]];
		Float pressure = (ni > n0) * (ni - n0) * A2 * mi;
		press[i] = pressure;
//		if(pressure != 0.0) PRINT(pressure);
	}
}

void simulation::calc_press_grad(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] != FLUID) continue;
		sksat::math::vector<Float> Acc;
		Acc.x = Acc.y = Acc.z = 0.0;
		Float pre_min = press[i];
		int ix = static_cast<int>((pos[i].x - min.x)*bucket::width_inv) + 1;
		int iy = static_cast<int>((pos[i].y - min.y)*bucket::width_inv) + 1;
		for(int jy=iy-1;jy<=iy+1;jy++){
			for(int jx=ix-1;jx<=ix+1;jx++){
				int jb = jy*bucket::nx + jx;
				int j = bucket::first[jb];
				if(j == -1) continue;
				for(;;){
					Float v0 = pos[j].x - pos[i].x;
					Float v1 = pos[j].y - pos[i].y;
					Float dist2 = v0*v0 + v1*v1;
					if(dist2 < r2){
						if(j!=i && type[j]!=GHOST){
							if(pre_min > press[j]) pre_min = press[j];
						}
					}
					j = bucket::next[j];
					if(j == -1) break;
				}
			}
		}
		for(int jy=iy-1;jy<=iy+1;jy++){
			for(int jx=ix-1;jx<=ix+1;jx++){
				int jb = jy * bucket::nx + jx;
				int j = bucket::first[jb];
				if(j == -1) continue;
				for(;;){
					Float v0 = pos[j].x - pos[i].x;
					Float v1 = pos[j].y - pos[i].y;
					Float dist2 = v0*v0 + v1*v1;
					if(dist2 < r2){
						if(j!=i && type[j]!=GHOST){
							Float dist = sqrt(dist2);
							Float w = weight(dist, r);
							w *= (press[j] - pre_min)/dist2;
							Acc.x += v0*w;
							Acc.y += v1*w;
						}
					}
					j = bucket::next[j];
					if(j == -1) break;
				}
			}
		}
		acc[i] = Acc * dens_inv[FLUID] * A3;
	}
}

void simulation::move_particle(){
	for(auto i=0;i<particle_number;i++){
		if(type[i] != FLUID) continue;
		vel[i] += acc[i] * dt;
		pos[i] += acc[i] * dt * dt;
		check_overflow(i);
//		if(check_overflow(i)) std::cout<<"move_particle: "<<acc[i].x<<std::endl;
		acc[i].x = acc[i].y = acc[i].z = 0.0;
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
		check_overflow(i);
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
		if(type[i] == GHOST) continue;
		f << type[i] << " "
			<< pos[i].x << " "
			<< pos[i].y << " "
			<< vel[i].x << " "
			<< vel[i].y << " "
			<< press[i] << endl;
	}

	file_number++;
}
