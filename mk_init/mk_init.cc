#include <iostream>
#include <fstream>
#include <vector>
#include "../common.hpp"

const char fname[] = "init.prof";
int dim = 2;

struct Particle {
	Float x,y;
	Type type;
};

template<typename T>
void ask(const char *str, T &v){
	std::cout<<str<<": ";
	std::cin>>v;
}

void error(const char *str){
	std::cerr<<str<<std::endl;
	exit(-1);
}

void error(){
	error("error");
}

void write_data(){
	using namespace rw;
	using std::endl;
	std::ofstream f(fname);
	f << 0 << " " << 0.0 << endl
		<< dim << endl
		<< r_in << " " << rw::r_out << endl
		<< theta <<" "<<rw::w<<endl;

	auto width = r_out - r_in;
	auto xmax = r_out+width;
	auto ymax = r_out+width;
	auto dx = xmax / 50;
	auto dy = ymax / 50;

	std::vector<Particle> v;

	for(auto y = -1*ymax;y<=ymax;y+=dy){
		for(auto x=-1*xmax;x<=xmax;x+=dx){
			Particle p;
			auto r2 = x*x + y*y;
			if(r2 > r_out*r_out){
				if(r2 > (r_out+width)*(r_out+width))
					continue;
				p.type = WALL;
			}else if(r2 > r_in*r_in){
				p.type = FLUID;
			}else if(width >= r_in){
				p.type = WALL;
			}else{
				if(r2 > (r_in-width)*(r_in-width))
					p.type = WALL;
				else
					continue;
			}
			p.x = x;
			p.y = y;
			v.push_back(p);
		}
	}

	std::cout<<"number of particles: "<<v.size()<<std::endl;
	f << v.size() << endl;

	for(auto it = v.begin();it!=v.end();it++){
//		if(it->type == WALL) continue;
		f << it->type << " "
			<< it->x << " "
			<< it->y << " "
			<< 0.0 << " "
			<< 0.0 << " "
			<< 0.0 << endl;
	}
}

int main(int argc, char **argv){
	using namespace rw;

	std::cout<<"make init situation."<<std::endl;

	ask("次元", dim);
	if(dim != 2 && dim != 3)
		error();
	ask("内半径", r_in);
	ask("外半径", r_out);
	if(r_in < 0.0 || r_out < 0.0)
		error();
	if(r_in >= r_out)
		error();
	ask("角速度", w);

	write_data();

	return 0;
}
