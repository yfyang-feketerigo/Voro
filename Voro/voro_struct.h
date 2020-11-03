#pragma once
#include"particle.h"
#include "configuration.h"
#include <voro++/voro++.hh>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
struct Voro_vertices
{
	size_t DIMENSION = 3;
	size_t center_pid = 0;
	std::vector<double> relative_p;
	std::vector<double> relative_r;
public:
	void compute_re_r();

	inline double max_re_r() const
	{
		return *std::max_element(relative_r.begin(), relative_r.end());
	}
};

struct Voro_volume
{
	size_t center_pid = 0;
	double volume = 0;
};
void voro_to_dump(const Configuration& config, string fname, const vector<Voro_vertices>& vertices, const vector<Voro_volume>& volume);
//void compute_voro(voro::container_periodic_poly& con);
const Voro_vertices& find_vertex(const vector<Voro_vertices>& vec_vertex, size_t center_pid);
vector<Voro_vertices> sort_voro_vertex(const vector<Particle>& vec_pa, const vector<Voro_vertices>& vec_vertex);
void config_to_vorocontainer(const Configuration& config, voro::container_periodic_poly& con);
void get_voro_info(voro::container_periodic_poly& con, size_t particle_num, vector<Voro_vertices>& vec_vertices, vector<Voro_volume>& vec_volume);


const Voro_volume& find_volume(const vector<Voro_volume>& vec_volume, size_t center_pid);
vector<Voro_volume> sort_voro_volume(const vector<Particle>& vec_pa, const vector<Voro_volume>& vec_volume);
//vector<Voro_volume> get_voro_volume(const Configuration& config, voro::container_periodic_poly& con);