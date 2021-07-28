#pragma once
#include"particle.h"
#include "configuration.h"
#include <voro++/voro++.hh>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>


// vertices of a voronoi cell
// DIMENSION is the dimension of space, center_pid is the id of center particle
// all vertices info stored as a one-dimension vector
// for each vertice the components store as three adjacent number
// e.g. in 3d space vector structure will be like {v1x, v1y, v1z, v2x, v2y, v2z ...}
struct Voro_vertices
{
	size_t DIMENSION = 3;
	size_t center_pid = 0;
	std::vector<double> relative_p;
	std::vector<double> relative_r;
public:
	void compute_relative_r();

	inline double max_re_r() const
	{
		return *std::max_element(relative_r.begin(), relative_r.end());
	}
};


// volume of a voronoi cell
// center_pid is the id of center particle
// volume is volume of voronoi cell, of course
struct Voro_volume
{
	size_t center_pid = 0;
	double volume = 0;
};


// coordination number based on voronoi cell
// which is the number of atoms forms the first voronoi cell, equal to number of faces of voronoi cell
struct Voro_coordination_num
{
	size_t center_pid = 0;
	size_t voro_coordination_num = 0;
};

// write given voronoi info into file as cols
// config: original configuration
// fname: output file name
// other parameters: information to output
void voro_to_dump(const Configuration& config, string fname, const vector<Voro_vertices>& vertices, const vector<Voro_volume>& volume, const vector<Voro_coordination_num>& cn);

// convert Configuration class to voro::container_periodic_poly class
// config: original configuration 
// con: voronoi container to be filled with configuration
void config_to_vorocontainer(const Configuration& config, voro::container_periodic_poly& con);

// compute voronoi & store voronoi info in vector, respectively
void get_voro_info(
	voro::container_periodic_poly& con,
	size_t particle_num,
	vector<Voro_vertices>& vec_vertices,
	vector<Voro_volume>& vec_volume,
	vector<Voro_coordination_num>& vec_cn
);

//x const Voro_vertices& find_vertex(const vector<Voro_vertices>& vec_vertex, size_t center_pid);
//x vector<Voro_vertices> sort_voro_vertex(const vector<Particle>& vec_pa, const vector<Voro_vertices>& vec_vertex);
//x const Voro_volume& find_volume(const vector<Voro_volume>& vec_volume, size_t center_pid);
//x vector<Voro_volume> sort_voro_volume(const vector<Particle>& vec_pa, const vector<Voro_volume>& vec_volume);
//x void compute_voro(voro::container_periodic_poly& con);

// find para in vector<T_para> for given center_pid particle 
template <typename T_para>
const T_para& find_para(const vector<T_para>& vec_para, size_t _center_pid)
{
	for (size_t i = 0; i < vec_para.size(); i++)
	{
		if (vec_para[i].center_pid == _center_pid)
		{
			return vec_para[i];
		}
	}
	throw std::exception(("vec_para not found: center pid: " + std::to_string(_center_pid)).c_str());
}

// sort para order, make para center_pid seq matches particle vector center_pid seq
template <typename T_para>
inline vector<T_para> match_vorocell_para(const vector<Particle>& vec_pa, const vector<T_para>& vec_para)
{
	if (vec_pa.size() != vec_para.size())
	{
		throw std::exception("ERROR: size of particle vector and para vector not match!");
	}

	bool flag_match = true;
	for (size_t i = 0; i < vec_pa.size(); i++)
	{
		flag_match = ((vec_pa[i].id == vec_para[i].center_pid) && flag_match);
	}
	if (flag_match)
	{
		return vec_para;
	}
	else
	{
		cerr << "match_vorocell_para WARNNING: id seq not match: vec_para & vec_pa" << endl;
		vector<T_para> match_vec_para(vec_para.size());
		for (size_t i = 0; i < vec_para.size(); i++)
		{
			match_vec_para[i] = find_para(vec_para, vec_pa[i].id);
		}
		cerr << "match_vorocell_para WARNNING: id seq re-matched: vec_para & vec_pa" << endl;
		return match_vec_para;
	}
}
//vector<Voro_volume> get_voro_volume(const Configuration& config, voro::container_periodic_poly& con);

//vector<Voro_coordination_num> get_voro_coordiantion_num(const voro::container_periodic_poly& con, size_t pid);