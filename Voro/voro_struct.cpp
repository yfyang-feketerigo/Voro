#include "voro_struct.h"

// compute length of relative vertices vector
void Voro_vertices::compute_relative_r()
{
	if (relative_p.size() % DIMENSION != 0)
	{
		throw "DIMENSION: " + to_string(DIMENSION) + " not match relative position vector size: " + to_string(relative_p.size());
	}
	relative_r.resize(relative_p.size() / DIMENSION);
	for (size_t i = 0; i < relative_r.size(); i++)
	{
		double re_r = 0;
		for (size_t j = 0; j < DIMENSION; j++)
		{
			re_r += std::pow(relative_p[i * DIMENSION + j], 2);
		}
		relative_r[i] = std::sqrt(re_r);
	}

}

// convert Configuration class to voro::container_periodic_poly class
// config: original configuration 
// con: voronoi container to be filled with configuration
void config_to_vorocontainer(const Configuration& config, voro::container_periodic_poly& con)
{
	con.clear();
	for (size_t i = 0; i < config.get_particle().size(); i++)
	{
		const Particle& pa = config.get_particle()[i];
		con.put(pa.id, pa.rx, pa.ry, pa.rz, radius(pa));
	}
}

// compute voronoi & store voronoi info in vector, respectively
void get_voro_info(voro::container_periodic_poly& con, size_t particle_num, vector<Voro_vertices>& vec_vertices, vector<Voro_volume>& vec_volume, vector<Voro_coordination_num>& vec_cn)
{
	// allocate variables & resize vector container 
	voro::voronoicell vc;
	vec_vertices.resize(particle_num);
	vec_volume.resize(particle_num);
	vec_cn.resize(particle_num);
	voro::c_loop_all_periodic con_looper(con); // NOTICE: c_loop_all_periodic DOES NOT loop container as original order
	// if exact order is need, could use:
	// voro::c_loop_order_periodic con_looper(con);

	size_t counter = 0;

	// start loop voronoi container 
	if (con_looper.start()) do
	{
		// excute vornoi compute
		// NOTICE: c_loop_all_periodic DOES NOT loop container as original order
		con.compute_cell(vc, con_looper);

		// pass out vertices relative position
		vec_vertices[counter].center_pid = con_looper.pid();
		vc.vertices(vec_vertices[counter].relative_p);

		// pass out voronoi volume 
		vec_volume[counter].center_pid = con_looper.pid();
		vec_volume[counter].volume = vc.volume();


		// pass out voronoi coordination number
		vec_cn[counter].center_pid = con_looper.pid();
		vec_cn[counter].voro_coordination_num = vc.number_of_faces();

		counter++;
	} while (con_looper.inc());
	if (counter != particle_num)
	{
		cerr << "WARNING: number of voro cells and particle number not match!" << endl;
	}
	return;
}

// write given voronoi info into file as cols
// config: original configuration
// fname: output file name
// other parameters: information to output
// NOTICE: since voro::c_loop_all_periodic is used to looper through voro_container, 
// the particle order will NOT be the same as orginal configuration (at least not always the same).
// Hence, <match_vorocell_para> is used to make particle order match each other.
// There will be some warnning when parameter of particle order reorganize, it's usually normal in this case
void voro_to_dump(const Configuration& config, string fname, const vector<Voro_vertices>& vertices, const vector<Voro_volume>& volume, const vector<Voro_coordination_num>& cn)
{
	// compute length of vertiveces 
	// and prepare output_vector max length of vertiveces
	vector<Voro_vertices> _vertices = match_vorocell_para(config.get_particle(), vertices);
	for (size_t i = 0; i < vertices.size(); i++)
	{
		_vertices[i].compute_relative_r();
	}
	vector<double> output_max_r(vertices.size());
	for (size_t i = 0; i < output_max_r.size(); i++)
	{
		if (config.get_particle()[i].id != _vertices[i].center_pid)
		{
			throw "voro cell center pid squence not match configuration when output max vertex";
		}
		output_max_r[i] = _vertices[i].max_re_r();
	}


	// prepare output_vector of voronoi volume
	vector<Voro_volume> _volume = match_vorocell_para(config.get_particle(), volume);
	vector<double> output_volume(_volume.size());
	for (size_t i = 0; i < output_volume.size(); i++)
	{
		if (config.get_particle()[i].id != _volume[i].center_pid)
		{
			throw "voro cell center pid squence not match configuration when output vec_volume";
		}
		output_volume[i] = _volume[i].volume;
	}

	// prepare output_vector of voronoi coordination number
	vector<Voro_coordination_num> _cn = match_vorocell_para(config.get_particle(), cn);
	vector<double> output_cn(_cn.size());
	for (size_t i = 0; i < cn.size(); i++)
	{
		if (config.get_particle()[i].id != _cn[i].center_pid)
		{
			throw std::exception("voro cell center pid squence not match configuration when output vec_cn");
		}
		output_cn[i] = _cn[i].voro_coordination_num;
	}

	config.para_to_dump(fname, { "max_r","vec_volume" ,"CN_voro" }, { output_max_r, output_volume, output_cn });
	return;
}

/*
const Voro_volume& find_volume(const vector<Voro_volume>& vec_volume, size_t center_pid)
{
	for (size_t i = 0; i < vec_volume.size(); i++)
	{
		if (vec_volume[i].center_pid == center_pid)
		{
			return vec_volume[i];
		}
	}
	throw "vec_volume not found: center pid: " + std::to_string(center_pid);
	return Voro_volume();

}

vector<Voro_volume> sort_voro_volume(const vector<Particle>& vec_pa, const vector<Voro_volume>& vec_volume)
{
	if (vec_pa.size() != vec_volume.size())
	{
		cerr << "WARNING: size of particle vector and vec_volume vector not match when sorting vec_volume" << endl;
	}
	vector<Voro_volume> sort_volume(vec_volume.size());
	for (size_t i = 0; i < vec_volume.size(); i++)
	{
		sort_volume[i] = find_volume(vec_volume, vec_pa[i].id);
	}
	return sort_volume;
}


vector<Voro_vertices> sort_voro_vertex(const vector<Particle>& vec_pa, const vector<Voro_vertices>& vec_vertex)
{
	if (vec_pa.size() != vec_vertex.size())
	{
		cerr << "WARNING: size of particle vector and vertex vector not match when sorting vec_volume" << endl;
	}
	vector<Voro_vertices> sort_vertex(vec_vertex.size());
	for (size_t i = 0; i < vec_vertex.size(); i++)
	{
		sort_vertex[i] = find_vertex(vec_vertex, vec_pa[i].id);
	}
	return sort_vertex;
}


vector<Voro_volume> get_voro_volume(const Configuration& config, voro::container_periodic_poly& con)
{
	voro::voronoicell vc;
	std::vector<Voro_vertices> volume;
	volume.resize(particle_num);
	voro::c_loop_all_periodic looper_equi(con);
	size_t counter = 0;
	if (looper_equi.start()) do
	{
		con.compute_cell(vc, looper_equi);
		volume[counter].center_pid = looper_equi.pid();
		vc.vertices(volume[counter].relative_p);
		counter++;
	} while (looper_equi.inc());
	if (counter != particle_num)
	{
		cerr << "WARNING: number of voro cells and particle number not match!" << endl;
	}
	return volume;
}

const Voro_vertices& find_vertex(const vector<Voro_vertices>& vec_vertex, size_t center_pid)
{
	for (size_t i = 0; i < vec_vertex.size(); i++)
	{
		if (vec_vertex[i].center_pid == center_pid)
		{
			return vec_vertex[i];
		}
	}
	throw "vertex not found: center pid: " + std::to_string(center_pid);
	return Voro_vertices();
}
*/

/*
void compute_voro(voro::container_periodic_poly& con)
{
	voro::voronoicell vc;
	voro::c_loop_all_periodic looper(con);
	if (looper.start()) do
	{
		con.compute_cell(vc, looper);

	} while (looper.inc());
}
*/