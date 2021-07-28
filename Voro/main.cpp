#include "configuration.h"
#include "voro_struct.h"
#include <voro++/voro++.hh>
#include <iostream>
#include <string>
#include <fstream>
#include <json/json.h>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

int main()
{
	try
	{
		//reading settings from Voro_settings.json
		boost::timer::auto_cpu_timer timer;
		Json::Value root;

		std::clog << "Reading Voro_settings.json..." << endl;

		ifstream f_settings;
		f_settings.open("Voro_settings.json", std::ios_base::binary);
		if (!f_settings.is_open())
		{
			std::cerr << "file \"Voro_settings.json\" open failed" << endl;
			throw "file \"Voro_settings.json\" open failed";
		}
		f_settings >> root;
		//

		//allocate settings variable
		string equi_fname = root["equi_config_fname"].asString();
		size_t end_step = root["end_step"].asLargestUInt();
		size_t delta_step = root["delta_step"].asLargestUInt();
		string fname_prefix = root["fname_prefix"].asString();
		string fname_postfix = root["fname_postfix"].asString();
		string data_fpath = root["data_fpath"].asString();
		string output_path = root["output_path"].asString();
		string equi_data_fpath = root["equi_data_fpath"].asString();
		RADIUS_CONST::type1_radius = root["type1_radius"].asDouble();
		RADIUS_CONST::type2_radius = root["type2_radius"].asDouble();
		//

		//file path check
		boost::filesystem::path boost_path_check(data_fpath);
		if (!(boost::filesystem::exists(boost_path_check) && boost::filesystem::is_directory(boost_path_check)))
		{
			cerr << "data file path: " << boost_path_check << " not exits" << endl;
			throw ("data file path: " + boost_path_check.string() + " not exits");
		}
		boost_path_check = boost::filesystem::path(equi_data_fpath);
		if (!(boost::filesystem::exists(boost_path_check) && boost::filesystem::is_directory(boost_path_check)))
		{
			cerr << "data file path: " << boost_path_check << " not exits" << endl;
			throw ("data file path: " + boost_path_check.string() + " not exits");
		}
		void mkdir(std::string path);
		mkdir(output_path);
		//


		//x string vertices_dir = output_path;
		//x mkdir(vertices_dir);
		//x string vertex_orders_dir = output_path + "vertex_orders_dir/";
		//x mkdir(vertex_orders_dir);
		//x string max_radius_squared = output_path + "max_radius_squared/";
		//x mkdir(max_radius_squared);
		//x string volume = output_path + "Voro_volume/";
		//x mkdir(volume);


		//processing equi(t0) configuration without shear
		Configuration config_equi
		(equi_data_fpath + equi_fname, Configuration::BoxType::orthogonal, Configuration::PairStyle::single);
		double xlo = config_equi.get_xlo();
		double xhi = config_equi.get_xhi();
		double ylo = config_equi.get_ylo();
		double yhi = config_equi.get_yhi();
		double zlo = config_equi.get_zlo();
		double zhi = config_equi.get_zhi();
		double xy = config_equi.get_xy();
		double xz = config_equi.get_xz();
		double yz = config_equi.get_yz();
		double lx = xhi - xlo;
		double ly = yhi - ylo;
		double lz = zhi - zlo;
		int voro_nx = root["voro_nx"].asInt();
		int voro_ny = root["voro_ny"].asInt();
		int voro_nz = root["voro_nz"].asInt();
		int voro_initmem = root["voro_initmem"].asInt();

		clog << "start processing Voro info..." << endl;
		clog << endl;
		voro::container_periodic_poly container_equi(lx, xy, ly, xz, yz, lz, voro_nx, voro_ny, voro_nz, voro_initmem);
		config_to_vorocontainer(config_equi, container_equi);
		std::vector<Voro_vertices> vec_vertices_equi;
		std::vector<Voro_volume> vec_volume_equi;
		std::vector<Voro_coordination_num> vec_cn_equi;
		get_voro_info(container_equi, config_equi.get_particle().size(), vec_vertices_equi, vec_volume_equi, vec_cn_equi);
		for (size_t i = 0; i < vec_vertices_equi.size(); i++)
		{
			vec_vertices_equi[i].compute_relative_r();
		}
		voro_to_dump(config_equi, output_path + "voro.0", vec_vertices_equi, vec_volume_equi, vec_cn_equi);
		container_equi.clear();


		// start time loop, processing shear(t) configuration 
		for (size_t i_step = delta_step; i_step <= end_step; i_step += delta_step)
		{
			string fname = fname_prefix + std::to_string(i_step) + fname_postfix;
			Configuration config_shear(data_fpath + fname, Configuration::BoxType::tilt);
			xlo = config_shear.get_xlo();
			xhi = config_shear.get_xhi();
			ylo = config_shear.get_ylo();
			yhi = config_shear.get_yhi();
			zlo = config_shear.get_zlo();
			zhi = config_shear.get_zhi();
			xy = config_shear.get_xy();
			yz = config_shear.get_yz();
			xz = config_shear.get_xz();
			lx = xhi - xlo;
			ly = yhi - ylo;
			lz = zhi - zlo;
			voro::container_periodic_poly container_shear(lx, xy, ly, xz, yz, lz, voro_nx, voro_ny, voro_nz, voro_initmem);
			config_to_vorocontainer(config_shear, container_shear);
			//container_shear.compute_all_cells();
			std::vector<Voro_vertices> vec_vertices_shear;
			std::vector<Voro_volume> vec_volume_shear;
			std::vector<Voro_coordination_num> vec_cn;
			get_voro_info(container_shear, config_shear.get_particle().size(), vec_vertices_shear, vec_volume_shear, vec_cn);
			string ofname = output_path + "voro." + to_string(i_step);
			voro_to_dump(config_shear, ofname, vec_vertices_shear, vec_volume_shear, vec_cn);
		}
	}
	catch (const string e)
	{
		cout << e << endl;
	}

}

void mkdir(std::string path)
{
	boost::filesystem::path bpath(path);
	if (!boost::filesystem::exists(path))
	{
		boost::filesystem::create_directories(bpath);
	}
}

