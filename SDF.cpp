
//#include "stdafx.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/igl_inline.h>
#include <igl/per_vertex_normals.h>
#include <igl/ply.h>
#include <igl/readPLY.h>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>

#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <array>

#include "tutorial_shared_path.h"
#include "octree.h"

using namespace Eigen;
using namespace std;

template <class T>
vector<array<T, 3>> Eigen2CPP(Matrix<T, Dynamic, 3> & mat) {
	/*
	This function serves as a bridge between Eigen and this C++ header project.
	Input: an Eigen matrix, that NEEDS to be dynamic rows and with 3 cols. The type is templated.
	Output: a C++ vector of arrays with 3 columns. The type is templated.
	*/
	T* T_ptr = mat.data();
	vector<T> col0(T_ptr, T_ptr + mat.rows());
	vector<T> col1(T_ptr + mat.rows(), T_ptr + 2 * mat.rows());
	vector<T> col2(T_ptr + 2 * mat.rows(), T_ptr + 3 * mat.rows());
	vector<array<T, 3>> out;
	for (size_t i = 0; i < mat.rows(); ++i) {
		out.push_back({ col0.at(i), col1.at(i), col2.at(i) });
	}
	return out;
}

int main(int argc, char *argv[])
{
	//Read data
	Matrix<double, Dynamic, 3> V;
	Matrix<int, Dynamic, 3> F;
	Matrix<double, Dynamic, 3> B;
	Matrix<double, Dynamic, 3> N;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/Alucy.obj", V, F);
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/buddha.obj", V, F);
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/bunny.obj", V, F);
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/bumpy-cube.obj", V, F);
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/arm.obj", V, F);
	//igl::readMESH(TUTORIAL_SHARED_PATH "/octopus-high.mesh", V, temp,F);
	//igl::readMESH(TUTORIAL_SHARED_PATH "/big-sigcat.mesh", V, temp, F);
	//igl::readMESH(TUTORIAL_SHARED_PATH "/hand.mesh", V, temp, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/lion.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);
	igl::barycenter(V, F, B);
	igl::per_vertex_normals(V, F, N);

	cout << "welcome to octreesdf, an octree implementation of the surface diameter function." << endl;
	cout << "this is still under development, todo: very big meshes, visualization." << endl;
	cout << "the mesh has v.rows(): " << V.rows() << ", f.rows(): " << F.rows() << ", b.rows():" << B.rows() << endl;
	cout << endl;

	vector<array<double, 3>> vecB = Eigen2CPP(B);
	vector<array<double, 3>> vecV = Eigen2CPP(V);
	vector<array<double, 3>> vecN = Eigen2CPP(N);
	vector<array<int, 3>> vecF = Eigen2CPP(F);
	for (size_t i = 0; i < vecN.size(); ++i) {
		vecN.at(i)[0] = -vecN.at(i)[0];
		vecN.at(i)[1] = -vecN.at(i)[1];
		vecN.at(i)[2] = -vecN.at(i)[2];
	}

	//Create the search tree
	double total1 = 0;
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	SDF tree(vecV, vecF, vecB);
	tree.build();
	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	cout << "Time taken by build(): " << time_span.count() << "s." << endl << endl;

	t1 = chrono::high_resolution_clock::now();
	vector<double> sdf = tree.query(vecN);
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	cout << "Time taken by query(): " << time_span.count() << "s." << endl << endl;


	//Adjust the values for visualization
	//Raw, the results of the sdf are too small for the visualiser to pick up on them
	double mymax = *max_element(sdf.begin(), sdf.end());
	double mean = 0;
	for (size_t i = 0; i < sdf.size(); ++i) {
		mean += sdf.at(i);
	}
	mean /= sdf.size();
	for (size_t i = 0; i < sdf.size(); ++i) {
		sdf.at(i) = log(sdf.at(i) + mean / mymax); 
	}

	//Creating the visualisation
	double* ptr = &sdf[0];
	Map<VectorXd> Z(ptr, sdf.size());
	MatrixXd C(V.rows(), 3);
	igl::jet(Z, true, C);

	//// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	//// Launch the viewer

	viewer.data().set_colors(C);
	viewer.data().show_lines = !viewer.data().show_lines;
	viewer.launch();
}