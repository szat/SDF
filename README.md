# SDF

This is an octree implementation of the Surface Distance Function, i.e. the penetration distance in the surface normal direction. 

This is a header project, i.e. the octree.h file is the only file needed. Note that it take in data points in the format vector<array<double,3>>, so that the Eigen MatriXd objects in Libigl are not directly usable. There is a small helper function "Eigen2CPP()" in SDF.cpp that converts MatrixXd into vector<array<double,3>>. 

Disclaimer: the following project will only compile if you have libigl to install. For that matter you will need to change the paths in the project sheets. This is only if you want to visualize your results. Otherwise, the octree.h along wih "Eigen2CPP" are enough. 

The code is able to handle arbitrarily big meshes, but I tested only up to 1M poligons. Multithreading should still be implemented, currently only 1 cpu is used. 
