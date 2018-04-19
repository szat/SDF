#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <array>
#include <cmath>
#include <memory>
#include <thread>

using namespace std;

class SDF {
private:
	vector<array<double, 3>> & V;
	vector<array<int, 3>> & F;
	vector<array<double, 3>> & bary;

	int max_depth;
	int depth;
	array<array<double, 3>, 2> box;
	vector<int> indices;
	array<unique_ptr<SDF>, 8> children;

	bool is_leaf() const;
	void init();
	void build_tree(); //makes the octree on the barycenters
	array<array<double, 3>, 2> build_boxes(); //returns the box of the current node
	array<array<double, 3>, 2> bb_leaf();
	vector<array<double, 3>> query(array<double, 3> & source, array<double, 3> & dir) const;
public:
	SDF(vector<array<double, 3>> & V, vector<array<int, 3>> & F, vector<array<double, 3>> & bary) : V(V), F(F), bary(bary) {};
	SDF(const SDF& other) = delete;
	SDF& operator=(const SDF& rhs) = delete;

	vector<double> query(vector<array<double, 3>> & v_normals) const;

	void build();
};

double dot(const array<double, 3> &A, const array<double, 3> &B) {
	//Returns the dot product, only for a 3 array
	double out = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
	return out;
}

array<double, 3> cross(const array<double, 3> &A, const array<double, 3> &B) {
	//Returns the cross product only for a 3 array
	array<double, 3> out;
	out[0] = A[1] * B[2] - A[2] * B[1];
	out[1] = A[2] * B[0] - A[0] * B[2];
	out[2] = A[0] * B[1] - A[1] * B[0];
	return out;
}

bool RayTriangle(const array<double, 3> & source, const array<double, 3> & dir, const array<array<double, 3>, 3> & tri, array<double, 3> & intersection) {
	/*
	Input: source point of ray, direction vector of ray, triangle defined by 3 points.
	Output: bool of whether there is intersection, intersection point by reference.
	Warning: This routine returns false and empty intersection in case the ray is on the same plane as the triangle and goes through.
	Link: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
	*/
	const double epsilon = 0.00000001;
	array<double, 3> v0 = tri[0];
	array<double, 3> v1 = tri[1];
	array<double, 3> v2 = tri[2];
	array<double, 3> h, s, q;
	double a, f, u, v;
	array<double, 3> edge1 = { v1[0] - v0[0], v1[1] - v0[1] , v1[2] - v0[2] };
	array<double, 3> edge2 = { v2[0] - v0[0], v2[1] - v0[1] , v2[2] - v0[2] };
	h = cross(dir, edge2);
	a = dot(edge1, h);
	if (a > -epsilon && a < epsilon) {
		//DOES NOT HANDLE CASE WHEN LINE IS PARALLEL TO TRIANGLE PLANE
		return false;
	}
	f = 1 / a;
	s = { source[0] - v0[0], source[1] - v0[1], source[2] - v0[2] };
	u = f * dot(s, h);
	if (u < 0.0 - epsilon || u > 1.0 + epsilon)
		return false;
	q = cross(s, edge1);
	v = f * dot(dir, q);
	if (v < 0.0 - epsilon || u + v > 1.0 + epsilon)
		return false;
	// At this stage we can compute t to find out where the intersection point is on the line.
	double t = f * dot(edge2, q);
	if (t > epsilon) { //ray intersection
		intersection[0] = source[0] + t * dir[0];
		intersection[1] = source[1] + t * dir[1];
		intersection[2] = source[2] + t * dir[2];
		return true;
	}
	else
		return false;
}

bool RayBox(const array<double, 3> & source, const array<double, 3> & dir, double low_t, double high_t, const array<array<double, 3>, 2> & box) {
	/*
	Input: source point of ray, direction vector of ray, interval for the scalar factor for direction vector, box defined by diagonal points low and high.
	Output: bool of whether there is intersection.
	Link: https://www.cs.utah.edu/~awilliam/box/box.pdf
	*/
	array<double, 3> inv_dir = { 1 / dir[0], 1 / dir[1], 1 / dir[2] };
	array<int, 3> sign;
	sign[0] = inv_dir[0] < 0;
	sign[1] = inv_dir[1] < 0;
	sign[2] = inv_dir[2] < 0;
	double tmin, tmax, tymin, tymax, tzmin, tzmax;
	tmin = (box[sign[0]][0] - source[0]) * inv_dir[0];
	tmax = (box[1 - sign[0]][0] - source[0]) * inv_dir[0];
	tymin = (box[sign[1]][1] - source[1]) * inv_dir[1];
	tymax = (box[1 - sign[1]][1] - source[1]) * inv_dir[1];
	if ((tmin > tymax) || (tymin > tmax))
		return false;
	if (tymin > tmin)
		tmin = tymin;
	if (tymax < tmax)
		tmax = tymax;
	tzmin = (box[sign[2]][2] - source[2]) * inv_dir[2];
	tzmax = (box[1 - sign[2]][2] - source[2]) * inv_dir[2];
	if ((tmin > tzmax) || (tzmin > tmax))
		return false;
	if (tzmin > tmin)
		tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;
	return ((tmin < high_t) && (tmax > low_t));
}

bool TriangleBox(const array<array<double, 3>, 3> & triangle, const array<array<double, 3>, 2> & box) {
	/*
	Input: triangle defined by 3 points, box defined by diagonal points low and high.
	Output: bool of whether there is intersection.
	*/
	//There is an intersection iff one of the points of the triangle is not in the box
	for (size_t i = 0; i < 3; ++i) {
		array<double, 3> pt = triangle[i];
		if (box[0][0] > pt[0] || pt[0] > box[1][0]) {
			return false;
		}
		if (box[0][1] > pt[1] || pt[1] > box[1][1]) {
			return false;
		}
		if (box[0][2] > pt[2] || pt[2] > box[1][2]) {
			return false;
		}
	}
	return true;
}

void SDF::build() {
	this->init();
	this->build_tree();
	this->build_boxes();
	cout << "idx size: " << this->indices.size() << endl;
}

void SDF::init() {
	for (size_t i = 0; i < this->bary.size(); ++i) {
		this->indices.push_back((int)i);
	}
	this->max_depth = 1 + log(this->bary.size()) / log(8);
	this->depth = 0;
	return;
}

void SDF::build_tree() {
	//Top down approach
	int nb = this->indices.size();
	//cout << "At depth " << this->depth << ", nb of pts is " << nb << endl;

	if (nb <= 1) {
		return;
	}
	if (this->depth == this->max_depth) {
		return;
	}

	array<double, 3> mean;
	mean[0] = bary.at(this->indices.at(0))[0];
	mean[1] = bary.at(this->indices.at(0))[1];
	mean[2] = bary.at(this->indices.at(0))[2];
	for (size_t i = 1; i < nb; ++i) {
		mean[0] = mean[0] + (bary.at(this->indices.at(i))[0] - mean[0]) / (i + 1);
		mean[1] = mean[1] + (bary.at(this->indices.at(i))[1] - mean[1]) / (i + 1);
		mean[2] = mean[2] + (bary.at(this->indices.at(i))[2] - mean[2]) / (i + 1);
	}

	array<vector<int>, 8> sub_indices;
	for (size_t i = 0; i < this->indices.size(); ++i) {
		array<double, 3> pt = bary.at(indices.at(i));
		if (pt[0] > mean[0]) {
			if (pt[1] > mean[1]) {
				if (pt[2] > mean[2]) {
					sub_indices[0].push_back(indices.at(i));
				}
				else {//pt[2] <= mean[2] 
					sub_indices[1].push_back(indices.at(i));
				}
			}
			else {   //pt[1] <= mean[1]
				if (pt[2] > mean[2]) {
					sub_indices[2].push_back(indices.at(i));
				}
				else {//pt[2] <= mean[2]
					sub_indices[3].push_back(indices.at(i));
				}
			}
		}
		else {		 //pt[0] <= mean[0]
			if (pt[1] > mean[1]) {
				if (pt[2] > mean[2]) {
					sub_indices[4].push_back(indices.at(i));
				}
				else {//pt[2] <= mean[2]
					sub_indices[5].push_back(indices.at(i));
				}
			}
			else {   //pt[1] <= mean[1]
				if (pt[2] > mean[2]) {
					sub_indices[6].push_back(indices.at(i));
				}
				else //pt[2] <= mean[2]
					sub_indices[7].push_back(indices.at(i));
			}
		}
	}
	indices.clear();

	for (size_t i = 0; i < 8; ++i) {
		if (sub_indices[i].size() > 0) {
			this->children[i] = unique_ptr<SDF>(new SDF(this->V, this->F, this->bary));
			this->children[i]->indices = sub_indices[i];
			sub_indices[i].clear();
			this->children[i]->max_depth = this->max_depth;
			this->children[i]->depth = this->depth + 1;
		}
		else {
			this->children[i] = nullptr;
			sub_indices[i].clear();
		}
	}

	//Recurse on children
	for (size_t i = 0; i < 8; ++i) {
		if (this->children[i] != nullptr) {
			this->children[i]->build_tree();
		}
	}
	return;
}

bool SDF::is_leaf() const {
	bool all_null = true;
	for (size_t i = 0; i < 8 && all_null; ++i) {
		if (this->children[i] != nullptr) {
			all_null = false;
		}
	}
	return all_null;
}

array<array<double, 3>, 2> SDF::bb_leaf() {
	array<array<double, 3>, 2> bb;
	array<int, 3> face = this->F.at(this->indices.at(0));
	bb[0] = this->V.at(face[0]);
	bb[1] = this->V.at(face[0]);
	for (size_t f = 0; f < this->indices.size(); ++f) {
		face = F.at(this->indices.at(f));
		for (size_t v = 0; v < 3; ++v) {
			array<double, 3> vertex = this->V.at(face[v]);
			for (size_t dim = 0; dim < 3; ++dim) {
				if (bb[0][dim] > vertex[dim]) { //low
					bb[0][dim] = vertex[dim];
				}
				if (bb[1][dim] < vertex[dim]) { //high
					bb[1][dim] = vertex[dim];
				}
			}
		}
	}
	return bb;
}

array<array<double, 3>, 2> SDF::build_boxes() {
	//Bottom up approach
	array<array<double, 3>, 2> bb;
	if (this->is_leaf() == true) {
		//trying to catch errors
		if (this->indices.size() == 0) { return bb; }
		else {
			bb = this->bb_leaf();
			this->box = bb;
			return bb;
		}
	}
	else {
		bool premier = true;
		for (size_t i = 0; i < 8; ++i) {
			if (this->children[i] != nullptr && premier) {
				premier = false;
				bb = this->children[i]->build_boxes();
			}
			if (this->children[i] != nullptr && !premier) {
				array<array<double, 3>, 2> temp = this->children[i]->build_boxes();
				for (size_t dim = 0; dim < 3; ++dim) {
					if (bb[0][dim] > temp[0][dim]) { //low
						bb[0][dim] = temp[0][dim];
					}
					if (bb[1][dim] < temp[1][dim]) { //high
						bb[1][dim] = temp[1][dim];
					}
				}
			}
		}
		this->box = bb;
		return bb;
	}
}

vector<array<double, 3>> SDF::query(array<double, 3> & source, array<double, 3> & dir) const {
	if (this->is_leaf() == true) { //now a leaf has a set of triangles
		array<double, 3> intersection;
		vector<array<double, 3>> out;
		array<array<double, 3>, 3> tri;

		for (size_t i = 0; i < this->indices.size(); ++i) {
			int tri_idx = (this->indices).at(i);
			array<int, 3> pts_idx = (this->F).at(tri_idx);
			tri[0] = (this->V).at(pts_idx[0]);
			tri[1] = (this->V).at(pts_idx[1]);
			tri[2] = (this->V).at(pts_idx[2]);
			bool is_hit = RayTriangle(source, dir, tri, intersection);
			if (is_hit == true) {
				out.push_back(intersection);
			}
		}
		return out;
	}
	else { //if this is not a leaf, query until you find the leaf and then concatenate up
		vector<array<double, 3>> out;
		for (size_t i = 0; i < 8; ++i) {
			if (this->children[i] != nullptr) {
				if (true == RayBox(source, dir, 0, DBL_MAX, this->children[i]->box)) { //leap of faith
					vector<array<double, 3>> temp = this->children[i]->query(source, dir);
					out.insert(end(out), begin(temp), end(temp));
				}
			}
		}
		return out; //may be empty
	}
}

vector<double> SDF::query(vector<array<double, 3>> & inv_normals) const {
	vector<double> sdf;
	for (size_t i = 0; i < inv_normals.size(); ++i) {
		vector<array<double, 3>> intersect = this->query((this->V).at(i), inv_normals.at(i));
		double my_min = DBL_MAX;
		for (size_t j = 0; j < intersect.size(); ++j) {
			array<double, 3> diff = { intersect.at(j)[0] - (this->V).at(i)[0], intersect.at(j)[1] - (this->V).at(i)[1],intersect.at(j)[2] - (this->V).at(i)[2] };
			double norm_sq = dot(diff, diff);
			if (my_min > norm_sq && norm_sq > 0.00000005) {
				my_min = norm_sq;
			}
		}
		if (my_min != DBL_MAX) {
			sdf.push_back(my_min);
		}
		else {
			sdf.push_back(0); //in case there was "no interection"
		}
	}
	return sdf;
};