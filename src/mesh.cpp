/***************************************************************************
 *   Copyright (C) 2009-2015 by Veselin Georgiev, Slavomir Kaslev et al    *
 *   admin@raytracing-bg.net                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/**
 * @File mesh.cpp
 * @Brief Contains implementation of the Mesh class
 */

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <SDL/SDL.h>
#include "mesh.h"
#include "constants.h"
#include "color.h"

using std::max;
using std::string;

Triple &Triple::operator+=(const Triple &other) {
	vertex += other.vertex;
	normal += other.normal;
	uv += other.uv;
	return *this;
}

Triple &Triple::operator*=(float scalar) {
	vertex *= scalar;
	normal *= scalar;
	uv *= scalar;
	return *this;
}

Triple operator+(Triple left, const Triple &right) {
	left += right;
	return left;
}

Triple operator*(double scalar, Triple current) {
	current *= scalar;
	return current;
}

inline bool compareEqual(const Vector &lhs, const Vector &rhs) {
	double eps = 1e-10;
	return fabs(lhs.x - rhs.x) < eps && fabs(lhs.y - rhs.y) < eps && fabs(lhs.z - rhs.z) < eps;
}
inline bool compareLess(const Vector &lhs, const Vector &rhs) {
	if (lhs.x != rhs.x) {
		return lhs.x < rhs.x;
	}
	if (lhs.y != rhs.y) {
		return lhs.y < rhs.y;
	}
	if (lhs.z != rhs.z) {
		return lhs.z < rhs.z;
	}
	return false;
}

inline void assignIndices(Triangle &t, const TripleIndices &indices0,
                          const TripleIndices &indices1,
                          const TripleIndices &indices2) {
	t.v[0] = indices0.v;
	t.n[0] = indices0.n;
	t.t[0] = indices0.uv;
	t.v[1] = indices1.v;
	t.n[1] = indices1.n;
	t.t[1] = indices1.uv;
	t.v[2] = indices2.v;
	t.n[2] = indices2.n;
	t.t[2] = indices2.uv;
}

bool operator<(const Triple& left, const Triple& right) {
	return compareLess(left.vertex, right.vertex);
}

bool operator==(const Triple& left, const Triple& right){
	return compareEqual(left.vertex, right.vertex);
}

/* This function checks if an element is a member of a vector and returns its
 * index if that is the case. If the element is not a member the function
 * returns -1 */
inline int memberOf(const vector<Vector>& vertices, const Vector& x){
	int verticesSIZE = vertices.size();
	for (int i = 0; i < verticesSIZE; i++){
		if (compareEqual(vertices[i], x)){
			return i;
		}
	}
	return -1;
}

/* The following method adds all elements from vector of triples to
 * three seperate vectors of Vectors. The function adds an element to
 * a vector only if it doesn't already exist in the vector. */
inline vector<TripleIndices> addTo(const vector<Triple> &source,
                                   vector<Vector> &targetV,
                                   vector<Vector> &targetN,
                                   vector<Vector> &targetUV) {
  TripleIndices targetSIZE(targetV.size(), targetN.size(), targetUV.size());
  size_t sourceSIZE = source.size();
  TripleIndices foundIndex;
  vector<TripleIndices> indices;
  for (size_t j = 0; j < sourceSIZE; j++) {
    foundIndex.v = memberOf(targetV, source[j].vertex);
	foundIndex.n = memberOf(targetN, source[j].normal);
	foundIndex.uv = memberOf(targetUV, source[j].uv);
    if (foundIndex.v == -1) {
	  foundIndex.v = targetSIZE.v;
	  targetV.push_back(source[j].vertex);
      targetSIZE.v++;
    }
	if (foundIndex.n == -1) {
		foundIndex.n = targetSIZE.n;
		targetN.push_back(source[j].normal);
		targetSIZE.n++;
	}
	if (foundIndex.uv == -1) {
		foundIndex.uv = targetSIZE.uv;
		targetUV.push_back(source[j].uv);
		targetSIZE.uv++;
	}
	indices.push_back(foundIndex);
  }
  return indices;
}

/* The following function interpolates with coefficents from the loop
 * subdivison scheme four points and returns a new one. Three of the
 * points are from the current triangle and one is from an adjacent one.
 * The function returns a new Triple representing a point on one of the
 * sides of the initial triangle */
inline Triple getEdgeTriple(const Triple &A, const Triple &B, const Triple &C,
	const Triple &D){
	double coef1 = 3.0 / 8, coef2 = 1.0 / 8;
	return coef1*(A + B) + coef2*(C + D);
}

/* This function takes an already existing triple representing a vertex
 * with normal and uv coordinate and a list of triples
 * representing all of its adjacent vertices. The function returns a new
 * shifted vertex of the current triangle. */
inline Triple getVertexTriple(const Triple &vertex,
                              vector<Triple> &adjacent) {
  /// First sort and remove the duplicate adjacent points
	std::sort(adjacent.begin(), adjacent.end());
	adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
  Triple sumAdj;
  size_t adjNum = adjacent.size();
  double coef;
  if (adjNum > 3) {
    coef = 1.0 / adjNum *
           (5.0 / 8 - pow((3.0 / 8 + 1.0 / 4 * cos(2.0 * PI / adjNum)), 2));
  } 
  else if (adjNum == 3){
	coef = 3.0 / 16;
  }
  for (auto &p : adjacent) {
    sumAdj += p;
  }

  return (1 - adjNum * coef) * vertex + coef * sumAdj;
}

/* The following method returns a vector of two lists of triangles by a given current
 * triangle. The first list contains the triangles without a common 
 * side with the current Triangle and the second the ones with. */
vector<vector<Triangle>> Mesh::getNeighbours(const Triangle &currTriangle) const {
	vector<vector<Triangle>> neighbours(2);
	size_t eqCount;
	for (auto &t : triangles) {
		eqCount = 0;
		for (auto &vertexNum : t.v) {
			if (vertexNum == currTriangle.v[0] || vertexNum == currTriangle.v[1] ||
				vertexNum == currTriangle.v[2]) {
				eqCount++;
			}
		}
		if (eqCount == 1 || eqCount == 2) {
			neighbours[eqCount - 1].push_back(t);
		}
	}
	return neighbours;
}

/* The following method finds the needed triples and passes them to the
 * getEdgeTriple function in order to calculate the new edge triples and
 * then return them as vector with size 3*/
vector<Triple> Mesh::computeEdgePoints(const Triangle &currentT,
                        vector<Triangle> &closeNeighbours) const {
  size_t closeSIZE = closeNeighbours.size();
  Triple current;
  vector<Triple> edgePoints(3);
  Triple A(vertices[currentT.v[0]], normals[currentT.n[0]], uvs[currentT.t[0]]),
      B(vertices[currentT.v[1]], normals[currentT.n[1]], uvs[currentT.t[1]]),
      C(vertices[currentT.v[2]], normals[currentT.n[2]], uvs[currentT.t[2]]);
  for(auto &t: closeNeighbours){
    bool found[3] = {false};
    for (auto &vIndex : t.v) {
      if (vIndex == currentT.v[0]) {
        vIndex = -1;
        found[0] = true;
      } 
	  else if (vIndex == currentT.v[1]) {
        vIndex = -1;
        found[1] = true;
      } 
	  else if (vIndex == currentT.v[2]) {
        vIndex = -1;
        found[2] = true;
      }
    }
	for (size_t i = 0; i < 3; i++){
      if (t.v[i] != -1) {
		  current.vertex = vertices[t.v[i]];
		  current.normal = normals[t.n[i]];
		  current.uv = uvs[t.t[i]];
      }
    }
    if (found[0] && found[1]) {
      edgePoints[0] = getEdgeTriple(A, B, C, current);
    } 
	else if (found[1] && found[2]) {
      edgePoints[1] = getEdgeTriple(B, C, A, current);
    } 
	else if (found[2] && found[0]) {
      edgePoints[2] = getEdgeTriple(A, C, B, current);
    }
  }
  return edgePoints;
}

/* Computes and returns the vertex triples which represent shifted vertices
 * with specially chosen coefficents from the loop subdivision scheme */
vector<Triple> Mesh::computeVertexTriples(const Triangle &currentT,
	const vector<vector<Triangle>> &neighbours) const {
	vector<vector<Triple>> adjacents(3);
	Triple current;
	/// Get all adjacent triples to the three vertices of the triangle
	for (size_t i = 0; i < 2; i++) {
		for (auto &t : neighbours[i]) {
			for (size_t j = 0; j < 3; j++) {
				size_t vIndex = currentT.v[j];
				for (auto &neighbIndex: t.v){
					if (vIndex == neighbIndex) {
						for (size_t k = 0; k < 3; k++) {
							if (vIndex != t.v[k]) {
								current.vertex = vertices[t.v[k]];
								current.normal = normals[t.n[k]];
								current.uv = uvs[t.t[k]];
								adjacents[j].push_back(current);
							}
						}
					}
				}
			}
		}
	}
	vector<Triple> vertexTriples(3);
	Triple currTVertex;
	/// Calculate and add the new shifted vertices to the vector
	for (size_t j = 0; j < 3; j++) {
		currTVertex.vertex = vertices[currentT.v[j]];
		currTVertex.normal = normals[currentT.n[j]];
		currTVertex.uv = uvs[currentT.t[j]];
		current = getVertexTriple(currTVertex, adjacents[j]);
		vertexTriples[j] = current;
	}
	return vertexTriples;
}

/// Construct and add the new triangles to the vector newTriangles
void Mesh::addNewTriangles(vector<Triangle> &newTriangles,
						   const vector<Vector> newVertices,
						   vector<TripleIndices> vertexIndices,
                           vector<TripleIndices> edgeIndices) const {
  Triangle newT;
  /// Add the triangle formed by the three edge points
  assignIndices(newT, edgeIndices[0], edgeIndices[1], edgeIndices[2]);
  newTriangles.push_back(newT);
  /// Add the remaining 3 triangles formed by 1 vertex Point and 2 Edge
  /// Points
  assignIndices(newT, edgeIndices[0], vertexIndices[1], edgeIndices[1]);
  newTriangles.push_back(newT);
  assignIndices(newT, edgeIndices[2], edgeIndices[1], vertexIndices[2]);
  newTriangles.push_back(newT);
  assignIndices(newT, vertexIndices[0], edgeIndices[0], edgeIndices[2]);
  newTriangles.push_back(newT);

  size_t tSIZE = newTriangles.size();
  for (size_t i = 1; i <= 4; i++) {
	  newTriangles[tSIZE - i].AB = newVertices[newTriangles[tSIZE - i].v[1]] - newVertices[newTriangles[tSIZE - i].v[0]];
	  newTriangles[tSIZE - i].AC = newVertices[newTriangles[tSIZE - i].v[2]] - newVertices[newTriangles[tSIZE - i].v[0]];
	  newTriangles[tSIZE - i].ABcrossAC = newTriangles[tSIZE - i].AB ^ newTriangles[tSIZE - i].AC;
	  newTriangles[tSIZE - i].gnormal = newTriangles[tSIZE - i].ABcrossAC;
	  newTriangles[tSIZE - i].gnormal.normalize();
  }
}
void Mesh::subdivide() {
  for (size_t step = 0; step < subdivSteps; step++) {
    vector<Triangle> newTriangles;
    vector<Vector> newVertices;
    vector<Vector> newNormals;
    vector<Vector> newUvs;
    printf("STEP %d: \n", step + 1);
    /// Iterate every triangle and generate a maximum 4 new for each of them
    for (auto &currentT : triangles) {
      /// get all the neighbours of the current triangle
      vector<vector<Triangle>> neighbours = getNeighbours(currentT);
      /// Compute the Vertex Points
      vector<Triple> vertexTriples = computeVertexTriples(currentT, neighbours);
      /// Compute the Edge Points
      vector<Triple> edgeTriples = computeEdgePoints(currentT, neighbours[1]);
      /// Add both vertex and edge triples to the appropriate vector
      /// and get their indices
      vector<TripleIndices> vertexIndices =
          addTo(vertexTriples, newVertices, newNormals, newUvs);
      vector<TripleIndices> edgeIndices =
          addTo(edgeTriples, newVertices, newNormals, newUvs);
      /// Add the new triangles
      addNewTriangles(newTriangles, newVertices, vertexIndices, edgeIndices);
    }
    /// Change the current mesh's Ts, Vs, uvs and normals with the new one
    triangles = newTriangles;
    normals = newNormals;
    uvs = newUvs;
    vertices = newVertices;
  }
}

void Mesh::beginRender()
{
	if (subdivSteps > 0){
		subdivide();
	}
	computeBoundingGeometry();
	kdroot = NULL;
	printf("Mesh loaded, %d triangles\n", int(triangles.size()));
	maxDepthSum = 0;
	numNodes = 0;
	if (triangles.size() > 50 && useKDTree) {
		Uint32 startBuild = SDL_GetTicks();
		kdroot = new KDTreeNode;
		vector<int> triangleList(triangles.size());
		std::iota(triangleList.begin(), triangleList.end(), 0);
		buildKD(kdroot, bbox, triangleList, 0);
		Uint32 endBuild = SDL_GetTicks();
		printf(" -> KDTree built in %.2lfs, avg depth = %.1lf\n", (endBuild - startBuild) / 1000.0, maxDepthSum / double(numNodes));
	}
	
	if (normals.size() <= 1 && autoSmooth) {
		normals.resize(vertices.size(), Vector(0, 0, 0)); // extend the normals[] array, and fill with zeros
		for (int i = 0; i < (int) triangles.size(); i++)
			for (int j = 0; j < 3; j++) {
				triangles[i].n[j] = triangles[i].v[j];
				normals[triangles[i].n[j]] += triangles[i].gnormal;
			}
		for (int i = 1; i < (int) normals.size(); i++)
			if (normals[i].lengthSqr() > 1e-9) normals[i].normalize();
	}
	// if the object is set to be smooth-shaded, but it lacks normals, we have to revert it to "faceted":
	if (normals.size() <= 1) faceted = true;
}

void Mesh::buildKD(KDTreeNode* node, BBox bbox, const vector<int>& triangleList, int depth)
{
	if (depth > MAX_TREE_DEPTH || int(triangleList.size()) < TRIANGLES_PER_LEAF) {
		maxDepthSum += depth;
		numNodes++;
		node->initLeaf(triangleList);
		return;
	}
	Axis axis = (Axis) (depth % 3); // TODO: could be better
	double leftLimit = bbox.vmin[axis];
	double rightLimit = bbox.vmax[axis];
	
	double optimalSplitPos = (leftLimit + rightLimit) * 0.5; // TODO: could be MUCH better!
	
	BBox bboxLeft, bboxRight;
	vector<int> trianglesLeft, trianglesRight;
	
	bbox.split(axis, optimalSplitPos, bboxLeft, bboxRight);
	for (auto triangleIdx: triangleList) {
		Triangle& T = this->triangles[triangleIdx];
		const Vector& A = this->vertices[T.v[0]];
		const Vector& B = this->vertices[T.v[1]];
		const Vector& C = this->vertices[T.v[2]];
		
		if (bboxLeft.intersectTriangle(A, B, C))
			trianglesLeft.push_back(triangleIdx);
		
		if (bboxRight.intersectTriangle(A, B, C))
			trianglesRight.push_back(triangleIdx);
	}
	node->initTreeNode(axis, optimalSplitPos);
	buildKD(&node->children[0],  bboxLeft,  trianglesLeft, depth + 1);
	buildKD(&node->children[1], bboxRight, trianglesRight, depth + 1);
}

void Mesh::computeBoundingGeometry()
{
	bbox.makeEmpty();
	
	for (auto& v: vertices) {
		bbox.add(v);
	}
}

Mesh::~Mesh()
{
	if (kdroot) delete kdroot;
}

inline double det(const Vector& a, const Vector& b, const Vector& c)
{
	return (a^b) * c;
}

bool intersectTriangleFast(const Ray& ray, const Vector& A, const Vector& B, const Vector& C, double& dist)
{
	Vector AB = B - A;
	Vector AC = C - A;
	Vector D = -ray.dir;
	//              0               A
	Vector H = ray.start - A;

	/* 2. Solve the equation:
	 *
	 * A + lambda2 * AB + lambda3 * AC = ray.start + gamma * ray.dir
	 *
	 * which can be rearranged as:
	 * lambda2 * AB + lambda3 * AC + gamma * D = ray.start - A
	 *
	 * Which is a linear system of three rows and three unknowns, which we solve using Carmer's rule
	 */

	// Find the determinant of the left part of the equation:
	Vector ABcrossAC = AB ^ AC;
	double Dcr = ABcrossAC * D;
	
	// are the ray and triangle parallel?
	if (fabs(Dcr) < 1e-12) return false;
	
	double lambda2 = ( ( H ^ AC) * D ) / Dcr;
	double lambda3 = ( (AB ^  H) * D ) / Dcr;
	double gamma   = ( ABcrossAC * H ) / Dcr;

	// is intersection behind us, or too far?
	if (gamma < 0 || gamma > dist) return false;

	// is the intersection outside the triangle?
	if (lambda2 < 0 || lambda2 > 1 || lambda3 < 0 || lambda3 > 1 || lambda2 + lambda3 > 1)
		return false;

	dist = gamma;
	
	
	return true;
}

bool Mesh::intersectTriangle(const RRay& ray, const Triangle& t, IntersectionInfo& info)
{
	if (backfaceCulling && dot(ray.dir, t.gnormal) > 0) return false;
	Vector A = vertices[t.v[0]];
	Vector H = ray.start - A;
	Vector D = ray.dir;

	double Dcr = - (t.ABcrossAC * D);

	if (fabs(Dcr) < 1e-12) return false;
	

	double rDcr = 1 / Dcr;
	double gamma = (t.ABcrossAC * H) * rDcr;
	if (gamma < 0 || gamma > info.distance) return false;

	Vector HcrossD = H^D;
	double lambda2 = (HcrossD * t.AC) * rDcr;
	if (lambda2 < 0 || lambda2 > 1) return false;
	
	double lambda3 = -(t.AB * HcrossD) * rDcr;
	if (lambda3 < 0 || lambda3 > 1) return false;
	
	if (lambda2 + lambda3 > 1) return false;
	
	info.distance = gamma;
	info.ip = ray.start + ray.dir * gamma;
	if (!faceted) {
		Vector nA = normals[t.n[0]];
		Vector nB = normals[t.n[1]];
		Vector nC = normals[t.n[2]];
		
		info.normal = nA + (nB - nA) * lambda2 + (nC - nA) * lambda3;
		info.normal.normalize();
	} else {
		info.normal = t.gnormal;
	}
	
	info.dNdx = t.dNdx;
	info.dNdy = t.dNdy;
			
	Vector uvA = uvs[t.t[0]];
	Vector uvB = uvs[t.t[1]];
	Vector uvC = uvs[t.t[2]];
	
	Vector uv = uvA + (uvB - uvA) * lambda2 + (uvC - uvA) * lambda3;
	info.u = uv.x;
	info.v = uv.y;
	info.geom = this;
	return true;
}

bool Mesh::intersectKD(KDTreeNode* node, const BBox& bbox, const RRay& ray, IntersectionInfo& info)
{
	if (node->axis == AXIS_NONE) {
		bool found = false;
		for (int& triIdx: (*node->triangles)) {
			if (intersectTriangle(ray, triangles[triIdx], info))
				found = true;
		}
		return (found && bbox.inside(info.ip));
	} else {
		BBox childBBox[2];
		bbox.split(node->axis, node->splitPos, childBBox[0], childBBox[1]);
		
		int childOrder[2] = { 0, 1 };
		if (ray.start[node->axis] > node->splitPos) {
			std::swap(childOrder[0], childOrder[1]);
		}
		
		BBox& firstBB = childBBox[childOrder[0]];
		BBox& secondBB = childBBox[childOrder[1]];
		KDTreeNode& firstChild = node->children[childOrder[0]];
		KDTreeNode& secondChild = node->children[childOrder[1]];
		// if the ray intersects the common wall between the two sub-boxes, then it invariably
		// intersects both boxes (we can skip the testIntersect() checks):
		// (see http://raytracing-bg.net/?q=node/68 )
		if (bbox.intersectWall(node->axis, node->splitPos, ray)) {
			if (intersectKD(&firstChild, firstBB, ray, info)) return true;
			return intersectKD(&secondChild, secondBB, ray, info);
		} else {
			// if the wall isn't hit, then we intersect exclusively one of the sub-boxes;
			// test one, if the test fails, then it's in the other:
			if (firstBB.testIntersect(ray))
				return intersectKD(&firstChild, firstBB, ray, info);
			else
				return intersectKD(&secondChild, secondBB, ray, info);
		}
		return false;
	}
}

bool Mesh::intersect(const Ray& _ray, IntersectionInfo& info)
{
	RRay ray(_ray);
	ray.prepareForTracing();
	if (!bbox.testIntersect(ray))
		return false;
	
	if (kdroot) {
		info.distance = INF;
		return intersectKD(kdroot, bbox, ray, info);
	} else {
		bool found = false;
		
		info.distance = INF;
		
		for (auto& T: triangles) {
			if (intersectTriangle(ray, T, info)) {
				found = true;
			}
		}
		
		return found;
	}
}

static int toInt(const string& s)
{
	if (s.empty()) return 0;
	int x;
	if (1 == sscanf(s.c_str(), "%d", &x)) return x;
	return 0;
}

static double toDouble(const string& s)
{
	if (s.empty()) return 0;
	double x;
	if (1 == sscanf(s.c_str(), "%lf", &x)) return x;
	return 0;
}

static void parseTrio(string s, int& vertex, int& uv, int& normal)
{
	vector<string> items = split(s, '/');
	// "4" -> {"4"} , "4//5" -> {"4", "", "5" }
	
	vertex = toInt(items[0]);
	uv = items.size() >= 2 ? toInt(items[1]) : 0;
	normal = items.size() >= 3 ? toInt(items[2]) : 0;
}

static Triangle parseTriangle(string s0, string s1, string s2)
{
	// "3", "3/4", "3//5", "3/4/5"  (v/uv/normal)
	Triangle T;
	parseTrio(s0, T.v[0], T.t[0], T.n[0]);
	parseTrio(s1, T.v[1], T.t[1], T.n[1]);
	parseTrio(s2, T.v[2], T.t[2], T.n[2]);
	return T;
}

static void solve2D(Vector A, Vector B, Vector C, double& x, double& y)
{
	// solve: x * A + y * B = C
	double mat[2][2] = { { A.x, B.x }, { A.y, B.y } };
	double h[2] = { C.x, C.y };
	
	double Dcr = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
	x =         (     h[0] * mat[1][1] -      h[1] * mat[0][1]) / Dcr;
	y =         (mat[0][0] *      h[1] - mat[1][0] *      h[0]) / Dcr;
}

bool Mesh::loadFromOBJ(const char* filename)
{
	FILE* f = fopen(filename, "rt");
	
	if (!f) return false;
	
	vertices.push_back(Vector(0, 0, 0));
	uvs.push_back(Vector(0, 0, 0));
	normals.push_back(Vector(0, 0, 0));
	
	char line[10000];
	
	while (fgets(line, sizeof(line), f)) {
		if (line[0] == '#') continue;
	
		vector<string> tokens = tokenize(line); 
		// "v 0 1    4" -> { "v", "0", "1", "4" }
		
		if (tokens.empty()) continue;
		
		if (tokens[0] == "v")
			vertices.push_back(
				Vector(
					toDouble(tokens[1]),
					toDouble(tokens[2]),
					toDouble(tokens[3])));
		
		if (tokens[0] == "vn")
			normals.push_back(
				Vector(
					toDouble(tokens[1]),
					toDouble(tokens[2]),
					toDouble(tokens[3])));

		if (tokens[0] == "vt")
			uvs.push_back(
				Vector(
					toDouble(tokens[1]),
					toDouble(tokens[2]),
					0));
		
		if (tokens[0] == "f") {
			for (int i = 0; i < int(tokens.size()) - 3; i++) {
				triangles.push_back(
					parseTriangle(tokens[1], tokens[2 + i], tokens[3 + i])
				);
			}
		}
	}
	
	fclose(f);
	
	for (auto& t: triangles) {
		Vector A = vertices[t.v[0]];
		Vector B = vertices[t.v[1]];
		Vector C = vertices[t.v[2]];
		Vector AB = B - A;
		Vector AC = C - A;
		t.AB = AB;
		t.AC = AC;
		t.gnormal = t.ABcrossAC = AB ^ AC;
		t.gnormal.normalize();
		
		// (1, 0) = px * texAB + qx * texAC; (1)
		// (0, 1) = py * texAB + qy * texAC; (2)
		
		Vector texA = uvs[t.t[0]];
		Vector texB = uvs[t.t[1]];
		Vector texC = uvs[t.t[2]];
		
		Vector texAB = texB - texA;
		Vector texAC = texC - texA;
		
		double px, py, qx, qy;
		solve2D(texAB, texAC, Vector(1, 0, 0), px, qx); // (1)
		solve2D(texAB, texAC, Vector(0, 1, 0), py, qy); // (2)
		
		t.dNdx = px * AB + qx * AC;
		t.dNdy = py * AB + qy * AC;
		t.dNdx.normalize();
		t.dNdy.normalize();
	}

	return true;
}
