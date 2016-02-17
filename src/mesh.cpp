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
using std::vector;

inline bool compareEqual(const Vector& lhs, const Vector& rhs) {
	double eps = 0.00000001;
	return fabs(lhs.x - rhs.x) < eps && fabs(lhs.y - rhs.y) < eps && fabs(lhs.z - rhs.z) < eps;
}
inline bool compareLess(const Vector& lhs, const Vector& rhs) {
	if (lhs.x != rhs.x){
		if (lhs.x < rhs.x){
			return true;
		}
		return false;
	}
	if (lhs.y != rhs.y){
		if (lhs.y < rhs.y){
			return true;
		}
		return false;
	}
	if (lhs.z != rhs.z){
		if (lhs.z < rhs.z){
			return true;
		}
		return false;
	}
	return false;
}

inline int memberOf(const vector<Vector>& vertices, const Vector& x){
	int verticesSIZE = vertices.size();
	for (int i = 0; i < verticesSIZE; i++){
		if (compareEqual(vertices[i], x)){
			return i;
		}
	}
	return -1;
}
/// Compute the barycentric coordinates of Points in T

void printVertices(const Triangle &t) {
	printf("(%d, %d, %d)\n", t.v[0], t.v[1], t.v[2]);
}

vector<Triangle> Mesh::getNeighbours(const Triangle &currTriangle) const {
	vector<Triangle> neighbours;
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
			neighbours.push_back(t);
		}
	}
	return neighbours;
}

Vector Mesh::getEdgePoint(const Vector &A, const Vector &B, const Vector &C,
	const Vector &D) const{
	double coeff1 = 3.0 / 8, coeff2 = 1.0 / 8;
	Vector shit = coeff1*(A + B) + coeff2*(C + D);
	///printf("shit: %f, %f, %f\n", shit.x, shit.y, shit.z);
	return coeff1*(A + B) + coeff2*(C + D);
}

Vector Mesh::getVertexPoint(const Vector &vertex,
                            vector<Vector> &adjascent) const {
  /// First sort and remove the duplicate adjacent points
	std::sort(adjascent.begin(), adjascent.end(), compareLess);
	adjascent.erase(std::unique(adjascent.begin(), adjascent.end(), compareEqual),
		adjascent.end());
  Vector sumAdj(0, 0, 0);
  size_t adjNum = adjascent.size();
  double coef;
  if (adjNum > 3) {
    coef = 1.0 / adjNum *
           (5.0 / 8 - pow((3.0 / 8 + 1.0 / 4 * cos(2.0 * PI / adjNum)), 2));
  } 
  else if (adjNum == 3){
	coef = 3.0 / 16;
  }
  for (auto &p : adjascent) {
    sumAdj += p;
  }

  return (1 - adjNum * coef) * vertex + coef * sumAdj;
}

void Mesh::computeEdgePoints(const Triangle &currentT,
							 vector<Triangle>& commonSideNeighbours,
                             Vector edgePoints[3]) const {
	size_t commonSIZE = commonSideNeighbours.size();
  for (size_t i = 0; i < commonSIZE; i++) {
    bool found[3] = {false};
	for (auto &vNum : commonSideNeighbours[i].v) {
      if (vNum == currentT.v[0]) {
        vNum = -1;
        found[0] = true;
      } 
	  else if (vNum == currentT.v[1]) {
        vNum = -1;
        found[1] = true;
      } 
	  else if (vNum == currentT.v[2]) {
        vNum = -1;
        found[2] = true;
      }
    }
    size_t vertexNum;
	for (auto &vNum : commonSideNeighbours[i].v) {
      if (vNum != -1) {
        vertexNum = vNum;
      }
    }
    if (found[0] && found[1]) {
      edgePoints[0] =
          getEdgePoint(vertices[currentT.v[0]], vertices[currentT.v[1]],
                       vertices[currentT.v[2]], vertices[vertexNum]);
    } 
	else if (found[1] && found[2]) {
      edgePoints[1] =
          getEdgePoint(vertices[currentT.v[1]], vertices[currentT.v[2]],
                       vertices[currentT.v[0]], vertices[vertexNum]);
    } 
	else if(found[2] && found[0]){
      edgePoints[2] =
          getEdgePoint(vertices[currentT.v[0]], vertices[currentT.v[2]],
                       vertices[currentT.v[1]], vertices[vertexNum]);
    }
  }
}

void Mesh::computeVertexPoints(const vector<Triangle>& neighbours,
                               const Triangle &currentT,
                               vector<Triangle>& commonSideNeighbours,
							   vector<vector<Vector>>& adjacents) const {
  size_t eqCount, commonSideIndex = 0;
  for (auto &neighbourT : neighbours) {
    eqCount = 0;
    /* Do the work of finding the adjacent points
    * for all three vertices of the triangle
    * At the same time take in an array of triangles
    * the triangles which got common side with our */
    for (auto &neighbourVNum : neighbourT.v) {
      for (size_t j = 0; j < 3; j++) {
        if (neighbourVNum == currentT.v[j]) {
          for (auto &neighbourVertex : neighbourT.v) {
            if (neighbourVertex != neighbourVNum) {
				adjacents[j].push_back(vertices[neighbourVertex]);
            }
          }
          eqCount++;
        }
      }
    }
    if (eqCount == 2) {
      commonSideNeighbours.push_back(neighbourT);
    }
  }
}

void Mesh::addNewTriangles(std::vector<Vector>& newVertices,
	std::vector<Triangle>& newTriangles,
	const Triangle& currentT,
	size_t vertexPointsIndices[3],
	size_t edgePointsIndices[3]) const {
    Triangle fourNewTriangles[4];
	/// Add the triangle formed by the three edge points
	for (size_t i = 0; i < 3; i++){
		fourNewTriangles[0].v[i] = edgePointsIndices[i];
	}
	/// Add the remaining 3 triangles formed by 1 vertex Point and 2 Edge
	/// Points
	fourNewTriangles[1].v[0] = edgePointsIndices[0];
	fourNewTriangles[1].v[1] = vertexPointsIndices[1];
	fourNewTriangles[1].v[2] = edgePointsIndices[1];

	fourNewTriangles[2].v[0] = edgePointsIndices[2];
	fourNewTriangles[2].v[1] = edgePointsIndices[1];
	fourNewTriangles[2].v[2] = vertexPointsIndices[2];

	fourNewTriangles[3].v[0] = vertexPointsIndices[0];
	fourNewTriangles[3].v[1] = edgePointsIndices[0];
	fourNewTriangles[3].v[2] = edgePointsIndices[2];
	
	/*for (auto& newTr : fourNewTriangles){
		for (size_t i = 0; i < 3; i++){
			printf("%f, %f, %f\n", newVertices[newTr.v[i]].x, newVertices[newTr.v[i]].y, newVertices[newTr.v[i]].z);
		}
		printf("\n");
	}*/

	for (auto &newT : fourNewTriangles) {
		newT.AB = newVertices[newT.v[1]] - newVertices[newT.v[0]];
		newT.AC = newVertices[newT.v[2]] - newVertices[newT.v[0]];
		newT.ABcrossAC = newT.AB ^ newT.AC;
		newT.gnormal = newT.ABcrossAC;
		newT.gnormal.normalize();
		if (uvs.size() > 0){
			newT.t[0] = currentT.t[0];
			newT.t[1] = currentT.t[1];
			newT.t[2] = currentT.t[2];
		}
		else{
			newT.t[0] = 0;
			newT.t[1] = 0;
			newT.t[2] = 0;
		}
		/*if (uvs.size() != 0 && normals.size() != 0) {
			for (size_t j = 0; j < 3; j++){
				barCoords = getBarCoords(newVertices[newT.v[j]], currentT);
				newUV = barCoords.x * uvs[currentT.t[0]] + barCoords.y * uvs[currentT.t[1]] + barCoords.z * uvs[currentT.t[2]];
				foundUVIndex = memberOf(newUVs, newUV);

				if (foundUVIndex == -1){
					newUVs.push_back(newUV);
					newT.t[j] = newUVsSIZE;
					newUVsSIZE++;
				}
				else{
					newT.t[j] = foundUVIndex;
				}

				newNormal = barCoords.x * normals[currentT.n[0]] + barCoords.y * normals[currentT.n[1]] + barCoords.z * normals[currentT.n[2]];
				newNormal.normalize();
				foundNormalIndex = memberOf(newNormals, newNormal);
				if (foundNormalIndex == -1){
					newNormals.push_back(newNormal);
					newT.n[j] = newNormalsSIZE;
					newNormalsSIZE++;
				}
				else{
					newT.n[j] = foundNormalIndex;
				}
			}
		}*/
		newTriangles.push_back(newT);
	}
}
void Mesh::subdivide() {

  for (size_t step = 0; step < subdivSteps; step++) {
    vector<Triangle> newTriangles;
    vector<Vector> newVertices;
    printf("STEP %d: \n", step + 1);
    /// Iterate every triangle and generate 4 new for each of them
    for (auto &currentT : triangles) {

      /// All the neighbours of the current triangle
      vector<Triangle> neighbours = getNeighbours(currentT);
      /// All the direct adjacent vertices of each point from the triangle
      vector<vector<Vector>> adjacents(3);
      /// All of the neighbours with a common side with the current triangle
      vector<Triangle> commonSideNeighbours;
      /// Collect all the edgepoints in this array
      Vector edgePoints[3];
      /// Compute the Vector Points and the Neighbours with a common side with the current triangle
      computeVertexPoints(neighbours, currentT, commonSideNeighbours,
                          adjacents);
      /// Compute the Edge Points
      computeEdgePoints(currentT, commonSideNeighbours, edgePoints);
      /// Add the new Points
	  size_t newVSIZE = newVertices.size();;
      size_t vertexPointsIndices[3];
	  size_t edgePointsIndices[3];
      Vector vertexPoint, edgePoint;
	  int foundIndex;
	  /// Add the new vertexPoints
      for (size_t j = 0; j < 3; j++) {
        vertexPoint = getVertexPoint(vertices[currentT.v[j]], adjacents[j]);
		foundIndex = memberOf(newVertices, vertexPoint);
        if (foundIndex == -1) {
			newVertices.push_back(vertexPoint);
			vertexPointsIndices[j] = newVSIZE;
			newVSIZE++;
        } 
		else {
			vertexPointsIndices[j] = foundIndex;
        }
      }
	  /// Add the new edge Points
      for (size_t j = 0; j < 3; j++) {
		  foundIndex = memberOf(newVertices, edgePoints[j]);
		  if (foundIndex == -1) {
			  newVertices.push_back(edgePoints[j]);
			  edgePointsIndices[j] = newVSIZE;
			  newVSIZE++;;
		  }
		  else {
			  edgePointsIndices[j] = foundIndex;
		  }
      }
      addNewTriangles(newVertices, newTriangles, currentT,
                      vertexPointsIndices, edgePointsIndices);
    }
    // Change the current mesh's Ts, Vs, uvs and normals with the new one
    triangles = newTriangles;
    vertices = newVertices;
  }
  if (!faceted){
	  autoSmooth = true;
	  normals.clear();
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
