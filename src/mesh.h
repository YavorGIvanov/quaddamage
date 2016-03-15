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
 * @File mesh.h
 * @Brief Contains the Mesh class.
 */
#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include "geometry.h"
#include "vector.h"
#include "bbox.h"

using std::vector;

struct KDTreeNode {
	Axis axis; // AXIS_NONE if this is a leaf node
	double splitPos;
	union {
		vector<int>* triangles;
		KDTreeNode* children;
	};
	//
	void initLeaf(const std::vector<int>& triangles)
	{
		axis = AXIS_NONE;
		this->triangles = new vector<int>(triangles);
	}
	
	void initTreeNode(Axis axis, double splitPos)
	{
		this->axis = axis;
		this->splitPos = splitPos;
		this->children = new KDTreeNode[2];
	}
	~KDTreeNode()
	{
		if (axis == AXIS_NONE)
			delete triangles;
		else
			delete [] children;
	}
};

struct Triple {
	Triple(){}
	Triple(Vector v, Vector n, Vector uv) : vertex(v), normal(n), uv(uv){}
	Vector vertex, normal, uv;
	Triple &operator+=(const Triple& other);
	Triple &operator*=(float scalar);
};
bool operator<(const Triple& left, const Triple& right);
bool operator==(const Triple& left, const Triple& right);
Triple operator+(Triple left, const Triple& right);
Triple operator*(double scalar, Triple current);

struct TripleIndices {
	size_t v, n, uv;
};

class Mesh: public Geometry {
	vector<Vector> vertices;
	vector<Vector> normals;
	vector<Vector> uvs;
	vector<Triangle> triangles;
	BBox bbox;
	
	KDTreeNode* kdroot;
	bool useKDTree;
	bool autoSmooth;
	int maxDepthSum;
	int numNodes;

	void computeBoundingGeometry();
	bool intersectTriangle(const RRay& ray, const Triangle& t, IntersectionInfo& info);
	void buildKD(KDTreeNode* node, BBox bbox, const vector<int>& triangleList, int depth);
	bool intersectKD(KDTreeNode* node, const BBox& bbox, const RRay& ray, IntersectionInfo& info);
	void subdivide();
	///Subdivision helper methods:
	friend Triple getEdgeTriple(const Triple& A, const Triple& B, const Triple& C,
							    const Triple& D);
	friend Triple getVertexTriple(const Triple& vertex, vector<Triple>& adjacent);
	vector<vector<Triangle>> getNeighbours(const Triangle &currTriangle) const;
	vector<Triple> computeEdgePoints(const Triangle& currentT,
									 const vector<Triangle>& closeNeighbours) const;
	vector<Triple> computeVertexTriples(const Triangle &currentT,
									    const vector<vector<Triangle>> &neighbours) const;
	void addNewTriangles(vector<Triangle> &newTriangles,
		const vector<Vector> newVertices,
		vector<TripleIndices> vertexIndices,
		vector<TripleIndices> edgeIndices) const;
	///Vector getBarCoords(const Vector& P, const Triangle& T) const;
        ///
public:
	
	bool faceted;
	bool backfaceCulling;
	int subdivSteps;
	Mesh() {
		faceted = false;
		useKDTree = true;
		backfaceCulling = true;
		autoSmooth = false;
		subdivSteps = 0;
		kdroot = NULL;
	}
	~Mesh();

	bool loadFromOBJ(const char* filename);
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getBoolProp("faceted", &faceted);
		pb.getBoolProp("backfaceCulling", &backfaceCulling);
		pb.getIntProp("subdivSteps", &subdivSteps, 0, 10);
		char fn[256];
		if (pb.getFilenameProp("file", fn)) {
			if (!loadFromOBJ(fn)) {
				pb.signalError("Could not parse OBJ file!");
			}
			
		} else {
			pb.requiredProp("file");
		}
		pb.getBoolProp("useKDTree", &useKDTree);
		pb.getBoolProp("autoSmooth", &autoSmooth);
	}
	
	void beginRender();
	
	bool intersect(const Ray& ray, IntersectionInfo& info);
};

#endif // __MESH_H__
