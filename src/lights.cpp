
#include "lights.h"
#include "util.h"
#include "random_generator.h"

void RectLight::beginFrame(void)
{
	center = T.point(Vector(0, 0, 0));
	a = T.point(Vector(-0.5, 0.0, -0.5));
	b = T.point(Vector( 0.5, 0.0, -0.5));
	c = T.point(Vector( 0.5, 0.0,  0.5));
	d = T.point(Vector(-0.5, 0.0,  0.5));
	float width = (float) (b - a).length();
	float height = (float) (b - c).length();
	area = width * height; // obtain the area of the light, in world space
}


int RectLight::getNumSamples()
{
	return xSubd * ySubd;
}

void RectLight::getNthSample(int sampleIdx, const Vector& shadePos,
						  Vector& samplePos, Color& color)
{
	Random& rnd = getRandomGen();
	double x = (sampleIdx % xSubd + rnd.randfloat()) / xSubd;
	double y = (sampleIdx / xSubd + rnd.randfloat()) / ySubd;
	
	samplePos = Vector(x - 0.5, 0, y - 0.5);
	
	Vector shadePos_LS = T.undoPoint(shadePos);
	
	if (shadePos_LS.y < 0) {
		float cosWeight = float(dot(Vector(0, -1, 0), shadePos_LS) / shadePos_LS.length());
		color = this->color * power * area * cosWeight;
	} else {
		color.makeZero();
	}
	
	samplePos = T.point(samplePos);
}
	
bool RectLight::intersect(const Ray& ray, double& intersectionDist)
{
	Ray ray_LS = T.undoRay(ray);
	// check if ray_LS (the incoming ray, transformed in local space) hits the oriented square 1x1, resting
	// at (0, 0, 0), pointing downwards:
	if (ray_LS.start.y >= 0) return false; // ray start is in the wrong subspace; no intersection is possible
	if (ray_LS.dir.y <= 0) return false; // ray direction points downwards; no intersection is possible
	double lengthToIntersection = -(ray_LS.start.y / ray_LS.dir.y); // intersect with XZ plane
	Vector p = ray_LS.start + ray_LS.dir * lengthToIntersection;
	if (fabs(p.x) < 0.5 && fabs(p.z) < 0.5) {
		// the hit point is inside the 1x1 square - calculate the length to the intersection:
		double distance = (T.point(p) - ray.start).length(); 
		
		if (distance < intersectionDist) {
			intersectionDist = distance;
			return true; // intersection found, and it improves the current closest dist
		}
	}
	return false;
}

float RectLight::solidAngle(const Vector& x)
{
	///Find the projections
	Vector prA, prB, prC, prD;
	prA = (a - x);
	prA.normalize();
	prB = (b - x);
	prB.normalize();
	prC = (c - x);
	prC.normalize();
	prD = (d - x);
	prD.normalize();

	///Find the normal vectors of the planes
	Vector n1, n2, n3, n4;
	n1 = prB^prA;
	n1.normalize();
	n2 = prC^prB;
	n2.normalize();
	n3 = prD^prC;
	n3.normalize();
	n4 = prA^prD;
	n4.normalize();

	///Find the dihedral angles
	float angle1, angle2, angle3, angle4;
	angle1 = acos(dot(n1, n2));
	angle2 = acos(dot(n2, n3));
	angle3 = acos(dot(n3, n4));
	angle4 = acos(dot(n4, n1));

	///Calculate the area
	return fabs(angle1 + angle2 + angle3 + angle4 - 2*PI);

	
}

