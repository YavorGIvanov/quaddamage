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
 * @File shading.cpp
 * @Brief Contains implementations of shader classes
 */
#include "shading.h"
#include "bitmap.h"
#include "constants.h"
#include <iostream>

bool visibilityCheck(const Vector &start, const Vector &end);

extern Vector lightPos;
extern double lightIntensity;
extern Color ambientLight;

void MandelbrotFract::setUpPalitre(){
	double divisor;
	for (unsigned i = 0; i < palitreSIZE; i++){
		divisor = (palitreSIZE - i)/3.6;
		palitre[i] = color2 / divisor;
	}
	isSetPalitre = true;
}

Color MandelbrotFract::sample(const IntersectionInfo &info) {
	
	double cReal = (-2.5 + (info.u * 3.5)); //< in the interval (-2.5,1)
	double cImag = (1 + (info.v * 2.0)); //< in the interval (-1, 1)
	double sc = log(scaling) / 1.5;
	cReal *= sc;
	cImag *= sc;
	unsigned MaxIterations = 50, iter = 0; //< 50 is a decent value
	double zReal = 0, zImag = 0, zRealSq = 0, zImagSq = 0;
	//Setting up the palitre for the color around the fractal
	if (!isSetPalitre){
		setUpPalitre();
	}

  while (iter < MaxIterations && zRealSq + zImagSq < 4.0) {
    zRealSq = zReal * zReal;
    zImagSq = zImag * zImag;
    zImag = 2.0 * zReal * zImag + cImag;
    zReal = zRealSq - zImagSq + cReal;
    iter++;
  }
  if (iter == MaxIterations){
	  return color1;
  }
  unsigned index = (unsigned)ceil(((double)iter / MaxIterations)*palitreSIZE) % palitreSIZE;
  return palitre[index];
}

Color CheckerTexture::sample(const IntersectionInfo& info)
{
	int x = (int) floor(info.u * scaling / 5.0);
	int y = (int) floor(info.v * scaling / 5.0);
	
	Color checkerColor = ((x + y) % 2 == 0) ? color1 : color2;
	return checkerColor;
}

double getLightContrib(const IntersectionInfo& info)
{
	double distanceToLightSqr = (info.ip - lightPos).lengthSqr();

	if (!visibilityCheck(info.ip + info.normal * 1e-6, lightPos)) {
		return 0;
	} else {
		return lightIntensity / distanceToLightSqr;
	}
}

Color Lambert::shade(const Ray& ray, const IntersectionInfo& info)
{
	Color diffuse = texture ? texture->sample(info) : this->color;
	
	Vector v1 = faceforward(ray.dir, info.normal);
	Vector v2 = lightPos - info.ip;
	v2.normalize();
	double lambertCoeff = dot(v1, v2);
	
	return ambientLight * diffuse
		+ diffuse * lambertCoeff * getLightContrib(info);
	
}

Color Phong::shade(const Ray& ray, const IntersectionInfo& info)
{
	Color diffuse = texture ? texture->sample(info) : this->color;
	
	Vector v1 = faceforward(ray.dir, info.normal);
	Vector v2 = lightPos - info.ip;
	v2.normalize();
	double lambertCoeff = dot(v1, v2);
	double fromLight = getLightContrib(info);
	
	Vector r = reflect(info.ip - lightPos, info.normal);
	Vector toCamera = -ray.dir;
	double cosGamma = dot(toCamera, r);
	double phongCoeff;
	if (cosGamma > 0)
		phongCoeff = pow(cosGamma, specularExponent);
	else
		phongCoeff = 0;
	
	return ambientLight * diffuse
		+ diffuse * lambertCoeff * fromLight
		+ Color(1, 1, 1) * (phongCoeff * specularMultiplier * fromLight);
}

BitmapTexture::BitmapTexture(const char* filename, double scaling)
{
	bitmap = new Bitmap;
	bitmap->loadImage(filename);
	this->scaling = 1/scaling;
}

BitmapTexture::~BitmapTexture() { delete bitmap; }

Color BitmapTexture::sample(const IntersectionInfo& info)
{
	int x = (int) floor(info.u * scaling * bitmap->getWidth());
	int y = (int) floor(info.v * scaling * bitmap->getHeight());
	// 0 <= x < bitmap.width
	// 0 <= y < bitmap.height
	x = (x % bitmap->getWidth());
	y = (y % bitmap->getHeight());
	if (x < 0) x += bitmap->getWidth();
	if (y < 0) y += bitmap->getHeight();
	
	return bitmap->getPixel(x, y);
}

extern Color raytrace(Ray ray);

Color Refl::shade(const Ray& ray, const IntersectionInfo& info)
{
	Vector n = faceforward(ray.dir, info.normal);

	if (glossiness == 1) {
		Ray newRay = ray;
		newRay.start = info.ip + n * 0.000001;
		newRay.dir = reflect(ray.dir, n);
		newRay.depth++; 
		
		return raytrace(newRay) * multiplier;
	} else {
		Color result(0, 0, 0);
		int count = numSamples;
		if (ray.depth > 0)
			count = 2;
		for (int i = 0; i < count; i++) {
			Vector a, b;
			orthonormalSystem(n, a, b);
			double angle = randomFloat() * 2 * PI;
			double radius = randomFloat() * 1;
			double x, y;
			x = cos(angle) * radius;
			y = sin(angle) * radius;
			//
			x *= tan((1 - glossiness) * PI/2);
			y *= tan((1 - glossiness) * PI/2);
			
			Vector modifiedNormal = n + a * x + b * y;

			Ray newRay = ray;
			newRay.start = info.ip + n * 0.000001;
			newRay.dir = reflect(ray.dir, modifiedNormal);
			newRay.depth++; 
			
			result += raytrace(newRay) * multiplier;
		}
		return result / count;
	}
}

inline Vector refract(const Vector& i, const Vector& n, double ior)
{
	double NdotI = (double) (i * n);
	double k = 1 - (ior * ior) * (1 - NdotI * NdotI);
	if (k < 0.0)		// Check for total inner reflection
		return Vector(0, 0, 0);
	return ior * i - (ior * NdotI + sqrt(k)) * n;
}

Color Refr::shade(const Ray& ray, const IntersectionInfo& info)
{
// ior = eta2 / eta1
	Vector refr;
	if (dot(ray.dir, info.normal) < 0) {
		// entering the geometry
		refr = refract(ray.dir, info.normal, 1 / ior_ratio);
	} else {
		// leaving the geometry
		refr = refract(ray.dir, -info.normal, ior_ratio);
	}
	if (refr.lengthSqr() == 0) return Color(1, 0, 0);
	Ray newRay = ray;
	newRay.start = info.ip - faceforward(ray.dir, info.normal) * 0.000001;
	newRay.dir = refr;
	newRay.depth++;
	return raytrace(newRay) * multiplier;
}

void Layered::addLayer(Shader* shader, Color blend, Texture* tex)
{
	if (numLayers == COUNT_OF(layers)) return;
	layers[numLayers++] = { shader, blend, tex };
}

Color Layered::shade(const Ray& ray, const IntersectionInfo& info)
{
	Color result(0, 0, 0);
	for (int i = 0; i < numLayers; i++) {
		Color fromLayer = layers[i].shader->shade(ray, info);
		Color blendAmount = layers[i].blend;
		if (layers[i].tex) blendAmount = blendAmount * layers[i].tex->sample(info);
		result = blendAmount * fromLayer + (Color(1, 1, 1) - blendAmount) * result;
	}
	return result;
}

inline double fresnel(const Vector& i, const Vector& n, double ior)
{
	// Schlick's approximation
	double f = sqr((1.0f - ior) / (1.0f + ior));
	double NdotI = (double) -dot(n, i);
	return f + (1.0f - f) * pow(1.0f - NdotI, 5.0f);
}

Color Fresnel::sample(const IntersectionInfo& info)
{
	double eta = ior;
	if (dot(info.normal, info.rayDir) > 0)
		eta = 1 / eta;
	Vector n = faceforward(info.rayDir, info.normal);
	double fr = fresnel(info.rayDir, n, eta);
	return Color(fr, fr, fr);
}
