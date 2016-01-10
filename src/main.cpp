#include <SDL/SDL.h>
#include <math.h>
#include <vector>
#include "vector.h"
#include "util.h"
#include "sdl.h"
#include "color.h"
#include "camera.h"
#include "geometry.h"
#include "shading.h"
#include "environment.h"
#include "mesh.h"
#include "random_generator.h"
#include "scene.h"
using std::vector;


Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE];

Color raytrace(Ray ray)
{
	if (ray.depth > scene.settings.maxTraceDepth) return Color(0, 0, 0);
	Node* closestNode = NULL;
	double closestDist = INF;
	IntersectionInfo closestInfo;
	for (auto& node: scene.nodes) {
		IntersectionInfo info;
		if (!node->intersect(ray, info)) continue;
		
		if (info.distance < closestDist) {
			closestDist = info.distance;
			closestNode = node;
			closestInfo = info;
		}
	}
	// check if we hit the sky:
	if (closestNode == NULL) {
		if (scene.environment) return scene.environment->getEnvironment(ray.dir);
		else return Color(0, 0, 0);
	} else {
		closestInfo.rayDir = ray.dir;
		if (closestNode->bump)
			closestNode->bump->modifyNormal(closestInfo);
		return closestNode->shader->shade(ray, closestInfo);
	}
}

bool visibilityCheck(const Vector& start, const Vector& end)
{
	Ray ray;
	ray.start = start;
	ray.dir = end - start;
	ray.dir.normalize();
	
	double targetDist = (end - start).length();
	
	for (auto& node: scene.nodes) {
		IntersectionInfo info;
		if (!node->intersect(ray, info)) continue;
		
		if (info.distance < targetDist) {
			return false;
		}
	}
	return true;
}

void debugRayTrace(int x, int y)
{
	Ray ray = scene.camera->getScreenRay(x, y);
	ray.debug = true;
	raytrace(ray);
}

void render()
{
	scene.beginFrame();
	const double kernel[5][2] = {
		{ 0.0, 0.0 },
		{ 0.6, 0.0 },
		{ 0.0, 0.6 },
		{ 0.3, 0.3 },
		{ 0.6, 0.6 },
	};
	vector<Rect> buckets = getBucketsList();
	for (Rect& r: buckets) {
		for (int y = r.y0; y < r.y1; y++)
			for (int x = r.x0; x < r.x1; x++) {
				if (scene.settings.wantAA) {
					Color sum(0, 0, 0);
					for (int i = 0; i < COUNT_OF(kernel); i++)
						sum += raytrace(scene.camera->getScreenRay(x + kernel[i][0], y + kernel[i][1]));
					vfb[y][x] = sum / double(COUNT_OF(kernel));
				} else {
					Ray ray = scene.camera->getScreenRay(x, y);
					vfb[y][x] = raytrace(ray);
				}
			}
		if (!displayVFBRect(r, vfb)) return;
	}
}

int renderSceneThread(void* /*unused*/)
{
	render();
	rendering = false;
	return 0;
}

const char* DEFAULT_SCENE = "data/hw9/dragon.qdmg";

int main ( int argc, char** argv )
{
	initRandom(42);
	Color::init_sRGB_cache();
	const char* sceneFile = argc == 2 ? argv[1] : DEFAULT_SCENE;
	if (!scene.parseScene(sceneFile)) {
		printf("Could not parse the scene!\n");
		return -1;
	}
	
	initGraphics(scene.settings.frameWidth, scene.settings.frameHeight);
	setWindowCaption("Quad Damage: preparing...");
	
	scene.beginRender();
	
	setWindowCaption("Quad Damage: rendering...");
	Uint32 startTicks = SDL_GetTicks();
	renderScene_threaded();
	Uint32 elapsedMs = SDL_GetTicks() - startTicks;
	printf("Render took %.2fs\n", elapsedMs / 1000.0f);
	setWindowCaption("Quad Damage: rendered in %.2fs\n", elapsedMs / 1000.0f);
	
	displayVFB(vfb);
	waitForUserExit();
	closeGraphics();
	printf("Exited cleanly\n");
	return 0;
}
