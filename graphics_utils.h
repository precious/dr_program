#ifndef GRAPHICS_UTILS_H
#define GRAPHICS_UTILS_H

#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cstdlib>
#include "types.h"
#include <utility>

namespace Graphics {
    extern GLboolean shouldRotate;
    extern Point viewerPosition;

    extern int width;
    extern int height;

    extern bool isLMousePressed;
    extern pair<double,Vector> rotationParams;

    void initGraphics(int,int,Point&);
    void draw(Object3D&,ParticlePolygon*,int);
    void quitGraphics(int);
}

#endif // GRAPHICS_UTILS_H
