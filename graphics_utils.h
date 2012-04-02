#ifndef GRAPHICS_UTILS_H
#define GRAPHICS_UTILS_H

#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cstdlib>
#include "types.h"

extern Point viewerPosition;

void setupGraphics(int,int,Point&);
void quit(int);

#endif // GRAPHICS_UTILS_H
