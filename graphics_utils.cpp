#include "graphics_utils.h"

Point viewerPosition(0,0,0);

void setupGraphics(int width, int height,Point &maxObjCoords) {
    if(SDL_Init(SDL_INIT_VIDEO) < 0) {
        cerr << "Video initialization failed: " << SDL_GetError() << endl;
        quit(1);
    }

    const SDL_VideoInfo* info = NULL;
    info = SDL_GetVideoInfo( );
    if(!info) {
        cerr << "Getting video info failed: " << SDL_GetError() << endl;
        quit(1);
    }

   SDL_GL_SetAttribute(SDL_GL_RED_SIZE,5);
   SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,5);
   SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,5);
   SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,16);
   SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);

   int bitsPerPixel = info->vfmt->BitsPerPixel;
   int flags = SDL_OPENGL | SDL_HWSURFACE | SDL_ASYNCBLIT | SDL_RESIZABLE;
   if(SDL_SetVideoMode(width,height,bitsPerPixel,flags) == 0) {
       cerr << "Setting video mode failed: " << SDL_GetError() << endl;
       quit(1);
   }

    float ratio = (float)width/(float)height;

    glShadeModel(GL_SMOOTH);

   glCullFace(GL_BACK);
   glFrontFace(GL_CCW);
   glEnable(GL_CULL_FACE);

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    /// TODO find max Z coord and set is as last argument in this func
    //gluPerspective(150.0, ratio, 1.0, 5000.0);

    /*double maxObjSize = 3000.0;
    double objHeight = (height > width)? maxObjSize : maxObjSize*height/width;
    double objWidth = (width > height)? maxObjSize : maxObjSize*width/height;*/
    /*glFrustum(-maxObjCoords.x,
              maxObjCoords.x,
              -maxObjCoords.y,
              maxObjCoords.y,
              3*maxObjCoords.z,
              1.0);*/
    gluPerspective(45,ratio,1.0,30*maxObjCoords.z); /////////////// 40
    viewerPosition.z = -20*maxObjCoords.z;
    ///cout << max(viewerPosition.z - maxObjCoords.z,1.0) << endl;
    ///cout << max(viewerPosition.z + maxObjCoords.z,2.0) << endl;
    //glOrtho(-objWidth/2,objWidth/2,-objHeight/2,objHeight/2,1.0,5000.0);

    //tempLine = //new Line(Point(-maxObjCoords.x,0,maxObjCoords.z),Point(maxObjCoords.x,0,maxObjCoords.z));
    //        new Line(viewerPosition,Point(100,12,33));
    //tempVector = new Vector(tempLine->b,tempLine->a);
}

void quit(int code) {
    SDL_Quit();
    exit(code);
}
