#include "graphics_utils.h"

Point Graphics::viewerPosition(0,0,0);
int Graphics::width = 0;
int Graphics::height = 0;
float Graphics::zoomFactor = 1.0;
bool Graphics::isLMousePressed = false;
double Graphics::rotationAngles[2] = {0,0};

void Graphics::initGraphics(int _width, int _height,Sphere &satelliteObj) {
    height = _height;
    width = _width;
    if(SDL_Init(SDL_INIT_VIDEO) < 0) {
        cerr << "Video initialization failed: " << SDL_GetError() << endl;
        quitGraphics(1);
    }

    const SDL_VideoInfo* info = NULL;
    info = SDL_GetVideoInfo();
    if(!info) {
        cerr << "getting video info failed: " << SDL_GetError() << endl;
        quitGraphics(1);
    }

    SDL_GL_SetAttribute(SDL_GL_RED_SIZE,5);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,5);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,5);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,16);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);

    int bitsPerPixel = info->vfmt->BitsPerPixel;
    int flags = SDL_OPENGL | SDL_HWSURFACE | SDL_ASYNCBLIT | SDL_RESIZABLE;
    if(SDL_SetVideoMode(width,height,bitsPerPixel,flags) == 0) {
        cerr << "setting video mode failed: " << SDL_GetError() << endl;
        quitGraphics(1);
    }

    glShadeModel(GL_SMOOTH);

    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);
    glEnable(GL_CULL_FACE);

    glViewport(0, 0, width, height);

    glPointSize(2);
}

void Graphics::quitGraphics(int code) {
    SDL_Quit();
    exit(code);
}

void Graphics::draw(Object3D &satelliteObj,Particle* particlesArray = NULL,int particlesNumber = 0)
{
    // colors
    static GLubyte purple[] = {255,150,255,0};
    static GLubyte grey[] = {100,100,100,0};
    static GLubyte red[] = {255,0,0,0};
    static GLubyte green[] = {0,255,0,0};
    static GLubyte blue[] = {0,0,255,0};

    glClearColor(255,255,255,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Projections matrix processing
    static float ratio = (float)width/(float)height;
    static double diameter = satelliteObj.radius*2;
    static GLdouble zNear = 0.0;
    static GLdouble zFar = zNear + 2*diameter;
    static GLdouble left = satelliteObj.center.x - diameter;
    static GLdouble right = satelliteObj.center.x + diameter;
    static GLdouble bottom = satelliteObj.center.y - diameter;
    static GLdouble top = satelliteObj.center.y + diameter;
    static bool firstTimeCall = true;

    if (firstTimeCall) {
        firstTimeCall = false;
        viewerPosition.z = 1.5*diameter;
        if (ratio < 1.0) { // width < height
            bottom /= ratio;
            top /= ratio;
        } else {
            left *= ratio;
            right *= ratio;
        }
    }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left*zoomFactor,right*zoomFactor,bottom*zoomFactor,top*zoomFactor,zNear,zFar);
    //glOrtho(left, right, bottom, top, zNear, zFar);

    // Modelview matrix processing
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    /* Move down the z-axis. */
    //glTranslatef(viewerPosition.x, viewerPosition.y, viewerPosition.z);
    gluLookAt(viewerPosition.x, viewerPosition.y, viewerPosition.z,
              satelliteObj.center.x, satelliteObj.center.y, satelliteObj.center.z,
              0.0, 1.0, 0.0);

    /* Rotation by mouse */
    if (rotationAngles[0])
        glRotated(rotationAngles[0],0,1,0);
    if (rotationAngles[1])
        glRotated(rotationAngles[1],1,0,0);

    // draw axes:
    glBegin(GL_LINES);
    glColor4ubv(red); glVertex3d(0,0,0); glVertex3d(1.5*diameter,0,0);
    glColor4ubv(green); glVertex3d(0,0,0); glVertex3d(0,1.5*diameter,0);
    glColor4ubv(blue); glVertex3d(0,0,0); glVertex3d(0,0,1.5*diameter);
    glEnd();

    // draw the particles
    if (particlesArray != NULL) {
        glBegin(GL_POINTS);
        for(int i = 0;i < particlesNumber;++i) {
            glColor4ubv((particlesArray[i].type == PTYPE_ELECTRON)? blue: red);
            glVertex3f(particlesArray[i].x,particlesArray[i].y,particlesArray[i].z);
        }
        glEnd();
    }

    // draw the object
    glColor4ubv(purple);
    vector<PlaneType> *coords = satelliteObj.polygons;
    for(vector<PlaneType>::iterator it = coords->begin();it != coords->end();it++) {
        glBegin(GL_LINE_LOOP);
        for(int i = 0;i < 3;i++)
            glVertex3d((*it).set[i].x,(*it).set[i].y,(*it).set[i].z);
        glEnd();
    }

    SDL_GL_SwapBuffers();
}
