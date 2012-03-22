#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <unistd.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <errno.h>

#include "types.h"
#include "read_file.h"
#include "geometry_utils.h"
#include "constants.h"

#define exitErr(msg) { cerr << msg << "\nerrno: " << errno << endl; exit(1); }
#define rand(max) rand()%max


static GLboolean should_rotate = GL_FALSE;
Line *tempLine;
Vector *tempVector;
int step = 30;

Point viewerPosition(0,0,0);

const char usage[] = "Usage:\n\t%s <filename>\n";

void quit(int code) {
    SDL_Quit();
    exit(code);
}

static void handleKeyDown(SDL_keysym* keysym)
{
    switch(keysym->sym) {
    case SDLK_ESCAPE:
        quit(0);
        break;
    case SDLK_SPACE:
        should_rotate = !should_rotate;
        break;

    case SDLK_KP_PLUS:
    case SDLK_PLUS:
        //viewerPosition.z += 20;
        tempLine->a().z += step;
        tempLine->b().z += step;
        break;
    case SDLK_KP_MINUS:
    case SDLK_MINUS:
        //viewerPosition.z -= 20;
        tempLine->a().z -= step;
        tempLine->b().z -= step;
        break;
    case SDLK_UP:
        cout << "up" << endl;
        tempLine->a().y += step;
        tempLine->b().y += step;
        break;
    case SDLK_DOWN:
        cout << "down" << endl;
        tempLine->a().y -= step;
        tempLine->b().y -= step;
        break;
    case SDLK_LEFT:
        cout << "left" << endl;
        tempLine->a().x -= step;
        tempLine->b().x -= step;
        break;
    case SDLK_RIGHT:
        cout << "right" << endl;
        tempLine->a().x += step;
        tempLine->b().x += step;
        break;
    default:
        break;
    }
}


void processEvents(void)
{
    /* Our SDL event placeholder. */
    SDL_Event event;

    /* Grab all the events off the queue. */
    while(SDL_PollEvent(&event)) {

        switch (event.type) {
        case SDL_KEYDOWN:
            handleKeyDown( &event.key.keysym );
            break;
        case SDL_QUIT:
            /* Handle quit requests (like Ctrl-c). */
            quit(0);
            break;
        }
    }
}


int processParticles(Object3D* satelliteObj) {
    real stepLength = Vector(satelliteObj->nearestPoint,
                             satelliteObj->furthermostPoint).length()/10.0;
    real timeInterval = stepLength/(2*ELECTRON_VELOCITY);
    cout << "stepLength: " << stepLength << endl;
    cout << "timeInterval: " << timeInterval << endl;

    // static bool firstTime = true;

    // static vector<Particle> particles(count);

    GenerativeSphere generativeSphere(satelliteObj->center(),
                                      10*Vector(satelliteObj->center(),satelliteObj->maxCoords).length(),
                                      satelliteObj->front.normalize()*satelliteObj->speed);
    Particle particle;
    int count = 10000;
    int numOfItersections = 0;
    for (int j = 0;j < count;j++) {
        particle = generativeSphere.generateParticle(PTYPE_ELECTRON);
        if (GeometryUtils::doesParticlesTrajectoryIntersectObject(particle,*satelliteObj)) {
            ++numOfItersections;
        }
    }

    cout << "coefficient: " << float(numOfItersections)/count << endl;

    /*if (true) { ////////////////////////////////////////////////////////////////////////////!!!!!!!!!!!!!!!
        firstTime = false;
        for(int i = 0;i < count;i++) {

            //cout << particles[i].step << " [" << i << "] " << endl;//////////////////////////////////////////
            //cout << particles[i] << " [" << i << "] " << endl;//////////////////////////////////////////
        }
    }/* else {
        for(int i = 0;i < count;i++)
            particles[i] = particles[i] + particles[i].step;
    }*/

    ////////////////////////////////////////////////////////////
    //cout << particles[0] << " [0] " << endl;//.step = Vector(particles[0],satelliteObj->center()).normalize()*particles[0].step.length();
    //cout << "normal: " << generativeSurface.normal << endl;/////////////////////////////
    ////////////////////////////////////////////////////////////

    /*Particle fastestParticle = GeometryUtils::getFastestParticle(particles,satelliteObj->center(),
                                [satelliteObj](Particle p) -> bool {
                                    return GeometryUtils::doesParticlesTrajectoryIntersectObject(p,*satelliteObj);
                                });
    if (fastestParticle == Point()) {
        return 0;
    }
    else
        return 1;
    cout << "fastest: " << fastestParticle << endl;
    cout << "center: " << satelliteObj->center() << endl;
    real time = Vector(fastestParticle,satelliteObj->center()).length() / fastestParticle.step.length();
    cout << "time: " << time << endl;*/

    return 0;
}


void draw(Object3D* satelliteObj)
{
    static float angle = 0.0f;
    static GLubyte purple[] = {255,   230, 255,   0 };
    static GLubyte grey[] = {200,200,200,0};
    static GLubyte black[] = {0,0,0,0};
    static GLubyte blue[] = {0,0,255,0};

    /*
    lengthGL lengthReal
            x = stepGL
*/


    glClearColor(255,255,255,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    /* Move down the z-axis. */
    glTranslatef(viewerPosition.x, viewerPosition.y, viewerPosition.z);

    /* Rotate. */
    glRotatef(angle, 0.0, 1.0, 0.0);

    if (should_rotate) {
        if (++angle > 360.0f) {
            angle = 0.0f;
        }
    }


    vector<PlaneType> *coords = satelliteObj->polygons;
    for(vector<PlaneType>::iterator it = coords->begin();it != coords->end();it++) {
        glBegin(GL_LINE_LOOP);
        glColor4ubv(purple);
        for(int i = 0;i < 3;i++)
            glVertex3d((*it).set[i].x,(*it).set[i].y,(*it).set[i].z);
    /*    if (tempVector->cos((*it).normal) < 0)
            continue;
        if (GeometryUtils::doesLineIntersectTriangle(*it,*tempLine)) {
            glBegin(GL_TRIANGLES);
            glColor4ubv(blue);
        } else {
            glBegin(GL_LINE_LOOP);
            glColor4ubv(purple);
        }
        for(int i = 0;i < 3;i++)
            glVertex3d((*it).set[i].x,(*it).set[i].y,(*it).set[i].z);
    */
        glEnd();
    }


    /*glColor4ubv(blue);
    glBegin(GL_POINTS);
    glVertex3f(satelliteObj->furthermostPoint.x,
               satelliteObj->furthermostPoint.y,
               satelliteObj->furthermostPoint.z);
    glVertex3f(satelliteObj->nearestPoint.x,
               satelliteObj->nearestPoint.y,
               satelliteObj->nearestPoint.z);
    for(int i = 0;i < count;i++) {
        glVertex3f(particles[i].x,particles[i].y,particles[i].z);
    }
    glEnd();*/

    glBegin(GL_LINES);
    glColor4ubv(grey);
    glVertex3d(tempLine->a().x,tempLine->a().y,tempLine->a().z);
    glVertex3d(tempLine->b().x,tempLine->b().y,tempLine->b().z);
    glEnd();

    SDL_GL_SwapBuffers();
}


void setupOpenGL(int width, int height,Point &maxObjCoords) {
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
    gluPerspective(40,ratio,1.0,6*maxObjCoords.z);
    viewerPosition.z = -5*maxObjCoords.z;
    ///cout << max(viewerPosition.z - maxObjCoords.z,1.0) << endl;
    ///cout << max(viewerPosition.z + maxObjCoords.z,2.0) << endl;
    //glOrtho(-objWidth/2,objWidth/2,-objHeight/2,objHeight/2,1.0,5000.0);
    tempLine = //new Line(Point(-maxObjCoords.x,0,maxObjCoords.z),Point(maxObjCoords.x,0,maxObjCoords.z));
            new Line(viewerPosition,Point(100,12,33));
    tempVector = new Vector(tempLine->b(),tempLine->a());
}


int main(int argc, char** argv) {
    if (argc == 1) {
        fprintf(stderr,usage,argv[0]);
        exit(1);
    }

    int sleepTime = 10000; //microsecond
    // getting coordinatates from file
    vector<PlaneType> *coordinatesList = getCoordinatesFromFile(argv[1]);

    // creating object using coordinates
    Object3D *satelliteObj = new Object3D(coordinatesList);

    srand(time(NULL));
    /*if(SDL_Init(SDL_INIT_VIDEO) < 0) {
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

    int width = 640;
    int height = 480;
    int bitsPerPixel = info->vfmt->BitsPerPixel;
    int flags = SDL_OPENGL | SDL_HWSURFACE | SDL_ASYNCBLIT | SDL_RESIZABLE;
    if(SDL_SetVideoMode(width,height,bitsPerPixel,flags) == 0) {
        cerr << "Setting video mode failed: " << SDL_GetError() << endl;
        quit(1);
    }
    // set appropriate OpenGL properties
    setupOpenGL(width,height,satelliteObj->maxCoords);*/

    cout << "size of Particle: " << sizeof(Particle) << endl;
    cout << "size of Point: " << sizeof(Point) << endl;
    cout << "size of Vector: " << sizeof(Vector) << endl;

    int count = 0;
    //while(!processParticles(satelliteObj)) cout << ++count << endl;////////////////////////////////////////
    //processParticles(satelliteObj);
    processParticles(satelliteObj);


    // main program loop
    /*while(true) {
        processEvents();
        draw(satelliteObj);
        usleep(sleepTime);
    }*/

    return 0;
}

