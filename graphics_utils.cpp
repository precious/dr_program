#include "graphics_utils.h"

GLboolean Graphics::shouldRotate = GL_FALSE;
Point Graphics::viewerPosition(0,0,0);

extern int Graphics::width = 0;
extern int Graphics::height = 0;

bool Graphics::isLMousePressed = false;
pair<double,Vector> Graphics::rotationParams(0,Vector());

void Graphics::initGraphics(int _width, int _height,Point &maxObjCoords) {
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
    double maxCoordinate = 2*max(maxObjCoords.x,max(maxObjCoords.y,maxObjCoords.z));
    glFrustum(-maxCoordinate,maxCoordinate,
              -maxCoordinate/ratio,maxCoordinate/ratio,
              10,20);
    ///gluPerspective(45,ratio,1.0,15*maxObjCoords.z); /////////////// 45 30
    viewerPosition.z = -20*maxObjCoords.z;
    ///cout << max(viewerPosition.z - maxObjCoords.z,1.0) << endl;
    ///cout << max(viewerPosition.z + maxObjCoords.z,2.0) << endl;
    //glOrtho(-objWidth/2,objWidth/2,-objHeight/2,objHeight/2,1.0,5000.0);

    //tempLine = //new Line(Point(-maxObjCoords.x,0,maxObjCoords.z),Point(maxObjCoords.x,0,maxObjCoords.z));
    //        new Line(viewerPosition,Point(100,12,33));
    //tempVector = new Vector(tempLine->b,tempLine->a);
}

void Graphics::quitGraphics(int code) {
    SDL_Quit();
    exit(code);
}

void Graphics::draw(Object3D &satelliteObj,ParticlePolygon* particlesArray = NULL,int particlesAmount = 0)
{
    static float angle = 0.0f;
    static GLubyte purple[] = {255,   150, 255,   0 };
    static GLubyte grey[] = {100,100,100,0};
    static GLubyte black[] = {0,0,0,0};
    static GLubyte blue[] = {0,0,255,0};

    /*
    lengthGL lengthReal
            x = stepGL
*/


    glClearColor(255,255,255,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    glMatrixMode(GL_MODELVIEW);
    static bool firstTime = true;
    if (firstTime) {
        firstTime = false;
        glLoadIdentity();
        /* Move down the z-axis. */
        glTranslatef(viewerPosition.x, viewerPosition.y, viewerPosition.z);
    }


    /* Rotation by mouse */
    if (rotationParams.first) {
        Vector result(0,0,0);
        GLdouble modelMatrix[16];
        glGetDoublev(GL_MODELVIEW, modelMatrix);
        //////////for(int k = 0;k< 16;++k,cout << modelMatrix[k] << ", ");cout << endl;////////////////////
        for(int i = 0;i < 3;++i)
            for(int j = 0;j < 3;++j)
                result[j] += rotationParams.second[i]*(/*modelMatrix[i*4 + j] + */(i == j)? 1: 0);
        glRotated(rotationParams.first,rotationParams.second.x,
                  rotationParams.second.y,rotationParams.second.z);
        ////////////////cout << result << endl;//////////////////////////////////////////
        //glRotated(rotationParams.first,result.x,result.y,result.z);
        rotationParams.first = 0;
    }
    //

    /* Rotate. */
    /////////////////glRotatef(angle, 0.0, 1.0, 0.0);
    if (shouldRotate) {
        if (++angle > 360.0f) {
            angle = 0.0f;
        }
    }

    //Vector tmpVector(Point(),viewerPosition);

    // draw the object
    glColor4ubv(purple);
    vector<PlaneType> *coords = satelliteObj.polygons;
    for(vector<PlaneType>::iterator it = coords->begin();it != coords->end();it++) {
        glBegin(GL_LINE_LOOP);
        //tmpVector.cos((*it).getNormal()) > 0? black: purple);
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

    // draw the particles
    if (particlesArray != NULL) {
        glBegin(GL_POINTS);
        glColor4ubv(grey);
        for(int i = 0;i < particlesAmount;++i) {
            glVertex3f(particlesArray[i].first.x,particlesArray[i].first.y,particlesArray[i].first.z);
        }
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

    /*glBegin(GL_LINES);
    glColor4ubv(grey);
    tempLine = new Line(Point(),viewerPosition);
    glVertex3d(tempLine->a.x,tempLine->a.y,tempLine->a.z);
    glVertex3d(tempLine->b.x,tempLine->b.y,tempLine->b.z);
    glEnd();
    */

    SDL_GL_SwapBuffers();
}
