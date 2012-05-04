#include <unistd.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <errno.h>
#include <unistd.h>
#include <assert.h>

#include "types.h"
#include "file_utils.h"
#include "geometry_utils.h"
#include "constants.h"
#include "data_utils.h"
#include "graphics_utils.h"

#define EXIT_ERR(msg) { cerr << msg << "\nerrno: " << errno << endl; Graphics::quitGraphics(1); }
#define rand(max) rand()%max
#define PRINTLN(arg) cout << arg << endl;
#define PRINT(arg) cout << arg && cout.flush();
#define COUT(args) cout << args << endl;

//Line *tempLine;
//Vector *tempVector;
//int step = 30;

const char usage[] = "Usage:\n\tprogram [-t NUMBER][-r RADIUS][-s TIME][-m][-v][-d][-l] <filename>\n\
        -t NUMBER - test probabilty with number of particles NUMBER\n\
        -r RADIUS - radius of generative sphere\n\
        -s TIME - time to sleep in microseconds\n\
        -m - model particles\n\
        -v - verbose mode\n\
        -d - draw (requires also -l option)\n\
        -l - run mainloop";

static void handleKeyDown(SDL_keysym* keysym)
{
    switch(keysym->sym) {
    case SDLK_ESCAPE:
        Graphics::quitGraphics(0);
        break;
    case SDLK_SPACE:
        Graphics::shouldRotate = !Graphics::shouldRotate;
        break;
    case SDLK_KP_PLUS:
    case SDLK_PLUS:
        //viewerPosition.z += 20;
        //tempLine->a.z += step;
        //tempLine->b.z += step;
        break;
    case SDLK_KP_MINUS:
    case SDLK_MINUS:
        //viewerPosition.z -= 20;
        //tempLine->a.z -= step;
        //tempLine->b.z -= step;
        break;
    case SDLK_UP:
        cout << "up" << endl;
        //tempLine->a.y += step;
        //tempLine->b.y += step;
        break;
    case SDLK_DOWN:
        cout << "down" << endl;
        //tempLine->a.y -= step;
        //tempLine->b.y -= step;
        break;
    case SDLK_LEFT:
        cout << "left" << endl;
        //tempLine->a.x -= step;
        //tempLine->b.x -= step;
        break;
    case SDLK_RIGHT:
        cout << "right" << endl;
        //tempLine->a.x += step;
        //tempLine->b.x += step;
        break;
    default:
        break;
    }
}


void processEvents(void)
{
    /* Our SDL event placeholder. */
    SDL_Event event;
    static Vector zAxis(0,0,10);

    /* Grab all the events off the queue. */
    while(SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_KEYDOWN:
            handleKeyDown( &event.key.keysym );
            break;
        case SDL_QUIT:
            /* Handle quit requests (like Ctrl-c). */
            Graphics::quitGraphics(0);
            break;
        case SDL_MOUSEBUTTONDOWN:
        case SDL_MOUSEBUTTONUP:
            Graphics::isLMousePressed = !Graphics::isLMousePressed;
            break;
        case SDL_MOUSEMOTION:
            if (Graphics::isLMousePressed) {
                Vector finish(event.motion.x,-event.motion.y,0);
                Vector start(finish.x - event.motion.xrel,finish.y + event.motion.yrel,0);
                double R = max(start.length(),finish.length()) + 2;
                finish.z = -sqrt(R*R - pow(finish.length(),2));

                ///cout << start.length() << ", " << R << ", " << finish.length() << endl;//////////////////////
                //cout << event.motion.x << "  " << event.motion.y << endl;

                start.z = -sqrt(R*R - pow(start.length(),2));

                Graphics::rotationParams.second = OrientedPlane(start,finish,POINT_OF_ORIGIN).getNormal();
                Graphics::rotationParams.first = GU::getDistanceBetweenPoints(start,finish)*0.25;
                // cout << motion.resized(10) << "   " << motion.vectorProduct(zAxis).resized(10) << endl;///////////////////////
            }
        }
    }
}


int processParticles(Object3D &satelliteObj) {
    real stepLength = GeometryUtils::getDistanceBetweenPoints(satelliteObj.nearestPoint,
                                                              satelliteObj.furthermostPoint)/10.0;
    real timeInterval = stepLength/(ELECTRON_VELOCITY); // average velocity
    cout << "stepLength: " << stepLength << endl;
    cout << "timeInterval: " << timeInterval << endl;

    // static bool firstTime = true;

    // static vector<Particle> particles(count);



    /*Particle fastestParticle = GeometryUtils::getFastestParticle(particles,satelliteObj->center(),
                                [satelliteObj](Particle p) -> bool {
                                    return GeometryUtils::doesParticlesTrajectoryIntersectObject(p,*satelliteObj);
                                });
    if (fastestParticle == Point()) {
        return 0;
    }
    else
        return 1;*/

    /*cout << "fastest: " << fastestParticle << endl;
    cout << "center: " << satelliteObj->center() << endl;
    real time = Vector(fastestParticle,satelliteObj->center()).length() / fastestParticle.step.length();
    cout << "time: " << time << endl;*/

    return 0;
}

void processParticlesWhichIntersectObject(ParticlePolygon *particlesArray,int &count,double timeStep) {
    for(int i = 0;i < count;++i) {
        particlesArray[i].first = particlesArray[i].first + particlesArray[i].first.step*timeStep;
        particlesArray[i].first.ttl -= timeStep;
    }
    // checking all particlees including he last one
    for(int i = 0;i < count - 1;) {
        if (particlesArray[i].first.ttl <= 0) {
            memcpy(particlesArray + i,particlesArray + count - 1,sizeof(ParticlePolygon));
            --count;
        } else {
            ++i;
        }
    }
    // checking the last particle
    if (count > 0 && particlesArray[count - 1].first.ttl <= 0)
        --count;
}

int initRandomParticles(Particle *particlesArray,int count,GenerativeSphere generativeSphere) {
    int n;
    for(n = 0;n < count;++n) {
        particlesArray[n] = generativeSphere.generateRandomParticle(PTYPE_ELECTRON);
    }
    return n;
}

int initParticlesWhichIntersectsObject(ParticlePolygon *particlesArray,int count,
                                       GenerativeSphere &generativeSphere,bool isParticleOnSphere) {
    int n;
    ParticlePolygon tmp(Particle(),NULL);
    for(n = 0;n < count;++n) {
        tmp = generativeSphere.generateParticleWhichIntersectsObject(PTYPE_ELECTRON,isParticleOnSphere);
        memcpy(particlesArray + n,&tmp,sizeof(ParticlePolygon));
    }
    return n;
}

int main(int argc, char** argv) {
    srand(time(NULL));
    cout.precision(16);
    cout.setf(ios::fixed, ios::floatfield);

    // process arguments
    int c;
    bool modelingFlag = false;
    bool verboseFlag = false;
    bool drawFlag = false;
    bool testModeFlag = false;
    bool mainloopFlag = false;
    char *filename = NULL;
    int testProbabilityCount = -1;
    int generativeSphereRadius = -1;
    int sleepTime = 0; //microsecond


    while ((c = getopt (argc, argv, ":mvdlxt:r:s:")) != -1) {
        switch(c) {
        case 't':
            testProbabilityCount = atoi(optarg);
            break;
        case 'r':
            generativeSphereRadius = atoi(optarg);
            break;
        case 's':
            sleepTime = atoi(optarg);
            break;
        case 'd':
            drawFlag = true;
            break;
        case 'v':
            verboseFlag = true;
            break;
        case 'm':
            modelingFlag = true;
            break;
        case 'x':
            testModeFlag = true;
            break;
        case 'l':
            mainloopFlag = true;
            break;
        case '?':
        default:
            EXIT_ERR(usage);
        }
    }
    if (optind == argc) {
        EXIT_ERR(usage);
    }
    filename = argv[optind];
    if (generativeSphereRadius < 0) generativeSphereRadius = DEFAULT_GENERATIVE_SPHERE_RADIUS;


    /*------------------------------------*/
    // getting coordinatates from file
    vector<PlaneType> *coordinatesList = getCoordinatesFromFile(filename);
    assert(coordinatesList != NULL);

    // creating object using coordinates
    Object3D satelliteObj(coordinatesList);

    GenerativeSphere generativeSphere(satelliteObj.center,
                                      generativeSphereRadius,//100*GeometryUtils::getDistanceBetweenPoints(satelliteObj->center(),satelliteObj->maxCoords),
                                      satelliteObj);


    if (testProbabilityCount > 0) {
        // allocating memory for particles array
        verboseFlag && PRINTLN("memory allocation");
        Particle *particlesArray = (Particle*)calloc(testProbabilityCount,sizeof(Particle));

        verboseFlag && COUT("memory usage: " << testProbabilityCount*sizeof(Particle)/(1024*1024.0) << " MB");
        verboseFlag && PRINTLN("particles generation");
        int n = initRandomParticles(particlesArray,testProbabilityCount,generativeSphere);
        assert(n == testProbabilityCount);

        int intersectionsCounter = 0;
        verboseFlag && PRINTLN("checking for intersections");
        for(int j = 0;j < n;++j) {
            if (GeometryUtils::doesParticlesTrajectoryIntersectObject(particlesArray[j],satelliteObj))
                ++intersectionsCounter;
            verboseFlag && (!(j%(n/20 + 1))) && PRINT('.');
        }
        verboseFlag && PRINTLN("");
        if (verboseFlag) {
            COUT("percentage: " << intersectionsCounter << "/" << n << " = " << intersectionsCounter/double(n)*100 << '%');
        } else {
            cout << intersectionsCounter/double(n) << endl;
        }
        free(particlesArray);
    }


    ParticlePolygon *particlesArray = NULL;
    int particlesAmount = 0;
    double timeStep = 0;
    function<unsigned long long()> particlesAmountGenerator;
    unsigned long long newParticlesAmount;
    if (modelingFlag) {
        // for generating count of particles that intersects object
        // parameters for generator obtained using script test_intersections.sh
        // see also results of this script in freq.ods
        double a = 0.000001;//0.0001608266666667; //0.0000874633333333;
        double sigma = 0.00000008;//0.0000080581883820; //0.0000049781511517;

        // constants
        const int density = 200000; // density is 0.2 cm^-3
        const float volume = 4.0/3*M_PI*pow(generativeSphereRadius,3);
        const unsigned long long totalAmount = density*volume;
        // see explanation in draft, page 1
        const int maxParticlesAmount = (a + 4*sigma)*totalAmount;
        GaussianDistributionGenerator *particlesAmountRateGenerator = getGaussianDistributionGenerator(a,sigma);
        particlesAmountGenerator = [particlesAmountRateGenerator,totalAmount,maxParticlesAmount]() -> unsigned long long
            {return min((*particlesAmountRateGenerator)()*totalAmount,double(maxParticlesAmount));};
        newParticlesAmount = particlesAmountGenerator();

        verboseFlag && PRINTLN("particles array initialization...");
        verboseFlag && COUT("(memory will be allocated: " << maxParticlesAmount*sizeof(Particle)/pow(1024.,2) << " MB)");
        particlesAmount = particlesAmountGenerator();
        particlesArray = (ParticlePolygon*)malloc(maxParticlesAmount*sizeof(ParticlePolygon));
        verboseFlag && PRINTLN(particlesAmount);
        assert(initParticlesWhichIntersectsObject(particlesArray,particlesAmount,generativeSphere,false) == particlesAmount);

        verboseFlag && PRINTLN("searching for fastest particle...");
        ParticlePolygon fastestPP = DU::reduce<ParticlePolygon*,ParticlePolygon>(
                    [](ParticlePolygon &p1,ParticlePolygon &p2) -> ParticlePolygon&
        {return (p2.first.step.length() > p1.first.step.length())? p2: p1;},
        particlesArray,particlesAmount);
        verboseFlag && COUT("fastest particle speed: " << fastestPP.first.step.length());

        double distanceStep = GeometryUtils::getDistanceBetweenPoints(satelliteObj.nearestPoint,
                                                                    satelliteObj.furthermostPoint)/2.;
        timeStep = distanceStep/ELECTRON_VELOCITY; // time to pass 1/10 of object for particle with average velocity

        verboseFlag && PRINTLN("decreasing distance to object for all particles");      
        // time during the fastest particle will reach object
        double distanceDelta = GU::getDistanceBetweenPoints(fastestPP.first,
                                                            GU::getPlaneAndLineIntersection(*fastestPP.second,
                                                                                            Line(fastestPP.first,fastestPP.first.step)));
        double timeDelta = (distanceDelta - 5*distanceStep)/fastestPP.first.step.length();

/*        DU::map<ParticlePolygon*,ParticlePolygon>(
                    [timeDelta](ParticlePolygon &pp) -> void {pp.first = pp.first + pp.first.step*timeDelta; pp.first.ttl -= timedelta;},
                    particlesArray,particlesAmount); /// TODO check this
*/
        verboseFlag && COUT("distanceStep: " << distanceStep << "; timeStep: " << timeStep);
    }

    // video mode initialization
    if (drawFlag) {
        verboseFlag && cout << "polygons: " << satelliteObj.polygons->size() << endl;
        verboseFlag && cout << "center: " << satelliteObj.center << endl;
        verboseFlag && cout << "radius: " << satelliteObj.radius << endl;

        // set appropriate OpenGL & properties SDL
        int width = 640;
        int height = 480;
        Graphics::initGraphics(width,height,satelliteObj.maxCoords);
    }

    timespec start, stop, *delta;
    int framesCount = 0;
    double seconds = 0;
    int frames = 0;

    // -------- main program loop --------
    if (mainloopFlag && (drawFlag || modelingFlag)) {
        while(true) {

            if (drawFlag) {
                processEvents();

                clock_gettime(CLOCK_ID,&start);
                Graphics::draw(satelliteObj,particlesArray,particlesAmount);
                clock_gettime(CLOCK_ID,&stop);

                delta = getTimespecDelta(&start,&stop);
                ++frames;
                seconds += delta->tv_sec + delta->tv_nsec/pow(10,9);
                if (seconds >= 1) {
                    framesCount += frames;
                    verboseFlag && COUT(frames/seconds << " fps; frames drawed: " << framesCount);
                    seconds = frames = 0;
                }

            }

            if (modelingFlag) {
                processParticlesWhichIntersectObject(particlesArray,particlesAmount,timeStep);
                if (particlesAmount < newParticlesAmount) {
                    initParticlesWhichIntersectsObject(particlesArray + particlesAmount,
                                                       newParticlesAmount - particlesAmount,generativeSphere,true);
                    particlesAmount = newParticlesAmount;
                    newParticlesAmount = particlesAmountGenerator();
                }
            }

            sleepTime && usleep(sleepTime);
            //cout << "count" << particlesAmount << endl;
        }
    }

    if (testModeFlag) {
            cout << GU::rotatePointAroundLine(Point(0,0,1),Line(Point(0,0,0),Point(1,0,0)),M_PI/2.) << endl;
    }

    if (particlesArray != NULL) {
        free(particlesArray);
    }

    Graphics::quitGraphics(0);

    return 0;
}

