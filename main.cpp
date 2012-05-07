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

const char usage[] = "Usage:\n\tprogram [-t NUMBER][-r RADIUS][-s TIME][-m][-v][-d][-l] <filename>\n\
    -t NUMBER - test probabilty with number of particles NUMBER\n\
    -r RADIUS - radius of generative sphere [not used]\n\
    -s TIME - time to sleep in microseconds\n\
    -m - model particles\n\
    -v - verbose mode\n\
    -d - draw (requires also -l option)\n\
    -l - run mainloop\n\
    -f SF - scale factor for coordinates in file to reduce them to SI\n\
        (default 0.001)\n\
    -n - total number of particles at time moment\n";

static void handleKeyDown(SDL_keysym* keysym)
{
    switch(keysym->sym) {
    case SDLK_ESCAPE:
        Graphics::quitGraphics(0);
        break;
    default:
        break;
    }
}

void processEvents(void)
{
    /* Our SDL event placeholder. */
    SDL_Event event;
    float zoomDelta = 0.01;
    float zoomDelta2 = 0.1;

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
            switch(event.button.button) {
            case SDL_BUTTON_LEFT:
                Graphics::isLMousePressed = true;
                break;
            case SDL_BUTTON_WHEELDOWN:
                if (Graphics::zoomFactor > 2*zoomDelta)
                    Graphics::zoomFactor -= (Graphics::zoomFactor > 2)? zoomDelta2: zoomDelta;
                break;
            case SDL_BUTTON_WHEELUP:
               Graphics::zoomFactor += (Graphics::zoomFactor > 2)? zoomDelta2: zoomDelta;
               break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            if (event.button.button == SDL_BUTTON_LEFT)
                Graphics::isLMousePressed = false;
            break;
        case SDL_MOUSEMOTION:
            if (Graphics::isLMousePressed) {
                double coef = 0.5;
                Graphics::rotationAngles[0] += -event.motion.xrel*coef;
                Graphics::rotationAngles[0] -= 360*int(Graphics::rotationAngles[0]/360);
                Graphics::rotationAngles[1] += -event.motion.yrel*coef;
                Graphics::rotationAngles[1] -= 360*int(Graphics::rotationAngles[1]/360);
            }
        }
    }
}

// for internal usge only
inline void finalizeParticle(Object3D &satelliteObj,Particle* particles,
                             unsigned long long &electronsNumber,unsigned long long &ionsNumber,int i) {
    particles[i].polygonIndex != -1 && COUT("Collision!" << satelliteObj.charge);///////////////////////////////////////////////////////
    switch(particles[i].type) {
    case PTYPE_ELECTRON:
        if (particles[i].polygonIndex != -1)
            satelliteObj.changeCharge(particles[i].polygonIndex,ELECTRON_ELECTRIC_CHARGE);
        electronsNumber--;
        break;
    case PTYPE_ION:
        if (particles[i].polygonIndex != -1)
            satelliteObj.changeCharge(particles[i].polygonIndex,ION_ELECTRIC_CHARGE);
        ionsNumber--;
        break;
    }
}

int processParticles(Object3D &satelliteObj,Particle* particles,
                     unsigned long long &electronsNumber,unsigned long long &ionsNumber,
                     double timeStep) {
    for(int i = 0;i < electronsNumber + ionsNumber;++i) {
        particles[i] = particles[i] + particles[i].step*timeStep;
        particles[i].ttl -= timeStep;
    }

    // checking all particlees excluding the last one
    for(int i = 0;i < electronsNumber + ionsNumber - 1;) {
        if (particles[i].ttl <= 0) {
            finalizeParticle(satelliteObj,particles,electronsNumber,ionsNumber,i);
            memcpy(particles + i,particles + electronsNumber + ionsNumber - 1,sizeof(Particle));
        } else {
            ++i;
        }
    }
    // checking the last particle
    int lastIndex = electronsNumber + ionsNumber - 1;
    if (lastIndex >= 0 && particles[lastIndex].ttl <= 0)
        finalizeParticle(satelliteObj,particles,electronsNumber,ionsNumber,lastIndex);
}

void processParticlesWhichIntersectObject(ParticlePolygon *particlesArray,int &count,double timeStep) {
    for(int i = 0;i < count;++i) {
        particlesArray[i].first = particlesArray[i].first + particlesArray[i].first.step*timeStep;
        particlesArray[i].first.ttl -= timeStep;
    }
    // checking all particlees excluding the last one
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
        particlesArray[n] = generativeSphere.generateParticleInSphere(PTYPE_ELECTRON);
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

int initParticlesOnSphere(Object3D& obj,
                          Particle *particlesArray,
                          int electronsNumber, int ionsNumber,
                          GenerativeSphere electronsGenerativeSphere,
                          GenerativeSphere ionsGenerativeSphere) {
    //COUT("generating..." << electronsNumber << " * " << ionsNumber << " * " << obj.charge);///////
    int n = 0;
    for(;n < electronsNumber;++n) {
        particlesArray[n] = electronsGenerativeSphere.generateParticleOnSphere(PTYPE_ELECTRON);

    }
    for(;n < electronsNumber + ionsNumber;++n) {
        particlesArray[n] = ionsGenerativeSphere.generateParticleOnSphere(PTYPE_ION);
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
    unsigned long long averageParticlesNumber = 1000000;

    while ((c = getopt (argc, argv, ":mvdlxt:r:s:f:t:n:")) != -1) {
        switch(c) {
        case 't':
            testProbabilityCount = atoi(optarg);
            break;
        case 'r':
            generativeSphereRadius = atoi(optarg);
            break;
        case 'f':
            File::scaleFactor = atof(optarg);
            break;
        case 's':
            sleepTime = atoi(optarg);
            break;
        case 'n':
            averageParticlesNumber = atoll(optarg);
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
    //if (generativeSphereRadius < 0) generativeSphereRadius = DEFAULT_GENERATIVE_SPHERE_RADIUS;


    /*------------------------------------*/
    // getting coordinatates from file
    vector<PlaneType> *coordinatesList = File::getCoordinatesFromFile(filename);
    assert(coordinatesList != NULL);

    // creating object using coordinates
    Object3D satelliteObj(coordinatesList);

    GenerativeSphere electronsGenerativeSphere(satelliteObj.center,
                                      ELECTRONS_GENERATIVE_SPHERE_RADIUS,
                                      satelliteObj);

    GenerativeSphere ionsGenerativeSphere(satelliteObj.center,
                                      IONS_GENERATIVE_SPHERE_RADIUS,
                                      satelliteObj);

    double electronsToIonsRatio = 1.*pow(ELECTRONS_GENERATIVE_SPHERE_RADIUS,3)*ELECTRONS_CONSISENCE/
            (pow(IONS_GENERATIVE_SPHERE_RADIUS,3)*IONS_CONSISENCE);
    unsigned long long  averageElectronsNumber = electronsToIonsRatio*averageParticlesNumber/(electronsToIonsRatio + 1);
    unsigned long long  averageIonsNumber = averageParticlesNumber/(electronsToIonsRatio + 1);


    if (testProbabilityCount > 0) {
        // allocating memory for particles array
        verboseFlag && PRINTLN("memory allocation");
        Particle *particlesArray = (Particle*)calloc(testProbabilityCount,sizeof(Particle));

        verboseFlag && COUT("memory usage: " << testProbabilityCount*sizeof(Particle)/(1024*1024.0) << " MB");
        verboseFlag && PRINTLN("particles generation");
        int n = initRandomParticles(particlesArray,testProbabilityCount,electronsGenerativeSphere);
        assert(n == testProbabilityCount);

        int intersectionsCounter = 0;
        verboseFlag && PRINTLN("checking for intersections");
        for(int j = 0;j < n;++j) {
            if (Geometry::doesParticlesTrajectoryIntersectObject(particlesArray[j],satelliteObj))
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

    Particle *particlesArray = NULL;
    double timeStep = 0;
    unsigned long long maxParticlesNumber = averageParticlesNumber*1.5;
    unsigned long long maxElectronsNumber = electronsToIonsRatio*maxParticlesNumber/(electronsToIonsRatio + 1);
    unsigned long long maxIonsNumber = maxParticlesNumber/(electronsToIonsRatio + 1);
    unsigned long long electronsNumber;
    unsigned long long ionsNumber;
    GaussianDistributionGenerator electronsNumberGenerator =
            Time::getGaussianDistributionGenerator(averageElectronsNumber,averageElectronsNumber*0.05);
    GaussianDistributionGenerator ionsNumberGenerator =
            Time::getGaussianDistributionGenerator(averageIonsNumber,averageIonsNumber*0.05);
    if (modelingFlag) {
        verboseFlag && PRINTLN("particles array initialization...");
        verboseFlag && COUT("(memory will be allocated: " << maxParticlesNumber*sizeof(Particle)/pow(1024.,2) << " MB)");
        electronsNumber = averageElectronsNumber;
        ionsNumber = averageIonsNumber;
        particlesArray = (Particle*)malloc(maxParticlesNumber*sizeof(Particle));
        verboseFlag && COUT("number of particles: " << electronsNumber + ionsNumber << endl << "initialization...");


        assert(initParticlesOnSphere(satelliteObj,particlesArray,electronsNumber,ionsNumber,
                                     electronsGenerativeSphere,ionsGenerativeSphere)
               == electronsNumber + ionsNumber);

        verboseFlag && PRINTLN("searching for fastest particle...");
        Object3D *satelliteObjPtr = &satelliteObj;
        Particle fastestParticle = Data::reduce<Particle*,Particle>([satelliteObjPtr](Particle &p1,Particle &p2) -> Particle& {
            if (p1.polygonIndex == -1) return p2;
            if (p2.polygonIndex == -1) return p1;
            return (p2.step.length() > p1.step.length())? p2: p1;
        }, particlesArray,ionsNumber + electronsNumber);
        verboseFlag && COUT("fastest particle speed: " << fastestParticle.step.length());

        double distanceStep = satelliteObj.radius/2.;
        timeStep = distanceStep/ELECTRON_VELOCITY_M; // time to do step for particle with average velocity

        verboseFlag && PRINTLN("decreasing distance to object for all particles");      
        // time during the fastest particle will reach object
        double distanceDelta = Geometry::getDistanceBetweenPointAndSphere(satelliteObj,fastestParticle);
        double timeDelta = (distanceDelta - 5*distanceStep)/fastestParticle.step.length();

        Data::map<Particle*,Particle>(
                    [timeDelta](Particle &pp) -> void {pp = pp + pp.step*timeDelta; pp.ttl -= timeDelta;},
                    particlesArray,ionsNumber + electronsNumber); /// TODO check this

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
        Graphics::initGraphics(width,height,satelliteObj);
    }

    timespec start, stop, *delta;
    int framesCount = 0;
    double seconds = 0;
    int frames = 0;

    // -------- main program loop --------
    unsigned long long newElectronsNumber = min<unsigned long long>(electronsNumberGenerator(),maxElectronsNumber);
    unsigned long long newIonsNumber = min<unsigned long long>(ionsNumberGenerator(),maxIonsNumber);

    if (mainloopFlag && (drawFlag || modelingFlag)) {
        while(true) {

            if (drawFlag) {
                processEvents();

                clock_gettime(CLOCK_ID,&start);
                Graphics::draw(satelliteObj,particlesArray,electronsNumber + ionsNumber);
                clock_gettime(CLOCK_ID,&stop);

                delta = Time::getTimespecDelta(&start,&stop);
                ++frames;
                seconds += delta->tv_sec + delta->tv_nsec/pow(10,9);
                if (seconds >= 1) {
                    framesCount += frames;
                    verboseFlag && COUT(frames/seconds << " fps; frames drawed: " << framesCount);
                    seconds = frames = 0;
                }

            }

            if (modelingFlag) {
                processParticles(satelliteObj,particlesArray,electronsNumber,ionsNumber,timeStep);
                for(int i = 0;i < electronsNumber + ionsNumber;++i)
                    if(particlesArray[i].polygonIndex != -1) {////////////////////////////////////////////////////////////////
                        COUT(particlesArray[i].step.length() << " - " << particlesArray[i].ttl << "  steps");
                    }
                // processing new particles if necessary
                if (electronsNumber < newElectronsNumber) {
                    initParticlesOnSphere(satelliteObj,particlesArray + electronsNumber + ionsNumber,
                                          newElectronsNumber - electronsNumber,0,
                                          electronsGenerativeSphere,ionsGenerativeSphere);
                    electronsNumber = newElectronsNumber;
                    newElectronsNumber = min<unsigned long long>(electronsNumberGenerator(),maxElectronsNumber);
                }
                if (ionsNumber < newIonsNumber) {
                    initParticlesOnSphere(satelliteObj,particlesArray + electronsNumber + ionsNumber,
                                          0,newIonsNumber - ionsNumber,
                                          electronsGenerativeSphere,ionsGenerativeSphere);
                    ionsNumber = newIonsNumber;
                    newIonsNumber = min<unsigned long long>(ionsNumberGenerator(),maxIonsNumber);
                }
            }

            sleepTime && usleep(sleepTime);
        }
    }

    if (testModeFlag) {
    }

    if (particlesArray != NULL) {
        free(particlesArray);
    }

    Graphics::quitGraphics(0);

    return 0;
}

