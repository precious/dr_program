#include <unistd.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <limits>

#include "types.h"
#include "file_utils.h"
#include "geometry_utils.h"
#include "constants.h"
#include "data_utils.h"
#include "graphics_utils.h"
#include "fortran_modules.h"

#define rand(max) rand()%max

const char usage[] = "Usage:\n\nprogram [-m][-v][-d][-x][-g][-t NUMBER]\
[-r RADIUS][-s TIME][-n N][-f SF][-i INTERVAL][-p STEP] <filename>\n\n\
    -t NUMBER - test probabilty with number of particles NUMBER\n\
    -r RADIUS - radius of generative sphere [not used]\n\
    -s TIME - time to sleep in microseconds\n\
    -m - model particles\n\
    -v - verbose mode\n\
    -d - draw\n\
    -a - draw axes\n\
    -f SF - scale factor for coordinates in file to reduce them to SI\n\
        (default 1)\n\
    -n N - total number of particles at time momentn\n\
    -i INTERVAL - interval to print measurings\n\
        (use 'i' prefix for number of steps or 's' for seconds)\n\
    -p STEP - step of particle measured in spacecrafts length\n\
        (default 0.25)\n\
    -x - file with complex data format\n";

namespace Globals {
    unsigned long long  realToModelNumber;
    double electricFieldToCharge = 0;
}

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

// for internal usage only
inline void finalizeParticle(Object3D &satelliteObj,Particle* particles,
                             unsigned long long &electronsNumber,unsigned long long &ionsNumber,int i) {
    switch(particles[i].type) {
    case PTYPE_ELECTRON:
        if (particles[i].polygonIndex != -1) {
            satelliteObj.changePlasmaCurrents(particles[i].polygonIndex,ELECTRON_CURRENT_DENSITY);
            satelliteObj.totalCharge += ELECTRON_ELECTRIC_CHARGE*Globals::realToModelNumber;
        }
        electronsNumber--;
        break;
    case PTYPE_ION:
        if (particles[i].polygonIndex != -1) {
            satelliteObj.changePlasmaCurrents(particles[i].polygonIndex,ION_CURRENT_DENSITY);
            satelliteObj.totalCharge += ION_ELECTRIC_CHARGE*Globals::realToModelNumber;
        }
        ionsNumber--;
        break;
    }
}

int processParticles(Object3D &satelliteObj,Particle* particles,
                     unsigned long long &electronsNumber,unsigned long long &ionsNumber,
                     double timeStep) {
    int finalizedNumber = 0;
    for(unsigned long long i = 0;i < electronsNumber + ionsNumber;++i) {
        particles[i] = particles[i] + particles[i].speed*timeStep;
        particles[i].ttl -= timeStep;
    }

    // checking all particles excluding the last one
    for(unsigned long long i = 0;i < electronsNumber + ionsNumber - 1;) {
        if (particles[i].ttl <= 0) {
            if (particles[i].polygonIndex != -1)
                ++finalizedNumber;
            finalizeParticle(satelliteObj,particles,electronsNumber,ionsNumber,i);
            memcpy(particles + i,particles + electronsNumber + ionsNumber - 1,sizeof(Particle));
        } else {
            ++i;
        }
    }
    // checking the last particle
    int lastIndex = electronsNumber + ionsNumber - 1;
    if (lastIndex >= 0 && particles[lastIndex].ttl <= 0) {
        if (particles[lastIndex].polygonIndex != -1)
            ++finalizedNumber;
        finalizeParticle(satelliteObj,particles,electronsNumber,ionsNumber,lastIndex);
    }
    return finalizedNumber;
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
    char *filename = NULL;
    int testProbabilityCount = -1;
    int generativeSphereRadius = -1;
    int sleepTime = 0; //microsecond
    double printInterval = 10000.0;
    double intervalInSteps = true;
    double distanceStepCoef = 0.25;
    unsigned long long averageParticlesNumber = 10000;
    float complexDataFileFlag = false;

    while ((c = getopt (argc, argv, ":vdamxgt:r:s:f:t:n:i:p:")) != -1) {
        switch(c) {
        case 'a':
            Graphics::drawAxes = true;
            break;
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
            complexDataFileFlag = true;
            break;
        case 'p':
            distanceStepCoef = atof(optarg);
            break;
        case 'i':
            if(optarg[0] == 'i')
            { printInterval = atof(optarg + 1); intervalInSteps = true; }
            else if (optarg[0] == 's')
            { printInterval = atof(optarg + 1); intervalInSteps = false; }
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
    vector<PlaneType> *coordinatesList;
    if (complexDataFileFlag) {
        coordinatesList = File::getCoordinatesFromSpecialFile(filename);
    } else {
        coordinatesList = File::getCoordinatesFromPlainFile(filename);
    }
    assert(coordinatesList != NULL);

    // creating object using coordinates
    Object3D satelliteObj(coordinatesList);

    GenerativeSphere electronsGenerativeSphere(satelliteObj.center,
                                      ELECTRONS_GENERATIVE_SPHERE_RADIUS,
                                      satelliteObj);
    GenerativeSphere ionsGenerativeSphere(satelliteObj.center,
                                      IONS_GENERATIVE_SPHERE_RADIUS,
                                      satelliteObj);    

    double electronsToIonsRatio = 1.*pow(ELECTRONS_GENERATIVE_SPHERE_RADIUS,3)*ELECTRONS_CONSISTENCE/
            (pow(IONS_GENERATIVE_SPHERE_RADIUS,3)*IONS_CONSISTENCE);
    unsigned long long  averageElectronsNumber = electronsToIonsRatio*averageParticlesNumber/(electronsToIonsRatio + 1);
    unsigned long long  averageIonsNumber = averageParticlesNumber/(electronsToIonsRatio + 1);
    Particle::electronTrajectoryCurrent =
            4*M_PI*pow(electronsGenerativeSphere.radius,2)*ELECTRON_CURRENT_DENSITY / averageElectronsNumber;
    Particle::ionTrajectoryCurrent =
            4*M_PI*pow(ionsGenerativeSphere.radius,2)*ION_CURRENT_DENSITY / averageIonsNumber;

    verboseFlag && COUT("electron trajectory current: " << Particle::electronTrajectoryCurrent);
    verboseFlag && COUT("ion trajectory current: " << Particle::ionTrajectoryCurrent);
    Globals::realToModelNumber = 4.0/3.0*M_PI*pow(ELECTRONS_GENERATIVE_SPHERE_RADIUS,3)
            *ELECTRONS_CONSISTENCE/averageElectronsNumber;

    verboseFlag && COUT("real number/modeln number: " << Globals::realToModelNumber);

    if (testProbabilityCount > 0) {
        // allocating memory for particles array
        verboseFlag && PRINTLN("memory allocation");
        Particle *particlesArray = (Particle*)calloc(testProbabilityCount,sizeof(Particle));

        verboseFlag && COUT("memory usage: " << testProbabilityCount*sizeof(Particle)/(1024*1024.0) << " MB");
        verboseFlag && PRINTLN("particles generation");
        electronsGenerativeSphere.populateArray(particlesArray,testProbabilityCount,PTYPE_ELECTRON,GEN_RANDOM);

        int intersectionsCounter = 0;
        verboseFlag && PRINTLN("checking for intersections");
        for(int j = 0;j < testProbabilityCount;++j) {
            if (Geometry::doesParticlesTrajectoryIntersectObject(particlesArray[j],satelliteObj))
                ++intersectionsCounter;
            verboseFlag && (!(j%(testProbabilityCount/20 + 1))) && PRINT('.');
        }
        verboseFlag && PRINTLN("");
        if (verboseFlag) {
            COUT("percentage: " << intersectionsCounter << "/" << testProbabilityCount
                 << " = " << intersectionsCounter/double(testProbabilityCount)*100 << '%');
        } else {
            cout << intersectionsCounter/double(testProbabilityCount) << endl;
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
        verboseFlag && COUT("average number of electrons: " << electronsNumber << ", ions: " << ionsNumber);
        verboseFlag && COUT("number of particles: " << electronsNumber + ionsNumber << endl << "initialization...");

        electronsGenerativeSphere.populateArray(particlesArray,electronsNumber,PTYPE_ELECTRON,GEN_ON_SPHERE);
        ionsGenerativeSphere.populateArray(particlesArray + electronsNumber,ionsNumber,PTYPE_ION,GEN_ON_SPHERE);

        verboseFlag && PRINTLN("searching for fastest particle...");
        Object3D *satelliteObjPtr = &satelliteObj;
        Particle fastestParticle = Data::reduce<Particle*,Particle>([satelliteObjPtr](Particle &p1,Particle &p2) -> Particle& {
            if (p1.polygonIndex == -1) return p2;
            if (p2.polygonIndex == -1) return p1;
            return (p2.speed.length() > p1.speed.length())? p2: p1;
        }, particlesArray,ionsNumber + electronsNumber);
        verboseFlag && COUT("fastest particle speed: " << fastestParticle.speed.length());

        double distanceStep = satelliteObj.radius*2.0*distanceStepCoef;
        timeStep = distanceStep/ELECTRON_VELOCITY_M; // time to do step for particle with average velocity

        verboseFlag && PRINTLN("decreasing distance to object for all particles");      
        // time during the fastest particle will reach object
        double distanceDelta = Geometry::getDistanceBetweenPointAndSphere(satelliteObj,fastestParticle);
        double timeDelta = (distanceDelta - 5*distanceStep)/fastestParticle.speed.length();

        for_each(particlesArray,particlesArray + ionsNumber + electronsNumber,
                    [timeDelta](Particle &pp) -> void {pp = pp + pp.speed*timeDelta; pp.ttl -= timeDelta;}); //TODO check this

        verboseFlag && COUT("distanceStep: " << distanceStep << "; timeStep: " << timeStep);
    }

    verboseFlag && COUT("polygons: " << satelliteObj.polygons->size());
    verboseFlag && COUT("center: " << satelliteObj.center);
    verboseFlag && COUT("radius: " << satelliteObj.radius);
    verboseFlag && COUT("capacitance: " << satelliteObj.capacitance());

    // video mode initialization
    if (drawFlag) {
        // set appropriate OpenGL & properties SDL
        int width = 640;
        int height = 480;
        Graphics::initGraphics(width,height);
    }

    timespec start, stop, *delta;
    int framesCount = 0;
    double seconds = 0;
    int frames = 0;

    // -------- main program loop --------
    unsigned long long newElectronsNumber = min<unsigned long long>(electronsNumberGenerator(),maxElectronsNumber);
    unsigned long long newIonsNumber = min<unsigned long long>(ionsNumberGenerator(),maxIonsNumber);
    double elapsedTime = 0.0;
    double timeToPrint = printInterval;
    double spacecraftCapacitance = satelliteObj.capacitance();
    double surfaceCharge;
    unsigned long long numberOfIntersections = 0;
    if (drawFlag || modelingFlag) {
        if (modelingFlag) {
            solveBoundaryProblem(coordinatesList,verboseFlag); // solve using fortran module
        }
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
                numberOfIntersections += processParticles(satelliteObj,particlesArray,electronsNumber,ionsNumber,timeStep);
                elapsedTime += timeStep;
                timeToPrint -= ((intervalInSteps)? 1: timeStep);
                surfaceCharge = satelliteObj.totalPlasmaCurrent*elapsedTime;
                if (timeToPrint <= 0) {
                    cout << satelliteObj.totalPlasmaCurrent << " " << surfaceCharge << " " << surfaceCharge/spacecraftCapacitance
                         << "   " << elapsedTime << "   " << numberOfIntersections*Globals::realToModelNumber << " " << satelliteObj.totalCharge << endl;
                    (timeToPrint = printInterval);
                }
                // processing new particles if necessary
                if (electronsNumber < newElectronsNumber) {
                    electronsGenerativeSphere.populateArray(particlesArray + electronsNumber + ionsNumber,
                                           newElectronsNumber - electronsNumber,PTYPE_ELECTRON,GEN_ON_SPHERE);
                    electronsNumber = newElectronsNumber;
                    newElectronsNumber = min<unsigned long long>(electronsNumberGenerator(),maxElectronsNumber);
                }
                if (ionsNumber < newIonsNumber) {
                    ionsGenerativeSphere.populateArray(particlesArray + electronsNumber + ionsNumber,
                                           newIonsNumber - ionsNumber,PTYPE_ION,GEN_ON_SPHERE);
                    ionsNumber = newIonsNumber;
                    newIonsNumber = min<unsigned long long>(ionsNumberGenerator(),maxIonsNumber);
                }
            }

            sleepTime && usleep(sleepTime);
        }
    }

    if (particlesArray != NULL) {
        free(particlesArray);
    }

    Graphics::quitGraphics(0);

    return 0;
}

