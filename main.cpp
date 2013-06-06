#include <unistd.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <limits>
#include <pthread.h>
#include <signal.h>

#include "types.h"
#include "file_utils.h"
#include "geometry_utils.h"
#include "constants.h"
#include "data_utils.h"
#include "graphics_utils.h"
#include "fortran_modules.h"

#define rand(max) rand()%max

const char usage[] = "Usage:\n\nprogram [OPTIONS] <filename>\n\n\
    -t NUMBER - test probabilty with number of particles NUMBER\n\
    -r RADIUS - radius of generative sphere [not used]\n\
    -s TIME - time to sleep in microseconds\n\
    -m - model particles\n\
    -o NUM - modeling type; 1, 2 or 3\n\
        1 - modeling without field affecting\n\
        2 - (default) most optimized modeling\n\
        3 - modeling best applicable for drawing\n\
    -v - verbose mode\n\
    -d - draw\n\
    -j - draw trajectories of particles\n\
    -a - draw axes\n\
    -c CHARGE - initial charge of spacecraft (default -0.0000005)\n\
    -f SF - scale factor for coordinates in file to reduce them to SI\n\
        (default 1)\n\
    -n N - total number of particles at time momentn\n\
    -i INTERVAL - interval to print measurings\n\
        (use 'i' prefix for number of steps or 's' for seconds)\n\
    -p STEP - step of particle measured in spacecrafts length\n\
        (default 0.25)\n\
    -x - file with complex data format\n\
    -h NUM - number of posix threads (default 1)\n";

namespace Globals {
    unsigned long long  realToModelNumber;
    long double initialCharge = -0.0000005; // -0.00000445;
    bool debug = true;
    bool pause = false;
    bool drawTrajectories = false;
    int modelingType = 2;
    unsigned int threadNum = 1;
}

static void handleKeyDown(SDL_keysym* keysym)
{
    switch(keysym->sym) {
    case SDLK_ESCAPE:
        Graphics::quitGraphics(0);
        break;
    case SDLK_p:
        Globals::pause = !Globals::pause;
        break;
    default:
        break;
    }
}

void processEvents(void)
{
    /* Our SDL event placeholder. */
    SDL_Event event;
    float zoomDelta = 0.05;
    float zoomDelta2 = 0.5;

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
    if (particles[i].behaviour == PARTICLE_WILL_INTERSECT_OBJ) {
        satelliteObj.totalCharge += PARTICLE_CHARGE(particles[i].type)*Globals::realToModelNumber;
        Globals::debug && COUT("charge delta = " << PARTICLE_CHARGE(particles[i].type)*Globals::realToModelNumber << ", totalCharge = " << satelliteObj.totalCharge);
        satelliteObj.changePlasmaCurrents((particles[i].type == PTYPE_ELECTRON)?
                                              Particle::electronTrajectoryCurrent:
                                              Particle::ionTrajectoryCurrent);
    }

    switch(particles[i].type) {
    case PTYPE_ELECTRON:
        electronsNumber--;
        break;
    case PTYPE_ION:
        ionsNumber--;
        break;
    }

    particles[i].finalize();
}

namespace Globals {
    struct _Env { // struct used in pthreads
        Object3D *satelliteObj;
        GenerativeSphere *electronsGenerativeSphere;
        GenerativeSphere *ionsGenerativeSphere;
        double timeStep;
    } env;
}

void* processParticlesArray(void *_args) { // function to call in pthread; Globals::env should be filled
    pair<Particle*,unsigned long long> *args = (pair<Particle*,unsigned long long>*)_args;
    Particle *particles = args->first;
    unsigned long long num = args->second;
//    COUT("num = " << num);
    Object3D &satelliteObj = *Globals::env.satelliteObj;
    double timeStep = Globals::env.timeStep;
    Vector fieldGrad;
    real fieldPot;

    for(unsigned long long i = 0;i < num;++i) {
        GenerativeSphere *gs = (particles[i].type == PTYPE_ELECTRON)? Globals::env.electronsGenerativeSphere: Globals::env.ionsGenerativeSphere;
        if (Geometry::getDistanceBetweenPoints(gs->center,particles[i]) > gs->radius) { // kick paricle if it has left the modeling sphere
            particles[i].ttl = -1;
            continue;
        }

        if (Globals::drawTrajectories)
            particles[i].addPreviousStates(particles[i]);

            // ---------------------------------------------------------------------------------------------------------------------
        if (Globals::modelingType == 1) { // modeling without affecting of field
            // ---------------------------------------------------------------------------------------------------------------------
            // if for some reason particle has left generative sphere - remove it
            if (particles[i].behaviour == PARTICLE_HAS_UNDEFINED_BEHAVIOUR) {
                real index = Geometry::getIndexOfPolygonThatParicleIntersects(satelliteObj,particles[i]);
                if (index == -1) {
                    particles[i].behaviour = PARTICLE_WILL_NOT_INTERSECT_OBJ;
                    particles[i].ttl = 100; // particle will be kicked
                } else {
                    particles[i].behaviour = PARTICLE_WILL_INTERSECT_OBJ;
                    particles[i].ttl = Geometry::getDistanceBetweenPointAndPlane(satelliteObj.polygons->at(index),particles[i]) /
                            particles[i].speed.length();
                }
            }
            particles[i] = particles[i] + particles[i].speed*timeStep;
            particles[i].ttl -= timeStep;
            // ---------------------------------------------------------------------------------------------------------------------
        } else if (Globals::modelingType == 2) { // most optimized modeling
            // ---------------------------------------------------------------------------------------------------------------------
            if (particles[i].behaviour == PARTICLE_WILL_INTERSECT_OBJ || particles[i].behaviour == PARTICLE_WILL_NOT_INTERSECT_OBJ) {
                particles[i] = particles[i] + particles[i].speed*timeStep;
                particles[i].ttl -= timeStep;
            } else { // particles[i].behaviour == PARTICLE_HAS_UNDEFINED_BEHAVIOUR
                real distanceToSatellite = Geometry::getDistanceBetweenPointAndSphere(satelliteObj,particles[i]);
                int index;

                if (distanceToSatellite == 0 || // if particle is inside satellite's sphere or too close to sphere and will be inside it soon
                        (distanceToSatellite < particles[i].speed.length()*timeStep &&
                        Geometry::doesLineIntersectSphere(Line(particles[i],particles[i].speed),satelliteObj)))  {
                    index = Geometry::getIndexOfPolygonThatParicleIntersects(satelliteObj,particles[i]);
                    if (index != -1) { // then particle will intersect object
                        particles[i].behaviour = PARTICLE_WILL_INTERSECT_OBJ;
                        particles[i].ttl = Geometry::getDistanceBetweenPointAndPlane(satelliteObj.polygons->at(index),particles[i]) /
                                particles[i].speed.length();
                        Globals::debug && COUT("particle will intersect object, ttl = " << particles[i].ttl <<", timestep = " << timeStep << ", steps = " << particles[i].ttl/timeStep <<", behaviour = " << particles[i].behaviour);
                        continue;
                    } else { // then particle will not intersect object
                        particles[i].behaviour = PARTICLE_WILL_NOT_INTERSECT_OBJ;
                        particles[i].ttl = 100; // particle will be kicked
                    }
                } else {
                    real distanceToCenterOfSatellite = Geometry::getDistanceBetweenPoints(satelliteObj.center,particles[i]);
                    resultf_(particles + i,&fieldPot,&fieldGrad); // get gradient of field in the current point
                    real electricField = satelliteObj.totalCharge/(4*M_PI*VACUUM_PERMITTIVITY*distanceToCenterOfSatellite*distanceToCenterOfSatellite);
                    fieldGrad.resize(electricField); // resize gradient vector according to current satellite charge by formula 1 in the draft
                    particles[i].affectField(fieldGrad,timeStep);
                    index = Geometry::getIndexOfPolygonThatParicleIntersects(satelliteObj,particles[i]);
                    if (index == -1 && sign(PARTICLE_CHARGE(particles[i].type)) == sign(satelliteObj.totalCharge)) {
                        // particle will be repeled by satellite
                        particles[i].behaviour = PARTICLE_WILL_NOT_INTERSECT_OBJ;
                        particles[i].ttl = 100; // particle will be kicked
                    } else if (index != -1 && sign(PARTICLE_CHARGE(particles[i].type)) == -sign(satelliteObj.totalCharge)) {
                        // particle will intersect satellite
                        particles[i].behaviour = PARTICLE_WILL_INTERSECT_OBJ;
                        particles[i].ttl = Geometry::getDistanceBetweenPointAndPlane(satelliteObj.polygons->at(index),particles[i]) /
                                particles[i].speed.length();
                        Globals::debug && COUT("particle will intersect satellite, ttl = " << particles[i].ttl <<", timestep = " << timeStep << ", steps = " << particles[i].ttl/timeStep <<", behaviour = " << particles[i].behaviour);
                    }
                }
            }
            // ---------------------------------------------------------------------------------------------------------------------
        } else { // modeling best applicable for drawing
            // ---------------------------------------------------------------------------------------------------------------------
            if (particles[i].behaviour == PARTICLE_WILL_INTERSECT_OBJ) {
                particles[i] = particles[i] + particles[i].speed*timeStep;
                particles[i].ttl -= timeStep;
            } else { // particles[i].behaviour == PARTICLE_HAS_UNDEFINED_BEHAVIOUR
                real distanceToSatellite = Geometry::getDistanceBetweenPointAndSphere(satelliteObj,particles[i]);
                int index;
                if (distanceToSatellite == 0 || // if particle is inside satellite's sphere or too close to sphere and will be inside it soon
                        (distanceToSatellite < particles[i].speed.length()*timeStep &&
                        Geometry::doesLineIntersectSphere(Line(particles[i],particles[i].speed),satelliteObj)))  {
                    index = Geometry::getIndexOfPolygonThatParicleIntersects(satelliteObj,particles[i]);
                    if (index != -1) { // then particle will intersect object
                        particles[i].behaviour = PARTICLE_WILL_INTERSECT_OBJ;
                        particles[i].ttl = Geometry::getDistanceBetweenPointAndPlane(satelliteObj.polygons->at(index),particles[i]) /
                                particles[i].speed.length();
                        Globals::debug && COUT("particle will intersect object, ttl = " << particles[i].ttl <<", timestep = " << timeStep << ", steps = " << particles[i].ttl/timeStep <<", behaviour = " << particles[i].behaviour);
                    }
                } else {
                    real distanceToCenterOfSatellite = Geometry::getDistanceBetweenPoints(satelliteObj.center,particles[i]);
                    resultf_(particles + i,&fieldPot,&fieldGrad); // get gradient of field in the current point
                    real electricField = satelliteObj.totalCharge/(4*M_PI*VACUUM_PERMITTIVITY*distanceToCenterOfSatellite*distanceToCenterOfSatellite);
                    fieldGrad.resize(electricField); // resize gradient vector according to current satellite charge by formula 1 in the draft
                    particles[i].affectField(fieldGrad,timeStep);
                }
            }
            // ---------------------------------------------------------------------------------------------------------------------
        } // end of modelling branch
    } // end of for loop
    return NULL;
}

int processParticles(Object3D &satelliteObj,Particle* particles,
                     unsigned long long &electronsNumber,unsigned long long &ionsNumber,
                     double timeStep,GenerativeSphere electronsGenerativeSphere,
                     GenerativeSphere ionsGenerativeSphere) {


    Globals::env.electronsGenerativeSphere = &electronsGenerativeSphere;
    Globals::env.ionsGenerativeSphere = &ionsGenerativeSphere;
    Globals::env.satelliteObj = &satelliteObj;
    Globals::env.timeStep = timeStep;

    if(Globals::threadNum == 1) {
        pair<Particle*,unsigned long long> args(particles, electronsNumber + ionsNumber);
        processParticlesArray(&args);
    } else {
        pthread_t *threads = new pthread_t[Globals::threadNum];
        int particlesPerThread = ceil(1.0*(electronsNumber + ionsNumber)/Globals::threadNum);
        int firtsParticleForCurrentThread = 0;
        int threadsStarted = 0;
//        COUT("---------------------------------------------------------------------------");
//        COUT("num = " << electronsNumber+ionsNumber);
        pair<Particle*,unsigned long long> **threadArgs = new pair<Particle*,unsigned long long>*[Globals::threadNum];
        for(;threadsStarted < Globals::threadNum;++threadsStarted) {
            if (firtsParticleForCurrentThread >= electronsNumber + ionsNumber)
                break;
            threadArgs[threadsStarted] =
                    new pair<Particle*,unsigned long long>(particles + firtsParticleForCurrentThread,
                                                           min(electronsNumber + ionsNumber - firtsParticleForCurrentThread,(unsigned long long)particlesPerThread));
            assert(pthread_create(threads + threadsStarted,NULL,processParticlesArray,threadArgs[threadsStarted]) == 0);
            firtsParticleForCurrentThread += particlesPerThread;
        }

        // wait for all threads
        for(int t = 0;t < threadsStarted;++t)
            pthread_join(threads[t], NULL);

        // clean threads args
        for(int t = 0;t < threadsStarted;++t)
            delete threadArgs[t];

        delete threadArgs;
        delete threads;
    }

    // checking all particles excluding the last one
    int finalizedNumber = 0;
    unsigned long long end = electronsNumber + ionsNumber;
    for(unsigned long long i = 0;i < end;) {
        if (particles[i].ttl <= 0) {
            if (particles[i].behaviour == PARTICLE_WILL_INTERSECT_OBJ)
                ++finalizedNumber;
            finalizeParticle(satelliteObj,particles,electronsNumber,ionsNumber,i);
            if (i != end - 1)
                memcpy(particles + i,particles + end - 1,sizeof(Particle));
            end--;
        } else {
            ++i;
        }
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

    while ((c = getopt (argc, argv, ":vdjamxgt:r:s:f:t:n:i:p:o:c:h:")) != -1) {
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
        case 'j':
            Globals::drawTrajectories = true;
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
        case 'o':
            Globals::modelingType = atoi(optarg);
            break;
        case 'c':
            Globals::initialCharge = strtold(optarg,NULL);
            break;
        case 'h':
            Globals::threadNum = atoi(optarg);
            assert(Globals::threadNum >= 1);
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
    satelliteObj.totalCharge = Globals::initialCharge;

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

        double distanceStep = satelliteObj.radius*2.0*distanceStepCoef;
        timeStep = distanceStep/ELECTRON_VELOCITY_M; // time to do step for particle with average velocity
        verboseFlag && COUT("distanceStep: " << distanceStep << "; timeStep: " << timeStep);
    }

    verboseFlag && COUT("polygons: " << satelliteObj.polygons->size());
    verboseFlag && COUT("center: " << satelliteObj.center);
    verboseFlag && COUT("radius: " << satelliteObj.radius);
    verboseFlag && COUT("capacitance: " << satelliteObj.capacitance());

    // video mode initialization
    if (drawFlag) {
        // set appropriate OpenGL & properties SDL
        int width = 1200;
        int height = 900;
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
        if (modelingFlag && Globals::modelingType != 1) {
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

            if (Globals::pause)
                continue;

            if (modelingFlag) {
                numberOfIntersections += processParticles(satelliteObj,particlesArray,electronsNumber,
                                                          ionsNumber,timeStep,electronsGenerativeSphere,ionsGenerativeSphere);
                elapsedTime += timeStep;
                timeToPrint -= ((intervalInSteps)? 1: timeStep);
                surfaceCharge = satelliteObj.totalPlasmaCurrent*elapsedTime;
                if (timeToPrint <= 0) {
                    cout << satelliteObj.totalPlasmaCurrent << " " << surfaceCharge << " " << surfaceCharge/spacecraftCapacitance
                         << "   " << elapsedTime << "   " << numberOfIntersections*Globals::realToModelNumber << " " << satelliteObj.totalCharge
                         << "   " << electronsNumber << "   " << ionsNumber << endl; // "   " << realEN << "   " << realIN << endl;
                    (timeToPrint = printInterval);
                }
                // producing new particles if necessary
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

