QT       += core
QT       -= gui

TARGET = program
CONFIG   += console
CONFIG   += DEBUG
CONFIG   -= app_bundle

TEMPLATE = app
SOURCES += main.cpp \
    geometry_utils.cpp \
    types.cpp \
    time_utils.cpp \
    constants.cpp \
    graphics_utils.cpp \
    file_utils.cpp

HEADERS += \
    geometry_utils.h \
    types.h \
    time_utils.h \
    constants.h \
    data_utils.h \
    graphics_utils.h \
    file_utils.h \
    fortran_modules.h

LIBS += -Wall -lGL -lGLU `sdl-config --cflags --libs` -lrt -lassimp -lgfortran fortran_modules/Kul.o
QMAKE_CXXFLAGS += -std=c++0x -Wno-write-strings
# -Werror=conversion

OTHER_FILES += \
    fortran_modules/gf \
    fortran_modules/t1Lapl3D.f90 \
    fortran_modules/ReadFile.f90 \
    fortran_modules/Kul.f90 \
    fortran_modules/t2Lapl3D.f90 \
    fortran_modules/t0Lapl3D.f90 \
    fortran_modules/Symmetry.f90 \
    fortran_modules/Lapl_3D.f90







