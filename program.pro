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
    file_utils.h

LIBS += -Wall -lGL -lGLU `sdl-config --cflags --libs` -lrt
QMAKE_CXXFLAGS += -std=c++0x -Wno-write-strings
# -Werror=conversion







