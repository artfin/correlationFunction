TEMPLATE = SUBDIRS

QT += core
QT -= gui

QMAKE_CXX = mpic++
QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -O2 -g -lstdc++

TARGET = correlationFunction
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    filereader.cpp \
    ar_he_pot.cpp \
    trajectory.cpp \
    awp.cpp \
    basis.cpp \
    gear.cpp \
    vmblock.cpp \
    matrix_he_ar.cpp \
    ar_he_pot_derivative.cpp \
    ar_he_dip_buryak_fit.cpp \
    fgauss.cpp \
    sampling/metropolis/metropolis.cpp

LIBS += -lmpi -pthread -lmpi_cxx
LIBS += -lgsl -lgslcblas

HEADERS += \
    awp.hpp \
    basis.hpp \
    gear.hpp \
    vmblock.hpp \
    matrix_he_ar.hpp \
    ar_he_pot.hpp \
    ar_he_pot_derivative.hpp \
    ar_he_dip_buryak_fit.hpp \
    constants.hpp \
    filereader.hpp \
    parameters.hpp \
    tags.hpp \
    trajectory.hpp
