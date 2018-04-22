QT += core
QT -= gui

TARGET = quasibound
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS = -std=c++11 -O2 -Wall -Wextra

LIBS = -lgsl -lgslcblas

SOURCES += main.cpp \
    ar_he_pot.cpp \
    ar_he_pot_derivative.cpp

HEADERS += \
    ar_he_pot.hpp \
    ar_he_pot_derivative.hpp

