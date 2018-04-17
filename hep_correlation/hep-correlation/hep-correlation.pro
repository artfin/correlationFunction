QT += core
QT -= gui

TARGET = hep-correlation
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -lm -Wall -O2

INCLUDEPATH += -I /home/artfin/Downloads/hep-mc-0.5/include

SOURCES += main.cpp \
    ar_he_dip_buryak_fit.cpp \
    ar_he_pot.cpp \
    ar_he_pot_derivative.cpp \
    awp.cpp \
    basis.cpp \
    gear.cpp \
    vmblock.cpp \
    fgauss.cpp

HEADERS += \
    ar_he_dip_buryak_fit.hpp \
    ar_he_pot.hpp \
    ar_he_pot_derivative.hpp \
    awp.hpp \
    basis.hpp \
    constants.hpp \
    gear.hpp \
    vmblock.hpp

