QT += core
QT -= gui

TARGET = dipole_integral
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -lm -O2

INCLUDEPATH += -I /home/artfin/Downloads/hep-mc-0.5/include

SOURCES += main.cpp \
    ar_he_pot.cpp \
    ar_he_dip_buryak_fit.cpp

HEADERS += \
    ar_he_dip_buryak_fit.hpp \
    ar_he_pot.hpp

