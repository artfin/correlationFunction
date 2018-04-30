QT += core
QT -= gui

TARGET = hep_fixed_r
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    ar_he_pes.cpp \
    ar_he_dip_buryak_fit.cpp \
    diatomic_equations_2d.cpp \
    ar_he_pes_der.cpp \
    gear.cpp \
    vmblock.cpp \
    awp.cpp \
    basis.cpp \
    fgauss.cpp \
    diatomic_equations_3d.cpp

QMAKE_CXXFLAGS += -std=c++11 -O2 -Wall -Wextra -g

INCLUDEPATH += -I /home/artfin/Downloads/hep-mc-0.5/include

INCLUDEPATH += -lgsl -lgscblas

HEADERS += \
    ar_he_pes.hpp \
    ar_he_dip_buryak_fit.hpp \
    diatomic_equations_2d.hpp \
    ar_he_pes_der.hpp \
    gear.hpp \
    vmblock.hpp \
    basis.hpp \
    awp.hpp \
    constants.hpp \
    diatomic_equations_3d.hpp
