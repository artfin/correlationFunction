QT += core
QT -= gui

TARGET = MH_sampling
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -O2 -Wall -Wextra

SOURCES += main.cpp \
    ar_he_pes.cpp \
    ar_he_dip_buryak_fit.cpp

HEADERS += \
    ar_he_pes.hpp \
    ar_he_dip_buryak_fit.hpp

