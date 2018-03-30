QT += core
QT -= gui

TARGET = sampling
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -O2 -g -lstdc++

SOURCES += main.cpp
