######################################################################
# Automatically generated by qmake (2.01a) Mon Sep 28 11:05:54 2009
######################################################################

TEMPLATE = app
TARGET = ex2
DESTDIR = ../local/bin
DEPENDPATH += .
INCLUDEPATH += . ..  /usr/local/include

LIBS += -L/usr/local/lib -lQatPlotWidgets -lQatPlotting -lQatDataModeling  -lQatDataAnalysis -lQatGenericFunctions -lgsl -lgslcblas -ldl

CONFIG += qt release c++11
QT     += widgets

# Input
SOURCES += *.cpp

QMAKE_DEL_FILE=rm -rf
QMAKE_DISTCLEAN += ../local
INCLUDEPATH += /usr/include/eigen3


mac {
  CONFIG -= app_bundle
  INCLUDEPATH -=  /usr/include/eigen3
  INCLUDEPATH += /usr/local/include/eigen3
}

