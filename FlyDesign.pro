#-------------------------------------------------
#
# Project created by QtCreator 2016-03-28T20:44:52
#
#-------------------------------------------------

QT       -= gui

TARGET = FlyDesign
TEMPLATE = lib

DEFINES += FLYDESIGN_LIBRARY

SOURCES += flydesign.cpp

HEADERS += flydesign.h\
        flydesign_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
