
QT           += xml svg
TEMPLATE      = lib
CONFIG       += qt plugin
TARGET        = $$qtLibraryTarget(fileio_json)

include(../../../common_options.qmake)
include(../fileio_plugin_options.pri)

HEADERS       =
SOURCES       = plugin.cpp

