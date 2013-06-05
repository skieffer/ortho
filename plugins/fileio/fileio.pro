

TEMPLATE = subdirs

SUBDIRS = builtinsvg \
            builtingml \
            builtinlayout \
	    builtintgf \
	    json

packagesExist(libcgraph) {
    SUBDIRS += graphviz
}

CONFIG += ordered
