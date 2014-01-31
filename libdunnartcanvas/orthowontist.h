/*
 * Dunnart - Constraint-based Diagram Editor
 *
 * Copyright (C) 2003-2007  Michael Wybrow
 * Copyright (C) 2006-2008  Monash University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301, USA.
 *
 *
 * Author(s): Steven Kieffer  <http://skieffer.info/>
*/

#include <string>
#include <sstream>
#include <QList>
#include <QMap>
#include <QPointF>
#include <QRectF>

#include "libogdf/ogdf/basic/Graph_d.h"
#include "libogdf/ogdf/tree/TreeLayout.h"
#include "libcola/compound_constraints.h"
#include "libcola/cola.h"
#include "libvpsc/rectangle.h"
#include "libdunnartcanvas/canvas.h"
//#include "libdunnartcanvas/BCLayout.h"

using namespace ogdf;
using namespace dunnart;

namespace dunnart {
class Canvas;
class ShapeObj;
class Connector;
}

namespace ogdf {
class GraphAttributes;
}

namespace ow {

typedef QMap<node, ShapeObj *> shapemap;
typedef QMap<edge, Connector *> connmap;


class BiComp;
typedef QList<BiComp*> bclist;

struct SeparatedAlignment {
    cola::SeparationConstraint* separation;
    cola::AlignmentConstraint* alignment;
    int rect1;
    int rect2;
    Canvas::AlignmentFlags af;
    ShapeObj *shape1;
    ShapeObj *shape2;
};

enum ConstraintType {
    ALIGNMENT, SEPARATION, DISTRIBUTION
};

struct DunnartConstraint {
    ConstraintType type;
    Canvas::Dimension dim;
    QList<CanvasItem*> items;
    double minSep;
};

class ExternalTree {
public:
    ExternalTree(node root, node rootInG, QList<node> nodes, QList<edge> edges,
                 shapemap nodeShapes, connmap edgeConns);
    QString listNodes(void);
    void colourShapes(void);
private:
    Graph *m_graph;
    GraphAttributes *m_ga;
    node m_root;
    node m_rootInG;
    shapemap m_dunnartShapes;
    connmap m_dunnartConns;
};

class BiComp {
public:
    BiComp(void);
    BiComp(QList<node> nodes, QList<edge> edges, QList<node> cutnodes,
           shapemap nodeShapes, connmap edgeConns);
    QList<node> cutnodes(void);
    QString listNodes(void);
    void colourShapes(void);
    ShapeObj *getShape(node m);
    void dfs(QMap<ShapeObj *, BiComp *> endpts, QList<BiComp *> &elements);
    BiComp *fuse(BiComp *other);
private:
    Graph *m_graph;
    GraphAttributes *m_ga;
    QList<node> m_cutnodes;
    shapemap m_dunnartShapes;
    connmap m_dunnartConns;
};

class Orthowontist {
public:
    Orthowontist(Canvas *canvas);
    void run1(QList<CanvasItem*> items);
private:
    void buildOGDFGraph(CanvasItemsList items,
            Graph& G, GraphAttributes& GA, shapemap &nodeShapes, connmap &edgeConns);
    void removeExternalTrees(QList<ExternalTree*> &EE, Graph &G,
                             shapemap &nodeShapes, connmap &edgeConns);
    void buildNBCs(QList<BiComp*> &BB, QSet<node> &cutnodes, Graph &G,
                             shapemap &nodeShapes, connmap &edgeConns);
    QList<BiComp*> fuseBCs(QList<BiComp*> bicomps);
    Canvas *m_canvas;
};

}














































