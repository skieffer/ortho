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

// TODO: Cleanup. Is this struct used anywhere?
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

enum ACAFlags {
    ACAHORIZ = 1,
    ACAVERT = 2,
    ACADELIB = 4,
    ACACONN = 8
};

struct ACASeparatedAlignment {
    cola::SeparationConstraint* separation;
    cola::AlignmentConstraint* alignment;
    int rect1;
    int rect2;
    ACAFlags af;
    ShapeObj *shape1;
    ShapeObj *shape2;
};

struct EdgeNode {
    EdgeNode(QRectF src, QRectF tgt)
        : srcRect(src),
          tgtRect(tgt)
    {
        bbox = srcRect.united(tgtRect);
    }
    void setIndices(int si, int ti) {
        srcIndex = si;
        tgtIndex = ti;
    }
    // If e.g. it is a vertical edge, then 'orientation' is VERTICAL,
    // while 'constraintDimension' is HORIZONTAL.
    vpsc::Dim orientation;
    vpsc::Dim constraintDimension;
    int srcIndex;
    int tgtIndex;
    QRectF srcRect;
    QRectF tgtRect;
    QRectF bbox;
};

class ExternalTree {
public:
    ExternalTree(node root, node rootInG, QList<node> nodes, QList<edge> edges,
                 shapemap nodeShapes, connmap edgeConns);
    QString listNodes(void);
    void colourShapes(void);
    ShapeObj *rootShape(void);
    void orientation(ogdf::Orientation orient) { m_orientation = orient; }
    void treeLayout(void);
    void updateShapePositions(void);
    QRectF rootlessBBox(void);
    void translate(QPointF p);
    void placeRootAt(QPointF p);
private:
    Graph *m_graph;
    GraphAttributes *m_ga;
    node m_root;
    node m_rootInG;
    shapemap m_dunnartShapes;
    connmap m_dunnartConns;
    ogdf::Orientation m_orientation;
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
    QList<ShapeObj*> allShapes(void);
    void dfs(QMap<ShapeObj *, BiComp *> endpts, QList<BiComp *> &elements);
    BiComp *fuse(BiComp *other);
    void addStubNodeForTree(ExternalTree *E, QSizeF size);
    void layout(void);
    void updateShapePositions(void);
    void addStubNodeShapesToCanvas(Canvas *canvas);
    cola::CompoundConstraints generateStubEdgeSepCos(vpsc::Dim dim,
        QList<EdgeNode> ens, QMap<node, int> nodeIndices, double gap);
    void translateTrees(void);
    void idealLength(double L) { m_idealLength = L; }
private:
    void postACACola(bool preventOverlaps, double idealLength,
                     QMap<node,int> nodeIndices, cola::CompoundConstraints sepcos);
    Graph *m_graph;
    GraphAttributes *m_ga;
    QList<node> m_cutnodes;
    shapemap m_dunnartShapes;
    connmap m_dunnartConns;
    QMap<node,ExternalTree*> m_stubnodesToTrees;
    QList<edge> m_stubedges;
    bool m_stubNodeShapesHaveBeenAddedToCanvas;
    double m_idealLength;
};

struct ConvTest1 : public cola::TestConvergence
{
    ConvTest1(const double d,const unsigned i)
        : TestConvergence(d,i),
          m_layout(NULL),
          m_count(1),
          name(""),
          minIterations(0)
    {}

    void setLayout(cola::ConstrainedFDLayout *layout) { m_layout = layout; }

    bool test1(const double new_stress, std::valarray<double> & X,
               std::valarray<double> & Y)
    {
        if (m_count < minIterations) {
            return false;
        } else {
            return cola::TestConvergence::operator()(new_stress, X, Y);
        }
    }

    bool operator()(const double new_stress, std::valarray<double> & X,
            std::valarray<double> & Y)
    {
        QString outFName = "Debug-"+name;
        bool converged = test1(new_stress, X, Y);
        cout << "stress="<<new_stress<<" iteration="<<m_count<<endl;
        std::stringstream ss;
        QString f = outFName+QString("-%1").arg(m_count++,4,10,QLatin1Char('0'));
        m_layout->outputInstanceToSVG(f.toStdString());
        return converged;
    }

    cola::ConstrainedFDLayout *m_layout;
    int m_count;
    QString name;
    int minIterations;
};

class ACALayout {
public:
    ACALayout(Graph &G, GraphAttributes &GA);
    QMap<node,int> nodeIndices(void) { return m_nodeIndices; }
    void idealLength(double L) { m_idealLength = L; }
    void preventOverlaps(bool b) { m_preventOverlaps = b; }
    void debugName(QString s) { m_debugName = s; }
    cola::CompoundConstraints ccs(void) { return m_ccs; }
    void run(void);
    void readPositions(Graph &G, GraphAttributes &GA);
    QMap<vpsc::Dim,EdgeNode> generateEdgeNodes(void);
private:
    void initialPositions(void);
    void moveCoincidentNodes(void);
    void initialLayout(void);
    void acaLoopOneByOne(void);
    void acaLoopAllAtOnce(void);
    void initAlignmentState(void);
    void updateAlignmentState(ACASeparatedAlignment *sa);
    ACASeparatedAlignment *chooseSA(void);
    bool createsCoincidence(int src, int tgt, ACAFlags af);
    double deflection(int src, int tgt, ACAFlags af);
    double bendPointPenalty(int src, int tgt, ACAFlags af);
    double leafPenalty(int src, int tgt);
    QMap<node,int> m_nodeIndices; // map node to index of rect representing it
    double m_idealLength;
    bool m_preventOverlaps;
    bool m_addBendPointPenalty;
    bool m_postponeLeaves;
    bool m_useNonLeafDegree;
    bool m_allAtOnce;
    QString m_debugName;
    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;
    cola::CompoundConstraints m_ccs;
    Matrix2d<int> alignmentState;
    QList<ACASeparatedAlignment*> sepAligns;
    QMap<int,int> deg2Nodes; // maps indices of degree-2 nodes to indices of their nbrs
    QList<int> leaves; // indices of rectangles that are leaves
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














































