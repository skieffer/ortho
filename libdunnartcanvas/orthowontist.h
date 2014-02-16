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
    void numberShapes(void);
    ShapeObj *rootShape(void);
    void orientation(ogdf::Orientation orient) { m_orientation = orient; }
    void treeLayout(void);
    void updateShapePositions(void);
    void orthogonalRouting(bool b);
    QRectF rootlessBBox(void);
    void translate(QPointF p);
    void placeRootAt(QPointF p);
    bool needsAlignmentOffset(void);
    double alignmentOffset(void);
private:
    Graph *m_graph;
    GraphAttributes *m_ga;
    node m_root;
    node m_rootInG;
    shapemap m_dunnartShapes;
    connmap m_dunnartConns;
    ogdf::Orientation m_orientation;
};

static int cmpEdgeEvent(const void *p1, const void *p2);

class Planarization {
public:
    Planarization(Graph &G, GraphAttributes &GA,
                  QMap<edge,int> alignments, QSizeF dummyNodeSize);
    void addDummyNodeShapesToCanvas(Canvas *canvas);

    enum EdgeType { HTYPE, VTYPE, DTYPE };

    struct Edge {
        Edge(EdgeType t, ogdf::edge e, GraphAttributes *ga) :
            m_etype(t), m_ogdfEdge(e), m_ga(ga) {
            double xs=m_ga->x(e->source()), ys=m_ga->y(e->source());
            double xt=m_ga->x(e->target()), yt=m_ga->y(e->target());
            x0 = min(xs,xt);
            x1 = max(xs,xt);
            y0 = min(ys,yt);
            y1 = max(ys,yt);
        }
        bool sharesAnEndptWith(Edge o) {
            node src = m_ogdfEdge->source();
            node tgt = m_ogdfEdge->target();
            node os = o.m_ogdfEdge->source();
            node ot = o.m_ogdfEdge->target();
            return src==os||src==ot||tgt==os||tgt==ot;
        }
        QPair<bool,QPointF> intersectDiagonals(Edge *d);
        double constCoord(void) {return m_etype==HTYPE ? y0 : x0;}
        double lowerBd(void) {return m_etype==HTYPE ? x0 : y0;}
        double upperBd(void) {return m_etype==HTYPE ? x1 : y1;}
        bool coversY(double y) {return y0 <= y && y <= y1;}
        bool coversX(double x) {return x0 <= x && x <= x1;}
        double x(double y) { return y1 == y0 ? 0 : x0+(x1-x0)*(y-y0)/(y1-y0); }
        double y(double x) { return x1 == x0 ? 0 : y0+(y1-y0)*(x-x0)/(x1-x0); }
        double slope(void) { return (x1-x0)/(y1-y0); }
        double yInt(void) { return y0 - slope()*x0; }
        GraphAttributes *m_ga;
        EdgeType m_etype;
        ogdf::edge m_ogdfEdge;
        double x0, x1, y0, y1;
        int m_openEdgeIndex;
    };

    enum EdgeEventType { HEDGE, VEDGE, DOPENX, DOPENY, DCLOSEX, DCLOSEY };

    struct EdgeEvent {
        EdgeEvent() {}
        EdgeEvent(EdgeEventType t, Edge *e) : m_eetype(t), m_edge(e) {}
        EdgeEventType m_eetype;
        Edge *m_edge;
    };

    struct DummyCross {
        DummyCross(node dn, edge de1s, edge de1t, edge de2s, edge de2t) :
            dNode(dn), dEdge1Src(de1s), dEdge1Tgt(de1t), dEdge2Src(de2s), dEdge2Tgt(de2t) {}
        node dNode;
        edge dEdge1Src;
        edge dEdge1Tgt;
        edge dEdge2Src;
        edge dEdge2Tgt;
    };

private:
    void planarizeHDCrossings(void);
    void planarizeVDCrossings(void);
    void planarizeDDCrossings(void);
    void addDummyCross(Edge *e1, Edge *e2, QPointF p);
    QSizeF m_dummyNodeSize;
    Graph *m_graph;
    GraphAttributes *m_ga;
    QMap<node,node> m_origNodes;
    QMap<edge,edge> m_origEdges;
    QMap<edge,int> m_alignments;
    QList<Edge*> mH;
    QList<Edge*> mV;
    QList<Edge*> mD;
    QList<node> m_dummyNodes;
    QList<edge> m_dummyEdges;
    QList<DummyCross*> m_dummyCrosses;
    shapemap m_dunnartShapes;
};

class BiComp {
public:
    BiComp(void);
    BiComp(QList<node> nodes, QList<edge> edges, QList<node> cutnodes,
           shapemap nodeShapes, connmap edgeConns);
    QList<node> cutnodes(void);
    QString listNodes(void);
    void colourShapes(void);
    void numberShapes(void);
    ShapeObj *getShape(node m);
    QList<ShapeObj*> allShapes(void);
    void dfs(QMap<ShapeObj *, BiComp *> endpts, QList<BiComp *> &elements);
    BiComp *fuse(BiComp *other);
    void addStubNodeForTree(ExternalTree *E, QSizeF size);
    void layout(void);
    void updateShapePositions(void);
    void orthogonalRouting(bool b);
    void addStubNodeShapesToCanvas(Canvas *canvas);
    void addDummyNodeShapesToCanvas(Canvas *canvas);
    cola::CompoundConstraints generateStubEdgeSepCos(vpsc::Dim dim,
        QList<EdgeNode> ens, QMap<node, int> nodeIndices, double gap);
    void translateTrees(void);
    void idealLength(double L) { m_idealLength = L; }
    void nodePadding(double P) { m_nodePadding = P; }
    void dummyNodeSize(QSizeF s) { m_dummyNodeSize = s; }
    double *edgeLengths(QMap<node,int> nodeIndices, std::vector<cola::Edge> colaEdges);
    QList<double> nodePadding(QMap<node,int> nodeIndices);
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
    bool m_dummyNodeShapesHaveBeenAddedToCanvas;
    double m_idealLength;
    double m_nodePadding;
    QSizeF m_dummyNodeSize;
    Planarization *m_planarization;
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
    void nodePadding(double P);
    void nodePadding(QList<double> P);
    void edgeLengths(double *eL) { m_edgeLengths = eL; }
    QMap<node,int> nodeIndices(void) { return m_nodeIndices; }
    void idealLength(double L) { m_idealLength = L; }
    void preventOverlaps(bool b) { m_preventOverlaps = b; }
    void debugName(QString s) { m_debugName = s; }
    std::vector<cola::Edge> colaEdges(void) { return es; }
    cola::CompoundConstraints ccs(void) { return m_ccs; }
    void run(void);
    void readPositions(Graph &G, GraphAttributes &GA);
    QMap<vpsc::Dim,EdgeNode> generateEdgeNodes(void);
    QList<int> delibAlignedWith(int i);
    bool delibAligned(int i, int j) { return alignmentState(i,j) & ACADELIB; }
    bool offsetAlignment(int l, int r, double offset);
    int alignment(edge e);
    QMap<edge,int> alignments(Graph &G);
private:
    void initialPositions(void);
    void moveCoincidentNodes(void);
    void initialLayout(void);
    void acaLoopOneByOne(void);
    void acaLoopAllAtOnce(void);
    void finalLayout(void);
    void initAlignmentState(void);
    void updateAlignmentState(ACASeparatedAlignment *sa);
    ACASeparatedAlignment *chooseSA(void);
    bool createsCoincidence(int src, int tgt, ACAFlags af);
    double deflection(int src, int tgt, ACAFlags af);
    double bendPointPenalty(int src, int tgt, ACAFlags af);
    double leafPenalty(int src, int tgt);
    QMap<node,int> m_nodeIndices; // map node to index of rect representing it
    double m_idealLength;
    double m_nodePadding;
    bool m_preventOverlaps;
    bool m_addBendPointPenalty;
    bool m_postponeLeaves;
    bool m_useNonLeafDegree;
    bool m_allAtOnce;
    QString m_debugName;
    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;
    double *m_edgeLengths;
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














































