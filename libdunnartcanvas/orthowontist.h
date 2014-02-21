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
#include "libogdf/ogdf/basic/CombinatorialEmbedding.h"
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
                  QMap<edge,int> alignments, QSizeF dummyNodeSize,
                  shapemap nodeShapes);
    void addDummyNodeShapesToCanvas(Canvas *canvas);
    void defineRootNodes(QList<node> roots);
    void chooseFDTreeFaces(void);
    void chooseCombTreeFaces(void);
    void idealLength(double L) { m_idealLength = L; }
    void delEdge(edge e) { m_graph->delEdge(e); }
    void addDummyEdge(node a, node b) {
        edge e = m_graph->newEdge(a,b);
        m_dummyEdges.append(e);
    }
    QString nodeIDString(node n) {
        QString id = "";
        if (m_dummyNodes.contains(n)) {
            id = "d" + QString::number(m_dummyNodes.indexOf(n));
        } else {
            int i = m_dunnartShapes.value(n)->internalId();
            id = QString::number(i);
        }
        return id;
    }

    struct Edge;

    struct Intersection {
        Intersection(Edge *o, QPointF p) : otherEdge(o), point(p) {}
        Edge *otherEdge;
        QPointF point;
    };

    enum EdgeType { HTYPE, VTYPE, DTYPE };

    struct Edge {
        Edge(EdgeType t, ogdf::edge e, GraphAttributes *ga, Planarization *p) :
            m_etype(t), m_ogdfEdge(e), m_ga(ga), m_plan(p) {
            x0=m_ga->x(e->source()), y0=m_ga->y(e->source());
            x1=m_ga->x(e->target()), y1=m_ga->y(e->target());
            xmin = min(x0,x1);
            xmax = max(x0,x1);
            ymin = min(y0,y1);
            ymax = max(y0,y1);
        }
        bool sharesAnEndptWith(Edge o) {
            node src = m_ogdfEdge->source();
            node tgt = m_ogdfEdge->target();
            node os = o.m_ogdfEdge->source();
            node ot = o.m_ogdfEdge->target();
            return src==os||src==ot||tgt==os||tgt==ot;
        }
        void connectCrossings(void);
        Edge *rejectHalf(ogdf::edge e1, ogdf::edge e2);
        QPair<bool,QPointF> intersect(Edge *f);
        QPair<bool,QPointF> intersectDiagonals(Edge *d);
        void addIntersection(Edge *o, QPointF p) { m_intersections.append(new Intersection(o,p)); }
        void addCrossing(node n);
        void processCrossing(Edge *e);
        double constCoord(void) {return m_etype==HTYPE ? y0 : x0;}
        double lowerBd(void) {return m_etype==HTYPE ? xmin : ymin;}
        double upperBd(void) {return m_etype==HTYPE ? xmax : ymax;}
        bool coversY(double y) {return ymin <= y && y <= ymax;}
        bool coversX(double x) {return xmin <= x && x <= xmax;}
        double x(double y) { return y1 == y0 ? 0 : x0+(x1-x0)*(y-y0)/(y1-y0); }
        double y(double x) { return x1 == x0 ? 0 : y0+(y1-y0)*(x-x0)/(x1-x0); }
        double slope(void) { return (x1-x0)/(y1-y0); }
        double yInt(void) { return y0 - slope()*x0; }
        Planarization *m_plan;
        GraphAttributes *m_ga;
        EdgeType m_etype;
        ogdf::edge m_ogdfEdge;
        double x0, y0, x1, y1; // coords of endpts
        double xmin, xmax, ymin, ymax; // intervals covered
        int m_openEdgeIndex;
        QList<Intersection*> m_intersections;
        QList<node> m_crossings;
    };

    enum EdgeEventType { HEDGE, VEDGE, DOPENX, DOPENY, DCLOSEX, DCLOSEY };

    struct EdgeEvent {
        EdgeEvent() {}
        EdgeEvent(EdgeEventType t, Edge *e) : m_eetype(t), m_edge(e) {}
        EdgeEventType m_eetype;
        Edge *m_edge;
    };

    struct DummyCross {
        // For each dummy node added, this struct records how the
        // four edges pair off, and thus how the neighbours of the
        // dummy node were originally connected.
        DummyCross(node dn, edge de1s, edge de1t, edge de2s, edge de2t) :
            dNode(dn), dEdge1Src(de1s), dEdge1Tgt(de1t), dEdge2Src(de2s), dEdge2Tgt(de2t) {}
        node dNode;
        edge dEdge1Src;
        edge dEdge1Tgt;
        edge dEdge2Src;
        edge dEdge2Tgt;
    };

    struct NodeCombStruct {
        // Allows to associate with a node the faces to which it
        // is adjacent, and the two adjEntries corresponding to
        // each face.
        NodeCombStruct(node n) : m_node(n) {}
        node m_node;
        QMap<face,adjEntry> aesPerFace;
        QList<face> faces(void) { return aesPerFace.uniqueKeys(); }
        QPair<node,node> nbrs(face f) {
            QList<adjEntry> aes = aesPerFace.values(f);
            QPair<node,node> p;
            // first one
            node n = aes.at(0)->theNode();
            node t = aes.at(0)->twinNode();
            p.first = m_node == n ? t : n;
            // second one
            n = aes.at(1)->theNode();
            t = aes.at(1)->twinNode();
            p.second = m_node == n ? t : n;
            return p;
        }
        /***
          * Computes a normal vector pointing into the face, and bisecting the
          * angle between the two adjEntries incident at this node.
          */
        QPointF normalIntoFace(face f, GraphAttributes &GA) {
            QList<adjEntry> aes = aesPerFace.values(f);
            QPointF p0 = normalIntoFace(aes.at(0), GA);
            QPointF p1 = normalIntoFace(aes.at(1), GA);
            QPointF n = p0 + p1;
            double nx = n.x(), ny = n.y();
            double nl = sqrt(nx*nx+ny*ny);
            return QPointF(nx/nl, ny/nl);
        }
    private:
        /***
          * Computes a normal direction vector pointing perpendicular to
          * the adjEntry ae and into the face it belongs to. This relies
          * on an adjEntry for a face always being oriented so that the
          * the face lies on the right.
          */
        QPointF normalIntoFace(adjEntry ae, GraphAttributes &GA) {
            //node src = ae->theNode();
            //node tgt = ae->twinNode();
            node src = ae->twinNode();
            node tgt = ae->theNode();
            double sx = GA.x(src), sy = GA.y(src);
            double tx = GA.x(tgt), ty = GA.y(tgt);
            // Let u = t - s.
            double ux = tx - sx, uy = ty - sy;
            // Rotate u 90 deg clockwise to get v.
            double vx = uy, vy = -ux;
            // Normalise.
            double vl = sqrt(vx*vx + vy*vy);
            vx /= vl; vy /= vl;
            return QPointF(vx,vy);
        }
    };

private:
    void findHDHVCrossings(void);
    void findVDCrossings(void);
    void findDDCrossings(void);
    void processCrossings(void);
    void simplePlanarize(void);
    void simplerPlanarize(void);
    void findExternalFace(void);
    void mapNodesToFaces(void);
    void computeMinNodeSep(void);
    QList<Edge*> addDummyCross(Edge *e1, Edge *e2, QPointF p);
    void addCrossing(Edge *e1, Edge *e2, QPointF p);
    QSizeF m_dummyNodeSize;
    double m_idealLength;
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
    QList<node> m_rootNodes;
    shapemap m_dunnartShapes;
    CombinatorialEmbedding *m_comb;
    QMap<node,int> m_nodeIndices; // map node to index of rect representing it
    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;
    face m_extFace;
    QMap<node,NodeCombStruct*> m_nodeComb;
    double m_minNodeSep;
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
    void layout2(void);
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

    // 2nd version
    void noteRoot(ExternalTree *E);
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

    // 2nd version
    QList<node> m2_rootNodes;
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
    void run2(QList<CanvasItem*> items);
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














































