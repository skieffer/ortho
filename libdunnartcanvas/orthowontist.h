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
    // inner structs

    struct TreeNode {
        TreeNode(node n, ExternalTree *owner) :
            parent(NULL), ogdfNode(n), isLeaf(false), rank(-1),
            isomNumber(-1), owningTree(owner) {}
        TreeNode *parent;
        node ogdfNode;
        bool isLeaf;
        int rank;
        // --------------------------------------------------
        // For use with AHU rooted tree isom. algorithm -----
        int isomNumber;
        QList<int> isomTuple;

        QString isomTupleString(void) const {
            QStringList L;
            foreach (int n, isomTuple) {
                L.append(QString::number(n));
            }
            return L.join(",");
        }

        //bool operator <(TreeNode *other);
        // --------------------------------------------------
        ExternalTree *owningTree;
        QList<TreeNode*> kids;
        void setCentre(QPointF p);
    };

    struct TreeIsomClass {
        TreeIsomClass(ExternalTree *r, int num) :
            numMembers(num), rep(r) {}
        int numMembers;
        ExternalTree *rep;
        bool operator <(const TreeIsomClass &other);
    };

    // public member functions
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

    void translate2(QPointF p);
    void rotate2(ogdf::Orientation ori);
private:
    Graph *m_graph;
    GraphAttributes *m_ga;
    node m_root;
    node m_rootInG;
    shapemap m_dunnartShapes;
    connmap m_dunnartConns;
    ogdf::Orientation m_orientation;

    QMap<node,TreeNode*> m2_nodesToTreeNodes;
    QMap<int,TreeNode*> m2_ranks; // multmap
    QList<TreeNode*> m2_leaves;

    QList<TreeNode*> leavesOfRank(int r) {
        QList<TreeNode*> L;
        foreach (TreeNode *tn, m2_ranks.values(r)) {
            if (tn->isLeaf) L.append(tn);
        }
        return L;
    }

    QList<TreeNode*> nonleavesOfRank(int r) {
        QList<TreeNode*> L;
        foreach (TreeNode *tn, m2_ranks.values(r)) {
            if (!tn->isLeaf) L.append(tn);
        }
        return L;
    }

    TreeNode* m2_root;
    QList<ExternalTree*> cTrees2(void);
    QList<TreeIsomClass*> getIsomClasses2(QList<ExternalTree*> trees);
    QList<TreeIsomClass*> getIsomClasses2(void);
public:
    QString computeIsomString2(void);
    bool symmetricLayout2(double g);
    int m2_depth; // number of ranks
    int m2_breadth; // max num nodes in any single rank
    int m2_numNodes; // total number of nodes in tree
    QString m2_isomString;
    bool m2_actuallySymmetric;
};

static int cmpEdgeEvent(const void *p1, const void *p2);

struct OrdAlign {
    OrdAlign(cola::SeparationConstraint *s,
             cola::AlignmentConstraint *a) :
        sep(s), algn(a) {}
    cola::SeparationConstraint *sep;
    cola::AlignmentConstraint *algn;
};

class Planarization;

struct SortableFace {
    SortableFace(face f0, node r0, Planarization *P0) :
        f(f0), r(r0), P(P0) {}
    face f;
    node r;
    Planarization *P;
};

static bool cmpTreesBySize(QPair<node,Planarization*> r1, QPair<node,Planarization*> r2);
static bool cmpFaces(SortableFace *s1, SortableFace *s2);

class Planarization {
public:
    Planarization(Graph &G, GraphAttributes &GA,
                  QMap<edge,int> alignments, QSizeF avgNodeSize,
                  shapemap nodeShapes, bool orthoDiagonals = false);
    void addDummyNodeShapesToCanvas(Canvas *canvas);
    bool cmpTreesBySize2(node r1, node r2);
    void defineRootNodes(QList<node> roots);
    void assignTrees(QMap<node,ExternalTree*> trees);
    void chooseFDTreeFaces(void);
    void chooseCombTreeFaces(void);
    void chooseGreedyTreeFaces(void);
    void addBendsForDiagonalEdges(void);
    double areaOfFace(face f);
    double areaOfFace2(face f);
    double areaOfFace3(face f);
    double areaOfPolygon(QList<QPointF> pts);
    double areaOfTriangle(QList<QPointF> pts);
    double areaOfTriangle(QPointF A, QPointF B, QPointF C);
    QRectF nodeRect(node n);
    QRectF bboxWithoutTrees(void);
    void setTreeSizes(QMap<node,QSizeF> sizes) { origRootToTreeSize = sizes; }
    void layoutTreeForRoot(ExternalTree *E, node root);
    void translateTree(ExternalTree *E, node root);
    void translateNodes(Graph &G, GraphAttributes &GA);
    void writeOutGraphWithStubs(QString fn);
    ogdf::Orientation treeOrientation(node root);
    void idealLength(double L) { m_idealLength = L; }
    void delEdge(edge e) { m_graph->delEdge(e); }
    void addDummyEdge(node a, node b, int alignment) {
        edge e = m_graph->newEdge(a,b);
        m_dummyEdges.append(e);
        m_alignments.insert(e,alignment);
    }
    QString nodeIDString(node n) {
        QString id = "";
        if (m_dummyNodes.contains(n)) {
            id = "d" + QString::number(m_dummyNodes.indexOf(n));
        } else if (m_stubNodes.contains(n)) {
            id = "s" + QString::number(m_stubNodes.indexOf(n));
        } else if (m_bendNodes.contains(n)) {
            id = "b" + QString::number(m_bendNodes.indexOf(n));
        } else {
            int i = m_dunnartShapes.value(n)->internalId();
            id = QString::number(i);
        }
        return id;
    }

    Avoid::Polygon nodeAvoidPolygon(node n);

    bool thereAreEdgesBetweenNodes(node s, node t);
    void expand(int steps);
    void removeOverlaps(void);
    void distribWithNbrStress(void);
    QString filename;
    QMap<node,QPointF> origRootToStubPos;
    QMap<node,QSizeF> origRootToTreeSize;

    struct AreaEdge;

    struct AreaPoint {
        AreaPoint(double u, double v) :
            x(u), y(v), e0(NULL), e1(NULL) {}
        double x;
        double y;
        AreaEdge *e0;
        AreaEdge *e1;
    };

    struct AreaEdge {
        AreaEdge(AreaPoint *u0, AreaPoint *u1) :
            v0(u0), v1(u1), open(false), top(NULL), twin(NULL) {}
        AreaPoint *v0;
        AreaPoint *v1;
        bool open;
        AreaPoint *top;
        AreaEdge *twin;
        double x(double y) {
            double x0 = v0->x, y0 = v0->y;
            double x1 = v1->x, y1 = v1->y;
            assert(y0!=y1);
            double m = (x1-x0)/(y1-y0);
            return x0 + m*(y-y0);
        }
        AreaPoint *p(double y) {
            double x0 = x(y);
            return new AreaPoint(x0,y);
        }
    };

    double areaOfPolygon(QList<AreaPoint*> pts);
    void rotateAwayHorizontals(QList<AreaPoint*>& pts);
    double triTrapArea(AreaPoint *a, AreaPoint *b, AreaPoint *c, AreaPoint *d);

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
        /***
          * Say which cardinal directions are free for this node in this face.
          * We use ints 0, 1, 2, 3, corresp. to powers of sqrt(-1).
          */
        QSet<int> freeCardinals(face f, GraphAttributes &GA) {
            double tolerance = pi/36.0;
            QList<adjEntry> aes = aesPerFace.values(f);
            adjEntry ae = aes.at(0);
            adjEntry ae0 = ae->theNode()==m_node ? aes.at(1) : ae;
            adjEntry ae1 = ae->theNode()==m_node ? ae : aes.at(1);
            // Now ae0 is the one that has m_node at its head, and ae1 has it at its tail.
            // So since the face is always on the (left-hand side?), ae0 will yield the lower
            // bound on the available angles, and ae1 the upper bound.
            node c = ae1->twinNode();
            node b = m_node;
            node a = ae0->theNode();
            double ax = GA.x(a), ay = GA.y(a);
            double bx = GA.x(b), by = GA.y(b);
            double cx = GA.x(c), cy = GA.y(c);
            double ub = atan2(by-ay,ax-bx), lb = atan2(by-cy,cx-bx); // the y's are negated since graphics y-axis goes downward
            if (ub <= lb) ub += 2*pi;
            // Now lb and ub are in the range from -pi to +3*pi.
            QSet<int> free;
            for (int k = -2; k <= 6; k++) {
                double a = k*pi/2;
                double aFlat = a - tolerance, aSharp = a + tolerance;
                int direc = (k+4) % 4;
                if (aFlat >= lb && aSharp <= ub) free.insert(direc);
            }
            return free;
        }

    private:
        /***
          * Computes a normal direction vector pointing perpendicular to
          * the adjEntry ae and into the face it belongs to. This relies
          * on an adjEntry for a face always being oriented so that the
          * the face lies on the right.
          */
        QPointF normalIntoFace(adjEntry ae, GraphAttributes &GA) {
            node src = ae->theNode();
            node tgt = ae->twinNode();
            double sx = GA.x(src), sy = GA.y(src);
            double tx = GA.x(tgt), ty = GA.y(tgt);
            // Let u = t - s.
            double ux = tx - sx, uy = ty - sy;
            // Rotate u 90 deg clockwise to get v.
            double vx = uy, vy = -ux;
            // Normalise.
            double vl = sqrt(vx*vx + vy*vy);

            if (vl==0) {
                qDebug() << "foo";
            }

            vx /= vl; vy /= vl;
            return QPointF(vx,vy);
        }
    };

//private:
    void findHDHVCrossings(void);
    void findVDCrossings(void);
    void findDDCrossings(void);
    void processCrossings(void);
    void simplePlanarize(void);
    void simplerPlanarize(void);
    void findExternalFace(void);
    void mapNodesToFaces(void);
    void computeMinNodeSep(void);
    void computeFreeSides(void);
    void offsetAlignment(node r, node s, double offset);
    bool areAligned(node r, node s);
    QList<Edge*> addDummyCross(Edge *e1, Edge *e2, QPointF p);
    void addCrossing(Edge *e1, Edge *e2, QPointF p);
    void indexNodesAndEdges(void);
    vpsc::Rectangle *vpscNodeRect(node n, bool doubleSize = false);
    double edgeLengthForNodes(node s, node t);
    cola::SeparationConstraint *sepCoForNodes(vpsc::Dim dim, node s, node t, double gap);
    OrdAlign *ordAlignForNodes(node s, node t, ACAFlags af, double offset = 0);
    cola::CompoundConstraints ordAlignsForEdges(void);
    QList<EdgeNode*> genEdgeNodesForFace(ACAFlags af0, face f);
    cola::CompoundConstraints genNodeEdgeSepCos(vpsc::Dim dim,
                                                QList<node> ns,
                                                QList<EdgeNode*> ens,
                                                double gap);
    cola::CompoundConstraints stubStubOP(void);
    cola::CompoundConstraints faceLiftForNode(face f0, node s0, double gap);
    bool isAligned(edge e) {
        ACAFlags af = (ACAFlags) m_alignments.value(e);
        return (af&ACAHORIZ) || (af&ACAVERT);
    }
    QSizeF m_avgNodeSize;
    double m_avgNodeDim;
    QSizeF m_dummyNodeSize;
    QSizeF m_stubNodeSize;
    double m_idealLength;
    Graph *m_graph;
    GraphAttributes *m_ga;
    QMap<node,node> m_origNodes;
    QMap<edge,edge> m_origEdges;
    QMap<node,QSet<int> > m_freeSides;
    QList<Edge*> mH;
    QList<Edge*> mV;
    QList<Edge*> mD;
    QList<node> m_normalNodes;
    QList<node> m_dummyNodes;
    QList<node> m_rootNodes;
    QList<node> m_stubNodes;
    QList<node> m_bendNodes;
    QList<edge> m_normalEdges;
    QList<edge> m_dummyEdges;
    QList<DummyCross*> m_dummyCrosses;
    shapemap m_dunnartShapes;
    CombinatorialEmbedding *m_comb;
    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;
    face m_extFace;
    QMap<node,NodeCombStruct*> m_nodeComb;
    double m_minNodeSep;
    QMap<edge,int> m_alignments;
    QMap<node,node> m_rootsToStubs;
    QMap<node,ExternalTree*> m_rootsToTrees;
    QList< QPair<node,node> > m_deletedDiagonals;

    // stub node to int in 0,1,2,3 to indicate cardinal
    // Domain consists of precisely those stub nodes
    // whose edges are to be aligned.
    QMap<node,int> m_stubAlignments;

    QMap<node,face> m_faceAssigns; // roots to faces
    QMap<node,int> m_nodeIndices; // map node to index of rect representing it
    QMap<edge,int> m_edgeIndices; // only for normal and dummy edges
};

static bool cmpAreaPointsByY(Planarization::AreaPoint *a,
                             Planarization::AreaPoint *b);

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
    QString filename;
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
    QMap<node,ExternalTree*> m2_rootsToTrees;
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
    QString filename;
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














































