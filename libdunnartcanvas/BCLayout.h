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

#include <QList>
#include <QMap>
#include <QPointF>
#include <QRectF>

#include "libogdf/ogdf/basic/Graph_d.h"
#include "libcola/compound_constraints.h"
#include "libvpsc/rectangle.h"
#include "libdunnartcanvas/canvas.h"

using namespace ogdf;

namespace ogdf {
class GraphAttributes;
}

namespace dunnart {

class Canvas;
class ShapeObj;
class Connector;

typedef QMap<node, ShapeObj *> shapemap;
typedef QMap<edge, Connector *> connmap;

class RootedTree;
typedef QList<RootedTree*> treelist;

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

class Chunk {
public:
    virtual void setRelPt(QPointF p) = 0;
    virtual void recursiveLayout(shapemap& origShapes,
                                 node origBaseNode, QPointF cardinal,
                                 QList<node> cutnodes) = 0;
    virtual void recursiveDraw(Canvas *canvas, QPointF p) = 0;
    virtual bool containsOriginalNode(node n) = 0;
    virtual QList<node> getCutNodes(void) = 0;
    virtual void setChildren(QList<Chunk*> children) = 0;
    virtual void setParent(BiComp *bc) = 0;
    virtual void setParentCutNode(node cn) = 0;
    virtual node getParentCutNode(void) = 0;
    QList<Chunk*> findNeighbours(node origCutNode, QList<Chunk*> allChunks);
    virtual QRectF bbox(void) = 0;
};

class RootedTree : public Chunk
{
public:
    RootedTree(QList<node>& nodes, QSet<node>& cutnodes);
    bool containsOriginalNode(node n);
    void setRelPt(QPointF p);
    void recursiveLayout(shapemap& origShapes,
                         node origBaseNode, QPointF cardinal,
                         QList<node> cutnodes);
    void recursiveDraw(Canvas *canvas, QPointF p);
    void constructDunnartGraph(shapemap& origShapes,
                               QPointF cardinal, QList<node> cutnodes);
    void colaTreeLayout(QPointF cardinal);
    void ogdfTreeLayout(QPointF cardinal);
    QList<node> getCutNodes(void);
    static QPointF nearestCardinal(QPointF v);
    void setChildren(QList<Chunk*> children);
    void setParent(BiComp *bc);
    void setParentCutNode(node cn);
    node getParentCutNode(void);
    QPointF baryCentre(void);
    QRectF bbox(void);
    void inferConstraints(void);
    void applyDunnartConstraints(void);
private:
    Graph *m_graph;
    node m_root;
    QMap<node,node> m_nodemap; // maps own nodes to orig. graph nodes
    QList<node> m_cutNodes; // subdomain of nodemap which are cutnodes
    QList<node> m_normalNodes; // subdomain of those which are not
    QMap<edge,edge> m_edgemap; // maps own edges to orig. graph edges

    QList<Chunk*> m_children;
    BiComp *m_parent;
    node m_parentCutNode; // belongs to orig. graph

    shapemap m_origShapeMap;
    shapemap m_ownShapeMap;
    connmap m_ownConnMap;
    QPointF m_basept;
    QPointF m_relpt;

    QList<DunnartConstraint*> m_dunnartConstraints;
    int m_orientation; // in {0,1,2,3} for tree layout
};

class ExternalTree
{
public:
    ExternalTree(QList<node> nodes, QList<edge> edges, node root, node taproot,
                 QMap<node,int> dunnartIDs, int taprootID);
    void treeLayout(void);
    QSizeF getBoundingBoxSize(void);
    node taproot(void);
    //Testing:
    QString listNodes(void);
private:
    void inferConstraints(Canvas::Dimension dim);
    QList<DunnartConstraint*> m_dunnartConstraints;

    Graph *m_graph;
    GraphAttributes *m_graphAttributes;
    QMap<node,int> m_dunnartIDs; // maps own nodes to IDs of corresp. Dunnart shapes
    node m_root;    // a node in this object's m_graph
    node m_tapRoot; // a node in the original graph G
    int m_tapRootID; // ID of Dunnart shape corresp. to m_tapRoot

    QMap<node,ShapeObj*> m_shapeMap; // map from OGDF nodes to Dunnart shapes
};

struct InternalTree
{
    InternalTree(QList<node> nodes, QSet<node> cutnodes,
                 QMap<node,BiComp*> pNodesToBCs, Graph G);
    // Unlike most of the classes, which have their own Graph objects, and
    // maintain mappings from their own nodes to the nodes of the original
    // graph (i.e. the Graph G declared first in BCLayout::ortholayout2),
    // the InternalTree struct just stores nodes of G itself.
    QList<node> nodes;          // all nodes
    QList<node> cutNodes;
    QList<node> nonCutNodes;
    QList<edge> edges;          // all edges
    QList<edge> cutEdges;       // one endpt is a cutnode
    QList<edge> nonCutEdges;    // neither endpt is a cutnode
    QMap<node,BiComp*> nodesToBCs;
};

enum MetaGraphNodeType {
    MetaIntern, MetaBicon, MetaExtern
};

class MetaGraph
{
public:
    MetaGraph(QList<ExternalTree *> XX, QList<BiComp *> BB,
              QList<InternalTree *> II, QMap<node,BiComp*> nodesToBCs);
    void acaLayout(void);
private:
    Graph *m_graph;
    GraphAttributes *m_graphAttributes;
    QMap<node,MetaGraphNodeType> m_nodeTypes;
    // Nodemaps are from own nodes to other objects or their nodes.
    QMap<node,node> m_internNodemap;
    QMap<node,BiComp*> m_biconNodemap;
    QMap<node,ExternalTree*> m_externNodemap;
};

class BiComp : public Chunk
{
public:
    BiComp(void);
    BiComp(QList<edge>& edges, QList<node>& nodes, QList<node>& cutNodes, Graph& G);
    void setRelPt(QPointF p);
    void removeSelf(Graph& G);
    size_t size(void);
    void constructDunnartGraph(shapemap& origShapes,
                               QPointF cardinal, QList<node> cutnodes);
    void colaLayout(void);

    // -----------
    // ACA methods -- for BiComp's own built in ACA layout -- deprecated
    Matrix2d<int> initACA(int N, QMap<node,int> nodeIndices);
    SeparatedAlignment *chooseSA(vpsc::Rectangles rs, Matrix2d<int> &alignmentState,
                                 QMap<node,int> nodeIndices);
    void updateAlignmentState(SeparatedAlignment *sa, Matrix2d<int> &alignmentState);
    bool createsCoincidence(int srcIndex, int tgtIndex, Canvas::AlignmentFlags af,
                            vpsc::Rectangles rs, Matrix2d<int> &alignmentState);
    double deflection(int srcIndex, int tgtIndex, Canvas::AlignmentFlags af,
                      vpsc::Rectangles rs);
    void applyDunnartSepAligns(Canvas *canvas);
    void removeDunnartSepAligns(Canvas *canvas);
    // -----------

    // For ACA layout performed by separate ACALayout object -- preferred method
    void acaLayout(void);
    QSizeF getOGDFBoundingBoxSize(void);

    void improveOrthogonalTopology(void);
    void recursiveLayout(shapemap& origShapes, node origBaseNode,
                         QPointF cardinal, QList<node> cutnodes);
    void recursiveDraw(Canvas *canvas, QPointF p);
    static QPointF nearestCardinal(QPointF v);
    bool containsOriginalNode(node n);
    bool containsOriginalEdge(edge e);
    QList<node> getCutNodes(void);
    QList<node> getAllOriginalNodes(void);
    void setChildren(QList<Chunk*> children);
    void setParent(BiComp *bc);
    void setParentCutNode(node cn);
    node getParentCutNode(void);
    QPointF baryCentre(void);
    QRectF bbox(void);
    bool coincidence(void);
    void jog(double scale);
    static int method;
    void dfs(QMap<node,BiComp*> endpts, QList<BiComp*> &elements);
    BiComp *fuse(BiComp *other);
    ShapeObj *getShapeForOriginalNode(node orig);
private:
    Graph *m_graph;
    GraphAttributes *m_graphAttributes;
    QMap<node,node> m_nodemap; // maps own nodes to orig. graph nodes
    QList<node> m_cutNodes; // subdomain of nodemap which are cutnodes
    QList<node> m_normalNodes; // subdomain of those which are not
    QMap<edge,edge> m_edgemap; // maps own edges to orig. graph edges

    Graph& copyGraph(QMap<node,node>& nodemap);

    // For Dunnart layout:
    shapemap m_origShapeMap;
    shapemap m_ownShapeMap;
    connmap m_ownConnMap;
    QPointF m_basept;
    QPointF m_relpt;
    QList<Chunk*> m_children;
    BiComp *m_parent;
    node m_parentCutNode;
    QList<Chunk*> findCutNodeNeighbours(node origCutNode, bclist& bcs, treelist& trees);
    QList<SeparatedAlignment*> m_sepAligns;
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

class ACALayout
{
public:
    ACALayout(Canvas *canvas, bool selection = false);
    ACALayout(Graph G, GraphAttributes GA);
    ACALayout(QList<ShapeObj*> shapes, QList<Connector*> connectors);
    void setIdealLength(double il);
    void run(void);

    void readLayout(Graph G, GraphAttributes &GA);
private:
    void initAlignmentState(void);
    void updateAlignmentState(ACASeparatedAlignment *sa);
    ACASeparatedAlignment *chooseSA(void);
    bool createsCoincidence(int src, int tgt, ACAFlags af);
    double deflection(int src, int tgt, ACAFlags af);

    void debugOutput(ACASeparatedAlignment *sa);

    QMap<node,int> m_ogdfNodeIndices;

    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;
    double idealLength;
    Matrix2d<int> alignmentState;
    QList<ACASeparatedAlignment*> sepAligns;
};

class BCLayout
{
public:
    BCLayout(Canvas *canvas);
    void ogdfGraph(Graph& G, shapemap& nodeShapes, connmap& edgeConns);

    static void extractSizes(shapemap& nodeShapes, GraphAttributes& GA);
    static void extractSizes(
            QMap<node,node>& newToOldNodes, shapemap& oldNodesToShapes,
            GraphAttributes& newNodesGA);
    static void extractPosAndSize(
            QMap<node,node>& newToOldNodes, shapemap& oldNodesToShapes,
            GraphAttributes& newNodesGA);
    static void extractPosAndSize(shapemap& newNodesToShapes, GraphAttributes& newNodesGA);
    static void injectPositions(shapemap& nodeShapes, GraphAttributes& GA);
    static void injectSizes(shapemap& nodeShapes, GraphAttributes& GA);

    static void buildBFSTree(QList<Chunk*> chunks, Chunk *root,
                             QList<node>& usedCutNodes);

    void applyFM3(void);

    QList<ExternalTree*> removeExternalTrees(Graph &G, shapemap nodeShapes);
    QList<BiComp*> getNontrivialBCs(Graph& G, QSet<node>& cutnodes);
    QMap<node,BiComp*> makeGnodeToBCmap(QList<BiComp*> BB);
    QList<BiComp*> fuseBCs(QList<BiComp*> bicomps);
    QMap<int,node> getConnComps(Graph& G);
    QMap<int,node> getConnComps2(Graph *G2, QMap<node,node>& nodeMapG2ToG);
    Graph *removeBiComps(Graph& G, bclist& bcs, QMap<node,node>& nodeMapNewToOld);
    void orthoLayout(int method);
    void ortholayout2(void);
    void layoutBCTrees(void);

private:
    Canvas *m_canvas;
};

}
