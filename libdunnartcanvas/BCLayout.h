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
    void setParentCutNode(node cn);
    node getParentCutNode(void);
    QPointF baryCentre(void);
    QRectF bbox(void);
private:
    Graph *m_graph;
    node m_root;
    QMap<node,node> m_nodemap; // maps own nodes to orig. graph nodes
    QList<node> m_cutNodes; // subdomain of nodemap which are cutnodes
    QList<node> m_normalNodes; // subdomain of those which are not
    QMap<edge,edge> m_edgemap; // maps own edges to orig. graph edges

    QList<Chunk*> m_children;
    node m_parentCutNode; // belongs to orig. graph

    shapemap m_origShapeMap;
    shapemap m_ownShapeMap;
    connmap m_ownConnMap;
    QPointF m_basept;
    QPointF m_relpt;
};

class BiComp : public Chunk
{
public:
    BiComp(QList<edge>& edges, QList<node>& nodes, QList<node>& cutNodes, Graph& G);
    void setRelPt(QPointF p);
    void removeSelf(Graph& G);
    size_t size(void);
    void constructDunnartGraph(shapemap& origShapes,
                               QPointF cardinal, QList<node> cutnodes);
    void colaLayout(void);
    void improveOrthogonalTopology(void);
    void recursiveLayout(shapemap& origShapes, node origBaseNode,
                         QPointF cardinal, QList<node> cutnodes);
    void recursiveDraw(Canvas *canvas, QPointF p);
    static QPointF nearestCardinal(QPointF v);
    bool containsOriginalNode(node n);
    bool containsOriginalEdge(edge e);
    QList<node> getCutNodes(void);
    void setChildren(QList<Chunk*> children);
    void setParentCutNode(node cn);
    node getParentCutNode(void);
    QPointF baryCentre(void);
    QRectF bbox(void);
    bool coincidence(void);
    void jog(double scale);
    static int method;
private:
    Graph *m_graph;
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
    node m_parentCutNode;
    QList<Chunk*> findCutNodeNeighbours(node origCutNode, bclist& bcs, treelist& trees);
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

    QList<BiComp*> getNontrivialBCs(Graph& G, QSet<node>& cutnodes);
    QMap<int,node> getConnComps(Graph& G);
    QMap<int,node> getConnComps2(Graph *G2, QMap<node,node>& nodeMapG2ToG);
    Graph *removeBiComps(Graph& G, bclist& bcs, QMap<node,node>& nodeMapNewToOld);
    void orthoLayout(int method);

    void layoutBCTrees(void);

private:
    Canvas *m_canvas;
};

}
