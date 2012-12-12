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
    virtual void recursiveLayout(
            shapemap origShapes, bclist bcs, treelist trees,
            node origBaseNode, QPointF cardinal) = 0;
};

class RootedTree : public Chunk
{
public:
    RootedTree(QList<node> nodes, Graph G);
    bool containsOriginalNode(node n);
    void setRelPt(QPointF p);
    void recursiveLayout(shapemap origShapes, bclist bcs, treelist trees,
                         node origBaseNode, QPointF cardinal);
    void constructDunnartGraph(shapemap origShapes);
private:
    Graph *m_graph;
    node m_root;
    QMap<node,node> m_nodemap; // maps own nodes to orig. graph nodes
    QMap<edge,edge> m_edgemap; // maps own edges to orig. graph edges

    shapemap m_origShapeMap;
    shapemap m_ownShapeMap;
    connmap m_ownConnMap;
    QPointF m_basept;
    QPointF m_relpt;
};

class BiComp : public Chunk
{
public:
    BiComp(QList<edge> edges, QList<node> nodes, QList<node> cutNodes);
    void removeSelf(Graph& G);
    size_t size(void);
    void constructDunnartGraph(shapemap origShapes);
    void improveOrthogonalTopology(void);
    void recursiveLayout(shapemap origShapes, bclist bcs, treelist trees,
                         node origBaseNode, QPointF cardinal);
    static QPointF nearestCardinal(QPointF v);
    bool containsOriginalNode(node n);
    void setRelPt(QPointF p);
private:
    Graph *m_graph;
    QMap<node,node> m_nodemap; // maps own nodes to orig. graph nodes
    QList<node> m_cutNodes; // subdomain of nodemap which are cutnodes
    QList<node> m_normalNodes; // subdomain of those which are not
    QMap<edge,edge> m_edgemap; // maps own edges to orig. graph edges

    // For Dunnart layout:
    shapemap m_origShapeMap;
    shapemap m_ownShapeMap;
    connmap m_ownConnMap;
    QPointF m_basept;
    QPointF m_relpt;
    QList<Chunk*> m_children;
    QPointF baryCentre(void);
    QList<Chunk*> findCutNodeNeighbours(node origCutNode, bclist bcs, treelist trees);
};

class BCLayout
{
public:
    BCLayout(Canvas *canvas);
    Graph ogdfGraph(shapemap& nodeShapes, connmap& edgeConns);
    static void extractSizes(shapemap nodeShapes, GraphAttributes& GA);
    static void injectPositions(shapemap nodeShapes, GraphAttributes GA);
    static void injectSizes(shapemap nodeShapes, GraphAttributes GA);
    void applyKM3(void);
    void layoutBCTrees(void);
    QList<BiComp*> getNontrivialBCs(Graph G);
    QMap<int,node> getConnComps(Graph G);
    void orthoLayout(void);

private:
    Canvas *m_canvas;
};

}