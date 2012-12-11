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

class RootedTree
{
    //...
};

class BiComp
{
public:
    BiComp(QList<edge> edges, QList<node> nodes, QList<node> cutNodes);
    void removeSelf(Graph& G);
private:
    Graph *m_graph;
    QMap<node,node> m_nodemap; // maps own nodes to orig. graph nodes
    QList<node> m_cutNodes; // subdomain of nodemap which are cutnodes
    QList<node> m_normalNodes; // subdomain of those which are not
    QMap<edge,edge> m_edgemap; // maps own edges to orig. graph edges
};

class BCLayout
{
public:
    BCLayout(Canvas *canvas);
    Graph ogdfGraph(shapemap& nodeShapes, connmap& edgeConns);
    void extractSizes(shapemap nodeShapes, GraphAttributes& GA);
    void injectPositions(shapemap nodeShapes, GraphAttributes GA);
    void applyKM3(void);
    void layoutBCTrees(void);
    QList<BiComp*> getNontrivialBCs(Graph G);
    QMap<int,node> getConnComps(Graph G);
    void orthoLayout(void);

private:
    Canvas *m_canvas;
};

}
