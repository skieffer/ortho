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

#include <cmath>

#include <QMap>
#include <QPair>
#include <QPointF>
#include <QSizeF>

#include "libdunnartcanvas/BCLayout.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/connector.h"
#include "libdunnartcanvas/pluginshapefactory.h"

#include "libogdf/ogdf/basic/EdgeArray.h"
#include "libogdf/ogdf/basic/Graph_d.h"
#include "libogdf/ogdf/basic/GraphAttributes.h"
#include "libogdf/ogdf/basic/SList.h"
#include "libogdf/ogdf/basic/simple_graph_alg.h"
#include "libogdf/ogdf/decomposition/BCTree.h"
#include "libogdf/ogdf/energybased/FMMMLayout.h"

using namespace ogdf;

namespace dunnart {

BCLayout::BCLayout(Canvas *canvas) :
    m_canvas(canvas)
{}

Graph BCLayout::ogdfGraph(shapemap &nodeShapes, connmap &edgeConns)
{
    Graph G;
    foreach (CanvasItem *item, m_canvas->items())
    {
        if (ShapeObj *shape = isShapeForLayout(item))
        {
            node n = G.newNode();
            nodeShapes.insert(n,shape);
        }
    }
    foreach (CanvasItem *item, m_canvas->items())
    {
        if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            QPair<ShapeObj*,ShapeObj*> endpts = conn->getAttachedShapes();
            node u = nodeShapes.key(endpts.first);
            node v = nodeShapes.key(endpts.second);
            edge e = G.newEdge(u,v);
            edgeConns.insert(e,conn);
        }
    }
    return G;
}

void BCLayout::extractSizes(shapemap nodeShapes, GraphAttributes &GA)
{
    foreach (node v, nodeShapes.keys())
    {
        ShapeObj *shape = nodeShapes.value(v);
        GA.width(v) = shape->size().width();
        GA.height(v) = shape->size().height();
    }
}

void BCLayout::injectPositions(shapemap nodeShapes, GraphAttributes GA)
{
    foreach (node v, nodeShapes.keys())
    {
        double cx = GA.x(v) + GA.width(v)/2.0;
        double cy = GA.y(v) + GA.height(v)/2.0;
        ShapeObj *shape = nodeShapes.value(v);
        shape->setCentrePos(QPointF(cx,cy));
    }
}

void BCLayout::applyKM3()
{
    shapemap nodeShapes;
    connmap edgeConns;
    Graph G = ogdfGraph(nodeShapes, edgeConns);
    GraphAttributes GA(G);
    extractSizes(nodeShapes, GA);
    FMMMLayout fm3;
    fm3.call(GA);
    injectPositions(nodeShapes, GA);
}

void BCLayout::orthoLayout()
{
    shapemap nodeShapes;
    connmap edgeConns;
    Graph G = ogdfGraph(nodeShapes, edgeConns);
    EdgeArray<int> bcs;
    biconnectedComponents(G, bcs);
    //...
}

#if 0
void BCLayout::applyKM3()
{
    // Construct GraphAttributes, and maps from ogdf nodes and edges
    // to Dunnart shapes and connectors.
    ogdf::Graph G;
    QMap<ogdf::node,ShapeObj*> nodeShapes;
    QMap<ogdf::edge,Connector*> edgeConns;
    foreach (CanvasItem *item, m_canvas->items())
    {
        if (ShapeObj *shape = isShapeForLayout(item))
        {
            ogdf::node n = G.newNode();
            nodeShapes.insert(n,shape);
        }
    }
    foreach (CanvasItem *item, m_canvas->items())
    {
        if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            QPair<ShapeObj*,ShapeObj*> endpts = conn->getAttachedShapes();
            ogdf::node u = nodeShapes.key(endpts.first);
            ogdf::node v = nodeShapes.key(endpts.second);
            ogdf::edge e = G.newEdge(u,v);
            edgeConns.insert(e,conn);
        }
    }
    ogdf::GraphAttributes GA(G);
    foreach (ogdf::node v, nodeShapes.keys())
    {
        ShapeObj *shape = nodeShapes.value(v);
        GA.width(v) = shape->size().width();
        GA.height(v) = shape->size().height();
    }

    // Run the FMMM algorithm.
    ogdf::FMMMLayout fm3;
    fm3.call(GA);

    // Apply the results to the canvas.
    foreach (ogdf::node v, nodeShapes.keys())
    {
        double cx = GA.x(v) + GA.width(v)/2.0;
        double cy = GA.y(v) + GA.height(v)/2.0;
        ShapeObj *shape = nodeShapes.value(v);
        shape->setCentrePos(QPointF(cx,cy));
    }
}
#endif

void BCLayout::layoutBCTrees()
{
    //Construct graph, and label nodes by ID.
    ogdf::Graph G;
    QMap<ogdf::node,ShapeObj*> nodeShapes;
    QMap<ogdf::edge,Connector*> edgeConns;
    foreach (CanvasItem *item, m_canvas->items())
    {
        if (ShapeObj *shape = isShapeForLayout(item))
        {
            // Construct OGDF node
            ogdf::node n = G.newNode();
            nodeShapes.insert(n,shape);
            // Label by ID.
            shape->setLabel(QString::number(shape->internalId()));
        }
    }
    foreach (CanvasItem *item, m_canvas->items())
    {
        if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            QPair<ShapeObj*,ShapeObj*> endpts = conn->getAttachedShapes();
            ogdf::node u = nodeShapes.key(endpts.first);
            ogdf::node v = nodeShapes.key(endpts.second);
            ogdf::edge e = G.newEdge(u,v);
            edgeConns.insert(e,conn);
        }
    }
    //Construct BC-tree.
    ogdf::BCTree bctree(G);
    // Draw it.
    ogdf::Graph B = bctree.bcTree();
    ogdf::node v;
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    ShapeObj *sh = NULL;
    QMap<ogdf::node,ShapeObj*> bcShapes;
    m_canvas->stop_graph_layout();
    // nodes
    forall_nodes(v,B)
    {
        sh = factory->createShape("org.dunnart.shapes.ellipse");
        ogdf::BCTree::BNodeType bnt = bctree.typeOfBNode(v);
        int n = bctree.numberOfNodes(v);
        double s = 25 + (150/3.14)*std::atan(double(n-1)*0.173);
        switch (bnt)
        {
        case ogdf::BCTree::BComp:
            sh->setPosAndSize(QPointF(0,0), QSizeF(s,s));
            sh->setFillColour(QColor(128,128,255));
            break;
        case ogdf::BCTree::CComp:
            sh->setPosAndSize(QPointF(0,0), QSizeF(25,25));
            sh->setFillColour(QColor(255,128,128));
            break;
        }
        m_canvas->addItem(sh);
        bcShapes.insert(v,sh);
    }
    ogdf::edge e;
    // edges
    forall_edges(e,B)
    {
        ogdf::node src = e->source();
        ogdf::node dst = e->target();
        ShapeObj *srcSh = bcShapes.value(src);
        ShapeObj *dstSh = bcShapes.value(dst);
        Connector *conn = new Connector();
        conn->initWithConnection(srcSh, dstSh);
        m_canvas->addItem(conn);
    }

    // Draw the "auxiliary graph"
    /*
    ogdf::Graph H = bctree.auxiliaryGraph();
    QMap<ogdf::node,ShapeObj*> auxShapes;
    // nodes
    forall_nodes(v,H)
    {
        sh = factory->createShape("org.dunnart.shapes.ellipse");
        sh->setPosAndSize(QPointF(0,0), QSizeF(25,25));
        sh->setFillColour(QColor(128,255,128));
        addItem(sh);
        auxShapes.insert(v,sh);
    }
    // edges
    forall_edges(e,H)
    {
        ogdf::node src = e->source();
        ogdf::node dst = e->target();
        ShapeObj *srcSh = auxShapes.value(src);
        ShapeObj *dstSh = auxShapes.value(dst);
        Connector *conn = new Connector();
        conn->initWithConnection(srcSh, dstSh);
        addItem(conn);
    }
    */

    // Draw the "auxiliary graph" -- draw only nodes and edges belonging to
    // biconnected components of 3 or more nodes
    ogdf::Graph H = bctree.auxiliaryGraph();

    QSet<ogdf::node> hNodes;
    QSet<ogdf::edge> hEdges;
    forall_nodes(v,B)
    {
        int n = bctree.numberOfNodes(v);
        if (n < 3) continue;
        ogdf::SList<ogdf::edge> edges = bctree.hEdges(v);
        for (ogdf::SListIterator<ogdf::edge> i = edges.begin(); i!=edges.end(); i++)
        {
            hNodes.insert( (*i)->source() );
            hNodes.insert( (*i)->target() );
            hEdges.insert( *i );
        }
    }
    QMap<ogdf::node,ShapeObj*> auxShapes;
    // nodes
    foreach (ogdf::node v, hNodes)
    {
        sh = factory->createShape("org.dunnart.shapes.ellipse");
        sh->setPosAndSize(QPointF(0,0), QSizeF(25,25));
        sh->setFillColour(QColor(128,255,128));
        m_canvas->addItem(sh);
        auxShapes.insert(v,sh);
    }
    // edges
    foreach (ogdf::edge e, hEdges)
    {
        ogdf::node src = e->source();
        ogdf::node dst = e->target();
        ShapeObj *srcSh = auxShapes.value(src);
        ShapeObj *dstSh = auxShapes.value(dst);
        Connector *conn = new Connector();
        conn->initWithConnection(srcSh, dstSh);
        m_canvas->addItem(conn);
    }


    // labels
    /*
    forall_nodes(v,G)
    {
        if (bctree.typeOfGNode(v)==ogdf::BCTree::CutVertex)
        {
            ShapeObj *origSh = nodeShapes.value(v);
            QString label = origSh->getLabel();
            //ogdf::node bcNode =
        }
    }
    */
    m_canvas->interrupt_graph_layout();
}


}
