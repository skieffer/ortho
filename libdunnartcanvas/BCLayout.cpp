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
#include <assert.h>

#include <QtGui>
#include <QList>
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
#include "libogdf/ogdf/basic/NodeArray.h"
#include "libogdf/ogdf/basic/Graph_d.h"
#include "libogdf/ogdf/basic/GraphAttributes.h"
#include "libogdf/ogdf/basic/SList.h"
#include "libogdf/ogdf/basic/simple_graph_alg.h"
#include "libogdf/ogdf/decomposition/BCTree.h"
#include "libogdf/ogdf/energybased/FMMMLayout.h"

using namespace ogdf;

namespace dunnart {

// ----------------------------------------------------------------------------
// RootedTree -----------------------------------------------------------------

RootedTree::RootedTree(QList<node> nodes, Graph G) :
    m_graph(NULL),
    m_root(NULL)
{
    // Construct own copy of graph, maintaining maps to original graph.
    m_graph = new Graph;
    QSet<edge> edges;
    foreach (node n, nodes)
    {
        // Create own node rep.
        node m = m_graph->newNode();
        m_nodemap.insert(m,n);
        // Note those edges that are internal to this node list.
        edge e = NULL;
        forall_adj_edges(e,n)
        {
            node q = NULL;
            if (n==e->source()) { q = e->target(); }
            else { q = e->source(); }
            assert(q);
            if (nodes.contains(q)) { edges.insert(e); }
        }
    }
    // Now process internal edges.
    foreach (edge e, edges)
    {
        node m1 = m_nodemap.key(e->source());
        node m2 = m_nodemap.key(e->target());
        edge f = m_graph->newEdge(m1,m2);
        m_edgemap.insert(f,e);
    }
}

void RootedTree::setRelPt(QPointF p)
{
    m_relpt = p;
}

bool RootedTree::containsOriginalNode(node n)
{
    return m_nodemap.values().contains(n);
}

void RootedTree::recursiveLayout(shapemap origShapes, bclist bcs, treelist trees,
                                 node origBaseNode, QPointF cardinal)
{
    constructDunnartGraph(origShapes);
    // TODO: Use origBaseNode and cardinal.
    // If origBaseNode is NULL, then it is to be ignored.
    // Otherwise, it is a node from the original graph, and we are meant to
    // layout this subgraph starting from its local copy of that node, and
    // moving in the direction given by cardinal.
}

void RootedTree::constructDunnartGraph(shapemap origShapes)
{
    // Store passed value.
    m_origShapeMap = origShapes;

    // Create Dunnart shape and connector objects for own graph.
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    foreach (node m, m_nodemap.keys())
    {
        ShapeObj *shape = factory->createShape("org.dunnart.shapes.ellipse");
        shape->setFillColour(QColor(192,0,0));
        m_ownShapeMap.insert(m,shape);
    }
    foreach (edge f, m_edgemap.keys())
    {
        node src = f->source();
        node dst = f->target();
        ShapeObj *srcSh = m_ownShapeMap.value(src);
        ShapeObj *dstSh = m_ownShapeMap.value(dst);
        Connector *conn = new Connector();
        conn->initWithConnection(srcSh, dstSh);
        conn->setRoutingType(Connector::orthogonal);
        m_ownConnMap.insert(f,conn);
    }

    // Compute initial layout by FMMM.
    GraphAttributes GA(*m_graph);
    // Use sizes from original graph.
    BCLayout::extractSizes(m_origShapeMap, GA);
    FMMMLayout fm3;
    fm3.call(GA);
    BCLayout::injectPositions(m_ownShapeMap, GA);
    BCLayout::injectSizes(m_ownShapeMap, GA);
}

// ----------------------------------------------------------------------------
// BiComp ---------------------------------------------------------------------

#ifdef GRAPHREF
BiComp::BiComp(QList<edge> edges, QList<node> nodes, QList<node> cutNodes)
{
    // Construct own copy of graph, maintaining maps to original graph.
    foreach (node n, nodes)
    {
        node m = m_graph.newNode();
        m_nodemap.insert(m,n);
        if (cutNodes.contains(n))
        {
            m_cutNodes.append(m);
        } else
        {
            m_normalNodes.append(m);
        }
    }
    foreach (edge e, edges)
    {
        node m1 = m_nodemap.key(e->source());
        node m2 = m_nodemap.key(e->target());
        edge f = m_graph.newEdge(m1,m2);
        m_edgemap.insert(f,e);
    }
}
#else
BiComp::BiComp(QList<edge> edges, QList<node> nodes, QList<node> cutNodes) :
    m_graph(NULL)
{
    // Construct own copy of graph, maintaining maps to original graph.
    m_graph = new Graph;
    foreach (node n, nodes)
    {
        node m = m_graph->newNode();
        m_nodemap.insert(m,n);
        if (cutNodes.contains(n))
        {
            m_cutNodes.append(m);
        } else
        {
            m_normalNodes.append(m);
        }
    }
    foreach (edge e, edges)
    {
        node m1 = m_nodemap.key(e->source());
        node m2 = m_nodemap.key(e->target());
        edge f = m_graph->newEdge(m1,m2);
        m_edgemap.insert(f,e);
    }
}
#endif

void BiComp::setRelPt(QPointF p)
{
    m_relpt = p;
}

void BiComp::removeSelf(Graph &G)
{
    static int i;
    // Remove all edges.
    foreach (edge e, m_edgemap.values())
    {
        i++;
        G.delEdge(e);
        assert(G.consistencyCheck());
        //G.hideEdge(e);
    }
    // Remove all nodes which are /not/ cut nodes....
    // This was causing an error in ogdf, for no clear reason.
    // But we actually don't need to do it. After removing the
    // edges, the nodes in m_normalNodes should all be isolated
    // nodes in G. So they will be disregarded when we compute
    // the nontrivial trees among the remaining connected components.
    /*
    foreach (node m, m_normalNodes)
    {
        node n = m_nodemap.value(m);
        G.delNode(n);
    }
    */
}

size_t BiComp::size()
{
    return m_nodemap.size();
}


// pass shapemap for the shapes in the original graph
void BiComp::constructDunnartGraph(shapemap origShapes)
{
    // Store passed value.
    m_origShapeMap = origShapes;

    // Create Dunnart shape and connector objects for own graph.
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    foreach (node m, m_nodemap.keys())
    {
        ShapeObj *shape = factory->createShape("org.dunnart.shapes.ellipse");
        shape->setFillColour(QColor(0,0,192));
        m_ownShapeMap.insert(m,shape);
    }
    foreach (edge f, m_edgemap.keys())
    {
        node src = f->source();
        node dst = f->target();
        ShapeObj *srcSh = m_ownShapeMap.value(src);
        ShapeObj *dstSh = m_ownShapeMap.value(dst);
        Connector *conn = new Connector();
        conn->initWithConnection(srcSh, dstSh);
        conn->setRoutingType(Connector::orthogonal);
        m_ownConnMap.insert(f,conn);
    }

    // Compute initial layout by FMMM.
    /*
#ifdef GRAPHREF
    GraphAttributes GA(m_graph);
#else
    GraphAttributes GA(*m_graph);
#endif
    */
    QMap<node,node> nodemapNewToOld;
    Graph G = copyGraph(nodemapNewToOld);
    GraphAttributes GA(G);
    // Use sizes from original graph.
    BCLayout::extractSizes(m_origShapeMap, GA);
    FMMMLayout fm3;
    fm3.call(GA);
    //BCLayout::injectPositions(m_ownShapeMap, GA);
    //BCLayout::injectSizes(m_ownShapeMap, GA);
    BCLayout::injectPositionsAndSizes(nodemapNewToOld, m_ownShapeMap, GA);
}

Graph& BiComp::copyGraph(QMap<node, node>& nodemap)
{
    Graph G;
    //QMap<node,node> nodemap;
    node v = NULL;
    forall_nodes(v,*m_graph)
    {
        node u = G.newNode();
        nodemap.insert(u,v);
    }
    edge e = NULL;
    forall_edges(e,*m_graph)
    {
        node s = e->source();
        node t = e->target();
        node m = nodemap.key(s);
        node n = nodemap.key(t);
        G.newEdge(m,n);
    }
    return G;
}

void BiComp::improveOrthogonalTopology()
{
    // MUST call constructDunnartGraph first!
    // Add graph to a fresh Dunnart canvas.
    Canvas canvas;
    foreach (ShapeObj *sh, m_ownShapeMap.values())
    {
        canvas.addItem(sh);
    }
    foreach (Connector *conn, m_ownConnMap.values())
    {
        canvas.addItem(conn);
    }
    // Turn on automatic layout.
    canvas.setOptAutomaticGraphLayout(true);
    // Improve orthogonal topology.
    canvas.improveOrthogonalTopology();
}

QPointF BiComp::baryCentre()
{
    // MUST call constructDunnartGraph first!
    QPointF b(0,0);
    int n = 0;
    foreach (ShapeObj *sh, m_ownShapeMap.values())
    {
        b += sh->centrePos();
        n++;
    }
    //return QPointF(b.x()/float(n), b.y()/float(n));
    return b*(1/float(n));
}

QPointF BiComp::nearestCardinal(QPointF v)
{
    double x = v.x(); double y = v.y();
    QPointF u = (abs(y) > abs(x) ?
                     ( y<0?QPointF(0,-1):QPointF(0,1) ) :
                     ( x<0?QPointF(-1,0):QPointF(1,0) ) );
    return u;
}

bool BiComp::containsOriginalNode(node n)
{
    return m_nodemap.values().contains(n);
}

bool BiComp::containsOriginalEdge(edge e)
{
    return m_edgemap.values().contains(e);
}

/**
 * Find all chunks that contain (an image of) the cut node from the original graph,
 * and which are different from the present chunk.
 */
QList<Chunk*> BiComp::findCutNodeNeighbours(node origCutNode, bclist bcs, treelist trees)
{
    QList<Chunk*> nbrs;
    foreach (BiComp *bc, bcs)
    {
        if (bc==this) continue;
        if (bc->containsOriginalNode(origCutNode)) nbrs.append(bc);
    }
    foreach (RootedTree *rt, trees)
    {
        if (rt->containsOriginalNode(origCutNode)) nbrs.append(rt);
    }
    return nbrs;
}

void BiComp::recursiveLayout(shapemap origShapes, bclist bcs, treelist trees,
                             node origBaseNode, QPointF cardinal)
{
    // TODO: Use origBaseNode and cardinal.
    // If origBaseNode is NULL, then it is to be ignored.
    // Otherwise, it is a node from the original graph, and we are meant to
    // layout this subgraph starting from its local copy of that node, and
    // moving in the direction given by cardinal.
    constructDunnartGraph(origShapes);
    QPointF bary = baryCentre();
    foreach (node cn, m_cutNodes)
    {
        // Compute cardinal direction from barycentre to cut node.
        ShapeObj *cnShape = m_ownShapeMap.value(cn);
        QPointF cnCentre = cnShape->centrePos();
        QPointF cardinal = nearestCardinal(cnCentre - bary);
        // Find neighbours connected through original cut node.
        node origCutNode = m_nodemap.value(cn);
        QList<Chunk*> nbrs = findCutNodeNeighbours(origCutNode, bcs, trees);
        // Add as children
        m_children.append(nbrs);
        // Layout each one, and set its relative point.
        int fixedSep = 50; // magic number
        foreach (Chunk *ch, nbrs)
        {
            QPointF p = cnCentre + fixedSep*cardinal;
            ch->setRelPt(p);
            ch->recursiveLayout(origShapes, bcs, trees, cn, cardinal);
        }
    }
}

// ----------------------------------------------------------------------------
// BCLayout -------------------------------------------------------------------

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

void BCLayout::injectPositions(shapemap nodeShapes, GraphAttributes& GA)
{
    foreach (node v, nodeShapes.keys())
    {
        double cx = GA.x(v) + GA.width(v)/2.0;
        double cy = GA.y(v) + GA.height(v)/2.0;
        ShapeObj *shape = nodeShapes.value(v);
        shape->setCentrePos(QPointF(cx,cy));
    }
}

void BCLayout::injectSizes(shapemap nodeShapes, GraphAttributes& GA)
{
    foreach (node v, nodeShapes.keys())
    {
        double w = GA.width(v);
        double h = GA.height(v);
        ShapeObj *shape = nodeShapes.value(v);
        shape->setSize(QSizeF(w,h));
    }
}

void BCLayout::injectPositionsAndSizes(
        QMap<node, node>& nodemap, shapemap nodeShapes, GraphAttributes &GA)
{
    foreach (node v, nodemap.keys())
    {
        double cx = GA.x(v) + GA.width(v)/2.0;
        double cy = GA.y(v) + GA.height(v)/2.0;
        double w = GA.width(v);
        double h = GA.height(v);
        node u = nodemap.value(v);
        ShapeObj *shape = nodeShapes.value(u);
        shape->setCentrePos(QPointF(cx,cy));
        shape->setSize(QSizeF(w,h));
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

QList<BiComp*> BCLayout::getNontrivialBCs(Graph G)
{
    QList<BiComp*> bicomps;
    BCTree bctree(G);
    Graph B = bctree.bcTree();
    node vB;
    forall_nodes(vB, B)
    {
        // Skip cutnodes.
        if (bctree.typeOfBNode(vB)!=BCTree::BComp) continue;

        // Skip the component if it is too small.
        int n = bctree.numberOfNodes(vB);
        if (n < 3) continue;

        // Otherwise get the edges of the original graph belonging to this component,
        // as well as the set of vertices, and subset of cut vertices.
        SList<edge> hEdges = bctree.hEdges(vB);
        QList<edge> gEdges;
        QSet<node> gNodes;
        QSet<node> gCutNodes;
        for (SListIterator<edge> i = hEdges.begin(); i!=hEdges.end(); i++)
        {
            edge orig = bctree.original(*i);
            gEdges.append(orig);
            node src = orig->source();
            gNodes.insert(src);
            if (bctree.typeOfGNode(src)==BCTree::CutVertex) gCutNodes.insert(src);
            node dst = orig->target();
            gNodes.insert(dst);
            if (bctree.typeOfGNode(dst)==BCTree::CutVertex) gCutNodes.insert(dst);
        }
        BiComp *bc = new BiComp(gEdges, gNodes.toList(), gCutNodes.toList());
        bicomps.append(bc);
    }
    return bicomps;
}

QMap<int,node> BCLayout::getConnComps(Graph G)
{
    NodeArray<int> ccomps(G);
    connectedComponents(G, ccomps);
    QMap<int,node> map;
    node v = NULL;
    forall_nodes(v,G)
    {
        int n = ccomps[v];
        map.insertMulti(n,v);
    }
    return map;
}

QMap<int,node> BCLayout::getConnComps2(Graph& G2, QMap<node, node> nodeMap)
{
    // G2 the second graph we construct; nodeMap the map that maps nodes
    // of G2 back to G (our first graph).

    NodeArray<int> ccomps(G2);
    connectedComponents(G2, ccomps);
    QMap<int,node> map;
    node v2 = NULL;
    forall_nodes(v2,G2)
    {
        int n = ccomps[v2];
        node v = nodeMap.value(v2);
        map.insertMulti(n,v);
    }
    return map;
}

Graph BCLayout::removeBiComps(Graph G, bclist bcs, QMap<node, node>& nodeMap)
{
    // We construct and return a new graph which is equivalent of G
    // with all edges in the BiComps in bcs deleted.
    // The node map will be filled in with a mapping from nodes of the new graph
    // back to the corresponding nodes in the old graph.

    Graph Gp;

    // Copy nodes.
    node v = NULL;
    forall_nodes(v,G)
    {
        node u = Gp.newNode();
        nodeMap.insert(u,v);
    }
    // Copy edges, except those in the BiComps.
    edge e = NULL;
    forall_edges(e,G)
    {
        bool keep = true;
        foreach (BiComp *bc, bcs)
        {
            if (bc->containsOriginalEdge(e))
            {
                keep = false;
                break;
            }
        }
        if (keep)
        {
            node s = e->source();
            node t = e->target();
            node m = nodeMap.key(s);
            node n = nodeMap.key(t);
            Gp.newEdge(m,n);
        }
    }

    return Gp;
}

void BCLayout::orthoLayout()
{
    shapemap nodeShapes;
    connmap edgeConns;
    Graph G = ogdfGraph(nodeShapes, edgeConns);

    // Get nontrivial biconnected components (size >= 3).
    QList<BiComp*> bicomps = getNontrivialBCs(G);

    // Make a copy of the graph.
    //Graph G2(G);

    // Remove them from the original graph, and get the remaining
    // connected components.
    /*
    foreach (BiComp *bc, bicomps)
    {
        bc->removeSelf(G);
    }
    assert(G.consistencyCheck());
    */

    QMap<node,node> nodeMapG2;
    Graph G2 = removeBiComps(G, bicomps, nodeMapG2);
    assert(G2.consistencyCheck());


    //QMap<int,node> ccomps = getConnComps(G);
    QMap<int,node> ccomps = getConnComps2(G2, nodeMapG2);


    // Form trees on remaining components, throwing away isolated
    // nodes (which must have been cutnodes shared only by nontrivial BCs).
    QList<RootedTree*> rtrees;
    foreach (int i, ccomps.keys().toSet())
    {
        QList<node> nodes = ccomps.values(i);
        if (nodes.size() < 2) continue;
        RootedTree *rtree = new RootedTree(nodes,G);
        rtrees.append(rtree);
    }

    // Now layout each component, and connect them together.
    // Begin with a largest biconnected component.
    size_t n = 0;
    BiComp *largest = NULL;
    foreach (BiComp *bc, bicomps)
    {
        size_t m = bc->size();
        if (m > n)
        {
            n = m;
            largest = bc;
        }
    }
    assert(largest);
    largest->recursiveLayout(nodeShapes, bicomps, rtrees, NULL, QPointF(0,0));
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
