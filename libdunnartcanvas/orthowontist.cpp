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
#include <vector>

#include <QtGui>
#include <QList>
#include <QMap>
#include <QPair>
#include <QPointF>
#include <QSizeF>
#include <QTime>

#include "libvpsc/rectangle.h"

#include "libcola/cola.h"
#include "libcola/compound_constraints.h"
#include "libcola/cc_nonoverlapconstraints.h"

#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/connector.h"
#include "libdunnartcanvas/pluginshapefactory.h"
#include "libdunnartcanvas/guideline.h"
#include "libdunnartcanvas/separation.h"
#include "libdunnartcanvas/distribution.h"

#include "libogdf/ogdf/basic/EdgeArray.h"
#include "libogdf/ogdf/basic/NodeArray.h"
#include "libogdf/ogdf/basic/geometry.h"
#include "libogdf/ogdf/basic/Graph_d.h"
#include "libogdf/ogdf/basic/GraphAttributes.h"
#include "libogdf/ogdf/basic/SList.h"
#include "libogdf/ogdf/basic/simple_graph_alg.h"

#include "libogdf/ogdf/decomposition/BCTree.h"

#include "libogdf/ogdf/energybased/FMMMLayout.h"
#include "libogdf/ogdf/energybased/SpringEmbedderFR.h"
#include "libogdf/ogdf/energybased/StressMajorizationSimple.h"

#include "libogdf/ogdf/planarlayout/PlanarDrawLayout.h"
#include "libogdf/ogdf/planarlayout/PlanarStraightLayout.h"

#include "libogdf/ogdf/planarity/EmbedderMaxFace.h"
#include "libogdf/ogdf/planarity/EmbedderMinDepth.h"
#include "libogdf/ogdf/planarity/EmbedderMinDepthMaxFace.h"
#include "libogdf/ogdf/planarity/EmbedderMaxFaceLayers.h"
#include "libogdf/ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h"
#include "libogdf/ogdf/planarity/EmbedderMinDepthPiTa.h"
#include "libogdf/ogdf/planarity/PlanarizationLayout.h"
#include "libogdf/ogdf/planarity/PlanarizationGridLayout.h"

#include "libogdf/ogdf/tree/TreeLayout.h"

#include "libdunnartcanvas/orthowontist.h"

using namespace ogdf;

namespace ow {

// ------------------------------------------------------------------
// BiComp -----------------------------------------------------------

BiComp::BiComp(void) {}

BiComp::BiComp(QList<node> nodes, QList<edge> edges, QList<node> cutnodes,
               shapemap nodeShapes, connmap edgeConns) {
    m_graph = new Graph();
    QMap<node,node> nodemap; // maps passed nodes to own nodes; for use in creating edges
    m_ga = new GraphAttributes(*m_graph);
    foreach (node n, nodes) {
        node m = m_graph->newNode();
        nodemap.insert(n,m);
        if (cutnodes.contains(n)) m_cutnodes.append(m);
        ShapeObj *sh = nodeShapes.value(n);
        m_dunnartShapes.insert(m,sh);
        // Initialize size and position from Dunnart shape.
        QSizeF size = sh->size();
        QPointF pos = sh->centrePos();
        m_ga->width(m) = size.width();
        m_ga->height(m) = size.height();
        m_ga->x(m) = pos.x();
        m_ga->y(m) = pos.y();
    }
    foreach (edge e, edges) {
        node m1 = nodemap.value(e->source());
        node m2 = nodemap.value(e->target());
        edge f = m_graph->newEdge(m1,m2);
        m_dunnartConns.insert(f, edgeConns.value(e));
    }
}

QList<node> BiComp::cutnodes(void) {
    return m_cutnodes;
}

QString BiComp::listNodes(void) {
    QString s = "";
    s += "BC nodes:\n  ";
    foreach (node m, m_dunnartShapes.keys()) {
        int id = m_dunnartShapes.value(m)->internalId();
        s += QString("%1, ").arg(id);
    }
    return s;
}

void BiComp::colourShapes(void) {
    foreach (node m, m_dunnartShapes.keys()) {
        ShapeObj *sh = m_dunnartShapes.value(m);
        QColor col = m_cutnodes.contains(m) ? QColor(0,192,0) : QColor(0,192,255);
        sh->setFillColour(col);
    }
}

ShapeObj *BiComp::getShape(node m) {
    return m_dunnartShapes.value(m);
}

/* Perform depth-first search through the network whose nodes are the
 * BiComps and whose hyperedges are the cutnodes.
 *
 * endpts maps cutnodes to the BCs that are their "end points";
 * elements is the list of BCs seen so far in the DFS.
 *
 * We simply add BCs to the elements list.
 * This is used in decomposing the network of BCs into its connected
 * components.
 */
void BiComp::dfs(QMap<ShapeObj*, BiComp*> endpts, QList<BiComp *> &elements)
{
    elements.append(this);
    foreach (node m, m_cutnodes) {
        ShapeObj *sh = m_dunnartShapes.value(m);
        QList<BiComp*> next = endpts.values(sh);
        foreach (BiComp *bc, next) {
            if (elements.contains(bc)) continue;
            bc->dfs(endpts,elements);
        }
    }
}

/* Create a new BiComp representing the fusion of this one and the
 * one passed, at all of their shared cutnodes.
 * Neither of the constituent BCs is altered.
 * If they do not share any cutnode, then the BC returned will simply
 * be a copy of this one.
 * Note: a BiComp object has a list of cutnodes. After fusion, the new
 * list will still include any cutnodes at which the constituent BCs were
 * fused, even if these do not connect to any other BCs or trees.
 */
BiComp *BiComp::fuse(BiComp *other)
{
    BiComp *fusion = new BiComp();
    // Give the fusion BiComp its Graph object.
    fusion->m_graph = new Graph;
    fusion->m_ga = new GraphAttributes(*(fusion->m_graph));
    // Below we'll use the following letters for nodes:
    // t: node in 't'his BiComp 'T'
    // o: node in 'o'ther BiComp 'O'
    // f: node in 'f'usion BiComp 'F'
    // g: node in original 'g'raph 'G'
    // Begin by giving F a copy of each non-cutnode in T.
    foreach (node t, this->m_dunnartShapes.keys()) {
        if (this->m_cutnodes.contains(t)) continue;
        node f = fusion->m_graph->newNode();
        fusion->m_ga->width(f) = this->m_ga->width(t);
        fusion->m_ga->height(f) = this->m_ga->height(t);
        fusion->m_ga->x(f) = this->m_ga->x(t);
        fusion->m_ga->y(f) = this->m_ga->y(t);
        ShapeObj *sh = this->m_dunnartShapes.value(t);
        fusion->m_dunnartShapes.insert(f,sh);
    }
    // Now give F a copy of each non-cutnode in O.
    foreach (node o, other->m_dunnartShapes.keys()) {
        if (other->m_cutnodes.contains(o)) continue;
        node f = fusion->m_graph->newNode();
        fusion->m_ga->width(f) = other->m_ga->width(o);
        fusion->m_ga->height(f) = other->m_ga->height(o);
        fusion->m_ga->x(f) = other->m_ga->x(o);
        fusion->m_ga->y(f) = other->m_ga->y(o);
        ShapeObj *sh = other->m_dunnartShapes.value(o);
        fusion->m_dunnartShapes.insert(f,sh);
    }
    // Now the cutnodes of T.
    foreach (node t, this->m_cutnodes) {
        node f = fusion->m_graph->newNode();
        fusion->m_cutnodes.append(f);
        fusion->m_ga->width(f) = this->m_ga->width(t);
        fusion->m_ga->height(f) = this->m_ga->height(t);
        fusion->m_ga->x(f) = this->m_ga->x(t);
        fusion->m_ga->y(f) = this->m_ga->y(t);
        ShapeObj *sh = this->m_dunnartShapes.value(t);
        fusion->m_dunnartShapes.insert(f,sh);
    }
    // And the cutnodes of O. This time check for duplicates
    // with T.
    foreach (node o, other->m_cutnodes) {
        ShapeObj *sh = other->m_dunnartShapes.value(o);
        if (this->m_dunnartShapes.values().contains(sh)) continue;
        node f = fusion->m_graph->newNode();
        fusion->m_cutnodes.append(f);
        fusion->m_ga->width(f) = other->m_ga->width(o);
        fusion->m_ga->height(f) = other->m_ga->height(o);
        fusion->m_ga->x(f) = other->m_ga->x(o);
        fusion->m_ga->y(f) = other->m_ga->y(o);
        fusion->m_dunnartShapes.insert(f,sh);
    }
    // Now the edges.
    // First the edges of T.
    foreach (edge e, this->m_dunnartConns.keys()) {
        Connector *c = this->m_dunnartConns.value(e);
        node e1 = e->source();
        node e2 = e->target();
        ShapeObj *sh1 = this->m_dunnartShapes.value(e1);
        ShapeObj *sh2 = this->m_dunnartShapes.value(e2);
        node f1 = fusion->m_dunnartShapes.key(sh1);
        node f2 = fusion->m_dunnartShapes.key(sh2);
        edge f = fusion->m_graph->newEdge(f1,f2);
        fusion->m_dunnartConns.insert(f,c);
    }
    // And finally the edges of O.
    foreach (edge e, other->m_dunnartConns.keys()) {
        Connector *c = other->m_dunnartConns.value(e);
        node e1 = e->source();
        node e2 = e->target();
        ShapeObj *sh1 = other->m_dunnartShapes.value(e1);
        ShapeObj *sh2 = other->m_dunnartShapes.value(e2);
        node f1 = fusion->m_dunnartShapes.key(sh1);
        node f2 = fusion->m_dunnartShapes.key(sh2);
        edge f = fusion->m_graph->newEdge(f1,f2);
        fusion->m_dunnartConns.insert(f,c);
    }
    return fusion;
}

// ------------------------------------------------------------------
// ExternalTree -----------------------------------------------------

ExternalTree::ExternalTree(node root, node rootInG, QList<node> nodes, QList<edge> edges,
                           shapemap nodeShapes, connmap edgeConns) :
    m_rootInG(rootInG)
{
    assert(nodes.size() > 1);
    m_graph = new Graph();
    QMap<node,node> nodemap; // maps passed nodes to own nodes; for use in creating edges
    m_ga = new GraphAttributes(*m_graph);
    foreach (node n, nodes) {
        node m = m_graph->newNode();
        if (n==root) m_root = m;
        nodemap.insert(n,m);
        ShapeObj *sh = nodeShapes.value(n);
        m_dunnartShapes.insert(m,sh);
        // Initialize size and position from Dunnart shape.
        QSizeF size = sh->size();
        QPointF pos = sh->centrePos();
        m_ga->width(m) = size.width();
        m_ga->height(m) = size.height();
        m_ga->x(m) = pos.x();
        m_ga->y(m) = pos.y();
    }
    foreach (edge e, edges) {
        node m1 = nodemap.value(e->source());
        node m2 = nodemap.value(e->target());
        edge f = m_graph->newEdge(m1,m2);
        m_dunnartConns.insert(f, edgeConns.value(e));
    }
}

QString ExternalTree::listNodes(void) {
    QString s = "";
    s += QString("Root: %1").arg(m_dunnartShapes.value(m_root)->internalId());
    s += QString("\n  Other nodes: ");
    foreach (node m, m_dunnartShapes.keys()) {
        if (m==m_root) continue;
        int id = m_dunnartShapes.value(m)->internalId();
        s += QString("%1, ").arg(id);
    }
    return s;
}

void ExternalTree::colourShapes(void) {
    foreach (node m, m_dunnartShapes.keys()) {
        ShapeObj *sh = m_dunnartShapes.value(m);
        QColor col = m==m_root ? QColor(255,192,0) : QColor(192,192,0);
        sh->setFillColour(col);
    }
}


// ------------------------------------------------------------------
// Orthowontist -----------------------------------------------------

Orthowontist::Orthowontist(Canvas *canvas) :
    m_canvas(canvas) {
}

void Orthowontist::run1(QList<CanvasItem*> items) {
    bool debug = true;
    shapemap nodeShapes;
    connmap edgeConns;
    Graph G;
    GraphAttributes GA(G);
    buildOGDFGraph(items,G,GA,nodeShapes,edgeConns);

    // 0. Compute average dimensions.
    double avgHeight=0, avgWidth=0, avgDim=0;
    node n;
    int N = 0;
    forall_nodes(n,G) {
        double w = GA.width(n), h = GA.height(n);
        avgWidth += w;
        avgHeight += h;
        N++;
        if (debug) {
            qDebug() << QString("node size: %1 x %2").arg(GA.width(n)).arg(GA.height(n));
        }
    }
    avgWidth /= N;
    avgHeight /= N;
    avgDim = (avgWidth+avgHeight)/2;

    // 1. Remove external trees from OGDF graph G.
    QList<ExternalTree*> EE;
    removeExternalTrees(EE,G,nodeShapes,edgeConns);
    // Debug:
    if (debug) {
        qDebug() << "External Trees:";
        foreach (ExternalTree *E, EE) {
            qDebug() << E->listNodes();
        }
    }

    // 2. Get nontrivial biconnected components (size >= 3), and get the
    // set of cutnodes in G.
    QList<BiComp*> BB;
    QSet<node> cutnodes;
    buildNBCs(BB, cutnodes, G, nodeShapes, edgeConns);
    // Fuse BCs that share cutnodes.
    BB = fuseBCs(BB);
    // Make a map to say to which BC each BC node belongs.
    // TODO: QMap<node,BiComp*> origNodesToBCs = makeGnodeToBCmap(BB);
    // Debug:
    if (debug) {
        qDebug() << "\nCompound nontrivial biconnected components:";
        foreach (BiComp *B, BB) {
            qDebug() << B->listNodes();
            B->colourShapes();
        }
        foreach (ExternalTree *E, EE) {
            E->colourShapes();
        }
    }

}

void Orthowontist::buildOGDFGraph(CanvasItemsList items,
        Graph &G, GraphAttributes &GA, shapemap &nodeShapes, connmap &edgeConns) {
    foreach (CanvasItem *item, items)
    {
        if (ShapeObj *shape = isShapeForLayout(item))
        {
            node n = G.newNode();
            QSizeF size = shape->size();
            GA.width(n) = size.width();
            GA.height(n) = size.height();
            nodeShapes.insert(n,shape);
            // Show IDs?
            bool showIDs = true;
            if (showIDs) {
                QString label = QString("%1").arg(shape->internalId());
                shape->setLabel(label);
            }
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
}

/** Expects a connected graph.
  * Prunes all external trees, and returns them in a list.
  * The graph itself IS altered.
  * Uses the idea from the TopoLayout paper.
  */
void Orthowontist::removeExternalTrees(QList<ExternalTree *> &EE, Graph &G,
                                       shapemap &nodeShapes, connmap &edgeConns) {
    // Make a map from degrees to lists of nodes that have that degree.
    // Initialize with one slot for each degree from 1 up to the maximum
    // degree occurring in the graph. This dense data structure might be
    // bad if there is a node of extremely high degree, but in our examples
    // we don't expect that.
    QMap< int,QSet<node> > nodesByDegree;
    // Find max degree D in graph.
    int D = 0;
    node v;
    forall_nodes(v,G) {
        int d = v->degree();
        if (d>D) D = d;
    }
    for (int i=1;i<=D;i++) {
        QSet<node> set;
        nodesByDegree.insert(i,set);
    }
    forall_nodes(v,G) {
        int d = v->degree();
        QSet<node> set = nodesByDegree.value(d);
        set.insert(v);
        nodesByDegree.insert(d,set);
    }

    // Create a new graph, where we'll make copies of the nodes we delete from G.
    Graph H;
    // Keep track of the mappings from the nodes and edges of H to the original
    // Dunnart graph.
    shapemap HNodesToDunnart;
    connmap HEdgesToDunnart;

    // For each node r of H, we will keep track of the node t in G which was its parent.
    // We imagine r as a "root" and t as a "taproot". If subsequently t gets removed from
    // G (and t' is its copy added to H), then we delete the mapping r-->t from the map.
    // This way the domain of the map always is precisely the set of root nodes in H.
    // (H will be such that each connected component is one of the external trees, and each
    // one will have a root.)
    QMap< node,node > rootsToTaproots; // Maps from H to G

    // Keep track of the nodes and edges belonging to each tree in H.
    QMap< node,QList<node> > rootsToNodes;
    QMap< node,QList<edge> > rootsToEdges;

    // The basic idea is simple: Continue deleting leaves until
    // there aren't any more.

    // Will need to map leaves in G to the edges of G that attached them, before removed.
    QMap<node,edge> leafEdges;

    while (nodesByDegree.value(1).size() > 0) {
        // Grab the degree-1 nodes, and replace with an empty set in the nodesByDegree map.
        QSet<node> degreeOneNodes = nodesByDegree.value(1);
        QSet<node> set;
        nodesByDegree.insert(1,set);
        foreach (node r, degreeOneNodes) {
            // Create a copy in H.
            node rH = H.newNode();
            // Map to original Dunnart shape:
            HNodesToDunnart.insert( rH, nodeShapes.value(r) );
            // Prepare lists of the nodes and edges for the tree of which rH is the new root.
            // Note: the root IS included in the list of nodes belonging to the tree.
            QList<node> treeNodes;
            QList<edge> treeEdges;
            treeNodes.append(rH);
            // Was r the taproot for any nodes already in H?
            // If so, it should now be connected to them.
            QList<node> children = rootsToTaproots.keys(r);
            // 'children' should be nonempty only on 2nd pass, and be subset of domain
            // of leafEdges by that iteration.
            foreach (node c, children) {
                edge f = H.newEdge(rH,c);
                HEdgesToDunnart.insert( f, edgeConns.value( leafEdges.value(c) ) );
                treeEdges.append(f);
                treeNodes.append(rootsToNodes.value(c));
                treeEdges.append(rootsToEdges.value(c));
                rootsToTaproots.remove(c);
                rootsToNodes.remove(c);
                rootsToEdges.remove(c);
            }
            rootsToNodes.insert(rH,treeNodes);
            rootsToEdges.insert(rH,treeEdges);
            // Get rH's taproot, i.e. r's parent in G.
            edge e = NULL;
            node t = NULL;
            forall_adj_edges(e,r) { // There should be only one!
                leafEdges.insert(r,e);
                t = r==e->source() ? e->target() : e->source();
                rootsToTaproots.insert(rH,t);
            }
            // Now we can delete r.
            G.delNode(r); // Removes node and all incident edges from the graph.
            // Let t bubble up in the lists by degree, since its degree just dropped.
            int n = t->degree();
            QSet<node> degNplus1 = nodesByDegree.value(n+1);
            degNplus1.remove(t);
            nodesByDegree.insert(n+1,degNplus1);
            QSet<node> degN = nodesByDegree.value(n);
            degN.insert(t);
            nodesByDegree.insert(n,degN);
        }
    }
    // Now combine trees T1,...,Tn sharing a common taproot t into a single tree T.
    // T will keep t as its taproot, and will get a copy of t in H as its root r.
    // So the interpretation will now be different: t will be G's copy of r, instead of
    // the node in G to which G's copy of r attaches.
    QMap<node,node> rToT;
    QSet<node> taproots = rootsToTaproots.values().toSet(); // use set since don't want duplicates
    foreach (node t, taproots) {
        // Create a copy in H.
        node tH = H.newNode();
        // Map to original Dunnart shape.
        HNodesToDunnart.insert( tH, nodeShapes.value(t) );
        // Prepare lists of the nodes and edges for the tree of which tH is the new root.
        QList<node> treeNodes;
        QList<edge> treeEdges;
        treeNodes.append(tH);
        // Get the roots of the trees that share the taproot t.
        QList<node> children = rootsToTaproots.keys(t);
        // All trees T1,...,Tn with roots in 'children' will now be combined.
        foreach (node c, children) {
            edge f = H.newEdge(tH,c);
            HEdgesToDunnart.insert( f, edgeConns.value( leafEdges.value(c) ) );
            treeEdges.append(f);
            treeNodes.append(rootsToNodes.value(c));
            treeEdges.append(rootsToEdges.value(c));
            rootsToTaproots.remove(c);
            rootsToNodes.remove(c);
            rootsToEdges.remove(c);
        }
        rootsToNodes.insert(tH,treeNodes);
        rootsToEdges.insert(tH,treeEdges);
        rToT.insert(tH,t);
    }
    foreach (node tH, rToT.keys()) {
        rootsToTaproots.insert(tH,rToT.value(tH));
    }
    // Now construct an ExternalTree object for each tree in H.
    foreach (node rH, rootsToTaproots.keys()) {
        node rG = rootsToTaproots.value(rH);
        QList<node> nodes = rootsToNodes.value(rH);
        QList<edge> edges = rootsToEdges.value(rH);
        ExternalTree *E = new ExternalTree(
                        rH, rG, nodes, edges, HNodesToDunnart, HEdgesToDunnart
                    );
        EE.append(E);
    }
}

void Orthowontist::buildNBCs(QList<BiComp *> &BB, QSet<node> &cutnodes, Graph &G,
                              shapemap &nodeShapes, connmap &edgeConns) {
    BCTree bctree(G);
    Graph B = bctree.bcTree();
    node vB;
    forall_nodes(vB, B){
        // Skip cutnodes.
        if (bctree.typeOfBNode(vB)!=BCTree::BComp) continue;

        // Skip the component if it is too small.
        int n = bctree.numberOfNodes(vB);
        if (n < 3) continue;

        // Otherwise get the edges of the first OGDF graph G belonging to this component,
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
            node tgt = orig->target();
            gNodes.insert(tgt);
            if (bctree.typeOfGNode(tgt)==BCTree::CutVertex) gCutNodes.insert(tgt);
        }
        QList<node> gNodeList = gNodes.toList();
        QList<node> gCutNodeList = gCutNodes.toList();
        BiComp *bc = new BiComp(gNodeList, gEdges, gCutNodeList, nodeShapes, edgeConns);
        BB.append(bc);
    }
    // Finally, compute the set of cutnodes.
    foreach (BiComp *bc, BB) {
        QList<node> nontrivialCutnodes = bc->cutnodes();
        cutnodes.unite(nontrivialCutnodes.toSet());
    }
}

/* Fuse BCs that share cutnodes. Return list of BCs thus obtained.
 */
QList<BiComp*> Orthowontist::fuseBCs(QList<BiComp*> bicomps)
{
    // Prepare multimap from cutnode shapes to BCs they lie in.
    QMap<ShapeObj*,BiComp*> endpts;
    foreach (BiComp *bc, bicomps) {
        QList<node> cn = bc->cutnodes();
        foreach (node c, cn) {
            ShapeObj *sh = bc->getShape(c);
            endpts.insertMulti(sh,bc);
        }
    }
    // Prepare list of "compound BCs".
    QList<BiComp*> cbc;
    while (!bicomps.empty()) {
        // Do depth-first search starting from first BC, and
        // following hyperedges as given by the endpts map.
        QList<BiComp*> elements;
        BiComp *b = bicomps.first();
        b->dfs(endpts,elements);
        // Remove elements found from the main list.
        bicomps = bicomps.toSet().subtract(elements.toSet()).toList();
        // Fuse the elements found, and store as compound BC.
        BiComp *compound = elements.takeFirst();
        while (!elements.empty()) compound = compound->fuse(elements.takeFirst());
        cbc.append(compound);
    }
    return cbc;
}

}












































