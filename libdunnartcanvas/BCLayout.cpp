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

#include "libdunnartcanvas/BCLayout.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/connector.h"
#include "libdunnartcanvas/pluginshapefactory.h"

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

using namespace ogdf;

namespace dunnart {

// ----------------------------------------------------------------------------
// Chunk ----------------------------------------------------------------------

QList<Chunk*> Chunk::findNeighbours(node origCutNode, QList<Chunk *> allChunks)
{
    QList<Chunk*> nbrs;
    foreach (Chunk *ch, allChunks)
    {
        if (ch==this) continue;
        if (ch->containsOriginalNode(origCutNode)) nbrs.append(ch);
    }
    return nbrs;
}

// ----------------------------------------------------------------------------
// RootedTree -----------------------------------------------------------------

RootedTree::RootedTree(QList<node>& nodes, QSet<node>& cutnodes) :
    m_graph(NULL),
    m_root(NULL),
    m_parentCutNode(NULL),
    m_basept(QPointF(0,0)),
    m_relpt(QPointF(0,0))
{
    // Construct own copy of graph, maintaining maps to original graph.
    m_graph = new Graph;
    QSet<edge> edges;
    foreach (node n, nodes)
    {
        // Create own node rep.
        node m = m_graph->newNode();
        m_nodemap.insert(m,n);
        if (cutnodes.contains(n)) { m_cutNodes.append(m); }
        else { m_normalNodes.append(m); }
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

QList<node> RootedTree::getCutNodes()
{
    QList<node> cns;
    foreach (node cn, m_cutNodes)
    {
        cns.append(m_nodemap.value(cn));
    }
    return cns;
}

void RootedTree::setRelPt(QPointF p)
{
    m_relpt = p;
}

bool RootedTree::containsOriginalNode(node n)
{
    return m_nodemap.values().contains(n);
}

void RootedTree::recursiveLayout(shapemap& origShapes, node origBaseNode, QPointF cardinal, QList<node> cutnodes)
{
    // TODO: Use origBaseNode and cardinal.
    // If origBaseNode is NULL, then it is to be ignored.
    // Otherwise, it is a node from the original graph, and we are meant to
    // layout this subgraph starting from its local copy of that node, and
    // moving in the direction given by cardinal.
    constructDunnartGraph(origShapes, cardinal, cutnodes);
    QPointF bary = baryCentre();
    int fixedSep = 300; // magic number -- not v. good magic...
    foreach (Chunk *ch, m_children)
    {
        node origCutNode = ch->getParentCutNode();
        node localCutNode = m_nodemap.key(origCutNode);
        // Compute cardinal direction from barycentre to cut node.
        ShapeObj *cnShape = m_ownShapeMap.value(localCutNode);
        QPointF cnCentre = cnShape->centrePos();
        QPointF cardinal = nearestCardinal(cnCentre - bary);
        // Set relative point of child, and recurse.
        QPointF p = cnCentre + fixedSep*cardinal;
        ch->setRelPt(p);
        ch->recursiveLayout(origShapes, origCutNode, cardinal, cutnodes);
        fixedSep += 300;
    }
}

void RootedTree::constructDunnartGraph(shapemap& origShapes, QPointF cardinal,
                                       QList<node> cutnodes)
{
    // Store passed value.
    m_origShapeMap = origShapes;

    // Create Dunnart shape and connector objects for own graph.
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    foreach (node m, m_nodemap.keys())
    {
        ShapeObj *shape = factory->createShape("org.dunnart.shapes.ellipse");
        shape->setPosAndSize(QPointF(0,0), QSizeF(30,30));
        node n = m_nodemap.value(m);
        if (cutnodes.contains(n))
        {
            shape->setFillColour(QColor(0,192,0));
            int i = cutnodes.indexOf(n);
            shape->setLabel(QString::number(i));
        } else {
            shape->setFillColour(QColor(192,0,0));
        }
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
    ogdfTreeLayout(cardinal);
    //colaTreeLayout(cardinal);
}

void RootedTree::colaTreeLayout(QPointF cardinal)
{
    QMap<node,int> nodeIndices;
    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;
    cola::CompoundConstraints ccs;
    // make rectangles
    node v = NULL;
    int i = 0;
    forall_nodes(v,*m_graph)
    {
        ShapeObj *sh = m_ownShapeMap.value(v);
        QRectF rect = sh->boundingRect();
        double x = rect.left(); double X = rect.right();
        double y = rect.top();  double Y = rect.bottom();
        rs.push_back( new vpsc::Rectangle(x,X,y,Y) );
        nodeIndices.insert(v,i);
        i++;
    }
    // make edges and sep co's
    edge e = NULL;
    forall_edges(e,*m_graph)
    {
        node src = e->source();
        node tgt = e->target();
        int srcIndex = nodeIndices.value(src);
        int tgtIndex = nodeIndices.value(tgt);
        // edge
        es.push_back( cola::Edge(srcIndex, tgtIndex) );
        // sep co
        vpsc::Dim d = (cardinal.x()==0?vpsc::YDIM:vpsc::XDIM);
        int l, r;
        if (cardinal.x() + cardinal.y() > 0) { l = srcIndex; r = tgtIndex; }
        else { l = tgtIndex; r = srcIndex; }
        double gap = 0;
        ccs.push_back( new cola::SeparationConstraint(d,l,r,gap,false) );
    }
    // create and run layout
    double idealLength = 100;
    cola::ConstrainedFDLayout *fdlayout =
            new cola::ConstrainedFDLayout(rs,es,idealLength,true);
    fdlayout->setConstraints(ccs);
    fdlayout->run(true,true);
    // extract the node positions
    forall_nodes(v,*m_graph)
    {
        ShapeObj *sh = m_ownShapeMap.value(v);
        int i = nodeIndices.value(v);
        vpsc::Rectangle *r = rs.at(i);
        double x = r->getCentreX();
        double y = r->getCentreY();
        sh->setCentrePos(QPointF(x,y));
    }
}

void RootedTree::ogdfTreeLayout(QPointF cardinal)
{
    TreeLayout treelayout;
    int x = cardinal.x(); int y = cardinal.y();
    int n = x*x - x + 2*y*y + y; // maps (x,y) to log[i](x+iy)
    Orientation orient;
    switch(n)
    {
    case 0:
        orient = leftToRight; break;
    case 1:
        orient = topToBottom; break;
    case 2:
        orient = rightToLeft; break;
    case 3:
        orient = bottomToTop; break;
    }
    treelayout.orientation(orient);
    treelayout.rootSelection(TreeLayout::rootByCoord);
    treelayout.orthogonalLayout(true);

    // Now we must give our local copy of the parentCutNode an extreme coordinate,
    // so that it will be chosen as root.
    node localParent = m_nodemap.key(m_parentCutNode);
    ShapeObj *rootShape = m_ownShapeMap.value(localParent);
    // When all initial positions are (0,0), this is easy.
    rootShape->setCentrePos(-100*cardinal);

    /* If initial positions are not all (0,0), we will need to compute max and mins.
    switch(n)
    {
    case 0:
        // left to right, so root needs smallest x coord
        //double minX = DBL_MAX;
        //foreach ()
        ;
    case 1:
        orient = bottomToTop; break;
    case 2:
        orient = rightToLeft; break;
    case 3:
        orient = topToBottom; break;
    }
    */

    // Call the layout.
    GraphAttributes newNodesGA(*m_graph);
    // Use positions and sizes from own shapes.
    BCLayout::extractPosAndSize(m_ownShapeMap, newNodesGA);
    treelayout.call(newNodesGA);
    BCLayout::injectPositions(m_ownShapeMap, newNodesGA);
}

void RootedTree::recursiveDraw(Canvas *canvas, QPointF p)
{
    // Compute base pt.
    m_basept = p + m_relpt;
    // Compute shift so that parent cut node (if exists) is centred at basept.
    QPointF shift = m_basept;
    if (m_parentCutNode)
    {
        QPointF cn = m_ownShapeMap.value(m_nodemap.key(m_parentCutNode))->centrePos();
        shift -= cn;
    }
    foreach (ShapeObj *sh, m_ownShapeMap.values())
    {
        QPointF c = sh->centrePos();
        //c += m_basept;
        c += shift;
        sh->setCentrePos(c);
        canvas->addItem(sh);
    }

    foreach (Connector *conn, m_ownConnMap.values())
    {
        canvas->addItem(conn);
    }

    foreach (Chunk *ch, m_children)
    {
        ch->recursiveDraw(canvas, m_basept);
    }
}

void RootedTree::setChildren(QList<Chunk *> children)
{
    m_children = children;
}

void RootedTree::setParentCutNode(node cn)
{
    m_parentCutNode = cn;
    // Now set edge directions so that all flow is forward from this root node.
    node root = m_nodemap.key(cn);
    QList<edge> seen;
    QList<node> queue;
    queue.append(root);
    while (!queue.empty())
    {
        node n = queue.takeFirst();
        edge e;
        forall_adj_edges(e,n)
        {
            if (seen.contains(e)) continue;
            // Reverse edge if n not source.
            if (e->target()==n) m_graph->reverseEdge(e);
            // Now mark edge as seen, and add target to queue
            seen.append(e);
            queue.append(e->target());
        }
    }
    assert(m_graph->consistencyCheck());
}

node RootedTree::getParentCutNode()
{
    return m_parentCutNode;
}

QPointF RootedTree::baryCentre()
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

QPointF RootedTree::nearestCardinal(QPointF v)
{
    QPointF N(0,-1);
    QPointF S(0, 1);
    QPointF E( 1,0);
    QPointF W(-1,0);
    double x = v.x(); double y = v.y();
    QPointF u = (abs(y) > abs(x) ?
                     ( y < 0 ? N : S ) :
                     ( x < 0 ? W : E ) );
    return u;
}

QRectF RootedTree::bbox()
{
    QList<ShapeObj*> shapes = m_ownShapeMap.values();
    ShapeObj *sh = shapes.first();
    QRectF box = sh->boundingRect();
    for (int i = 1; i < shapes.size(); i++)
    {
        sh = shapes.at(i);
        QRectF b = sh->boundingRect();
        box = box.united(b);
    }
    return box;
}

// ----------------------------------------------------------------------------
// BiComp ---------------------------------------------------------------------

int BiComp::method = -1;

BiComp::BiComp(QList<edge>& edges, QList<node>& nodes, QList<node>& cutNodes, Graph& G) :
    m_graph(NULL),
    m_basept(QPointF(0,0)),
    m_relpt(QPointF(0,0)),
    m_parentCutNode(NULL)
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
    //diag:
    /*
    node vG;
    forall_nodes(vG,G)
    {
        bool okay = false;
        if (m_nodemap.values().contains(vG)) {
            okay = true;
        }
        qDebug() << okay;
    }
    */
    // Result: as expected, some are in the map, others are not.
    //
}

QList<node> BiComp::getCutNodes()
{
    QList<node> cns;
    foreach (node cn, m_cutNodes)
    {
        cns.append(m_nodemap.value(cn));
    }
    return cns;
}

void BiComp::setChildren(QList<Chunk *> children)
{
    m_children = children;
}

void BiComp::setRelPt(QPointF p)
{
    m_relpt = p;
}

void BiComp::removeSelf(Graph &G)
{
    // Remove all edges.
    foreach (edge e, m_edgemap.values())
    {
        G.delEdge(e);
        //G.hideEdge(e);
        assert(G.consistencyCheck());
    }
}

size_t BiComp::size()
{
    return m_nodemap.size();
}


// pass shapemap for the shapes in the original graph
void BiComp::constructDunnartGraph(shapemap& origShapes,
                                   QPointF cardinal, QList<node> cutnodes)
{
    // Create Dunnart shape and connector objects for own graph.
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    foreach (node m, m_nodemap.keys())
    {
        ShapeObj *shape = factory->createShape("org.dunnart.shapes.ellipse");
        //shape->setPosAndSize(QPointF(0,0), QSizeF(30,30));
        shape->setCentrePos(QPointF(0,0));
        shape->setSize(QSizeF(30,30));
        node n = m_nodemap.value(m);
        if (cutnodes.contains(n))
        {
            shape->setFillColour(QColor(0,192,0));
            int i = cutnodes.indexOf(n);
            shape->setLabel(QString::number(i));
        } else {
            shape->setFillColour(QColor(0,0,192));
        }
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

    // Compute initial layout
    GraphAttributes newNodesGA(*m_graph);
    //BCLayout::extractSizes(m_nodemap, origShapes, newNodesGA);

    // Get distinct initial positions
    // FIXME: just use spiral initial positions, as usual;
    // don't need random, and it takes longer!
    while (coincidence()) { jog(100.0); }
    // Displace cut node in opposite of cardinal direction /if/ this is not
    // the root component of the layout.
    //
    // FIXME:
    // Really, we should be using cola for this initial layout, and we should
    // use separation constraints to say that the cut node be on the side
    // where we want it.
    //
    if (m_parentCutNode)
    {
        node root = m_nodemap.key(m_parentCutNode);
        ShapeObj *sh = m_ownShapeMap.value(root);
        QPointF c = sh->centrePos();
        c -= 100.0*cardinal;
        sh->setCentrePos(c);
    }
    // Load positions and sizes into GraphAttributes object
    BCLayout::extractPosAndSize(m_ownShapeMap, newNodesGA);

    FMMMLayout fm3;
    SpringEmbedderFR sefr;
    StressMajorization stmj;

    PlanarizationLayout planar;
    //planar.setEmbedder(new EmbedderMaxFace);
    //planar.setEmbedder(new EmbedderMinDepth);
    //planar.setEmbedder(new EmbedderMinDepthPiTa);
    //planar.setEmbedder(new EmbedderMinDepthMaxFace);
    //planar.setEmbedder(new EmbedderMaxFaceLayers);
    planar.setEmbedder(new EmbedderMinDepthMaxFaceLayers);

    PlanarizationGridLayout planarGrid;
    PlanarDrawLayout planarDraw;
    PlanarStraightLayout planarStraight;

    //int layoutAlg = 4;
    //switch(layoutAlg)
    switch(BiComp::method)
    {
    case 0:
        // Cola FD layout
        colaLayout();
        break;
    case 1:
        // FM3
        fm3.call(newNodesGA);
        break;
    case 2:
        // Fruchterman-Reingold spring embedder
        sefr.call(newNodesGA);
        break;
    case 3:
        // Kamada-Kawai
        stmj.call(newNodesGA);
        break;
    case 4:
        // OGDF planarization layout
        planar.call(newNodesGA);
        break;
    case 5:
        // OGDF planarization grid layout
        planarGrid.call(newNodesGA);
        break;
    case 6:
        // other OGDF planar layout methods
        //planarDraw.call(newNodesGA);
        //planarStraight.call(newNodesGA);
        break;
    default:
        break;
    }
    if (BiComp::method != 0) {
        BCLayout::injectPositions(m_ownShapeMap, newNodesGA);
    }
}

void BiComp::colaLayout()
{
    QMap<node,int> nodeIndices;
    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;
    cola::CompoundConstraints ccs;
    // make rectangles
    node v = NULL;
    int i = 0;
    forall_nodes(v,*m_graph)
    {
        ShapeObj *sh = m_ownShapeMap.value(v);
        QRectF rect = sh->boundingRect();
        QPointF c = sh->centrePos();
        double x = rect.left() + c.x(); double X = rect.right() + c.x();
        double y = rect.top() + c.y();  double Y = rect.bottom() + c.y();
        rs.push_back( new vpsc::Rectangle(x,X,y,Y) );
        nodeIndices.insert(v,i);
        i++;
    }
    // make edges
    edge e = NULL;
    forall_edges(e,*m_graph)
    {
        node src = e->source();
        node tgt = e->target();
        int srcIndex = nodeIndices.value(src);
        int tgtIndex = nodeIndices.value(tgt);
        // edge
        es.push_back( cola::Edge(srcIndex, tgtIndex) );
    }
    // create and run layout
    double idealLength = 100;
    cola::ConstrainedFDLayout *fdlayout =
            new cola::ConstrainedFDLayout(rs,es,idealLength,false);

    // As an experiment, can we add a constraint to align all
    // the rectangles vertically?
    /*
    cola::AlignmentConstraint *ac = new cola::AlignmentConstraint(vpsc::XDIM);
    for (int j = 0; j < rs.size(); j++) {
        ac->addShape(j,0);
    }
    ccs.push_back(ac);
    */
    // Answer: Yes, this works perfectly.

    // Do an initial FD layout.
    fdlayout->setConstraints(ccs);
    fdlayout->run(true,true);

#define doACA
#ifdef doACA
    // _ACA_
    // Let N be the number of nodes.
    int N = rs.size();
    // Prepare the alignment state matrix.
    Matrix2d<int> alignmentState = initACA(N, nodeIndices);
    // Start main loop.
    SeparatedAlignment *sa = chooseSA(rs, alignmentState);
    while (sa) {
        // Add the new separated alignment constraints.
        ccs.push_back(sa->separation);
        ccs.push_back(sa->alignment);
        // Redo the layout, with the new constraints.
        fdlayout->setConstraints(ccs);
        fdlayout->run(true,true);
        // Update state.
        updateAlignmentState(sa, alignmentState);
        // Choose next SA.
        sa = chooseSA(rs, alignmentState);
    }
#endif

    // extract the node positions
    forall_nodes(v,*m_graph)
    {
        ShapeObj *sh = m_ownShapeMap.value(v);
        int i = nodeIndices.value(v);
        vpsc::Rectangle *r = rs.at(i);
        double x = r->getCentreX();
        double y = r->getCentreY();
        sh->setCentrePos(QPointF(x,y));
    }
}

Matrix2d<int> BiComp::initACA(int N, QMap<node,int> nodeIndices)
{
    Matrix2d<int> alignmentState = Matrix2d<int>(N,N);
    // Initialize with zeros.
    for (uint i = 0; i < N; i++) {
        for (uint j = 0; j < N; j++) {
            alignmentState(i,j) = 0;
        }
    }
    // Note connected nodes.
    edge ed = NULL;
    forall_edges(ed,*m_graph) {
        node src = ed->source();
        node tgt = ed->target();
        int srcIndex = nodeIndices.value(src);
        int tgtIndex = nodeIndices.value(tgt);
        alignmentState(srcIndex,tgtIndex) = Canvas::Connected;
        alignmentState(tgtIndex,srcIndex) = Canvas::Connected;
    }
    return alignmentState;
}

SeparatedAlignment *BiComp::chooseSA(vpsc::Rectangles rs, Matrix2d<int> &alignmentState)
{
    SeparatedAlignment *sa = NULL;
    edge ed = NULL;
    forall_edges(ed,*m_graph) {
        node src = ed->source();
        node tgt = ed->target();
        // oops, need nodeIndices again as argument... TODO
        int srcIndex = nodeIndices.value(src);
        int tgtIndex = nodeIndices.value(tgt);
        // If already aligned, continue to next edge.
        if (alignmentState(srcIndex,tgtIndex) & (Canvas::Horizontal|Canvas::Vertical)) continue;
        // Otherwise consider the rectangles.
        vpsc::Rectangle rsrc = rs.at(srcIndex);
        vpsc::Rectangle rtgt = rs.at(tgtIndex);
        //...
    }
    return sa;
}

void BiComp::updateAlignmentState(SeparatedAlignment *sa, Matrix2d<int> &alignmentState)
{
    // TODO
}

// FIXME: Do this right, in n log n time, instead of n^2, if we're really doing this.
bool BiComp::coincidence()
{
    QList<ShapeObj*> shapes = m_ownShapeMap.values();
    for (int i = 0; i < shapes.size(); i++) {
        QPointF p1 = shapes.at(i)->centrePos();
        for (int j = i+1; j < shapes.size(); j++) {
            QPointF p2 = shapes.at(j)->centrePos();
            if (p1==p2) return true;
        }
    }
    return false;
}

void BiComp::jog(double scale)
{
    QTime t = QTime::currentTime();
    int seed = t.msecsTo(QTime(0,0,0,0));
    srand(seed);
    foreach (ShapeObj *sh, m_ownShapeMap.values())
    {
        QPointF c = sh->centrePos();
        double dx = (rand() / static_cast<double>( RAND_MAX ) - 0.5)*scale;
        double dy = (rand() / static_cast<double>( RAND_MAX ) - 0.5)*scale;
        c += QPointF(dx,dy);
        sh->setCentrePos(c);
    }
}

Graph& BiComp::copyGraph(QMap<node, node>& nodemap)
{
    Graph *G = new Graph();
    node v = NULL;
    forall_nodes(v,*m_graph)
    {
        node u = G->newNode();
        nodemap.insert(u,v);
    }
    edge e = NULL;
    forall_edges(e,*m_graph)
    {
        node s = e->source();
        node t = e->target();
        node m = nodemap.key(s);
        node n = nodemap.key(t);
        G->newEdge(m,n);
    }
    return *G;
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
    QPointF N(0,-1);
    QPointF S(0, 1);
    QPointF E( 1,0);
    QPointF W(-1,0);
    double x = v.x(); double y = v.y();
    QPointF u = (abs(y) > abs(x) ?
                     ( y < 0 ? N : S ) :
                     ( x < 0 ? W : E ) );
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
 * Find all chunks that contain (an image of) a cut node from the original graph,
 * and which are different from the present chunk.
 */
QList<Chunk*> BiComp::findCutNodeNeighbours(node origCutNode, bclist& bcs, treelist& trees)
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

void BiComp::setParentCutNode(node cn)
{
    m_parentCutNode = cn;
}

node BiComp::getParentCutNode()
{
    return m_parentCutNode;
}

QRectF BiComp::bbox()
{
    QList<ShapeObj*> shapes = m_ownShapeMap.values();
    ShapeObj *sh = shapes.first();
    QRectF box = sh->boundingRect();
    for (int i = 1; i < shapes.size(); i++)
    {
        sh = shapes.at(i);
        QRectF b = sh->boundingRect();
        box = box.united(b);
    }
    return box;
}

void BiComp::recursiveLayout(shapemap& origShapes, node origBaseNode,
                             QPointF cardinal, QList<node> cutnodes)
{
    // TODO: Use origBaseNode and cardinal.
    // If origBaseNode is NULL, then it is to be ignored.
    // Otherwise, it is a node from the original graph, and we are meant to
    // layout this subgraph starting from its local copy of that node, and
    // moving in the direction given by cardinal.
    constructDunnartGraph(origShapes, cardinal, cutnodes);
    QPointF bary = baryCentre();
    int fixedSep = 300; // magic number
    foreach (Chunk *ch, m_children)
    {
        node origCutNode = ch->getParentCutNode();
        node localCutNode = m_nodemap.key(origCutNode);
        // Compute cardinal direction from barycentre to cut node.
        ShapeObj *cnShape = m_ownShapeMap.value(localCutNode);
        QPointF cnCentre = cnShape->centrePos();
        QPointF cardinal = nearestCardinal(cnCentre - bary);
        // Set relative point of child, and recurse.
        QPointF p = cnCentre + fixedSep*cardinal;
        ch->setRelPt(p);
        ch->recursiveLayout(origShapes, origCutNode, cardinal, cutnodes);
        fixedSep += 300;
    }
}

void BiComp::recursiveDraw(Canvas *canvas, QPointF p)
{
    // Compute base pt.
    m_basept = p + m_relpt;
    // Compute shift so that parent cut node (if exists) is centred at basept.
    QPointF shift = m_basept;
    if (m_parentCutNode)
    {
        QPointF cn = m_ownShapeMap.value(m_nodemap.key(m_parentCutNode))->centrePos();
        shift -= cn;
    }
    foreach (ShapeObj *sh, m_ownShapeMap.values())
    {
        QPointF c = sh->centrePos();
        //c += m_basept;
        c += shift;
        sh->setCentrePos(c);
        canvas->addItem(sh);
    }

    foreach (Connector *conn, m_ownConnMap.values())
    {
        canvas->addItem(conn);
    }

    foreach (Chunk *ch, m_children)
    {
        ch->recursiveDraw(canvas, m_basept);
    }
}

// ----------------------------------------------------------------------------
// BCLayout -------------------------------------------------------------------

BCLayout::BCLayout(Canvas *canvas) :
    m_canvas(canvas)
{}

/* Builds the tree structure on the Chunks simply by setting their own internal
 * records of what their children are, and what their "parent cut nodes" are.
 * Also records the list of cut nodes that were actually used.
 */
void BCLayout::buildBFSTree(QList<Chunk *> chunks, Chunk *root, QList<node>& usedCutNodes)
{
    QList<Chunk*> queue;
    QSet<Chunk*> seen;
    queue.push_back(root);
    seen.insert(root);
    while (!queue.empty())
    {
        Chunk *ch = queue.takeFirst();
        QList<node> cutnodes = ch->getCutNodes();
        QList<Chunk*> children;
        foreach (node cn, cutnodes)
        {
            QList<Chunk*> nbrs = ch->findNeighbours(cn, chunks);
            QSet<Chunk*> nbrSet = nbrs.toSet();
            nbrSet.subtract(seen);
            QList<Chunk*> nbrList = nbrSet.toList();
            foreach (Chunk *ch, nbrList) {
                ch->setParentCutNode(cn);
                usedCutNodes.append(cn);
            }
            children.append(nbrList);
            seen.unite(nbrSet);
            queue.append(nbrList);
        }
        ch->setChildren(children);
    }
}

void BCLayout::ogdfGraph(Graph& G, shapemap &nodeShapes, connmap &edgeConns)
{
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
}

void BCLayout::extractSizes(shapemap& nodeShapes, GraphAttributes& GA)
{
    foreach (node v, nodeShapes.keys())
    {
        ShapeObj *shape = nodeShapes.value(v);
        GA.width()[v] = shape->size().width();
        GA.height()[v] = shape->size().height();
    }
}

void BCLayout::extractSizes(
        QMap<node, node> &newToOldNodes, shapemap &oldNodesToShapes, GraphAttributes &newNodesGA)
{
    foreach (node u, newToOldNodes.keys())
    {
        node v = newToOldNodes.value(u);
        ShapeObj *shape = oldNodesToShapes.value(v);
        newNodesGA.width()[u] = shape->size().width();
        newNodesGA.height()[u] = shape->size().height();
    }
}

void BCLayout::extractPosAndSize(
        QMap<node, node> &newToOldNodes, shapemap &oldNodesToShapes, GraphAttributes &newNodesGA)
{
    foreach (node u, newToOldNodes.keys())
    {
        node v = newToOldNodes.value(u);
        ShapeObj *shape = oldNodesToShapes.value(v);
        newNodesGA.width()[u] = shape->size().width();
        newNodesGA.height()[u] = shape->size().height();
        newNodesGA.x(u) = shape->centrePos().x();
        newNodesGA.y(u) = shape->centrePos().y();
    }
}

void BCLayout::extractPosAndSize(shapemap &newNodesToShapes, GraphAttributes &newNodesGA)
{
    foreach (node u, newNodesToShapes.keys())
    {
        ShapeObj *shape = newNodesToShapes.value(u);
        newNodesGA.width()[u] = shape->size().width();
        newNodesGA.height()[u] = shape->size().height();
        newNodesGA.x(u) = shape->centrePos().x();
        newNodesGA.y(u) = shape->centrePos().y();
    }
}

void BCLayout::injectPositions(shapemap& nodeShapes, GraphAttributes& GA)
{
    foreach (node v, nodeShapes.keys())
    {
        double cx = GA.x(v) + GA.width(v)/2.0;
        double cy = GA.y(v) + GA.height(v)/2.0;
        ShapeObj *shape = nodeShapes.value(v);
        shape->setCentrePos(QPointF(cx,cy));
    }
}

void BCLayout::injectSizes(shapemap& nodeShapes, GraphAttributes& GA)
{
    foreach (node v, nodeShapes.keys())
    {
        double w = GA.width(v);
        double h = GA.height(v);
        ShapeObj *shape = nodeShapes.value(v);
        shape->setSize(QSizeF(w,h));
    }
}

QList<BiComp*> BCLayout::getNontrivialBCs(Graph& G, QSet<node>& cutnodes)
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
            node tgt = orig->target();
            gNodes.insert(tgt);
            if (bctree.typeOfGNode(tgt)==BCTree::CutVertex) gCutNodes.insert(tgt);
            //diag:
            /*
            bool srcokay = false;
            bool tgtokay = false;
            node vG;
            forall_nodes(vG,G)
            {
                if (vG==src) {srcokay = true;}
                if (vG==tgt) {tgtokay = true;}
            }
            qDebug() << srcokay << tgtokay;
            */
            // Result: Yes, the src and tgt nodes do always belong to G.
            //
        }
        QList<node> gNodeList = gNodes.toList();
        QList<node> gCutNodeList = gCutNodes.toList();
        BiComp *bc = new BiComp(gEdges, gNodeList, gCutNodeList, G);
        bicomps.append(bc);
    }
    node vG;
    forall_nodes(vG, G)
    {
        if (bctree.typeOfGNode(vG)==BCTree::CutVertex) cutnodes.insert(vG);
    }
    return bicomps;
}

QMap<int,node> BCLayout::getConnComps(Graph& G)
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

QMap<int,node> BCLayout::getConnComps2(Graph *G2, QMap<node, node>& nodeMapG2ToG)
{
    // G2: the second graph we construct
    // nodeMapG2ToG: map from nodes of G2 back to nodes of G, the graph
    // whose connected components we actually will return.
    // We return the connected components in the form of a multimap from
    // integers to nodes. The components are assigned integer labels, and
    // the map sends integer n to all those nodes belonging to the ccomp
    // whose label is n.

    NodeArray<int> ccomps(*G2);
    connectedComponents(*G2, ccomps);
    QMap<int,node> map;
    node v2 = NULL;
    forall_nodes(v2,*G2)
    {
        int n = ccomps[v2];
        node v = nodeMapG2ToG.value(v2);
        map.insertMulti(n,v);
    }
    return map;
}

Graph *BCLayout::removeBiComps(Graph& G, bclist& bcs, QMap<node, node>& nodeMapNewToOld)
{
    // We construct and return a new graph which is equivalent to what you get
    // by starting with G and then deleting all edges in the BiComps listed in
    // the list bcs.
    //
    // The node map will be filled in with a mapping from nodes of the new graph
    // back to the corresponding nodes in the old graph.

    Graph *Gp = new Graph;

    // Copy nodes.
    node v = NULL;
    forall_nodes(v,G)
    {
        node u = Gp->newNode();
        nodeMapNewToOld.insert(u,v);
    }
    // Copy edges, except those in the BiComps listed in bcs.
    edge e = NULL;
    // FIXME: Might it be more efficient to first construct the set S
    // of all original edges from the BCs, and /then/ pass through all
    // the edges of G and keep only those that are not in S?
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
            node m = nodeMapNewToOld.key(s);
            node n = nodeMapNewToOld.key(t);
            Gp->newEdge(m,n);
        }
    }

    return Gp;
}

void BCLayout::orthoLayout(int method)
{
    BiComp::method = method;
    shapemap nodeShapes;
    connmap edgeConns;
    Graph *Gp = new Graph;
    Graph G = *Gp;
    ogdfGraph(G, nodeShapes, edgeConns);
    //diag:
    /*
    node v;
    forall_nodes(v,G)
    {
        bool okay = false;
        if (nodeShapes.keys().contains(v))
        {
            okay = true;
        }
        qDebug() << okay;
    }
    */
    // Result: Yes, every node in G is a key in nodeShapes.
    //

    // Get nontrivial biconnected components (size >= 3), and get the set of cutnodes in G.
    QSet<node> cutnodes;
    QList<BiComp*> bicomps = getNontrivialBCs(G, cutnodes);

    // Get a new graph isomorphic to the result of removing the BC edges from G.
    QMap<node,node> nodeMapG2ToG;
    Graph *G2 = removeBiComps(G, bicomps, nodeMapG2ToG);
    assert(G2->consistencyCheck());

    // Get connected components of original graph G corresponding to those of G2.
    QMap<int,node> ccomps = getConnComps2(G2, nodeMapG2ToG);

    // Form trees on those components, throwing away isolated
    // nodes (which must have been nodes belonging to nontrivial BCs).
    QList<RootedTree*> rtrees;
    foreach (int i, ccomps.keys().toSet())
    {
        QList<node> nodes = ccomps.values(i);
        if (nodes.size() < 2) continue;
        RootedTree *rtree = new RootedTree(nodes, cutnodes);
        rtrees.append(rtree);
    }

    // Choose a largest biconnected component to be root.
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

    /* TODO
     * ====
     *
     * Here is where I will start to change the algorithm, according to
     * our latest ideas on how reassembly should work.
     *
     * For now all we have is a BFS turning all the Chunks into a tree,
     * in fact a rooted tree where we have chosen a largest Chunk as the
     * root.
     *
     */

    // Build tree of components.
    QList<Chunk*> allChunks;
    foreach (BiComp *bc, bicomps) allChunks.append(dynamic_cast<Chunk*>(bc));
    foreach (RootedTree *rt, rtrees) allChunks.append(dynamic_cast<Chunk*>(rt));
    QList<node> usedCutNodes;
    BCLayout::buildBFSTree(allChunks, largest, usedCutNodes);

    // Layout and draw
    largest->recursiveLayout(nodeShapes, NULL, QPointF(0,0), usedCutNodes);
    m_canvas->stop_graph_layout();
    largest->recursiveDraw(m_canvas, QPointF(0,0));
    m_canvas->interrupt_graph_layout();
}

void BCLayout::applyFM3()
{
    shapemap nodeShapes;
    connmap edgeConns;
    Graph *G = new Graph;
    ogdfGraph(*G, nodeShapes, edgeConns);
    GraphAttributes GA(*G);
    extractSizes(nodeShapes, GA);
    FMMMLayout fm3;
    fm3.call(GA);
    injectPositions(nodeShapes, GA);
}

void BCLayout::layoutBCTrees()
{
    //Construct graph, and label nodes by ID.
    ogdf::Graph G;
    QMap<ogdf::node,ShapeObj*> nodeShapes;
    QMap<ogdf::edge,Connector*> edgeConns;
    // First get all the nodes.
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
    // Now do the edges on a separate, second pass, so that all their
    // endpoints already exist.
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
    int p = 0;
    forall_nodes(v,B)
    {
        sh = factory->createShape("org.dunnart.shapes.ellipse");
        ogdf::BCTree::BNodeType bnt = bctree.typeOfBNode(v);
        int n = bctree.numberOfNodes(v);
        // Set the radius of the shape based on the number of nodes
        // in the BC represented by v, if it is a BComp.
        double s = 25 + (150/3.14)*std::atan(double(n-1)*0.173);
        // Choose an initial point.
        double x = 0, y = 0;
        switch (p%4)
        {
        case 0:
            x = (double)p; y = 0; break;
        case 1:
            x = 0; y = (double)p; break;
        case 2:
            x = -(double)p; y = 0; break;
        case 3:
            x = 0; y = -(double)p; break;
        }
        p++;
        // Construct the shape.
        switch (bnt)
        {
        case ogdf::BCTree::BComp:
            sh->setPosAndSize(QPointF(x,y), QSizeF(s,s));
            sh->setFillColour(QColor(128,128,255));
            break;
        case ogdf::BCTree::CComp:
            sh->setPosAndSize(QPointF(x,y), QSizeF(25,25));
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

    // Looks like we tried to use the auxiliaryGraph function, but
    // didn't like what we were getting?
    //ogdf::Graph H = bctree.auxiliaryGraph();

    QSet<ogdf::node> hNodes;
    QSet<ogdf::edge> hEdges;
    forall_nodes(v,B)
    {
        int n = bctree.numberOfNodes(v);
        if (n < 3) continue;
        // If we get this far, then v represents a BComp containing at least
        // 3 nodes of the original graph.
        // We now use the hEdges function which returns a list of all edges
        // of the original graph belonging to this BComp.
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
    p = 0;
    foreach (ogdf::node v, hNodes)
    {
        // Choose an initial point.
        double x = 0, y = 0;
        switch (p%4)
        {
        case 0:
            x = (double)p; y = 0; break;
        case 1:
            x = 0; y = (double)p; break;
        case 2:
            x = -(double)p; y = 0; break;
        case 3:
            x = 0; y = -(double)p; break;
        }
        p++;
        // Construct the shape.
        sh = factory->createShape("org.dunnart.shapes.ellipse");
        sh->setPosAndSize(QPointF(x,y), QSizeF(25,25));
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
