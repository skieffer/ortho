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

void RootedTree::recursiveLayout(shapemap& origShapes, node origBaseNode,
                                 QPointF cardinal, QList<node> cutnodes)
{
    // TODO: Use origBaseNode and cardinal.
    // If origBaseNode is NULL, then it is to be ignored.
    // Otherwise, it is a node from the original graph, and we are meant to
    // layout this subgraph starting from its local copy of that node, and
    // moving in the direction given by cardinal.
    constructDunnartGraph(origShapes, cardinal, cutnodes);
    QPointF bary = baryCentre();
    int magicNumber = 0;
    int fixedSep = magicNumber;
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
        fixedSep += magicNumber;
    }
}

void RootedTree::constructDunnartGraph(shapemap& origShapes, QPointF cardinal,
                                       QList<node> cutnodes)
{
//#define multicoloured
#ifdef multicoloured
    QString shapeName = "org.dunnart.shapes.ellipse";
    QSizeF shapeSize = QSizeF(30,30);
    QColor cutnodeColor = QColor(0,192,0);
    QColor normalnodeColor = QColor(192,0,0);
#else
    QString shapeName = "org.dunnart.shapes.rectangle";
    QSizeF shapeSize = QSizeF(30,30);
    QColor cutnodeColor = QColor(255,192,0);
    QColor normalnodeColor = QColor(255,192,0);
#endif
    // Store passed value.
    m_origShapeMap = origShapes;

    // Create Dunnart shape and connector objects for own graph.
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    foreach (node m, m_nodemap.keys())
    {
        ShapeObj *shape = factory->createShape(shapeName);
        shape->setPosAndSize(QPointF(0,0), shapeSize);
        node n = m_nodemap.value(m);
        if (cutnodes.contains(n))
        {
            if (n==m_parentCutNode) {
                // Replace shape by the shape that m_parent associated with
                // this node.
                shape = m_parent->getShapeForOriginalNode(n);
            } else {
                shape->setFillColour(cutnodeColor);
                int i = cutnodes.indexOf(n);
                shape->setLabel(QString::number(i));
            }
        } else {
            shape->setFillColour(normalnodeColor);
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
    m_orientation = n;
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
    // Infer constraints.
    inferConstraints();
}

/* After calling ogdfTreeLayout, we use this method to infer apt constraints
 * to preserve the layout that our shapes have been given.
 * Store them in m_dunnartConstraints.
 */
void RootedTree::inferConstraints()
{
    // We use a top-down search through the tree. This works because
    // the tree already has a layout. (If you were starting from scratch
    // you would have to work bottom-up.)
    node localRoot = m_nodemap.key(m_parentCutNode);
    QList<node> queue;
    QList<node> parentQueue; // keeps track of parents
    queue.push_back(localRoot);
    parentQueue.push_back(NULL);
    while (!queue.empty()) {
        node A = queue.takeFirst();
        node parent = parentQueue.takeFirst();
        // Make list of the children of A.
        QList<node> children;
        edge e = NULL;
        forall_adj_edges(e,A) {
            node B = A==e->source() ? e->target() : e->source();
            if (B==parent) continue;
            children.append(B);
        }
        int numChil = children.size();
        // If no children, continue.
        if (numChil==0) continue;
        // Else add children to the queue, and add as many copies of
        // node A to the parent queue.
        queue.append(children);
        for(int i=0;i<numChil;i++){ parentQueue.append(A); }

        // Define constraints -------------------------------------------
        DunnartConstraint *dc;

        // Align children.
        dc = new DunnartConstraint();
        dc->type = ALIGNMENT;
        dc->dim = m_orientation%2==0 ? Canvas::VERT : Canvas::HORIZ;
        for(int i=0;i<numChil;i++){ dc->items.append(m_ownShapeMap.value(children.at(i))); }
        m_dunnartConstraints.append(dc);

        // Let's try skipping the rest if the parent is the root node.
        if (A==localRoot) continue;

        // Separate parent and children.
        dc = new DunnartConstraint();
        dc->type = SEPARATION;
        dc->dim = m_orientation%2==0 ? Canvas::HORIZ : Canvas::VERT;
        dc->items.append(m_ownShapeMap.value(A));
        dc->items.append(m_ownShapeMap.value(children.at(0)));
        // Determine the min sep.
        QPointF p = m_ownShapeMap.value(A)->centrePos();
        QPointF q = m_ownShapeMap.value(children.at(0))->centrePos();
        double z = m_orientation%2==0 ? p.x() : p.y();
        double w = m_orientation%2==0 ? q.x() : q.y();
        dc->minSep = fabs(z-w);
        m_dunnartConstraints.append(dc);

        // Centre parent above children, and distribute children.
        if (numChil==1) {
            // In this case need just a single alignment.
            dc = new DunnartConstraint();
            dc->type = ALIGNMENT;
            dc->dim = m_orientation%2==0 ? Canvas::HORIZ : Canvas::VERT;
            dc->items.append(m_ownShapeMap.value(A));
            dc->items.append(m_ownShapeMap.value(children.at(0)));
            m_dunnartConstraints.append(dc);
        } else {
            // In this case need two distributions.
            // First distribute all the children:
            dc = new DunnartConstraint();
            dc->type = DISTRIBUTION;
            dc->dim = m_orientation%2==0 ? Canvas::VERT : Canvas::HORIZ;
            for(int i=0;i<numChil;i++){ dc->items.append(m_ownShapeMap.value(children.at(i))); }
            m_dunnartConstraints.append(dc);
            // Now the parent and the two outer children:
            dc = new DunnartConstraint();
            dc->type = DISTRIBUTION;
            dc->dim = m_orientation%2==0 ? Canvas::VERT : Canvas::HORIZ;
            dc->items.append(m_ownShapeMap.value(A));
            // Find the two outermost children.
            node minChild = NULL, maxChild = NULL;
            double minZ = DBL_MAX, maxZ = DBL_MIN;
            foreach (node c, children) {
                QPointF p = m_ownShapeMap.value(c)->centrePos();
                double z = m_orientation%2==0 ? p.y() : p.x();
                if (z<minZ) {minZ=z; minChild=c;}
                if (z>maxZ) {maxZ=z; maxChild=c;}
            }
            dc->items.append(m_ownShapeMap.value(minChild));
            dc->items.append(m_ownShapeMap.value(maxChild));
            m_dunnartConstraints.append(dc);
        }
        // End constraint definitions -----------------------------------
    }
}

/* To be called by recursiveDraw. Applies the Dunnart constraints that
 * were inferred by inferConstraints.
 */
void RootedTree::applyDunnartConstraints()
{
    foreach (DunnartConstraint *dc, m_dunnartConstraints) {
        switch (dc->type) {
        case ALIGNMENT:
        {
            atypes at = dc->dim==Canvas::VERT ? ALIGN_CENTER : ALIGN_MIDDLE;
            createAlignment(at,dc->items);
            break;
        }
        case SEPARATION:
        {
            dtype dt = dc->dim==Canvas::VERT ? SEP_VERTICAL : SEP_HORIZONTAL;
            bool sort = true;
            createSeparation(NULL,dt,dc->items,dc->minSep,sort);
            break;
        }
        case DISTRIBUTION:
        {
            dtype dt = dc->dim==Canvas::VERT ? DIST_MIDDLE : DIST_CENTER;
            createDistribution(NULL,dt,dc->items);
            break;
        }
        }
    }
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

    // Apply any constraints inferred from ogdfTreeLayout.
    applyDunnartConstraints();
}

void RootedTree::setChildren(QList<Chunk *> children)
{
    m_children = children;
}

void RootedTree::setParent(BiComp *bc)
{
    m_parent = bc;
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
// ExternalTree ---------------------------------------------------------------

ExternalTree::ExternalTree(QList<node> nodes, QList<edge> edges, node root,
                           node taproot, QMap<node,int> dunnartIDs, int taprootID) :
    m_tapRoot(taproot),
    m_tapRootID(taprootID)
{
    m_graph = new Graph;
    QMap<node,node> nodemap; // maps passed nodes (from H) to own nodes
    foreach (node n, nodes) { // n belongs to H
        node m = m_graph->newNode(); // m belongs to own graph
        nodemap.insert(n,m);
        m_dunnartIDs.insert(m,dunnartIDs.value(n));
        if (n==root) m_root = m;
    }
    foreach (edge e, edges) {
        node m1 = nodemap.value(e->source());
        node m2 = nodemap.value(e->target());
        m_graph->newEdge(m1,m2);
    }
}

QString ExternalTree::listNodes()
{
    QString s = "";
    s += QString("Root: %1").arg(m_dunnartIDs.value(m_root));
    s += QString("\n  Taproot: %1").arg(m_tapRootID);
    s += QString("\n  Other nodes: ");
    foreach (node m, m_dunnartIDs.keys()) {
        if (m==m_root) continue;
        int id = m_dunnartIDs.value(m);
        s += QString("%1, ").arg(id);
    }
    return s;
}

void ExternalTree::treeLayout(void)
{
    // Set up an OGDF TreeLayout object.
    TreeLayout tL;
    tL.orientation(leftToRight);
    tL.rootSelection(TreeLayout::rootByCoord);
    tL.orthogonalLayout(true);
    // Set GraphAttributes.
    // For now we use a default size of 30x30 for the nodes.
    m_graphAttributes = new GraphAttributes(*m_graph);
    node n = NULL;
    forall_nodes(n,*m_graph) {
        m_graphAttributes->width(n) = 30;
        m_graphAttributes->height(n) = 30;
        if (n==m_root) {
            // Since we lay out left to right, we can make m_root
            // be selected by the TreeLayout as the root by setting it off
            // to the left. All other nodes will get an initial position of
            // (0,0), and the root will get (-100,0).
            m_graphAttributes->x(n) = -100;
        }
    }
    // Do the layout.
    tL.call(*m_graphAttributes);
}

/** Assuming the tree layout has already been performed, infer constraints
  * to maintain it.
  *
  * The Canvas::Dimension dim describes the dimension in which this tree
  * is oriented. E.g. if the root is at the top and the leaves at the bottom,
  * then dim == Canvas::VERT.
  */
void ExternalTree::inferConstraints(Canvas::Dimension dim)
{
    Canvas::Dimension counterDim = dim==Canvas::VERT ? Canvas::HORIZ : Canvas::VERT;
    // We use a top-down search through the tree. This works because
    // the tree already has a layout. (If you were starting from scratch
    // you would have to work bottom-up.)
    QList<node> queue;
    QList<node> parentQueue; // keeps track of parents
    queue.push_back(m_root);
    parentQueue.push_back(NULL);
    while (!queue.empty()) {
        node A = queue.takeFirst();
        node parent = parentQueue.takeFirst();
        // Make list of the children of A.
        QList<node> children;
        edge e = NULL;
        forall_adj_edges(e,A) {
            node B = A==e->source() ? e->target() : e->source();
            if (B==parent) continue;
            children.append(B);
        }
        int numChil = children.size();
        // If no children, continue.
        if (numChil==0) continue;
        // Else add children to the queue, and add as many copies of
        // node A to the parent queue.
        queue.append(children);
        for(int i=0;i<numChil;i++){ parentQueue.append(A); }

        // Define constraints -------------------------------------------
        DunnartConstraint *dc;

        // TODO (Copy from RootedTree, and modify as necessary.)
        dc = new DunnartConstraint();
        dc->type = ALIGNMENT;
        dc->dim = counterDim;
        for(int i=0;i<numChil;i++){ dc->items.append(m_shapeMap.value(children.at(i))); }
        m_dunnartConstraints.append(dc);

        // For now we skip the rest if the parent is the root node.
        if (A==m_root) continue;

        // Separate parent and children.
        dc = new DunnartConstraint();
        dc->type = SEPARATION;
        dc->dim = dim;
        dc->items.append(m_shapeMap.value(A));
        dc->items.append(m_shapeMap.value(children.at(0)));
        // Determine the min sep.
        QPointF p = m_shapeMap.value(A)->centrePos();
        QPointF q = m_shapeMap.value(children.at(0))->centrePos();
        double z = dim==Canvas::HORIZ ? p.x() : p.y();
        double w = dim==Canvas::HORIZ ? q.x() : q.y();
        dc->minSep = fabs(z-w);
        m_dunnartConstraints.append(dc);

        // Centre parent above children, and distribute children.
        if (numChil==1) {
            // In this case need just a single alignment.
            dc = new DunnartConstraint();
            dc->type = ALIGNMENT;
            dc->dim = dim;
            dc->items.append(m_shapeMap.value(A));
            dc->items.append(m_shapeMap.value(children.at(0)));
            m_dunnartConstraints.append(dc);
        } else {
            // In this case need two distributions.
            // First distribute all the children:
            dc = new DunnartConstraint();
            dc->type = DISTRIBUTION;
            dc->dim = counterDim;
            for(int i=0;i<numChil;i++){ dc->items.append(m_shapeMap.value(children.at(i))); }
            m_dunnartConstraints.append(dc);
            // Now the parent and the two outer children:
            dc = new DunnartConstraint();
            dc->type = DISTRIBUTION;
            dc->dim = counterDim;
            dc->items.append(m_shapeMap.value(A));
            // Find the two outermost children.
            node minChild = NULL, maxChild = NULL;
            double minZ = DBL_MAX, maxZ = DBL_MIN;
            foreach (node c, children) {
                QPointF p = m_shapeMap.value(c)->centrePos();
                double z = dim==Canvas::HORIZ ? p.y() : p.x();
                if (z<minZ) {minZ=z; minChild=c;}
                if (z>maxZ) {maxZ=z; maxChild=c;}
            }
            dc->items.append(m_shapeMap.value(minChild));
            dc->items.append(m_shapeMap.value(maxChild));
            m_dunnartConstraints.append(dc);
        }
        // End constraint definitions -----------------------------------
    }
}

// ----------------------------------------------------------------------------
// InternalTree ---------------------------------------------------------------

InternalTree::InternalTree(QList<node> nodes, QSet<node> cutnodes)
{
    // TODO
}


// ----------------------------------------------------------------------------
// BiComp ---------------------------------------------------------------------

int BiComp::method = -1;

BiComp::BiComp() :
    m_graph(NULL),
    m_basept(QPointF(0,0)),
    m_relpt(QPointF(0,0)),
    m_parentCutNode(NULL)
{}

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

/* Return a list of the cutnodes /of the original graph/ that
 * belong to this BC.
 */
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

void BiComp::setParent(BiComp *bc)
{
    m_parent = bc;
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

/* Return ShapeObj associated in course of constructDunnartGraph with
 * the passed node orig, which must be a member of the /original/ graph.
 * So constructDunnartGraph must have already been executed by now!
 */
ShapeObj *BiComp::getShapeForOriginalNode(node orig)
{
    // Get own node corresp. to orig.
    node m = m_nodemap.key(orig);
    // Get shape.
    ShapeObj *sh = m_ownShapeMap.value(m);
    return sh;
}

// pass shapemap for the shapes in the original graph
void BiComp::constructDunnartGraph(shapemap& origShapes,
                                   QPointF cardinal, QList<node> cutnodes)
{
#ifdef multicoloured
    QString shapeName = "org.dunnart.shapes.ellipse";
    QSizeF shapeSize = QSizeF(30,30);
    QColor cutnodeColor = QColor(0,192,0);
    QColor normalnodeColor = QColor(0,0,192);
#else
    QString shapeName = "org.dunnart.shapes.rectangle";
    QSizeF shapeSize = QSizeF(30,30);
    QColor cutnodeColor = QColor(255,192,0);
    QColor normalnodeColor = QColor(255,192,0);
#endif
    // Create Dunnart shape and connector objects for own graph.
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    foreach (node m, m_nodemap.keys())
    {
        ShapeObj *shape = factory->createShape(shapeName);
        //shape->setPosAndSize(QPointF(0,0), QSizeF(30,30));
        shape->setCentrePos(QPointF(0,0));
        shape->setSize(shapeSize);
        node n = m_nodemap.value(m);
        if (cutnodes.contains(n))
        {
            shape->setFillColour(cutnodeColor);
            int i = cutnodes.indexOf(n);
            shape->setLabel(QString::number(i));
        } else {
            shape->setFillColour(normalnodeColor);
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

void BiComp::acaLayout(void)
{
    m_graphAttributes = new GraphAttributes(*m_graph);
    // For now we just use a default size of 30x30 for the nodes.
    node n = NULL;
    forall_nodes(n,*m_graph) {
        m_graphAttributes->width(n) = 30;
        m_graphAttributes->height(n) = 30;
    }
    ACALayout aca = ACALayout(*m_graph, *m_graphAttributes);
    aca.run();
    aca.readLayout(*m_graph, *m_graphAttributes);
}

void BiComp::colaLayout(void)
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
    //qDebug() << "run once";
    fdlayout->run(true,true);

#define doACA
#ifdef doACA
    // _ACA_
    // Let N be the number of nodes.
    int N = rs.size();
    // Prepare the alignment state matrix.
    Matrix2d<int> alignmentState = initACA(N, nodeIndices);
    // Start main loop.
    SeparatedAlignment *sa = chooseSA(rs, alignmentState, nodeIndices);
    int n = 2; // for debugging purpose only
    while (sa) {
        // Add the new separated alignment constraints.
        ccs.push_back(sa->separation);
        ccs.push_back(sa->alignment);
        // Redo the layout, with the new constraints.
        //qDebug() << "run" << n; n++;

        QString rects = "";
        rects += QString("%1: ").arg(n,3);
        for (int i = 0; i < rs.size(); i++) {
            vpsc::Rectangle *r = rs.at(i);
            double x = r->getCentreX();
            double y = r->getCentreY();
            rects += QString("(%1,%2) ").arg(x,5,'f',2).arg(y,5,'f',2);
        }
        //qDebug() << rects;
        QString sepalgn = "";
        sepalgn += QString("%1").arg(sa->rect1,2);
        sepalgn += sa->af == Canvas::Horizontal ? "-" : "|";
        sepalgn += QString("%1").arg(sa->rect2,2);
        //qDebug() << sepalgn;
        n++;

        delete fdlayout;
        fdlayout = new cola::ConstrainedFDLayout(rs,es,idealLength,false);
        fdlayout->setConstraints(ccs);
        fdlayout->run(true,true);
        // Update state.
        updateAlignmentState(sa, alignmentState);
        // Store SA and choose next one.
        m_sepAligns.append(sa);
        sa = chooseSA(rs, alignmentState, nodeIndices);
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
        delete r;
        sh->setCentrePos(QPointF(x,y));
    }
    // Clean up.
    delete fdlayout;
}

/* Apply all the separated alignments specified in m_sepAligns.
 */
void BiComp::applyDunnartSepAligns(Canvas *canvas)
{
    foreach (SeparatedAlignment *sa, m_sepAligns) {
        ShapeObj *s = sa->shape1, *t=sa->shape2;
        CanvasItemList items;
        items.append(s); items.append(t);
        // Separation
        dtype dt = sa->af==Canvas::Vertical ? SEP_VERTICAL : SEP_HORIZONTAL;
        double minDist = (s->width()+t->width())/2.0;
        bool sort = true;
        createSeparation(NULL,dt,items,minDist,sort);
        // Alignment
        atypes at = sa->af==Canvas::Vertical ? ALIGN_CENTER : ALIGN_MIDDLE;
        createAlignment(at,items);
    }
}

/* Remove all the separated alignments specified in m_sepAligns.
 */
void BiComp::removeDunnartSepAligns(Canvas *canvas)
{
    // TODO
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
    //qDebug() << "Building alignment state table.";
    edge ed = NULL;
    forall_edges(ed,*m_graph) {
        node src = ed->source();
        node tgt = ed->target();
        int srcIndex = nodeIndices.value(src);
        int tgtIndex = nodeIndices.value(tgt);

        //qDebug() << QString("endpts: %1,%2").arg(srcIndex).arg(tgtIndex);

        alignmentState(srcIndex,tgtIndex) = Canvas::Connected;
        alignmentState(tgtIndex,srcIndex) = Canvas::Connected;
    }

    //-----------------------------------
    // Debug: show alignment state table.
    /*
    for (uint i = 0; i < N; i++) {
        QString row = "";
        for (uint j = 0; j < N; j++) {
            row += QString("%1 ").arg(alignmentState(i,j),2);
        }
        qDebug() << row;
    }
    */
    //-----------------------------------

    return alignmentState;
}

SeparatedAlignment *BiComp::chooseSA(vpsc::Rectangles rs, Matrix2d<int> &alignmentState,
                                     QMap<node, int> nodeIndices)
{
    SeparatedAlignment *sa = NULL;
    double minDeflection = 1.0;
    edge ed = NULL;
    forall_edges(ed,*m_graph) {
        node src = ed->source();
        node tgt = ed->target();
        int srcIndex = nodeIndices.value(src);
        int tgtIndex = nodeIndices.value(tgt);

        //qDebug() << QString("chooseSA considering endpts: %1,%2").arg(srcIndex).arg(tgtIndex);

        // If already aligned, skip this edge.
        int astate = alignmentState(srcIndex,tgtIndex);
        if (astate & (Canvas::Horizontal|Canvas::Vertical)) {
            //qDebug() << "    already aligned -- skip";
            continue;
        }
        // Consider horizontal alignment.
        if (!createsCoincidence(srcIndex,tgtIndex,Canvas::Horizontal,rs,alignmentState)) {
            double dH = deflection(srcIndex,tgtIndex,Canvas::Horizontal,rs);
            //qDebug() << QString("    deflection: %1").arg(dH);
            if (dH < minDeflection) {
                minDeflection = dH;
                if (!sa) sa = new SeparatedAlignment;
                sa->af = Canvas::Horizontal;
                sa->rect1 = srcIndex;
                sa->rect2 = tgtIndex;
            }
        }
        // Consider vertical alignment.
        if (!createsCoincidence(srcIndex,tgtIndex,Canvas::Vertical,rs,alignmentState)) {
            double dV = deflection(srcIndex,tgtIndex,Canvas::Vertical,rs);
            //qDebug() << QString("    deflection: %1").arg(dV);
            if (dV < minDeflection) {
                minDeflection = dV;
                if (!sa) sa = new SeparatedAlignment;
                sa->af = Canvas::Vertical;
                sa->rect1 = srcIndex;
                sa->rect2 = tgtIndex;
            }
        }
    }
    // Did we find an alignment?
    if (sa) {
        // If so, then complete the SeparatedAlignment object.
        int srcIndex = sa->rect1;
        int tgtIndex = sa->rect2;

        /*
        QString state = QString("Nodes %1,%2 had state %3.").arg(srcIndex).arg(tgtIndex)
                .arg(alignmentState(srcIndex,tgtIndex));
        qDebug() << state;
        */

        vpsc::Rectangle *rsrc = rs.at(srcIndex);
        vpsc::Rectangle *rtgt = rs.at(tgtIndex);
        // Separation Constraint
        vpsc::Dim sepDim = sa->af == Canvas::Horizontal ? vpsc::XDIM : vpsc::YDIM;
        vpsc::Dim algnDim = sa->af == Canvas::Horizontal ? vpsc::YDIM : vpsc::XDIM;
        double sep = sa->af == Canvas::Horizontal ?
                    (rsrc->width()+rtgt->width())/2.0 : (rsrc->height()+rtgt->height())/2.0;
        int l = sa->af == Canvas::Horizontal ?
                    (rsrc->getCentreX() < rtgt->getCentreX() ? srcIndex : tgtIndex) :
                    (rsrc->getCentreY() < rtgt->getCentreY() ? srcIndex : tgtIndex);
        int r = l == srcIndex ? tgtIndex : srcIndex;
        sa->separation = new cola::SeparationConstraint(sepDim,l,r,sep);
        // Alignment Constraint
        sa->alignment = new cola::AlignmentConstraint(algnDim);
        sa->alignment->addShape(l,0);
        sa->alignment->addShape(r,0);
        // Store the shapes too, for use later in applying a Dunnart constraint.
        node src = nodeIndices.key(srcIndex);
        node tgt = nodeIndices.key(tgtIndex);
        sa->shape1 = m_ownShapeMap.value(src);
        sa->shape2 = m_ownShapeMap.value(tgt);
    }
    return sa;
}

/* Given that the SeparatedAlignment sa has just been applied, update the
 * alignmentState matrix to reflect the new alignments that now exist.
 */
void BiComp::updateAlignmentState(SeparatedAlignment *sa, Matrix2d<int> &alignmentState)
{
    // Which dimension?
    Canvas::AlignmentFlags af = sa->af;
    // Get the indices of the two rectangles r1 and r2.
    int i1 = sa->rect1, i2 = sa->rect2;
    // Get the sets of indices of nodes already aligned with r1 and r2.
    QList<int> A1, A2;
    for (int j = 0; j < alignmentState.cols; j++) {
        if (alignmentState(i1,j) & af) A1.append(j);
        if (alignmentState(i2,j) & af) A2.append(j);
    }
    // r1 and r2 are aligned deliberately.
    alignmentState(i1,i2) |= af | Canvas::Deliberate;
    alignmentState(i2,i1) |= af | Canvas::Deliberate;
    // r1 is now aligned with each element of A2, and vice versa.
    foreach (int k, A2) {
        alignmentState(i1,k) |= af;
        alignmentState(k,i1) |= af;
    }
    foreach (int k, A1) {
        alignmentState(i2,k) |= af;
        alignmentState(k,i2) |= af;
    }
}

/* Say whether the proposed alignment would create an edge coincidence.
 */
bool BiComp::createsCoincidence(int srcIndex, int tgtIndex, Canvas::AlignmentFlags af,
                                vpsc::Rectangles rs, Matrix2d<int> &alignmentState)
{
    vpsc::Rectangle *r1 = rs.at(srcIndex), *r2 = rs.at(tgtIndex);
    double z1 = af == Canvas::Horizontal ? r1->getCentreX() : r1->getCentreY();
    double z2 = af == Canvas::Horizontal ? r2->getCentreX() : r2->getCentreY();
    double lowCoord = z1<z2 ? z1 : z2;
    double highCoord = z1<z2 ? z2 : z1;
    int lowIndex = z1<z2 ? srcIndex : tgtIndex;
    int highIndex = z1<z2 ? tgtIndex : srcIndex;
    // Let L and H be the low and high shapes respectively.
    // We consider each node U which is already aligned with either L or H.
    // Any such node must have lower coord than L if it is connected to L, and
    // higher coord than H if it is connected to H. If either of those conditions
    // fails, then we predict coincidence.
    bool coincidence = false;
    for (int j = 0; j < alignmentState.cols; j++) {
        if (j==lowIndex || j==highIndex) continue;
        vpsc::Rectangle *r = rs.at(j);
        int lj = alignmentState(lowIndex, j);
        int hj = alignmentState(highIndex, j);
        if (lj&af || hj&af) {
            double z = af==Canvas::Horizontal ? r->getCentreX() : r->getCentreY();
            // low shape
            if ( lj&Canvas::Connected && lowCoord < z ) {
                coincidence = true; break;
            }
            // high shape
            if ( hj&Canvas::Connected && z < highCoord ) {
                coincidence = true; break;
            }
        }
    }
    return coincidence;
}

/* Compute a score in [0.0, 1.0] measuring how far the "edge" in question E is
 * deflected from horizontal or vertical, depending on the passed alignment flag.
 * The "edge" E is the straight line from the centre of the src rectangle to that
 * of the tgt rectangle, as given by the src and tgt indices, and the vector of
 * Rectangles.
 *
 * If t is the angle that E makes with the positive x-axis, then the score we return
 * is sin^2(t) if af == Horizontal, and cos^2(t) if af == Vertical.
 *
 * So smaller deflection scores mean edges that are closer to axis-aligned.
 */
double BiComp::deflection(int srcIndex, int tgtIndex, Canvas::AlignmentFlags af, vpsc::Rectangles rs)
{
    vpsc::Rectangle *s = rs.at(srcIndex), *t = rs.at(tgtIndex);
    double sx=s->getCentreX(), sy=s->getCentreY(), tx=t->getCentreX(), ty=t->getCentreY();
    double a;
    if (ty>sy) {
        a = atan2(ty-sy,tx-sx);
    } else if (ty<sy) {
        a = atan2(sy-ty,sx-tx);
    } else {
        a = 0;
    }
    double sn = sin(a);
    double sn2 = sn*sn;
    double dfl = af==Canvas::Horizontal ? sn2 : 1 - sn2;
    return dfl;
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
    int magicNumber = 0;
    int fixedSep = magicNumber;
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
        fixedSep += magicNumber;
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

    // Apply any separated alignments created if we used ACA.
    applyDunnartSepAligns(canvas);
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
void BiComp::dfs(QMap<node, BiComp *> endpts, QList<BiComp *> &elements)
{
    elements.append(this);
    foreach (node m, m_cutNodes) {
        node n = m_nodemap.value(m);
        QList<BiComp*> next = endpts.values(n);
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
    // Below we'll use the following letters for nodes:
    // t: node in this BiComp T
    // o: node in other BiComp O
    // f: node in fusion BiComp F
    // g: node in original graph G
    // Begin by giving F a copy of each "normal node" in T.
    foreach (node t, this->m_normalNodes) {
        node f = fusion->m_graph->newNode();
        fusion->m_normalNodes.append(f);
        node g = this->m_nodemap.value(t);
        fusion->m_nodemap.insert(f,g);
    }
    // Now give F a copy of each "normal node" in O.
    foreach (node o, other->m_normalNodes) {
        node f = fusion->m_graph->newNode();
        fusion->m_normalNodes.append(f);
        node g = other->m_nodemap.value(o);
        fusion->m_nodemap.insert(f,g);
    }
    // Now the cutnodes of T.
    foreach (node t, this->m_cutNodes) {
        node f = fusion->m_graph->newNode();
        fusion->m_cutNodes.append(f);
        node g = this->m_nodemap.value(t);
        fusion->m_nodemap.insert(f,g);
    }
    // And the cutnodes of O. This time check for duplicates
    // with T.
    foreach (node o, other->m_cutNodes) {
        node g = other->m_nodemap.value(o);
        if (this->m_nodemap.values().contains(g)) continue;
        node f = fusion->m_graph->newNode();
        fusion->m_cutNodes.append(f);
        fusion->m_nodemap.insert(f,g);
    }
    // Now the edges.
    // First the edges of G corresponding to those in T.
    foreach (edge e, this->m_edgemap.values()) {
        node f1 = fusion->m_nodemap.key(e->source());
        node f2 = fusion->m_nodemap.key(e->target());
        edge f = fusion->m_graph->newEdge(f1,f2);
        fusion->m_edgemap.insert(f,e);
    }
    // And finally the edges of G corresp. those in O.
    foreach (edge e, other->m_edgemap.values()) {
        node f1 = fusion->m_nodemap.key(e->source());
        node f2 = fusion->m_nodemap.key(e->target());
        edge f = fusion->m_graph->newEdge(f1,f2);
        fusion->m_edgemap.insert(f,e);
    }
    return fusion;
}

// ----------------------------------------------------------------------------
// ACALayout ------------------------------------------------------------------

/** Build an ACA layout object for some nodes on the canvas.
  * If selection = false then apply to all the nodes on the canvas,
  * otherwise just the nodes in the current selection.
  */
ACALayout::ACALayout(Canvas *canvas, bool selection)
{
    // TODO
}

/** Build an ACA layout object for the graph defined by the passed shapes and connectors.
  */
ACALayout::ACALayout(QList<ShapeObj *> shapes, QList<Connector *> connectors)
{
    // TODO
}

/** Build an ACA layout object for the passed OGDF graph.
  */
ACALayout::ACALayout(Graph G, GraphAttributes GA) :
    idealLength(100)
{
    node n = NULL;
    int i = 0;
    forall_nodes(n,G) {
        double x = GA.x(n), y = GA.y(n);
        double X = x + GA.width(n), Y = y + GA.height(n);
        rs.push_back( new vpsc::Rectangle(x,X,y,Y) );
        m_ogdfNodeIndices.insert(n,i);
        i++;
    }
    edge e = NULL;
    forall_edges(e,G) {
        node src = e->source();
        node tgt = e->target();
        int srcIndex = m_ogdfNodeIndices.value(src);
        int tgtIndex = m_ogdfNodeIndices.value(tgt);
        es.push_back( cola::Edge(srcIndex, tgtIndex) );
    }
}

/** You must call this after running the layout in order to have the
  * resulting positions put into the GraphAttributes object.
  */
void ACALayout::readLayout(Graph G, GraphAttributes &GA)
{
    node n = NULL;
    forall_nodes(n,G) {
        vpsc::Rectangle *r = rs.at(m_ogdfNodeIndices.value(n));
        GA.x(n) = r->getMinX();
        GA.y(n) = r->getMinY();
    }
}

void ACALayout::setIdealLength(double il)
{
    idealLength = il;
}

/** Run the layout algorithm.
  */
void ACALayout::run(void)
{
    // Do an initial unconstrained FD layout.
    cola::ConstrainedFDLayout *fdlayout =
            new cola::ConstrainedFDLayout(rs,es,idealLength,false);
    fdlayout->run(true,true);
    // Prepare the alignment state matrix.
    initAlignmentState();
    // Start main loop.
    cola::CompoundConstraints ccs;
    ACASeparatedAlignment *sa = chooseSA();
    while (sa) {
        //debugOutput(sa);
        // Add the new separated alignment constraints.
        ccs.push_back(sa->separation);
        ccs.push_back(sa->alignment);
        // Redo the layout, with the new constraints.
        delete fdlayout;
        fdlayout = new cola::ConstrainedFDLayout(rs,es,idealLength,false);
        fdlayout->setConstraints(ccs);
        fdlayout->run(true,true);
        // Update state.
        updateAlignmentState(sa);
        // Store SA and choose next one.
        sepAligns.append(sa);
        sa = chooseSA();
    }
}

void ACALayout::debugOutput(ACASeparatedAlignment *sa)
{
    QString rects = "";
    for (int i = 0; i < rs.size(); i++) {
        vpsc::Rectangle *r = rs.at(i);
        double x = r->getCentreX();
        double y = r->getCentreY();
        rects += QString("(%1,%2) ").arg(x,5,'f',2).arg(y,5,'f',2);
    }
    qDebug() << rects;
    QString sepalgn = "";
    sepalgn += QString("%1").arg(sa->rect1,2);
    sepalgn += sa->af == ACAHORIZ ? "-" : "|";
    sepalgn += QString("%1").arg(sa->rect2,2);
    qDebug() << sepalgn;
}

void ACALayout::initAlignmentState(void)
{
    int N = rs.size();
    alignmentState = Matrix2d<int>(N,N);
    // Initialize with zeros.
    for (uint i = 0; i < N; i++) {
        for (uint j = 0; j < N; j++) {
            alignmentState(i,j) = 0;
        }
    }
    // Note connections.
    foreach (cola::Edge e, es) {
        int src = e.first, tgt = e.second;
        alignmentState(src,tgt) = ACACONN;
        alignmentState(tgt,src) = ACACONN;
    }
}

void ACALayout::updateAlignmentState(ACASeparatedAlignment *sa)
{
    // Which dimension?
    ACAFlags af = sa->af;
    // Get the indices of the two rectangles r1 and r2.
    int i1 = sa->rect1, i2 = sa->rect2;
    // Get the sets of indices of nodes already aligned with r1 and r2.
    QList<int> A1, A2;
    for (int j = 0; j < alignmentState.cols; j++) {
        if (alignmentState(i1,j) & af) A1.append(j);
        if (alignmentState(i2,j) & af) A2.append(j);
    }
    // r1 and r2 are aligned deliberately.
    alignmentState(i1,i2) |= af | ACADELIB;
    alignmentState(i2,i1) |= af | ACADELIB;
    // r1 is now aligned with each element of A2, and vice versa.
    foreach (int k, A2) {
        alignmentState(i1,k) |= af;
        alignmentState(k,i1) |= af;
    }
    foreach (int k, A1) {
        alignmentState(i2,k) |= af;
        alignmentState(k,i2) |= af;
    }
}

ACASeparatedAlignment *ACALayout::chooseSA(void)
{
    ACASeparatedAlignment *sa = NULL;
    double minDeflection = 1.0;
    // Consisder each edge for potential alignment.
    foreach (cola::Edge e, es) {
        int src = e.first, tgt = e.second;
        // If already aligned, skip this edge.
        int astate = alignmentState(src,tgt);
        if (astate & (ACAHORIZ|ACAVERT)) {
            //qDebug() << "    already aligned -- skip";
            continue;
        }
        // Consider horizontal alignment.
        if (!createsCoincidence(src,tgt,ACAHORIZ)) {
            double dH = deflection(src,tgt,ACAHORIZ);
            //qDebug() << QString("    deflection: %1").arg(dH);
            if (dH < minDeflection) {
                minDeflection = dH;
                if (!sa) sa = new ACASeparatedAlignment;
                sa->af = ACAHORIZ;
                sa->rect1 = src;
                sa->rect2 = tgt;
            }
        }
        // Consider vertical alignment.
        if (!createsCoincidence(src,tgt,ACAVERT)) {
            double dV = deflection(src,tgt,ACAVERT);
            //qDebug() << QString("    deflection: %1").arg(dV);
            if (dV < minDeflection) {
                minDeflection = dV;
                if (!sa) sa = new ACASeparatedAlignment;
                sa->af = ACAVERT;
                sa->rect1 = src;
                sa->rect2 = tgt;
            }
        }
    }
    // Did we find an alignment?
    if (sa) {
        // If so, then complete the ACASeparatedAlignment object.
        int src = sa->rect1;
        int tgt = sa->rect2;
        //qDebug() << QString("Nodes %1,%2 had state %3.").arg(src).arg(tgt).arg(alignmentState(src,tgt));
        vpsc::Rectangle *rsrc = rs.at(src);
        vpsc::Rectangle *rtgt = rs.at(tgt);
        // Separation Constraint
        vpsc::Dim sepDim = sa->af == ACAHORIZ ? vpsc::XDIM : vpsc::YDIM;
        vpsc::Dim algnDim = sa->af == ACAHORIZ ? vpsc::YDIM : vpsc::XDIM;
        double sep = sa->af == ACAHORIZ ?
                    (rsrc->width()+rtgt->width())/2.0 : (rsrc->height()+rtgt->height())/2.0;
        int l = sa->af == Canvas::Horizontal ?
                    (rsrc->getCentreX() < rtgt->getCentreX() ? src : tgt) :
                    (rsrc->getCentreY() < rtgt->getCentreY() ? src : tgt);
        int r = l == src ? tgt : src;
        sa->separation = new cola::SeparationConstraint(sepDim,l,r,sep);
        // Alignment Constraint
        sa->alignment = new cola::AlignmentConstraint(algnDim);
        sa->alignment->addShape(l,0);
        sa->alignment->addShape(r,0);
    }
    return sa;
}

/* Say whether the proposed alignment would create an edge coincidence.
 */
bool ACALayout::createsCoincidence(int src, int tgt, ACAFlags af)
{
    vpsc::Rectangle *r1 = rs.at(src), *r2 = rs.at(tgt);
    double z1 = af == ACAHORIZ ? r1->getCentreX() : r1->getCentreY();
    double z2 = af == ACAHORIZ ? r2->getCentreX() : r2->getCentreY();
    double lowCoord = z1<z2 ? z1 : z2;
    double highCoord = z1<z2 ? z2 : z1;
    int lowIndex = z1<z2 ? src : tgt;
    int highIndex = z1<z2 ? tgt : src;
    // Let L and H be the low and high shapes respectively.
    // We consider each node U which is already aligned with either L or H.
    // Any such node must have lower coord than L if it is connected to L, and
    // higher coord than H if it is connected to H. If either of those conditions
    // fails, then we predict coincidence.
    bool coincidence = false;
    for (int j = 0; j < alignmentState.cols; j++) {
        if (j==lowIndex || j==highIndex) continue;
        vpsc::Rectangle *r = rs.at(j);
        int lj = alignmentState(lowIndex, j);
        int hj = alignmentState(highIndex, j);
        if (lj&af || hj&af) {
            double z = af==ACAHORIZ ? r->getCentreX() : r->getCentreY();
            // low shape
            if ( lj&ACACONN && lowCoord < z ) {
                coincidence = true; break;
            }
            // high shape
            if ( hj&ACACONN && z < highCoord ) {
                coincidence = true; break;
            }
        }
    }
    return coincidence;
}

/* Compute a score in [0.0, 1.0] measuring how far the "edge" in question E is
 * deflected from horizontal or vertical, depending on the passed alignment flag.
 * The "edge" E is the straight line from the centre of the src rectangle to that
 * of the tgt rectangle, as given by the src and tgt indices, and the vector of
 * Rectangles.
 *
 * If t is the angle that E makes with the positive x-axis, then the score we return
 * is sin^2(t) if af == Horizontal, and cos^2(t) if af == Vertical.
 *
 * So smaller deflection scores mean edges that are closer to axis-aligned.
 */
double ACALayout::deflection(int src, int tgt, ACAFlags af)
{
    vpsc::Rectangle *s = rs.at(src), *t = rs.at(tgt);
    double sx=s->getCentreX(), sy=s->getCentreY(), tx=t->getCentreX(), ty=t->getCentreY();
    double a;
    if (ty>sy) {
        a = atan2(ty-sy,tx-sx);
    } else if (ty<sy) {
        a = atan2(sy-ty,sx-tx);
    } else {
        a = 0;
    }
    double sn = sin(a);
    double sn2 = sn*sn;
    double dfl = af==ACAHORIZ ? sn2 : 1 - sn2;
    return dfl;
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
            foreach (Chunk *ch2, nbrList) {
                ch2->setParentCutNode(cn);
                BiComp *bc = dynamic_cast<BiComp*>(ch);
                ch2->setParent(bc);
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
#define showIDs
#ifdef showIDs
            QString label = QString("%1").arg(shape->internalId());
            shape->setLabel(label);
#endif
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

/* Fuse BCs that share cutnodes. Return list of BCs thus obtained.
 */
QList<BiComp*> BCLayout::fuseBCs(QList<BiComp *> bicomps)
{
    // Prepare multimap from cutnodes to BCs they lie in.
    QMap<node,BiComp*> endpts;
    foreach (BiComp *bc, bicomps) {
        QList<node> cn = bc->getCutNodes();
        foreach (node nd, cn) {
            endpts.insertMulti(nd,bc);
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
    qDebug() << QString("Found %1 compound BCs.").arg(cbc.size());
    return cbc;
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
    //Graph *Gp = new Graph;
    //Graph G = *Gp;
    Graph G;
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
    // Fuse BCs that share cutnodes.
    bicomps = fuseBCs(bicomps);

    // Get a new graph isomorphic to the result of removing the BC edges from G.
    QMap<node,node> nodeMapG2ToG;
    Graph *G2 = removeBiComps(G, bicomps, nodeMapG2ToG);
    assert(G2->consistencyCheck());

    // Get connected components of original graph G corresponding to those of G2.
    QMap<int,node> ccomps = getConnComps2(G2, nodeMapG2ToG);

    // Form trees on those components, throwing away isolated
    // nodes (which must have been the noncutnodes belonging to nontrivial BCs).
    QList<RootedTree*> rtrees;
    foreach (int i, ccomps.keys().toSet())
    {
        QList<node> nodes = ccomps.values(i);
        if (nodes.size() < 2) continue;
        RootedTree *rtree = new RootedTree(nodes, cutnodes);
        rtrees.append(rtree);
    }

    /* TODO
     * ====
     *
     * Change the algorithm starting here, according to
     * our latest ideas on how reassembly should work.
     *
     * For now all we have is a BFS turning all the Chunks into a tree,
     * in fact a rooted tree where we have chosen a largest Chunk as the
     * root.
     *
     */

    // The following code assumes there is at least one BC.
    // So this does not handle the case in which the graph itself is at tree!

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

/* Latest attempt (28 Nov 2013) at the Orthowontist layout algorithm.
  * Operates on all the nodes and edges currently on the Dunnart canvas.
  * Expects this to be a connected graph.
  */
void BCLayout::ortholayout2(void)
{
    BiComp::method = 0;
    shapemap nodeShapes;
    connmap edgeConns;
    Graph G;
    ogdfGraph(G, nodeShapes, edgeConns);

    // 1. Compute external trees.
    QList<ExternalTree*> XX = removeExternalTrees(G, nodeShapes);
    // Test:
    /*
    qDebug() << "External Trees:";
    foreach (ExternalTree *X, XX) {
        qDebug() << X->listNodes();
    }
    */

    // 2. Get nontrivial biconnected components (size >= 3), and get the set of cutnodes in G.
    QSet<node> cutnodes;
    QList<BiComp*> BB = getNontrivialBCs(G, cutnodes);
    // Fuse BCs that share cutnodes.
    BB = fuseBCs(BB);

    // 3. Compute internal trees.
    // Get a new graph isomorphic to the result of removing the BC edges from G.
    QMap<node,node> nodeMapG2ToG;
    Graph *G2 = removeBiComps(G, BB, nodeMapG2ToG);
    // Get connected components of original graph G corresponding to those of G2.
    QMap<int,node> ccomps = getConnComps2(G2, nodeMapG2ToG);
    // Form trees on those components, throwing away isolated
    // nodes (which must have been the noncutnodes belonging to nontrivial BCs).
    QList<InternalTree*> II;
    foreach (int i, ccomps.keys().toSet())
    {
        QList<node> nodes = ccomps.values(i);
        if (nodes.size() < 2) continue;
        InternalTree *I = new InternalTree(nodes, cutnodes);
        II.append(I);
    }

    // 4. Lay out each B in BB by FD+ACA.
    foreach (BiComp *B, BB) {
        B->acaLayout();
    }

    // 5. Lay out each X in XX by OGDF TreeLayout.
    foreach (ExternalTree *X, XX) {
        X->treeLayout();
    }

    // 6. Build the metagraph M, reading Bbar and Xbar sizes from layouts
    //    of corresponding B and X.
    Graph *M = new Graph();
    GraphAttributes *MA = new GraphAttributes(*M);
    buildMetagraph(*M,*MA,XX,BB,II,cutnodes);

    // 7. Lay out M using FD+ACA.
    ACALayout aca = ACALayout(*M, *MA);
    aca.run();
    aca.readLayout(*M, *MA);

    // 8. Place each B where Bbar lies, orienting to minimize stress.
    // TODO

    // 9. Place each X where Xbar lies, orienting according to alignment
    //    of the edge ( t(X), r(X) ).
    // TODO

    // 10. Draw it.
    // TODO

}

Graph *BCLayout::buildMetagraph(Graph &M, GraphAttributes &MA,
                     QList<ExternalTree *> XX, QList<BiComp *> BB, QList<InternalTree *> II,
                     QSet<node> cutnodes)
{
    // TODO
    return NULL;
}

/* Expects a connected graph.
  * Prunes all external trees, and returns them in a list.
  * The graph itself IS altered.
  * Uses the idea from the TopoLayout paper.
  */
QList<ExternalTree*> BCLayout::removeExternalTrees(Graph &G, shapemap nodeShapes)
{
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
    // For debugging purposes, let's keep track of the id numbers of the Dunnart ShapeObj
    // objects to which the nodes in H correspond.
    QMap<node, int> dunnartIDs; // domain is H
    // For each node r of H, we will keep track of the node t in G which was its parent.
    // We imagine r as a "root" and t as a "taproot". If subsequently t gets removed from
    // G (and t' is its copy added to H), then we delete the mapping r-->t from the map.
    // This way the domain of the map always is precisely the set of root nodes in H.
    // (H will be such that each connected component is one of the external trees, and each
    // one will have a root.)
    QMap< node,node > rootsToTaproots; // Maps from H to G!
    // We also keep track of the nodes and edges belonging to each tree in H.
    QMap< node,QList<node> > rootsToNodes;
    QMap< node,QList<edge> > rootsToEdges;
    // The basic idea is simple: Continue deleting nodes of degree 1, until
    // there aren't any more.
    while (nodesByDegree.value(1).size() > 0) {
        // Grab the degree-1 nodes, and replace with an empty set in the nodesByDegree map.
        QSet<node> degreeOneNodes = nodesByDegree.value(1);
        QSet<node> set;
        nodesByDegree.insert(1,set);
        foreach (node r, degreeOneNodes) {
            // Create a copy in H.
            node rH = H.newNode();
            // Note ID of corresp. Dunnart shape obj.
            ShapeObj *sh = nodeShapes.value(r);
            dunnartIDs.insert(rH,sh->internalId());
            // Prepare lists of the nodes and edges for the tree of which rH is the new root.
            // Note: the root IS included in the list of nodes belonging to the tree.
            QList<node> treeNodes;
            treeNodes.append(rH);
            QList<edge> treeEdges;
            // Was r the taproot for any nodes already in H?
            // If so, it should now be connected to them.
            QList<node> children = rootsToTaproots.keys(r);
            foreach (node c, children) {
                edge f = H.newEdge(rH,c);
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
                t = r==e->source() ? e->target() : e->source();
                rootsToTaproots.insert(rH,t);
            }
            // Now we can delete r.
            G.delNode(r);
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
    /*
    // Compute the connected components of H, and node the root node of each.
    QMap<int,node> CC;
    QMap<int,node> root;
    NodeArray<int> cc(H);
    ogdf::connectedComponents(H,cc);
    node a = NULL;
    forall_nodes(a,H) {
        int n = cc[a];
        CC.insertMulti(n,a);
        if (rootsToTaproots.keys().contains(a)) root.insert(n,a);
    }
    */
    // Now construct an ExternalTree object for each connected component in H,
    // and return the list of these.
    QList<ExternalTree*> XX;
    /*
    foreach (int n, CC.keys()) {
        QList<node> nodes = CC.values(n);   // list of nodes in H
        node r = root.value(n);             // node in H
        node t = rootsToTaproots.value(r);  // node in G
        ExternalTree *X = new ExternalTree(nodes,r,t);
        XX.append(X);
    }
    */
    foreach (node r, rootsToTaproots.keys()) {
        node t = rootsToTaproots.value(r);
        QList<node> nodes = rootsToNodes.value(r);
        QList<edge> edges = rootsToEdges.value(r);
        ShapeObj *tapRootShape = nodeShapes.value(t);
        ExternalTree *X = new ExternalTree(nodes,edges,r,t,dunnartIDs,tapRootShape->internalId());
        XX.append(X);
    }
    return XX;
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
