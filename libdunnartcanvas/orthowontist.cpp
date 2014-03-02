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
#include <QtAlgorithms>

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
#include "libogdf/ogdf/basic/EdgeComparer.h"
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

//#define CGAL
#ifdef CGAL
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Polygon_2.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Polygon_2<K> Polygon_2;
#endif

using namespace ogdf;

namespace ow {


// ------------------------------------------------------------------
// BiComp -----------------------------------------------------------

BiComp::BiComp(void) :
    m_stubNodeShapesHaveBeenAddedToCanvas(false),
    m_dummyNodeShapesHaveBeenAddedToCanvas(false),
    m_idealLength(100),
    m_nodePadding(50)
{}

BiComp::BiComp(QList<node> nodes, QList<edge> edges, QList<node> cutnodes,
               shapemap nodeShapes, connmap edgeConns) :
    m_stubNodeShapesHaveBeenAddedToCanvas(false),
    m_dummyNodeShapesHaveBeenAddedToCanvas(false),
    m_idealLength(100),
    m_nodePadding(50)
{
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

void BiComp::numberShapes(void) {
    foreach (node m, m_dunnartShapes.keys()) {
        ShapeObj *sh = m_dunnartShapes.value(m);
        sh->setLabel(QString::number(sh->internalId()));
    }
}

ShapeObj *BiComp::getShape(node m) {
    return m_dunnartShapes.value(m);
}

QList<ShapeObj*> BiComp::allShapes(void) {
    return m_dunnartShapes.values();
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

void BiComp::noteRoot(ExternalTree *E) {
    // Get the Dunnart shape that is the root of the tree.
    ShapeObj *Eshape = E->rootShape();
    // Get own node representing that shape.
    node ownRootNode = m_dunnartShapes.key(Eshape);
    // Note that this node is a root.
    m2_rootsToTrees.insert(ownRootNode,E);
}

void BiComp::addStubNodeForTree(ExternalTree *E, QSizeF size) {
    // Get the Dunnart shape that is the root of the tree.
    ShapeObj *Eshape = E->rootShape();
    // Get own node representing that shape.
    node ownRootNode = m_dunnartShapes.key(Eshape);
    // Create a stub node to represent the tree.
    node stub = m_graph->newNode();
    // Connect it to ownRootNode.
    edge e = m_graph->newEdge(ownRootNode,stub);
    // Register stub edge.
    m_stubedges.append(e);
    // And a shape for the stub node.
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    ShapeObj *stubShape = factory->createShape("org.dunnart.shapes.rect");
    stubShape->setFillColour(QColor(0,192,0));
    //stubShape->setLabel(QString::number(stubShape->internalId()));
    // Map the stub node to the shape.
    m_dunnartShapes.insert(stub,stubShape);
    // Map stub to tree.
    m_stubnodesToTrees.insert(stub,E);
    // Set size of node and shape
    m_ga->width(stub) = size.width();
    m_ga->height(stub) = size.height();
    stubShape->setSize(size);
    // Set position to be at a displacement from root of length
    // equal to the avg dimension of the passed 'size', and
    // random angle.
    //
    // TODO (maybe):
    // Later we might choose a good face for the tree, instead of
    // using a random angle.
    double d = (size.width()+size.height())/2;
    double t = (double)(qrand());
    double s = sin(t), c = cos(t);
    double dx = d*c, dy = d*s;
    double x0 = m_ga->x(ownRootNode), y0 = m_ga->y(ownRootNode);
    double x1 = x0 + dx, y1 = y0 + dy;
    m_ga->x(stub) = x1;
    m_ga->y(stub) = y1;
    stubShape->setCentrePos(QPointF(x1,y1));
}

void BiComp::layout(void) {
    bool debug = true;

    // 1. Build and run ACA Layout.
    ACALayout *aca = new ACALayout(*m_graph, *m_ga);

    QMap<node,int> nodeIndices = aca->nodeIndices();
    std::vector<cola::Edge> colaEdges = aca->colaEdges();

    aca->debugName("BC");

    //aca->idealLength(m_idealLength);
    //aca->nodePadding(m_nodePadding);

    aca->idealLength(1);
    aca->edgeLengths(edgeLengths(nodeIndices,colaEdges));
    aca->nodePadding(nodePadding(nodeIndices));

    aca->run();
    aca->readPositions(*m_graph, *m_ga);

    // 1.5. Build planarization.
    //m_planarization = new Planarization(*m_graph, *m_ga,
    //                                          aca->alignments(*m_graph), m_dummyNodeSize);
    // ... TODO ...

    // 2. Lay out external trees.
    //    Set stubnode sizes according to bounding boxes of laid out trees.
    foreach (node stub, m_stubnodesToTrees.keys()) {
        ExternalTree *E = m_stubnodesToTrees.value(stub);
        ShapeObj *Eshape = E->rootShape();
        node root = m_dunnartShapes.key(Eshape);
        // Determine closest cardinal direction of line from tap to stub.
        double rx=m_ga->x(root),ry=m_ga->y(root);
        double sx=m_ga->x(stub),sy=m_ga->y(stub);
        double vx=sx-rx, vy=sy-ry;
        ogdf::Orientation orient = vx >= vy ?
                    (vx >= -vy ? leftToRight : topToBottom) :
                    (vx >= -vy ? bottomToTop : rightToLeft) ;
        // Lay out tree.
        E->orientation(orient);
        E->treeLayout();
        // Read size.
        QSizeF size = E->rootlessBBox().size();
        // If stub edge was aligned by ACA, offset that alignment if necessary.
        int si = nodeIndices.value(stub), ri = nodeIndices.value(root);
        if (aca->delibAligned(si,ri) && E->needsAlignmentOffset()) {
            double offset = E->alignmentOffset();
            assert(aca->offsetAlignment(ri,si,offset));
        }
        // Set size in stub node and corresponding shape object.
        m_ga->width(stub) = size.width();
        m_ga->height(stub) = size.height();
        ShapeObj *stubShape = m_dunnartShapes.value(stub);
        stubShape->setSize(size);
    }

    // 3. Build EdgeNodes, and generate sep-co's between them and the stubnodes.
    QMap<vpsc::Dim,EdgeNode> ens = aca->generateEdgeNodes();
    double padding = m_nodePadding;
    cola::CompoundConstraints hSepcos =
            generateStubEdgeSepCos(vpsc::HORIZONTAL, ens.values(vpsc::HORIZONTAL),
                                   nodeIndices, padding);
    cola::CompoundConstraints vSepcos =
            generateStubEdgeSepCos(vpsc::VERTICAL, ens.values(vpsc::VERTICAL),
                                   nodeIndices, padding);
    cola::CompoundConstraints sepcos;
    foreach (cola::CompoundConstraint *cc, hSepcos) sepcos.push_back(cc);
    foreach (cola::CompoundConstraint *cc, vSepcos) sepcos.push_back(cc);

    if (debug) {
        foreach (cola::CompoundConstraint *cc, sepcos) {
            cola::SeparationConstraint *sc = dynamic_cast<cola::SeparationConstraint*>(cc);
            QList<QString> names;
            int l = sc->left(), r = sc->right();
            node nl = nodeIndices.key(l), nr = nodeIndices.key(r);
            ShapeObj *shl = m_dunnartShapes.value(nl);
            ShapeObj *shr = m_dunnartShapes.value(nr);
            names.append((QString("%1").arg(shl->internalId())));
            names.append((QString("%1").arg(shr->internalId())));

            QString rel = sc->dimension() == vpsc::HORIZONTAL ? "left of" : "above";
            qDebug() << QString("sepco %1 %2 %3 by %4")
                        .arg(names.at(0))
                        .arg(rel)
                        .arg(names.at(1))
                        .arg(sc->gap);
        }
    }

    // 4. Run FD layout again.
    bool preventOverlaps = true;
    // Keep the ACA sep-cos.
    foreach (cola::CompoundConstraint *cc, aca->ccs()) sepcos.push_back(cc);
    postACACola(preventOverlaps, m_idealLength, nodeIndices, sepcos);

    // 5. Translate trees.
    translateTrees();

    m_planarization = new Planarization(*m_graph, *m_ga,
                                              aca->alignments(*m_graph), m_dummyNodeSize, m_dunnartShapes);
    // ...


}

void BiComp::layout2(void) {
    bool debug = true;

    // 1. Build and run ACA Layout.
    ACALayout *aca = new ACALayout(*m_graph, *m_ga);

    QMap<node,int> nodeIndices = aca->nodeIndices();
    std::vector<cola::Edge> colaEdges = aca->colaEdges();

    aca->debugName("BC");

    //aca->idealLength(m_idealLength);
    //aca->nodePadding(m_nodePadding);

    aca->idealLength(1);
    aca->edgeLengths(edgeLengths(nodeIndices,colaEdges));
    //aca->nodePadding(nodePadding(nodeIndices));

    aca->run();
    aca->readPositions(*m_graph, *m_ga);

    // 1.25. Give trees an initial layout, just so we know their sizes.
    QMap<node,QSizeF> treeSizes;
    foreach (node root, m2_rootsToTrees.keys()) {
        ExternalTree *E = m2_rootsToTrees.value(root);
        E->treeLayout();
        QSizeF size = E->rootlessBBox().size();
        treeSizes.insert(root,size);
    }

    // 1.5. Build planarization.
    m_planarization = new Planarization(*m_graph, *m_ga,
        aca->alignments(*m_graph), m_dummyNodeSize, m_dunnartShapes);
    m_planarization->filename = filename;

    m_planarization->removeOverlaps();

    m_planarization->setTreeSizes(treeSizes);
    m_planarization->defineRootNodes(m2_rootsToTrees.keys());
    m_planarization->idealLength(m_idealLength);
    //m_planarization->chooseFDTreeFaces();
    //m_planarization->chooseCombTreeFaces();
    m_planarization->chooseGreedyTreeFaces();

    // Now redo tree layouts, with the chosen orientations.
    foreach (node root, m2_rootsToTrees.keys()) {
        ogdf::Orientation ori = m_planarization->treeOrientation(root);
        ExternalTree *E = m2_rootsToTrees.value(root);
        E->orientation(ori);
        E->treeLayout();
        QSizeF size = E->rootlessBBox().size();
        treeSizes.insert(root,size);
    }
    // Inform of the new sizes.
    m_planarization->setTreeSizes(treeSizes);
    // Expand to make room for trees.
    m_planarization->expand(10);

    // Lay out with "neighbour stress", in order to evenly
    // distribute nodes.
    m_planarization->distribWithNbrStress();

    // Update node positions.
    m_planarization->translateNodes(*m_graph, *m_ga);

    // Translate trees.
    foreach (node root, m2_rootsToTrees.keys()) {
        ExternalTree *E = m2_rootsToTrees.value(root);
        m_planarization->translateTree(E,root);
    }


    /*
    // 2. Lay out external trees.
    //    Set stubnode sizes according to bounding boxes of laid out trees.
    foreach (node root, m2_rootsToTrees.keys()) {
        ExternalTree *E = m2_rootsToTrees.value(root);
        m_planarization->layoutTreeForRoot(E,root);
        //QPointF pos = m_planarization->origRootToStubPos.value(root);
    }
    */

    /*
    foreach (node stub, m_stubnodesToTrees.keys()) {
        ExternalTree *E = m_stubnodesToTrees.value(stub);
        ShapeObj *Eshape = E->rootShape();
        node root = m_dunnartShapes.key(Eshape);
        // Determine closest cardinal direction of line from tap to stub.
        double rx=m_ga->x(root),ry=m_ga->y(root);
        double sx=m_ga->x(stub),sy=m_ga->y(stub);
        double vx=sx-rx, vy=sy-ry;
        ogdf::Orientation orient = vx >= vy ?
                    (vx >= -vy ? leftToRight : topToBottom) :
                    (vx >= -vy ? bottomToTop : rightToLeft) ;
        // Lay out tree.
        E->orientation(orient);
        E->treeLayout();
        // Read size.
        QSizeF size = E->rootlessBBox().size();
        // If stub edge was aligned by ACA, offset that alignment if necessary.
        int si = nodeIndices.value(stub), ri = nodeIndices.value(root);
        if (aca->delibAligned(si,ri) && E->needsAlignmentOffset()) {
            double offset = E->alignmentOffset();
            assert(aca->offsetAlignment(ri,si,offset));
        }
        // Set size in stub node and corresponding shape object.
        m_ga->width(stub) = size.width();
        m_ga->height(stub) = size.height();
        ShapeObj *stubShape = m_dunnartShapes.value(stub);
        stubShape->setSize(size);
    }

    // 3. Build EdgeNodes, and generate sep-co's between them and the stubnodes.
    QMap<vpsc::Dim,EdgeNode> ens = aca->generateEdgeNodes();
    double padding = m_nodePadding;
    cola::CompoundConstraints hSepcos =
            generateStubEdgeSepCos(vpsc::HORIZONTAL, ens.values(vpsc::HORIZONTAL),
                                   nodeIndices, padding);
    cola::CompoundConstraints vSepcos =
            generateStubEdgeSepCos(vpsc::VERTICAL, ens.values(vpsc::VERTICAL),
                                   nodeIndices, padding);
    cola::CompoundConstraints sepcos;
    foreach (cola::CompoundConstraint *cc, hSepcos) sepcos.push_back(cc);
    foreach (cola::CompoundConstraint *cc, vSepcos) sepcos.push_back(cc);

    if (debug) {
        foreach (cola::CompoundConstraint *cc, sepcos) {
            cola::SeparationConstraint *sc = dynamic_cast<cola::SeparationConstraint*>(cc);
            QList<QString> names;
            int l = sc->left(), r = sc->right();
            node nl = nodeIndices.key(l), nr = nodeIndices.key(r);
            ShapeObj *shl = m_dunnartShapes.value(nl);
            ShapeObj *shr = m_dunnartShapes.value(nr);
            names.append((QString("%1").arg(shl->internalId())));
            names.append((QString("%1").arg(shr->internalId())));

            QString rel = sc->dimension() == vpsc::HORIZONTAL ? "left of" : "above";
            qDebug() << QString("sepco %1 %2 %3 by %4")
                        .arg(names.at(0))
                        .arg(rel)
                        .arg(names.at(1))
                        .arg(sc->gap);
        }
    }

    // 4. Run FD layout again.
    bool preventOverlaps = true;
    // Keep the ACA sep-cos.
    foreach (cola::CompoundConstraint *cc, aca->ccs()) sepcos.push_back(cc);
    postACACola(preventOverlaps, m_idealLength, nodeIndices, sepcos);
    */

    // 5. Translate trees.
    translateTrees();



    // ...

}

void BiComp::updateShapePositions(void) {
    node n;
    forall_nodes(n, *m_graph) {
        if (m_stubnodesToTrees.contains(n) && !m_stubNodeShapesHaveBeenAddedToCanvas) continue;
        ShapeObj *sh = m_dunnartShapes.value(n);
        double x = m_ga->x(n), y = m_ga->y(n);
        sh->setCentrePos(QPointF(x,y));
    }
}

void BiComp::orthogonalRouting(bool b) {
    dunnart::Connector::RoutingType t = b ? Connector::orthogonal : Connector::polyline;
    foreach (Connector *conn, m_dunnartConns.values()) {
        if (conn==NULL) continue;
        conn->setRoutingType(t);
    }
}

void BiComp::addStubNodeShapesToCanvas(Canvas *canvas) {
    //canvas->stop_graph_layout();
    foreach (node stub, m_stubnodesToTrees.keys()) {
        ShapeObj *sh = m_dunnartShapes.value(stub);
        canvas->addItem(sh);
    }
    //canvas->restart_graph_layout();
    m_stubNodeShapesHaveBeenAddedToCanvas = true;
}

void BiComp::addDummyNodeShapesToCanvas(Canvas *canvas) {
    m_planarization->addDummyNodeShapesToCanvas(canvas);
    m_dummyNodeShapesHaveBeenAddedToCanvas = true;
}

cola::CompoundConstraints BiComp::generateStubEdgeSepCos(
        vpsc::Dim dim, QList<EdgeNode> ens, QMap<node, int> nodeIndices, double gap)
{
    double pad = gap/2;
    cola::CompoundConstraints ccs;
    unsigned int priority = cola::PRIORITY_NONOVERLAP - 1; //copied from colafd.cpp. What does it do?
    cola::NonOverlapConstraints *noc =
            new cola::NonOverlapConstraints(priority);
    QList<node> sns = m_stubnodesToTrees.keys();
    int Ne = ens.size();
    int Ns = sns.size();
    // Variables and Rectangles.
    vpsc::Variables vs;
    vpsc::Rectangles rs;
    // Add edgenodes first.
    for (int i = 0; i < Ne; i++) {
        EdgeNode E = ens.at(i);
        noc->addShape(i, E.bbox.width()/2 + pad, E.bbox.height()/2 + pad);
        double p = dim==vpsc::HORIZONTAL ? E.bbox.center().x() : E.bbox.center().y();
        vpsc::Variable *v = new vpsc::Variable(i,p);
        vs.push_back(v);
        vpsc::Rectangle *r = new vpsc::Rectangle(E.bbox.left(),E.bbox.right(),E.bbox.top(),E.bbox.bottom());
        rs.push_back(r);
    }
    // Add stubnodes.
    for (int j = 0; j < Ns; j++) {
        node s = sns.at(j);
        double w=m_ga->width(s), h=m_ga->height(s);
        double cx=m_ga->x(s), cy=m_ga->y(s);
        double p = dim == vpsc::HORIZONTAL ? cx : cy;
        noc->addShape(Ne + j, w/2 + pad, h/2 + pad);
        vpsc::Variable *v = new vpsc::Variable(Ne + j,p);
        vs.push_back(v);
        vpsc::Rectangle *r = new vpsc::Rectangle(cx-w/2,cx+w/2,cy-h/2,cy+h/2);
        rs.push_back(r);
    }
    // Generate vpsc constraints.
    vpsc::Constraints cs;
    noc->generateSeparationConstraints(dim,vs,cs,rs);
    // Convert into cola constraints.
    foreach (vpsc::Constraint *c, cs) {
        int l = c->left->id, r = c->right->id;

        // Reject constraints that hold between two edgenodes or two stubnodes.
        if ( (l < Ne && r < Ne) || (l >= Ne && r >= Ne) ) continue;

        // Convert IDs to indices assigned by original ACA layout.
        l = l < Ne ? ens.at(l).srcIndex : nodeIndices.value(sns.at(l-Ne));
        r = r < Ne ? ens.at(r).srcIndex : nodeIndices.value(sns.at(r-Ne));

        cola::SeparationConstraint *s = new cola::SeparationConstraint(dim,l,r,c->gap);
        ccs.push_back(s);
    }
    // Clean up
    // FIXME: why does this cause segfaults?
    //foreach (vpsc::Variable *v, vs) delete v;
    return ccs;
}

void BiComp::translateTrees(void) {
    foreach (node stub, m_stubnodesToTrees.keys()) {
        // Translate all but root node of tree into stubnode.
        ExternalTree *E = m_stubnodesToTrees.value(stub);
        double cx=m_ga->x(stub), cy=m_ga->y(stub);
        double w=m_ga->width(stub), h=m_ga->height(stub);
        double x = cx-w/2, y = cy-h/2;
        QPointF p = E->rootlessBBox().topLeft();
        double dx=x-p.x(), dy=y-p.y();
        E->translate(QPointF(dx,dy));
        // Place root node of tree over own copy of root.
        ShapeObj *sh = E->rootShape();
        node root = m_dunnartShapes.key(sh);
        cx = m_ga->x(root);
        cy = m_ga->y(root);
        E->placeRootAt(QPointF(cx,cy));
    }
}

double *BiComp::edgeLengths(QMap<node, int> nodeIndices, std::vector<cola::Edge> colaEdges) {
    int m = colaEdges.size();
    double *eL = new double[m];
    int i = 0;
    foreach (cola::Edge e, colaEdges) {
        int srcIdx = e.first, tgtIdx = e.second;
        node src = nodeIndices.key(srcIdx), tgt = nodeIndices.key(tgtIdx);
        double ws=m_ga->width(src), hs=m_ga->height(src);
        double wt=m_ga->width(tgt), ht=m_ga->height(tgt);
        eL[i++] = (ws+hs+wt+ht)/2;
        // Skip the rest for now.
        continue;

        // If neither node is a stubnode, then use own ideal length as length.
        if (!m_stubnodesToTrees.contains(src) && !m_stubnodesToTrees.contains(tgt)) {
            eL[i++] = m_idealLength;
        } else {
            // Otherwise at least one is a stubnode.
            // Should be exactly one, except in a degenerate case for which,
            // for now, we throw an error.
            assert(!m_stubnodesToTrees.contains(src) || !m_stubnodesToTrees.contains(tgt));
            // Use their average dimension.
            double ws=m_ga->width(src), hs=m_ga->height(src);
            double wt=m_ga->width(tgt), ht=m_ga->height(tgt);
            eL[i++] = (ws+hs+wt+ht)/2;
        }
    }
    return eL;
}

QList<double> BiComp::nodePadding(QMap<node, int> nodeIndices) {
    QList<double> P;
    QList<node> nodes = nodeIndices.keys();
    int N = nodes.size();
    for (int i = 0; i < N; i++) {
        node n = nodeIndices.key(i);
        double p = m_stubnodesToTrees.contains(n) ? 0 : m_nodePadding;
        P.append(p);
    }
    return P;
}

void BiComp::postACACola(bool preventOverlaps, double idealLength,
                         QMap<node, int> nodeIndices, cola::CompoundConstraints sepcos) {

    vpsc::Rectangles rs;
    std::vector<cola::Edge> es;

    // Build rectangles and edges.
    int N = nodeIndices.values().size();
    for (int i = 0; i < N; i++) {
        node n = nodeIndices.key(i);
        double w = m_ga->width(n), h = m_ga->height(n);
        w += m_nodePadding;
        h += m_nodePadding;
        double x = m_ga->x(n) - w/2, y = m_ga->y(n) - h/2;
        double X = x + w, Y = y + h;
        rs.push_back( new vpsc::Rectangle(x,X,y,Y) );
    }

    edge e = NULL;
    forall_edges(e,*m_graph) {
        node src = e->source();
        node tgt = e->target();
        int srcIndex = nodeIndices.value(src);
        int tgtIndex = nodeIndices.value(tgt);
        es.push_back(cola::Edge(srcIndex, tgtIndex) );
    }

    // Do layout.
    cola::ConstrainedFDLayout *fdlayout =
            new cola::ConstrainedFDLayout(rs,es,idealLength,preventOverlaps);
    fdlayout->setConstraints(sepcos);

    ConvTest1 *test = new ConvTest1(1e-3,100);
    test->setLayout(fdlayout);
    test->name = QString("BC-S4-postACA");
    //test->minIterations = 50;
    fdlayout->setConvergenceTest(test);
    fdlayout->run(true,true);

    // Read positions.
    node n;
    forall_nodes(n,*m_graph) {
        vpsc::Rectangle *r = rs.at(nodeIndices.value(n));
        m_ga->x(n) = r->getCentreX();
        m_ga->y(n) = r->getCentreY();
    }
}

// ------------------------------------------------------------------
// ACALayout --------------------------------------------------------

/** Build an ACA layout object for the passed OGDF graph.
  */
ACALayout::ACALayout(Graph &G, GraphAttributes &GA) :
    m_idealLength(100),
    m_nodePadding(50),
    m_preventOverlaps(false),
    m_addBendPointPenalty(true),
    m_postponeLeaves(true),
    m_useNonLeafDegree(true),
    m_allAtOnce(false),
    m_debugName(""),
    m_edgeLengths(NULL)
{
    QMap<int,int> deg; // multimap from rect indices to indices of neighbours
    node n = NULL;
    int i = 0;
    forall_nodes(n,G) {
        double w = GA.width(n), h = GA.height(n);
        double x = GA.x(n) - w/2, y = GA.y(n) - h/2;
        double X = x + w, Y = y + h;
        rs.push_back( new vpsc::Rectangle(x,X,y,Y) );
        m_nodeIndices.insert(n,i);
        i++;
    }
    edge e = NULL;
    forall_edges(e,G) {
        node src = e->source();
        node tgt = e->target();
        int srcIndex = m_nodeIndices.value(src);
        int tgtIndex = m_nodeIndices.value(tgt);
        es.push_back( cola::Edge(srcIndex, tgtIndex) );
        deg.insertMulti(srcIndex,tgtIndex);
        deg.insertMulti(tgtIndex,srcIndex);
    }
    // Note leaves.
    foreach (int i, deg.uniqueKeys()) {
        QList<int> J = deg.values(i);
        if (J.size()==1) {
            leaves.append(i);
        }
    }
    if (m_useNonLeafDegree) {
        // Use non-leaf degree.
        QMap<int,int> nld;
        foreach (int i, deg.uniqueKeys()) {
            QList<int> J = deg.values(i);
            foreach (int j, J) {
                if (!leaves.contains(j)) {
                    nld.insertMulti(i,j);
                }
            }
        }
        deg = nld;
    }
    // Build deg2Nodes structure.
    foreach (int i, deg.uniqueKeys()) {
        QList<int> J = deg.values(i);
        qDebug() << QString("Node %1 degree %2").arg(i).arg(J.size());
        if (J.size()==2) {
            deg2Nodes.insertMulti(i,J.at(0));
            deg2Nodes.insertMulti(i,J.at(1));
        }
    }
}

void ACALayout::nodePadding(double P) {
    m_nodePadding = P;
    double h = P/2;
    foreach (vpsc::Rectangle *r, rs) {
        double x = r->getMinX(), X = r->getMaxX();
        double y = r->getMinY(), Y = r->getMaxY();
        x -= h; X += h; y -= h; Y += h;
        r->setMinD(0,x);
        r->setMaxD(0,X);
        r->setMinD(1,y);
        r->setMaxD(1,Y);
    }
}

void ACALayout::nodePadding(QList<double> P) {
    int i = 0;
    foreach (vpsc::Rectangle *r, rs) {
        double h = P.at(i++)/2;
        double x = r->getMinX(), X = r->getMaxX();
        double y = r->getMinY(), Y = r->getMaxY();
        x -= h; X += h; y -= h; Y += h;
        r->setMinD(0,x);
        r->setMaxD(0,X);
        r->setMinD(1,y);
        r->setMaxD(1,Y);
    }
}

/** Run the layout algorithm.
  */
void ACALayout::run(void)
{
    initialLayout();
    if (m_allAtOnce) {
        acaLoopAllAtOnce();
    } else {
        acaLoopOneByOne();
    }
    //finalLayout();
}

void ACALayout::readPositions(Graph &G, GraphAttributes &GA)
{
    node n = NULL;
    forall_nodes(n,G) {
        vpsc::Rectangle *r = rs.at(m_nodeIndices.value(n));
        GA.x(n) = r->getCentreX();
        GA.y(n) = r->getCentreY();
    }
}

QMap<vpsc::Dim,EdgeNode> ACALayout::generateEdgeNodes(void) {
    QMap<vpsc::Dim,EdgeNode> ens;
    int N = rs.size();
    for (int i = 0; i < N; i++) {
        if (leaves.contains(i)) continue;
        for (int j = i+1; j < N; j++) {
            if (leaves.contains(j)) continue;
            int as = alignmentState(i,j);
            if ( (as & ACACONN) && (as & (ACAHORIZ | ACAVERT)) ) {
                vpsc::Rectangle *ri = rs.at(i);
                vpsc::Rectangle *rj = rs.at(j);
                double xi = ri->getMinX(), Xi = ri->getMaxX();
                double yi = ri->getMinY(), Yi = ri->getMaxY();
                double xj = rj->getMinX(), Xj = rj->getMaxX();
                double yj = rj->getMinY(), Yj = rj->getMaxY();
                QRectF Ri(QPointF(xi,yi),QPointF(Xi,Yi));
                QRectF Rj(QPointF(xj,yj),QPointF(Xj,Yj));
                EdgeNode E(Ri,Rj);
                E.setIndices(i,j);
                E.orientation = as & ACAHORIZ ? vpsc::HORIZONTAL : vpsc::VERTICAL;
                E.constraintDimension = as & ACAHORIZ ? vpsc::VERTICAL : vpsc::HORIZONTAL;
                ens.insertMulti(E.constraintDimension,E);
            }
        }
    }
    return ens;
}

/** Return list of all indices j such that rectangle i was deliberately
  * aligned with rectangle j in the ACA process.
  */
QList<int> ACALayout::delibAlignedWith(int i) {
    QList<int> J;
    for (int j = 0; j < alignmentState.cols; j++) {
        if (alignmentState(i,j) & ACADELIB) J.append(j);
    }
    return J;
}

/** Search for an alignment between rectangles of indices l and r.
  * If find such an alignment, add the passed offset to r and return true.
  * The change persists in the stored sepAligns.
  * Else return false.
  */
bool ACALayout::offsetAlignment(int l, int r, double offset) {
    bool found = false;
    cola::CompoundConstraints ccs;
    foreach (ACASeparatedAlignment *sa, sepAligns) {
        if ( (sa->rect1==l && sa->rect2==r) || (sa->rect1==r && sa->rect2==l) ) {
            vpsc::Dim algnDim = sa->af == ACAHORIZ ? vpsc::YDIM : vpsc::XDIM;
            sa->alignment = new cola::AlignmentConstraint(algnDim);
            sa->alignment->addShape(l,0);
            sa->alignment->addShape(r,offset);
            found = true;
        }
        ccs.push_back(sa->alignment);
        ccs.push_back(sa->separation);
    }
    m_ccs = ccs;
    return found;
}

int ACALayout::alignment(edge e) {
    node src = e->source(), tgt = e->target();
    int srcIndex = m_nodeIndices.value(src);
    int tgtIndex = m_nodeIndices.value(tgt);
    return alignmentState(srcIndex,tgtIndex);
}

QMap<edge,int> ACALayout::alignments(Graph &G) {
    QMap<edge,int> A;
    edge e;
    forall_edges(e,G) {
        A.insert( e, alignment(e) );
    }
    return A;
}

void ACALayout::initialPositions(void) {
    // TODO: In the future we might choose a good outer face, and positions
    // for all other nodes inside this.
    //
    // For now we just make sure no two nodes are coincident.
    moveCoincidentNodes();
}

void ACALayout::moveCoincidentNodes(void) {
    // TODO

    // Plan:
    // 1. Sort by x-coord.
    // 2. If any nodes share an x-coord, make sure their y-coords differ.
    //      A. Sort by y-coord.
    //      B. If yn < yn+1 = yn+2 = ... = yn+k < yn+k+1,
    //         then let dy = (yn+k+1 - yn)/(k+1),
    //         and for i = 1, 2, ..., k set yn+i = yn + i*dy.
}

void ACALayout::initialLayout(void) {
    // Do an initial unconstrained FD layout
    // with no overlap prevention, and long ideal edge length.
    bool preventOverlaps = false;
    double iL = m_idealLength;
    cola::ConstrainedFDLayout *fdlayout = m_edgeLengths == NULL ?
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps) :
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps,
                                          false,10.0,m_edgeLengths);
    ConvTest1 *test = new ConvTest1(1e-3,100);
    //test->minIterations = 100;
    test->setLayout(fdlayout);
    test->name = m_debugName+QString("-S1-noOP");
    fdlayout->setConvergenceTest(test);
    fdlayout->run(true,true);
    delete fdlayout;

    bool justFirst = false;
    if (justFirst) return;

    // Do another FD layout this time with
    // overlap prevention.
    preventOverlaps = true;
    // FIXME Deleting 'test' is causing a segfault.
    // Why? Accept memory leak for now.
    // delete test;

    //fdlayout = new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps);
    fdlayout = m_edgeLengths == NULL ?
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps) :
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps,
                                          false,10.0,m_edgeLengths);

    test = new ConvTest1(1e-4,200);
    //test->minIterations = 100;
    test->setLayout(fdlayout);
    test->name = m_debugName+QString("-S2-yesOP");
    fdlayout->setConvergenceTest(test);
    fdlayout->run(true,true);
    delete fdlayout;
}

void ACALayout::finalLayout(void) {
    // Keep the ACA constraints, but turn off OP.
    // (This allows stub nodes to retract toward their owners if they
    //  got separated. So this is specific to the NoNo layout, and should
    //  eventually be refactored out of the ACALayout object.) (FIXME)
    bool preventOverlaps = false;
    double iL = m_idealLength;
    cola::ConstrainedFDLayout *fdlayout = m_edgeLengths == NULL ?
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps) :
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps,
                                          false,10.0,m_edgeLengths);
    fdlayout->setConstraints(m_ccs);
    ConvTest1 *test = new ConvTest1(1e-3,100);
    //test->minIterations = 100;
    test->setLayout(fdlayout);
    test->name = m_debugName+QString("-S3-noOP-again");
    fdlayout->setConvergenceTest(test);
    fdlayout->skip_attractive_forces = false;
    fdlayout->run(true,true);
    delete fdlayout;

    bool justFirst = false;
    if (justFirst) return;

    // Do another FD layout this time with
    // overlap prevention.
    preventOverlaps = true;
    // FIXME Deleting 'test' is causing a segfault.
    // Why? Accept memory leak for now.
    // delete test;

    //fdlayout = new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps);
    fdlayout = m_edgeLengths == NULL ?
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps) :
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps,
                                          false,10.0,m_edgeLengths);
    fdlayout->setConstraints(m_ccs);
    test = new ConvTest1(1e-4,200);
    //test->minIterations = 100;
    test->setLayout(fdlayout);
    test->name = m_debugName+QString("-S3-yesOP-again");
    fdlayout->setConvergenceTest(test);
    fdlayout->run(true,true);
    delete fdlayout;
}

void ACALayout::acaLoopOneByOne(void) {
    double iL = m_idealLength;
    //double iL = 300;
    // Prepare the alignment state matrix.
    initAlignmentState();
    // Start main loop.
    cola::CompoundConstraints ccs;
    ACASeparatedAlignment *sa = chooseSA();
    bool preventOverlaps = true;
    int N = 0;
    cola::ConstrainedFDLayout *fdlayout = NULL;
    while (sa) {
        // Add the new separated alignment constraints.
        ccs.push_back(sa->separation);
        ccs.push_back(sa->alignment);
        // Redo the layout, with the new constraints.
        if (fdlayout!=NULL) delete fdlayout;

        //fdlayout = new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps);
        fdlayout = m_edgeLengths == NULL ?
                new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps) :
                new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps,
                                              false,10.0,m_edgeLengths);

        ConvTest1 *test = new ConvTest1(1e-3,3);
        test->setLayout(fdlayout);
        test->name = m_debugName+QString("-S3-ACA-%1").arg(++N,4,10,QLatin1Char('0'));
        fdlayout->setConvergenceTest(test);
        fdlayout->setConstraints(ccs);
        fdlayout->run(true,true);
        // Update state.
        updateAlignmentState(sa);
        // Store SA and choose next one.
        sepAligns.append(sa);
        sa = chooseSA();
    }
    if (fdlayout!=NULL) delete fdlayout;
    // Save the compound constraints.
    m_ccs = ccs;
}

void ACALayout::acaLoopAllAtOnce(void) {
    double iL = m_idealLength;
    // Prepare the alignment state matrix.
    initAlignmentState();
    // Start main loop.
    cola::CompoundConstraints ccs;
    ACASeparatedAlignment *sa = chooseSA();
    bool preventOverlaps = true;
    int N = 0;
    while (sa) {
        // Add the new separated alignment constraints.
        ccs.push_back(sa->separation);
        ccs.push_back(sa->alignment);
        // Update state.
        updateAlignmentState(sa);
        // Store SA and choose next one.
        sepAligns.append(sa);
        sa = chooseSA();
    }
    // Do the layout with the new constraints.

    //cola::ConstrainedFDLayout *fdlayout = new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps);
    cola::ConstrainedFDLayout *fdlayout = m_edgeLengths == NULL ?
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps) :
            new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps,
                                          false,10.0,m_edgeLengths);

    ConvTest1 *test = new ConvTest1(1e-3,3);
    test->setLayout(fdlayout);
    test->name = m_debugName+QString("-S3-ACA-%1").arg(++N,4,10,QLatin1Char('0'));
    fdlayout->setConvergenceTest(test);
    fdlayout->setConstraints(ccs);
    fdlayout->run(true,true);
    delete fdlayout;
    // FIXME: can we delete 'test' without segfault?
    // Save the compound constraints.
    m_ccs = ccs;
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
    double minDeflection = 10.0;
    // Consider each edge for potential alignment.
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
            if (m_addBendPointPenalty) dH += bendPointPenalty(src,tgt,ACAHORIZ);
            if (m_postponeLeaves) dH += leafPenalty(src,tgt);
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
            if (m_addBendPointPenalty) dV += bendPointPenalty(src,tgt,ACAVERT);
            if (m_postponeLeaves) dV += leafPenalty(src,tgt);
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
        // FIXME: Canvas::Horizontal? Shouldn't it be ACAHORIZ?
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
    double dx = tx-sx, dy = ty-sy;
    double dx2 = dx*dx, dy2 = dy*dy;
    double l = dx2 + dy2;
    double dfl = af==ACAHORIZ ? dy2/l : dx2/l;
    return dfl;
}

double ACALayout::bendPointPenalty(int src, int tgt, ACAFlags af)
{
    double penalty = 2;
    ACAFlags op = af==ACAHORIZ ? ACAVERT : ACAHORIZ;
    if (deg2Nodes.contains(src)) {
        QList<int> J = deg2Nodes.values(src);
        int j = J.at(0) == tgt ? J.at(1) : J.at(0);
        int as = alignmentState(src,j);
        if (as & op) return penalty;
    }
    if (deg2Nodes.contains(tgt)) {
        QList<int> J = deg2Nodes.values(tgt);
        int j = J.at(0) == src ? J.at(1) : J.at(0);
        int as = alignmentState(tgt,j);
        if (as & op) return penalty;
    }
    return 0;
}

double ACALayout::leafPenalty(int src, int tgt)
{
    double penalty = 3;
    return leaves.contains(src) || leaves.contains(tgt) ? penalty : 0;
}


// ------------------------------------------------------------------
// ExternalTree -----------------------------------------------------

ExternalTree::ExternalTree(node root, node rootInG, QList<node> nodes, QList<edge> edges,
                           shapemap nodeShapes, connmap edgeConns) :
    m_rootInG(rootInG),
    m_orientation(ogdf::leftToRight)
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

void ExternalTree::numberShapes(void) {
    foreach (node m, m_dunnartShapes.keys()) {
        ShapeObj *sh = m_dunnartShapes.value(m);
        sh->setLabel(QString::number(sh->internalId()));
    }
}

ShapeObj *ExternalTree::rootShape(void) {
    return m_dunnartShapes.value(m_root);
}

void ExternalTree::treeLayout(void) {
    TreeLayout tL;
    tL.orientation(m_orientation);
    tL.rootSelection(TreeLayout::rootIsSource);
    tL.orthogonalLayout(true);
    tL.call(*m_ga);
}

void ExternalTree::updateShapePositions(void) {
    node n;
    forall_nodes(n,*m_graph) {
        ShapeObj *sh = m_dunnartShapes.value(n);
        double x = m_ga->x(n), y = m_ga->y(n);
        sh->setCentrePos(QPointF(x,y));
    }
}

void ExternalTree::orthogonalRouting(bool b) {
    dunnart::Connector::RoutingType t = b ? Connector::orthogonal : Connector::polyline;
    foreach (Connector *conn, m_dunnartConns.values()) {
        if (conn==NULL) continue;
        conn->setRoutingType(t);
    }
}

QRectF ExternalTree::rootlessBBox(void) {
    double x = DBL_MAX, X = DBL_MIN;
    double y = DBL_MAX, Y = DBL_MIN;
    node n = NULL;
    forall_nodes(n,*m_graph) {
        // We leave the root node out of the bounding box.
        if (n == m_root) continue;
        //
        double cx = m_ga->x(n);
        double cy = m_ga->y(n);
        double w = m_ga->width(n);
        double h = m_ga->height(n);
        double u = cx - w/2.0;
        double v = cy - h/2.0;
        double U = u + w;
        double V = v + h;
        if (u<x) x = u;
        if (U>X) X = U;
        if (v<y) y = v;
        if (V>Y) Y = V;
    }
    return QRectF(QPointF(x,y),QPointF(X,Y));
}

void ExternalTree::translate(QPointF p) {
    node n = NULL;
    forall_nodes(n,*m_graph) {
        m_ga->x(n) += p.x();
        m_ga->y(n) += p.y();
    }
}

void ExternalTree::placeRootAt(QPointF p) {
    m_ga->x(m_root) = p.x();
    m_ga->y(m_root) = p.y();
}

/** A tree might need an alignment offset only if its
  * root is of degree 1.
  */
bool ExternalTree::needsAlignmentOffset(void) {
    return m_root->degree() == 1;
}

/** Returns an offset (which may be negative) such that
  * the constraint root_z + offset = stub_z will align
  * the tree and the root properly, where z is x if the
  * stub edge is aligned vertically, y if horizontally.
  */
double ExternalTree::alignmentOffset(void) {
    double offset = 0;
    if (needsAlignmentOffset()) {
        // In this case the root has just one child. Get it.
        node c; edge e;
        forall_adj_edges(e,m_root) { // should be only one
            c = m_root == e->source() ? e->target() : e->source();
        }
        QRectF bbox = rootlessBBox();
        if (m_orientation==leftToRight || m_orientation==rightToLeft) {
            offset = bbox.center().y() - m_ga->y(c);
        } else {
            offset = bbox.center().x() - m_ga->x(c);
        }
    }
    return offset;
}

// ------------------------------------------------------------------
// Orthowontist -----------------------------------------------------

Orthowontist::Orthowontist(Canvas *canvas) :
    m_canvas(canvas) {
}

void Orthowontist::run1(QList<CanvasItem*> items) {
    bool debug = true;
    bool useColours = false;
    bool showNumbers = false;
    bool drawStubnodes = false;
    bool drawDummynodes = true;
    shapemap nodeShapes;
    connmap edgeConns;
    Graph G;
    GraphAttributes GA(G);
    buildOGDFGraph(items,G,GA,nodeShapes,edgeConns);

    // 0. Compute average and min dimensions.
    double avgHeight=0, avgWidth=0, avgDim=0;
    double minHeight=DBL_MAX, minWidth=DBL_MAX, minDim=DBL_MAX;
    node n;
    int N = 0;
    forall_nodes(n,G) {
        double w = GA.width(n), h = GA.height(n);
        avgWidth += w;
        avgHeight += h;
        minHeight = h < minHeight ? h : minHeight;
        minWidth = w < minWidth ? w : minWidth;
        N++;
        if (debug) {
            qDebug() << QString("node size: %1 x %2").arg(GA.width(n)).arg(GA.height(n));
        }
    }
    avgWidth /= N;
    avgHeight /= N;
    avgDim = (avgWidth+avgHeight)/2;
    minDim = minHeight < minWidth ? minHeight : minWidth;

    // Ideal length and node padding
    double idealLength = 2*avgDim;
    double nodePadding = avgDim;

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
    // Debug:
    if (debug) {
        qDebug() << "\nCompound nontrivial biconnected components:";
        foreach (BiComp *B, BB) {
            qDebug() << B->listNodes();
            if (useColours) B->colourShapes();
        }
        foreach (ExternalTree *E, EE) {
            if (useColours) E->colourShapes();
        }
    }

    // 3. Make a map to say to which BC each BC shape belongs.
    QMap<ShapeObj*,BiComp*> shapesToBCs;
    foreach (BiComp *B, BB) {
        QList<ShapeObj*> shapes = B->allShapes();
        foreach (ShapeObj *sh, shapes) {
            shapesToBCs.insert(sh,B);
        }
    }

    // 4. TODO: compute the internal trees.

    // 5. Tack onto each B a "stub node" representing each of
    // the external trees that are attached to it.

    // Set stub size to be either half the minimum dimension in the graph,
    // or else the average dimension in the graph.
    bool useTinyStubs = true;
    QSizeF stubsize = useTinyStubs ? QSizeF(minDim/2,minDim/2) : QSizeF(avgDim,avgDim);

    foreach (ExternalTree *E, EE) {
        ShapeObj *sh = E->rootShape();
        BiComp *B = shapesToBCs.value(sh);
        B->addStubNodeForTree(E,stubsize);
    }

    if (showNumbers) {
        foreach (BiComp *B, BB) {
            B->numberShapes();
        }
        foreach (ExternalTree *E, EE) {
            E->numberShapes();
        }
    }

    // 6. Lay out each B in BB.
    foreach (BiComp *B, BB) {
        B->idealLength(idealLength);
        B->nodePadding(nodePadding);
        B->dummyNodeSize(QSizeF(avgDim,avgDim));
        B->layout();
    }

    if (debug) {

        m_canvas->stop_graph_layout();

        if (drawStubnodes) {
            foreach (BiComp *B, BB) {
                B->addStubNodeShapesToCanvas(m_canvas);
            }
        }

        if (drawDummynodes) {
            foreach (BiComp *B, BB) {
                B->addDummyNodeShapesToCanvas(m_canvas);
            }
        }

        foreach (ExternalTree *E, EE) {
            E->updateShapePositions();
            E->orthogonalRouting(true);
        }
        foreach (BiComp *B, BB) {
            B->updateShapePositions();
            //B->orthogonalRouting(true);
        }
        m_canvas->restart_graph_layout();

    }

    // ...

}

void Orthowontist::run2(QList<CanvasItem*> items) {
    bool debug = true;
    bool useColours = false;
    bool showNumbers = false;
    bool drawStubnodes = false;
    bool drawDummynodes = false;
    shapemap nodeShapes;
    connmap edgeConns;
    Graph G;
    GraphAttributes GA(G);
    buildOGDFGraph(items,G,GA,nodeShapes,edgeConns);

    // 0. Compute average and min dimensions.
    double avgHeight=0, avgWidth=0, avgDim=0;
    double minHeight=DBL_MAX, minWidth=DBL_MAX, minDim=DBL_MAX;
    node n;
    int N = 0;
    forall_nodes(n,G) {
        double w = GA.width(n), h = GA.height(n);
        avgWidth += w;
        avgHeight += h;
        minHeight = h < minHeight ? h : minHeight;
        minWidth = w < minWidth ? w : minWidth;
        N++;
        if (debug) {
            qDebug() << QString("node size: %1 x %2").arg(GA.width(n)).arg(GA.height(n));
        }
    }
    avgWidth /= N;
    avgHeight /= N;
    avgDim = (avgWidth+avgHeight)/2;
    minDim = minHeight < minWidth ? minHeight : minWidth;

    // Ideal length and node padding
    double idealLength = 2*avgDim;
    double nodePadding = avgDim;

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
    // Debug:
    if (debug) {
        qDebug() << "\nCompound nontrivial biconnected components:";
        foreach (BiComp *B, BB) {
            qDebug() << B->listNodes();
            if (useColours) B->colourShapes();
        }
        foreach (ExternalTree *E, EE) {
            if (useColours) E->colourShapes();
        }
    }

    // 3. Make a map to say to which BC each BC shape belongs.
    QMap<ShapeObj*,BiComp*> shapesToBCs;
    foreach (BiComp *B, BB) {
        QList<ShapeObj*> shapes = B->allShapes();
        foreach (ShapeObj *sh, shapes) {
            shapesToBCs.insert(sh,B);
        }
    }

    // 4. TODO: compute the internal trees.

    if (showNumbers) {
        foreach (BiComp *B, BB) {
            B->numberShapes();
        }
    }

    // Tell the BCs which are their root nodes.
    foreach (ExternalTree *E, EE) {
        ShapeObj *sh = E->rootShape();
        BiComp *B = shapesToBCs.value(sh);
        B->noteRoot(E);
    }

    // 6. Lay out each B in BB.
    foreach (BiComp *B, BB) {
        B->idealLength(idealLength);
        B->nodePadding(nodePadding);
        B->dummyNodeSize(QSizeF(avgDim,avgDim));
        B->filename = filename;
        B->layout2();
    }

    /*
    // Tack onto each B a "stub node" representing each of
    // the external trees that are attached to it.
    // Set stub size to be either half the minimum dimension in the graph,
    // or else the average dimension in the graph.
    bool useTinyStubs = true;
    QSizeF stubsize = useTinyStubs ? QSizeF(minDim/2,minDim/2) : QSizeF(avgDim,avgDim);
    foreach (ExternalTree *E, EE) {
        ShapeObj *sh = E->rootShape();
        BiComp *B = shapesToBCs.value(sh);
        B->addStubNodeForTree(E,stubsize);
    }
    if (showNumbers) {
        foreach (BiComp *B, BB) {
            B->numberShapes();
        }
        foreach (ExternalTree *E, EE) {
            E->numberShapes();
        }
    }
    */

    if (debug) {

        m_canvas->stop_graph_layout();

        if (drawStubnodes) {
            foreach (BiComp *B, BB) {
                B->addStubNodeShapesToCanvas(m_canvas);
            }
        }

        if (drawDummynodes) {
            foreach (BiComp *B, BB) {
                B->addDummyNodeShapesToCanvas(m_canvas);
            }
        }


        foreach (ExternalTree *E, EE) {
            E->updateShapePositions();
            E->orthogonalRouting(true);
        }


        foreach (BiComp *B, BB) {
            B->updateShapePositions();
            //B->orthogonalRouting(true);
        }
        m_canvas->restart_graph_layout();

    }


    // ...

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
            /*
            // Show IDs?
            bool showIDs = true;
            if (showIDs) {
                QString label = QString("%1").arg(shape->internalId());
                shape->setLabel(label);
            }
            */
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












































