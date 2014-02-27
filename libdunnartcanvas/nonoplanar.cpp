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

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Polygon_2.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Polygon_2<K> Polygon_2;

using namespace ogdf;

namespace ow {

// ------------------------------------------------------------------
// Planarization ----------------------------------------------------

Planarization::Planarization(Graph &G, GraphAttributes &GA,
                             QMap<edge,int> alignments,
                             QSizeF avgNodeSize, shapemap nodeShapes) :
    m_avgNodeSize(avgNodeSize)
{
    m_avgNodeDim = (avgNodeSize.width()+avgNodeSize.height())/2;
    m_dummyNodeSize = m_avgNodeSize;
    m_stubNodeSize = QSizeF(m_avgNodeSize.width()/5,m_avgNodeSize.height()/5);
    // Build local copy of original graph.
    m_graph = new Graph();
    m_ga = new GraphAttributes(*m_graph,
                   GraphAttributes::nodeGraphics |
                   GraphAttributes::nodeStyle |      // <-- Enables node stroke and filling
                   GraphAttributes::edgeGraphics
    );
    node n;
    forall_nodes(n, G) {
        node m = m_graph->newNode();
        m_origNodes.insert(m,n);
        ShapeObj *sh = nodeShapes.value(n);
        m_dunnartShapes.insert(m,sh);
        m_ga->x(m) = GA.x(n);
        m_ga->y(m) = GA.y(n);
        m_ga->width(m) = GA.width(n);
        m_ga->height(m) = GA.height(n);
    }
    edge e;
    forall_edges(e, G) {
        node m1 = m_origNodes.key(e->source());
        node m2 = m_origNodes.key(e->target());
        edge f = m_graph->newEdge(m1,m2);
        m_origEdges.insert(f,e);
        // Check alignment
        int af = alignments.value(e);
        m_alignments.insert(f,af);
        if (af & ACAHORIZ) {
            mH.append(new Edge(HTYPE,f,m_ga,this));
        } else if (af & ACAVERT) {
            mV.append(new Edge(VTYPE,f,m_ga,this));
        } else {
            mD.append(new Edge(DTYPE,f,m_ga,this));
        }
    }

    // Planarize by replacing each crossing with a dummy node.
    //
    // TODO: Should perturb graph if necessary so that no edge crosses over
    // the centre of any node, and so that never more than 2 edges intersect
    // at any given point. (For now we just hope that this is rare!)
    //

    simplerPlanarize();

    /*
    findHDHVCrossings();
    processCrossings();
    findVDCrossings();
    processCrossings();
    findDDCrossings();
    processCrossings();
    */

    m_ga->writeGML("planarization.gml");

    // Order the edges around each node corresponding to
    // the embedding given by the GraphAttributes object.
    NodeArray<SListPure<adjEntry> > adjList(*m_graph);
    EdgeComparer* ec = new EdgeComparer(*m_ga);
    node v;
    adjEntry ae;
    forall_nodes(v, *m_graph) {
        forall_adj(ae, v) {
            adjList[v].pushBack(ae);
        }
        adjList[v].quicksort(*ec);
        m_graph->sort(v, adjList[v]);
    }
    delete ec;

    m_comb = new CombinatorialEmbedding(*m_graph);
    mapNodesToFaces();
    findExternalFace();
    computeMinNodeSep();
    computeFreeSides();

    //DEBUG------------------------------------------------------------------
    /*
    qDebug() << QString("Number of faces: %1").arg(m_comb->numberOfFaces());
    face f = m_comb->firstFace();
    int I = m_comb->maxFaceIndex();
    for (int i = 0; i <= I; i++) {
        QString announce = QString("Face %1:").arg(i);
        qDebug() << announce;
        adjEntry ae = f->firstAdj();
        int J = f->size();
        QString nodeList = "    ";
        for (int j = 0; j < J; j++) {
            node n = ae->theNode();
            QString id = "";
            if (m_dummyNodes.contains(n)) {
                id = "d" + QString::number(m_dummyNodes.indexOf(n));
            } else {
                int i = m_dunnartShapes.value(n)->internalId();
                id = QString::number(i);
            }
            nodeList += id + " " ;
            // Next ae
            ae = f->nextFaceEdge(ae);
        }
        qDebug() << nodeList;
        // Next face.
        f = f->succ();
    }
    */
    //END DEBUG--------------------------------------------------------------

}

void Planarization::findExternalFace(void) {
    bool debug = true;
    // Find a node mx of minimal x-coordinate.
    node mx = NULL, n = NULL;
    forall_nodes(n,*m_graph) {
        if (mx==NULL || m_ga->x(n) < m_ga->x(mx)) mx = n;
    }
    if (debug) qDebug() << QString("Leftmost node ")+nodeIDString(mx);
    // Since there are no nodes to the left of mx, every neighbour a is
    // such that the edge (mx,a) makes an angle from 0 to pi (inclusive)
    // with the vertical vector (0,-1).
    // This allows us to identify a "first" edge in the clockwise ordering.
    // Namely, find a neighbouring node a0 for which the normalised dot
    // product of a0-mx with the vector (0,-1) is largest.
    // Then the edge (mx,a0) will come first in this ordering.
    double mxx = m_ga->x(mx), mxy = m_ga->y(mx);
    adjEntry ae = NULL;
    adjEntry ae0 = NULL;
    node a0 = NULL;
    double dp0 = -2; // All normalized dp's will be in [-1,1], so -2 == -infty
    forall_adj(ae,mx) {
        node a = ae->twinNode();
        double ax = m_ga->x(a), ay = m_ga->y(a);
        double dx = mxx-ax, dy = mxy-ay;
        double l = sqrt(dx*dx+dy*dy);
        double dp = dy/l;
        if (dp > dp0) {
            ae0 = ae;
            a0 = a;
            dp0 = dp;
        }
    }
    assert(ae0!=NULL);
    if (debug) qDebug() << QString("Next node ")+nodeIDString(a0);
    // Now form the list of adjacency elements around the external face.
    // The list begins with ae0, and always proceeds to the successor in
    // the clockwise order around the next node.
    // For some reason this is the cyclic "predecessor" in OGDF, not "successor".
    QList<adjEntry> aes;
    aes.append(ae0);
    ae = ae0->twin()->cyclicPred();
    while (ae != ae0) {
        aes.append(ae);
        ae = ae->twin()->cyclicPred();
    }
    // Debug output:
    if (debug) {
        QString s = "External face nodes: ";
        foreach (adjEntry z, aes) {
            node n = z->theNode();
            s += nodeIDString(n) + " " ;
        }
        qDebug() << s;
    }
    // Finally determine which of the faces is defined by these adjEntries.
    int N = aes.size();
    m_extFace = NULL;
    face f = m_comb->firstFace();
    int I = m_comb->maxFaceIndex();
    for (int i = 0; i <= I; i++) {
        if (debug) {
            qDebug() << "=======================================================";
            qDebug() << "Considering face "+QString::number(i);
        }
        // If the face has the wrong number of edges, then it is not a match.
        if (f->size()!=N) {
            if (debug) qDebug() << "Wrong number of edges.";
            f = f->succ();
            continue;
        }
        // Consider the first ae in face f.
        adjEntry ae = f->firstAdj();
        if (debug) qDebug() << "Begins with node "+nodeIDString(ae->theNode());
        // If it is not in aes at all, then f is not the external face.
        if (!aes.contains(ae)) {
            f = f->succ();
            continue;
        }
        // Consider the next ae in the face f.
        // Is this the next ae, on either side, in the external face that we found?
        int j0 = aes.indexOf(ae);
        ae = f->nextFaceEdge(ae);
        if (debug) {
            qDebug() << QString("Agrees at index %1 in our list.").arg(j0);
            qDebug() << "Next node is "+nodeIDString(ae->theNode());
        }
        int dir = 0;
        if (ae==aes.at( (j0+1) % N )) {
            dir = 1;
        } else if (ae==aes.at( (j0-1) % N )) {
            dir = -1;
        }
        // If not, then give up on face f.
        if (dir==0) {
            f = f->succ();
            continue;
        }
        // If so, then continue in that same direction, and see if
        // f has all the same ae's as aes.
        if (debug) qDebug() << "This agrees in direction: "+QString::number(dir);
        int j = 2;
        while (j < N) {
            ae = f->nextFaceEdge(ae);
            if (debug) qDebug() << "Next node is "+nodeIDString(ae->theNode());
            int k = (j0+dir*j) % N;
            if (ae!=aes.at(k)) {
                if (debug) qDebug() << QString("Disagrees at index %1.").arg(k);
                break;
            } else {
                if (debug) qDebug() << QString("Agrees at index %1.").arg(k);
            }
            j++;
        }
        if (j == N) {
            // In this case we have found the external face.
            m_extFace = f;
            break;
        }
        // Otherwise, try the next face.
        f = f->succ();
    }
    assert(m_extFace!=NULL);
    qDebug() << QString("External face is number %1").arg(m_extFace->index());
}

/***
  * Actually this approach with 'addCrossing' and 'connectCrossings' would also
  * work with the better approach which sorts the vert and horiz events.
  * TODO: Implement this.
  */
void Planarization::simplerPlanarize(void) {
    QList<Edge*> allEdges;
    allEdges.append(mH);
    allEdges.append(mV);
    allEdges.append(mD);
    int N = allEdges.size();
    for (int i = 0; i < N; i++) {
        Edge *e = allEdges.at(i);
        for (int j = i+1; j < N; j++) {
            Edge *f = allEdges.at(j);
            QPair<bool,QPointF> X = e->intersect(f);
            if (X.first) {
                QPointF p = X.second;
                addCrossing(e,f,p);
            }
        }
        e->connectCrossings();
    }
}

/***
  * This has a stupid O(m^2) runtime, but the faster approaches
  * seem to be vastly more complex. ... FIXME
  */
void Planarization::simplePlanarize(void) {
    QList<Edge*> allEdges;
    allEdges.append(mH);
    allEdges.append(mV);
    allEdges.append(mD);
    while (!allEdges.empty()) {
        Edge *e = allEdges.takeFirst();
        QList<Edge*> theRest(allEdges);
        foreach (Edge *f, theRest) {
            QPair<bool,QPointF> X = e->intersect(f);

            //debug
            int n = allEdges.size();
            int m = theRest.size();
            qDebug() << QString::number(n);
            //

            if (X.first) {
                QPointF p = X.second;
                QList<Edge*> newEdges = addDummyCross(e,f,p);
                allEdges.prepend(newEdges.at(1));
                allEdges.prepend(newEdges.at(0));
            }
        }
    }
}

/***
  * The edge makes itself now a representative of edge e1, and returns
  * a new edge representing e2.
  * This is meant to be used when the edge has been split into e1 and e2,
  * so it is as though the edge shortens to e1 and "rejects" the other half
  * of itself as a new edge, e2.
  */
Planarization::Edge* Planarization::Edge::rejectHalf(edge e1, edge e2) {
    m_ogdfEdge = e1;
    x0=m_ga->x(e1->source()), y0=m_ga->y(e1->source());
    x1=m_ga->x(e1->target()), y1=m_ga->y(e1->target());
    xmin = min(x0,x1);
    xmax = max(x0,x1);
    ymin = min(y0,y1);
    ymax = max(y0,y1);
    return new Edge(m_etype,e2,m_ga,m_plan);
}

void Planarization::addDummyNodeShapesToCanvas(Canvas *canvas) {
    //canvas->stop_graph_layout();
    foreach (ShapeObj *sh, m_dunnartShapes.values()) {
        canvas->addItem(sh);
    }
    //canvas->restart_graph_layout();
}

void Planarization::defineRootNodes(QList<node> roots) {
    foreach (node n, roots) {
        node m = m_origNodes.key(n);
        m_rootNodes.append(m);
    }
}

bool cmpTreesBySize(QPair<node,Planarization*> r1, QPair<node,Planarization*> r2) {
    Planarization *P = r1.second;
    node n1 = P->m_origNodes.value(r1.first);
    node n2 = P->m_origNodes.value(r2.first);
    QSizeF s1 = P->origRootToTreeSize.value(n1);
    QSizeF s2 = P->origRootToTreeSize.value(n2);
    double A1 = s1.width() * s1.height();
    double A2 = s2.width() * s2.height();
    return A2 < A1;
}

bool cmpFaces(SortableFace *s1, SortableFace *s2) {
    Planarization *P = s1->P;
    GraphAttributes *GA = P->m_ga;
    node r = s1->r; // should equal s2->r

    //for debugging:
    QString foo = P->nodeIDString(r);
    if (foo=="10") {
        qDebug() << foo;
    }
    //

    face f1 = s1->f, f2 = s2->f;
    // Prefer exterior face, primarily.
    face ext = P->m_extFace;
    if (f1==ext && f2!=ext) return true;
    if (f2==ext && f1!=ext) return false;
    // Secondarily, prefer a face in which the stub edge can be aligned horiz. or vert.
    QSet<int> C1 = P->m_nodeComb.value(r)->freeCardinals(f1,*GA);
    QSet<int> C2 = P->m_nodeComb.value(r)->freeCardinals(f2,*GA);
    if (!C1.empty() && C2.empty()) return true;
    if (C1.empty() && !C2.empty()) return false;
    // Thirdly and finally, prefer a face of greater area.
    return P->areaOfFace(f1) > P->areaOfFace(f2);
}

void Planarization::computeFreeSides(void) {
    // TODO
}

double Planarization::areaOfFace(face f) {
    Polygon_2 p;
    adjEntry ae = f->firstAdj();
    int J = f->size();
    for (int j = 0; j < J; j++) {
        node n = ae->theNode();
        double x = m_ga->x(n), y = m_ga->y(n);
        p.push_back(Point(x,y));
        ae = f->nextFaceEdge(ae);
    }
    double A = abs(p.area());
    return A;
}

QRectF Planarization::bboxWithoutTrees(void) {
    double x = DBL_MAX, X = DBL_MIN;
    double y = DBL_MAX, Y = DBL_MIN;
    node n = NULL;
    forall_nodes(n,*m_graph) {
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

QRectF Planarization::nodeRect(node n) {
    double cx = m_ga->x(n);
    double cy = m_ga->y(n);
    double w = m_ga->width(n);
    double h = m_ga->height(n);
    double x = cx - w/2.0;
    double y = cy - h/2.0;
    double X = x + w;
    double Y = y + h;
    return QRectF(QPointF(x,y),QPointF(X,Y));
}

vpsc::Rectangle *Planarization::vpscNodeRect(node n) {
    double cx = m_ga->x(n);
    double cy = m_ga->y(n);
    double w = m_ga->width(n);
    double h = m_ga->height(n);
    double x = cx - w/2.0;
    double y = cy - h/2.0;
    double X = x + w;
    double Y = y + h;
    return new vpsc::Rectangle(x,X,y,Y);
}

QList<EdgeNode*> Planarization::genEdgeNodesForFace(face f) {
    QList<EdgeNode*> ens;
    adjEntry ae = f->firstAdj();
    int J = f->size();
    for (int j = 0; j < J; j++) {
        node n = ae->theNode(), t = ae->twinNode();
        QRectF Rn = nodeRect(n), Rt = nodeRect(t);
        EdgeNode *E = new EdgeNode(Rn,Rt);
        int ni = m_nodeIndices.value(n), ti = m_nodeIndices.value(t);
        E->setIndices(ni,ti);
        ens.append(E);
        // Next ae
        ae = f->nextFaceEdge(ae);
    }
    return ens;
}

/***
  * Generate separation constraints to keep the nodes a distance of gap away from
  * the EdgeNodes, in the dimension requested.
  */
cola::CompoundConstraints Planarization::genNodeEdgeSepCos(
        vpsc::Dim dim, QList<node> ns, QList<EdgeNode *> ens, double gap) {
    cola::CompoundConstraints ccs;
    double pad = gap/2;
    unsigned int priority = cola::PRIORITY_NONOVERLAP - 1; //copied from colafd.cpp. What does it do?
    cola::NonOverlapConstraints *noc =
            new cola::NonOverlapConstraints(priority);
    int Ne = ens.size();
    int Nn = ns.size();
    // Variables and Rectangles.
    vpsc::Variables vs;
    vpsc::Rectangles rs;
    // Add edgenodes first.
    for (int i = 0; i < Ne; i++) {
        EdgeNode E = *ens.at(i);
        noc->addShape(i, E.bbox.width()/2 + pad, E.bbox.height()/2 + pad);
        double p = dim==vpsc::HORIZONTAL ? E.bbox.center().x() : E.bbox.center().y();
        vpsc::Variable *v = new vpsc::Variable(i,p);
        vs.push_back(v);
        vpsc::Rectangle *r = new vpsc::Rectangle(E.bbox.left(),E.bbox.right(),E.bbox.top(),E.bbox.bottom());
        rs.push_back(r);
    }
    // Add nodes.
    for (int j = 0; j < Nn; j++) {
        node s = ns.at(j);
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

        // Reject constraints that hold between two edgenodes or two nodes.
        if ( (l < Ne && r < Ne) || (l >= Ne && r >= Ne) ) continue;

        // Convert IDs to global indices.
        l = l < Ne ? ens.at(l)->srcIndex : m_nodeIndices.value(ns.at(l-Ne));
        r = r < Ne ? ens.at(r)->srcIndex : m_nodeIndices.value(ns.at(r-Ne));

        cola::SeparationConstraint *s = new cola::SeparationConstraint(dim,l,r,c->gap);
        ccs.push_back(s);
    }
    // Clean up
    // FIXME: why does this cause segfaults?
    //foreach (vpsc::Variable *v, vs) delete v;
    return ccs;
}

/***
  * Return separation constraints to create space for node s0 in face f0.
  */
cola::CompoundConstraints Planarization::faceLiftForNode(face f0, node s0, double gap) {
    cola::CompoundConstraints ccs;
    QList<EdgeNode*> ens = genEdgeNodesForFace(f0);
    QList<node> ns;
    ns.append(s0);
    cola::CompoundConstraints ccx = genNodeEdgeSepCos(vpsc::HORIZONTAL,ns,ens,gap);
    cola::CompoundConstraints ccy = genNodeEdgeSepCos(vpsc::VERTICAL,  ns,ens,gap);
    foreach (cola::CompoundConstraint *cc, ccx) ccs.push_back(cc);
    foreach (cola::CompoundConstraint *cc, ccy) ccs.push_back(cc);
    return ccs;
}

OrdAlign *Planarization::ordAlignForNodes(node s, node t, ACAFlags af) {
    double sx = m_ga->x(s), sy = m_ga->y(s);
    double sw = m_ga->width(s), sh = m_ga->height(s);
    double tx = m_ga->x(t), ty = m_ga->y(t);
    double tw = m_ga->width(t), th = m_ga->height(t);
    vpsc::Dim sepDim = af == ACAHORIZ ? vpsc::XDIM : vpsc::YDIM;
    vpsc::Dim algnDim = af == ACAHORIZ ? vpsc::YDIM : vpsc::XDIM;
    double sep = af == ACAHORIZ ? (sw+tw)/2.0 : (sh+th)/2.0;
    int i = m_nodeIndices.value(s), j = m_nodeIndices.value(t);
    int l = af == ACAHORIZ ? (sx<tx ? i : j) : (sy<ty ? i : j);
    int r = l == i ? j : i;
    cola::SeparationConstraint *sepc = new cola::SeparationConstraint(sepDim,l,r,sep);
    cola::AlignmentConstraint *algn = new cola::AlignmentConstraint(algnDim);
    algn->addShape(l,0);
    algn->addShape(r,0);
    OrdAlign *oa = new OrdAlign(sepc,algn);
    return oa;
}

/***
  * Create all "ordered alignment" constraints for aligned edges.
  * (I.e. one alignment and one separation constraint for each.)
  */
cola::CompoundConstraints Planarization::ordAlignsForEdges(void) {
    cola::CompoundConstraints ccs;
    // Edges in graph
    edge e;
    forall_edges(e,*m_graph) {
        int a = m_alignments.value(e);
        if (a!=0) {
            node s = e->source(), t = e->target();
            ACAFlags af = (ACAFlags) a;
            OrdAlign *oa = ordAlignForNodes(s,t,af);
            ccs.push_back(oa->sep);
            ccs.push_back(oa->algn);
        }
    }
    // Stub edges
    foreach (node r, m_rootNodes) {
        node s = m_rootsToStubs.value(r);
        if (m_stubAlignments.keys().contains(s)) {
            int a = m_stubAlignments.value(s);
            ACAFlags af = (ACAFlags) a;
            OrdAlign *oa = ordAlignForNodes(r,s,af);
            ccs.push_back(oa->sep);
            ccs.push_back(oa->algn);
        }
    }
    return ccs;
}

double Planarization::edgeLengthForNodes(node s, node t) {
    double sw = m_ga->width(s), sh = m_ga->height(s);
    double tw = m_ga->width(t), th = m_ga->height(t);
    double d = m_avgNodeDim;
    return (sw+sh)/4 + d + (tw+th)/4;
}

void Planarization::expand(int steps) {
    for (int n = 1; n <= steps; n++) {
        // 1. Define rs, es, ccs, eLengths.

        // rs
        vpsc::Rectangles rs;
        foreach (node u, m_normalNodes) rs.push_back(vpscNodeRect(u));
        foreach (node u, m_dummyNodes)  rs.push_back(vpscNodeRect(u));
        foreach (node s, m_stubNodes) {
            // Set size.
            node r = m_rootsToStubs.key(s);
            node r0 = m_origNodes.value(r);
            QSizeF S1 = origRootToTreeSize.value(r0);
            QSizeF S0 = m_stubNodeSize;
            double w0 = S0.width(), h0 = S0.height();
            double w1 = S1.width(), h1 = S1.height();
            double t = ((double)n) / ((double)steps);
            double w = w0 + (w1 - w0)*t;
            double h = h0 + (h1 - h0)*t;
            m_ga->width(s) = w;
            m_ga->height(s) = h;
            // Make rectangle.
            rs.push_back(vpscNodeRect(s));
        }

        // es and eLengths
        std::vector<cola::Edge> es;
        int Ne = m_normalEdges.size() + m_dummyEdges.size() + m_rootNodes.size();
        double *eLengths = new double[Ne];
        int j = 0;
        // edges actually in the graph
        foreach (edge e, m_normalEdges) {
            node s = e->source(), t = e->target();
            int si = m_nodeIndices.value(s), ti = m_nodeIndices.value(t);
            es.push_back( cola::Edge(si,ti) );
            eLengths[j] = edgeLengthForNodes(s,t);
            j++;
        }
        foreach (edge e, m_dummyEdges) {
            node s = e->source(), t = e->target();
            int si = m_nodeIndices.value(s), ti = m_nodeIndices.value(t);
            es.push_back( cola::Edge(si,ti) );
            eLengths[j] = edgeLengthForNodes(s,t);
            j++;
        }
        // stub edges
        foreach (node r, m_rootNodes) {
            node s = m_rootsToStubs.value(r);
            int ri = m_nodeIndices.value(r), si = m_nodeIndices.value(s);
            es.push_back( cola::Edge(ri,si) );
            eLengths[j] = edgeLengthForNodes(r,s);
            j++;
        }

        // ccs
        cola::CompoundConstraints ccs;
        // A. ordered alignments for aligned edges
        cola::CompoundConstraints oa = ordAlignsForEdges();
        foreach (cola::CompoundConstraint *cc, oa) ccs.push_back(cc);
        // B. face-lift constraints for trees
        foreach (node r, m_rootNodes) {
            face f = m_faceAssigns.value(r);
            node s = m_rootsToStubs.value(r);
            double gap = m_avgNodeDim;
            cola::CompoundConstraints fl = faceLiftForNode(f,s,gap);
            foreach (cola::CompoundConstraint *cc, fl) ccs.push_back(cc);
        }

        // 2. Run the FD layout.
        bool preventOverlaps = true;
        double iL = 1.0;
        cola::ConstrainedFDLayout *fdlayout =
                new cola::ConstrainedFDLayout(rs,es,iL,preventOverlaps,
                                              false,10.0,eLengths);
        fdlayout->setConstraints(ccs);
        ConvTest1 *test = new ConvTest1(1e-4,100);
        //test->minIterations = 10;
        test->setLayout(fdlayout);
        test->name = QString("Expand-%1").arg(n);
        fdlayout->setConvergenceTest(test);
        fdlayout->run(true,true);
        delete fdlayout;

        // 3. Update m_ga.
        node n;
        forall_nodes(n,*m_graph) {
            vpsc::Rectangle *r = rs.at(m_nodeIndices.value(n));
            m_ga->x(n) = r->getCentreX();
            m_ga->y(n) = r->getCentreY();
        }
    }
}

void Planarization::chooseGreedyTreeFaces(void) {
    bool debug = true;
    // Sort trees by descending size
    QList< QPair<node,Planarization*> > treez;
    foreach (node r, m_rootNodes) {
        treez.append(QPair<node,Planarization*>(r,this));
    }
    qSort(treez.begin(), treez.end(), cmpTreesBySize);
    QList<node> trees;
    for (int i = 0; i < treez.size(); i++) {
        QPair<node,Planarization*> foo = treez.at(i);
        trees.append(foo.first);
    }
    // Place the trees, in order of descending size.
    // Keep track of stub nodes in case we want to remove them.
    QList<node> stubs;
    // Get the initial bounding box, before we start to place any trees.
    QRectF bbox = bboxWithoutTrees();
    foreach (node r0, trees) {
        // debug
        if (debug) {
            qDebug() << QString("Placing tree %1.").arg(nodeIDString(r0));
            node o = m_origNodes.value(r0);
            QSizeF s = origRootToTreeSize.value(o);
            double A = s.width() * s.height();
            qDebug() << QString("  Tree area: %1").arg(A);
        }
        // Sort the faces for this tree.
        QList<face> F = m_nodeComb.value(r0)->faces();
        QList<SortableFace*> SF;
        foreach (face f, F) {
            SortableFace *sf = new SortableFace(f,r0,this);
            SF.append(sf);
        }
        qSort(SF.begin(), SF.end(), cmpFaces);
        QList<face> faces;
        for (int i = 0; i < SF.size(); i++) {
            SortableFace *sf = SF.at(i);
            faces.append(sf->f);
        }
        // debug
        if (debug) {
            QString s = "  Face order:\n";
            foreach (face f, faces) {
                s += QString::number(f->index()) + " ";
            }
            qDebug() << s;
        }
        // Choose the first face.
        face f0 = faces.first();
        // Place the tree at r0 in face f0
        m_faceAssigns.insert(r0,f0);

        // Choose a direction for the stub edge.
        NodeCombStruct *ncs = m_nodeComb.value(r0);
        QSet<int> cards = ncs->freeCardinals(f0,*m_ga);
        QPointF n;
        int cardinalDirec = -1;
        if (cards.empty()) {
            // In this case just use the normal into the face.
            n = m_nodeComb.value(r0)->normalIntoFace(f0,*m_ga);
        } else {
            // Choose to align the stub edge in one of the available cardinal directions.
            // TODO: use the bounding box to make an educated choice.
            // For now, we just choose the first one.
            cardinalDirec = cards.toList().first();
            // Map cardinal direction numbers 0,1,2,3 to corresp. unit vectors
            // (1,0), (0,-1), (-1,0), (0,1), in excessively clever way.
            // (Recall, y-axis increases downward.)
            int d = cardinalDirec;
            int d2 = d*d;
            int d3 = d2*d;
            double x = (d3-3*d2-d+3)/3;
            double y = (-d3+6*d2-8*d)/3;
            n = QPointF(x,y);
        }

        //DEBUG
        if (isnan(n.x())) {
            QString a = nodeIDString(r0);
            QList<adjEntry> aes = m_nodeComb.value(r0)->aesPerFace.values(f0);
            QString b1 = nodeIDString(aes.at(0)->theNode());
            QString b2 = nodeIDString(aes.at(0)->twinNode());
            QString c1 = nodeIDString(aes.at(1)->theNode());
            QString c2 = nodeIDString(aes.at(1)->twinNode());
            qDebug() << QString("ERROR: Singular normal: %1, %2, %3, %4, %5.").arg(a).arg(b1).arg(b2).arg(c1).arg(c2);
        }
        //END DEBUG

        // For now just do a simple scaling.
        double rw = m_ga->width(r0), rh = m_ga->height(r0);
        double sw = m_stubNodeSize.width(), sh = m_stubNodeSize.height();
        double d = ( sqrt(rw*rw+rh*rh) + sqrt(sw*sw+sh*sh) ) / 2;
        // Create stub node and place it.
        node stub = m_graph->newNode();
        m_stubNodes.append(stub);
        stubs.append(stub);
        m_rootsToStubs.insert(r0,stub);
        double x0 = m_ga->x(r0), y0 = m_ga->y(r0);
        double x1 = x0 + d*n.x(), y1 = y0 + d*n.y();
        m_ga->x(stub) = x1;
        m_ga->y(stub) = y1;
        m_ga->width(stub) = sw;
        m_ga->height(stub) = sh;
        // Record alignment if there was one.
        if (cardinalDirec >= 0) m_stubAlignments.insert(stub,cardinalDirec);
        // Record chosen position for retrieval by BC.
        // FIXME: Do we actually need this?
        node rOrig = m_origNodes.value(r0);
        QPointF pos(x1,y1);
        origRootToStubPos.insert(rOrig,pos);
    }

    // Diagnostic output:
    writeOutGraphWithStubs("greedy_FC_"+filename+".gml");

    // Cleanup
    bool removeStubs = false;
    if (removeStubs) {
        foreach (node s, stubs) {
            m_graph->delNode(s);
        }
    }

    // Index the nodes.
    indexNodesAndEdges();
}

void Planarization::indexNodesAndEdges(void) {
    // NODES
    // First figure out which are the normal nodes.
    // (Dummy nodes and stub nodes are already known.)
    node n;
    forall_nodes(n,*m_graph) {
        if (!m_dummyNodes.contains(n)) m_normalNodes.append(n);
    }
    // Now make the mapping from nodes to indices, which saves us
    // from having to search through lists for nodes later.
    m_nodeIndices.clear();
    int i = 0;
    // 1. Normal nodes
    foreach (node n, m_normalNodes) {
        m_nodeIndices.insert(n,i); i++;
    }
    // 2. Dummy nodes
    foreach (node n, m_dummyNodes) {
        m_nodeIndices.insert(n,i); i++;
    }
    // 3. Stub nodes
    foreach (node n, m_stubNodes) {
        m_nodeIndices.insert(n,i); i++;
    }

    // EDGES
    // Figure out which are the normal ones
    edge e;
    forall_edges(e,*m_graph) {
        if (!m_dummyEdges.contains(e)) m_normalEdges.append(e);
    }
    // Now make mapping.
    m_edgeIndices.clear();
    int j = 0;
    foreach (edge e, m_normalEdges) {
        m_edgeIndices.insert(e,j); j++;
    }
    foreach (edge e, m_dummyEdges) {
        m_edgeIndices.insert(e,j); j++;
    }
}

void Planarization::writeOutGraphWithStubs(QString fn) {
    Graph G;
    GraphAttributes GA(G);
    QMap<node,node> nodemap;
    node n;
    forall_nodes(n,*m_graph) {
        node m = G.newNode();
        nodemap.insert(n,m);
        GA.x(m) = m_ga->x(n);
        GA.y(m) = m_ga->y(n);
        GA.width(m) = m_ga->width(n);
        GA.height(m) = m_ga->height(n);
    }
    edge e;
    forall_edges(e,*m_graph) {
        node m1 = nodemap.value(e->source());
        node m2 = nodemap.value(e->target());
        G.newEdge(m1,m2);
    }
    foreach (node r, m_rootsToStubs.keys()) {
        node s = m_rootsToStubs.value(r);
        node m1 = nodemap.value(r);
        node m2 = nodemap.value(s);
        G.newEdge(m1,m2);
    }
    char *c = new char[fn.size()+1];
    for (int i = 0; i < fn.size(); i++) c[i] = fn.at(i).toAscii();
    c[fn.size()] = (char)(0);
    GA.writeGML(c);
}

void Planarization::chooseCombTreeFaces(void) {
    bool debug = true;
    // Assign weights to the faces.
    // Initialize to 0.
    QMap<face,double> weight;
    face f = m_comb->firstFace();
    int I = m_comb->maxFaceIndex();
    for (int i = 0; i <= I; i++) {
        weight.insert(f,0);
        f = f->succ();
    }
    // Add tree weights to the totals for the faces.
    foreach (node r, m_rootNodes) {
        QList<face> faces = m_nodeComb.value(r)->faces();
        int n = faces.size();
        // Each face adjacent to r gets 1/n added to its total weight.
        foreach (face f, faces) {
            if (f==m_extFace) continue; // external face is infinite sink
            double u = weight.value(f);
            u += 1/(double)n;
            weight.insert(f,u);
        }
    }
    // debug
    if (debug) {
        qDebug() << "==============================================================";
        qDebug() << "COMBINATORIAL TREE FACE ASSIGNMENT";
        face f = m_comb->firstFace();
        int I = m_comb->maxFaceIndex();
        for (int i = 0; i <= I; i++) {
            qDebug() << QString("Face %1 has weight %2.").arg(i).arg(weight.value(f));
            f = f->succ();
        }
    }
    // Prepare list of trees to be placed.
    QList<node> trees = QList<node>(m_rootNodes);
    int numTrees = m_rootNodes.size();
    // Keep track of stub nodes in case we want to remove them.
    QList<node> stubs;
    // Keep track of stub edges so they can be added at the end.
    // It is important that they not be added during the process of tree placement, or
    // else they will mess up the adjacency structure in the graph.
    QList< QPair<node,node> > stubEdges;
    // Loop until no trees remain.
    while (!trees.empty()) {
        // debug
        if (debug) {
            qDebug() << "---------------------------------------------------------------";
            QString s = "Trees remain at: ";
            foreach (node t, trees) s += nodeIDString(t) + " ";
            qDebug() << s;
        }
        // Choose a tree having an adjacent face of least weight.
        double u0 = numTrees+1; // effectively +infty since no face has more weight than numTrees
        node r0 = NULL;
        face f0 = NULL;
        int i0 = 0;
        double n0 = 1;
        QList<face> F0;
        for (int i = 0; i < trees.size(); i++) {
            node r = trees.at(i);
            QList<face> F = m_nodeComb.value(r)->faces();
            int n = F.size();
            foreach (face f, F) {
                double u = weight.value(f);
                if (u < u0) {
                    r0 = r; f0 = f; i0 = i; n0 = n; F0 = F; u0= u;
                }
            }
        }
        // List all faces of this tree that are of minimal weight.
        QList<face> M;
        foreach (face f, F0) {
            if (weight.value(f)==u0) M.append(f);
        }
        // debug
        if (debug) {
            qDebug() << "Chose tree at " + nodeIDString(r0);
            qDebug() << QString("Faces of minimal weight %1 are:").arg(u0);
            QString s = "";
            foreach (face f, M) s += QString::number(f->index()) + " ";
            qDebug() << s;
        }
        // Choose a face of minimal weight.
        //if (M.size() >= 0) {
        if (M.size() < 2) {
            f0 = M.first();
        } else {
            // If there are two or more faces of minimal weight, choose one with greatest area.
            double A0 = 0;
            f0 = NULL;
            foreach (face f, M) {
                Polygon_2 p;
                adjEntry ae = f->firstAdj();
                int J = f->size();
                for (int j = 0; j < J; j++) {
                    node n = ae->theNode();
                    double x = m_ga->x(n), y = m_ga->y(n);
                    //if (debug) qDebug() << QString("Node %1 at %2, %3.").arg(nodeIDString(n)).arg(x).arg(y);
                    if (debug) qDebug() << QString("Node at %1, %2.").arg(x).arg(y);
                    p.push_back(Point(x,y));
                    // Next ae
                    ae = f->nextFaceEdge(ae);
                }
                double A = abs(p.area());
                if (debug) qDebug() << QString("Face %1 has area %2.").arg(f->index()).arg(A);
                if (A > A0) {
                    A0 = A;
                    f0 = f;
                }
            }
            if (debug) qDebug() << QString("Chose face %1.").arg(f0->index());
        }
        assert(f0!=NULL);
        // Remove r0 from list of trees to be placed.
        trees.removeAt(i0);
        // Move all weight for r0 into f0.
        foreach (face f, F0) {
            if (f==m_extFace) continue; // external face is infinite sink
            double u = weight.value(f);
            u += f==f0 ? (n0-1)/n0 : -1/n0;
            weight.insert(f,u);
        }
        // debug
        if (debug) {
            qDebug() << "New weights after assigning tree to face:";
            face f = m_comb->firstFace();
            int I = m_comb->maxFaceIndex();
            for (int i = 0; i <= I; i++) {
                qDebug() << QString("    Face %1 has weight %2.").arg(i).arg(weight.value(f));
                f = f->succ();
            }
        }
        // Place the tree at r0 in face f0
        QPointF n = m_nodeComb.value(r0)->normalIntoFace(f0,*m_ga);

        //DEBUG
        if (isnan(n.x())) {
            QString a = nodeIDString(r0);
            QList<adjEntry> aes = m_nodeComb.value(r0)->aesPerFace.values(f0);
            QString b1 = nodeIDString(aes.at(0)->theNode());
            QString b2 = nodeIDString(aes.at(0)->twinNode());
            QString c1 = nodeIDString(aes.at(1)->theNode());
            QString c2 = nodeIDString(aes.at(1)->twinNode());
            qDebug() << QString("ERROR: Singular normal: %1, %2, %3, %4, %5.").arg(a).arg(b1).arg(b2).arg(c1).arg(c2);
        }
        //END DEBUG

        // For now just do a simple scaling.
        double rw = m_ga->width(r0), rh = m_ga->height(r0);
        double sw = m_stubNodeSize.width(), sh = m_stubNodeSize.height();
        double d = ( sqrt(rw*rw+rh*rh) + sqrt(sw*sw+sh*sh) ) / 2;
        // Create stub node and place it.
        node stub = m_graph->newNode();
        stubs.append(stub);
        m_rootsToStubs.insert(r0,stub);
        double x0 = m_ga->x(r0), y0 = m_ga->y(r0);
        double x1 = x0 + d*n.x(), y1 = y0 + d*n.y();
        m_ga->x(stub) = x1;
        m_ga->y(stub) = y1;
        m_ga->width(stub) = sw;
        m_ga->height(stub) = sh;
        // Record chosen position for retrieval by BC.
        node rOrig = m_origNodes.value(r0);
        QPointF pos(x1,y1);
        origRootToStubPos.insert(rOrig,pos);
        // Note edge, to be added later.
        stubEdges.append(QPair<node,node>(r0,stub));
    }
    // Now can add the edges.
    for (int i = 0; i < stubEdges.size(); i++) {
        QPair<node,node> e = stubEdges.at(i);
        m_graph->newEdge(e.first,e.second);
    }
    // Output
    bool writeout = true;
    if (writeout) {
        // Write out the results.
        QString fn = "comb_FC_"+filename+".gml";
        char *c = new char[fn.size()+1];
        for (int i = 0; i < fn.size(); i++) c[i] = fn.at(i).toAscii();
        c[fn.size()] = (char)(0);
        m_ga->writeGML(c);
    }
    // Cleanup
    bool removeStubs = false;
    if (removeStubs) {
        foreach (node s, stubs) {
            m_graph->delNode(s);
        }
    }
}

void Planarization::layoutTreeForRoot(ExternalTree *E, node root) {
    node r = m_origNodes.key(root);
    node s = m_rootsToStubs.value(r);
    // Determine closest cardinal direction of line from tap to stub.
    double rx=m_ga->x(r),ry=m_ga->y(r);
    double sx=m_ga->x(s),sy=m_ga->y(s);
    double vx=sx-rx, vy=sy-ry;
    ogdf::Orientation orient = vx >= vy ?
                (vx >= -vy ? leftToRight : topToBottom) :
                (vx >= -vy ? bottomToTop : rightToLeft) ;
    // Lay out tree.
    E->orientation(orient);
    E->treeLayout();
    // Read size.
    QSizeF size = E->rootlessBBox().size();
    // If stub edge was aligned, offset that alignment if necessary.
    if (areAligned(r,s) && E->needsAlignmentOffset()) {
        double offset = E->alignmentOffset();
        offsetAlignment(r,s,offset);
    }
    // Set size in stub node.
    m_ga->width(s) = size.width();
    m_ga->height(s) = size.height();
}

void Planarization::offsetAlignment(node r, node s, double offset) {
    // TODO
}

bool Planarization::areAligned(node r, node s) {
    // TODO
    return false;
}

void Planarization::computeMinNodeSep(void) {
    // There is an O(n log n) algorithm for this by Shamos and Hoey
    // (and even an O(n) expected time probabilistic version),
    // but for now we just do the dumb O(n^2) comparison of all pairs of nodes.
    double d = DBL_MAX;
    QList<node> nodes;
    node n;
    forall_nodes(n,*m_graph) {
        nodes.append(n);
    }
    int N = nodes.size();
    for (int i = 0; i + 1 < N; i++) {
        node n = nodes.at(i);
        double x1 = m_ga->x(n), y1 = m_ga->y(n);
        for (int j = i + 1; j < N; j++) {
            node m = nodes.at(j);
            double x2 = m_ga->x(m), y2 = m_ga->y(m);
            double dx = x2-x1, dy = y2-y1;
            double l = sqrt(dx*dx+dy*dy);
            if (l < d) d = l;
        }
    }
    m_minNodeSep = d;
}

void Planarization::mapNodesToFaces(void) {
    bool debug = true;
    if (debug) qDebug() << QString("Number of faces: %1").arg(m_comb->numberOfFaces());
    // Initialize node combinatorial info map
    node n;
    forall_nodes(n,*m_graph) {
        m_nodeComb.insert(n, new NodeCombStruct(n));
    }
    // Iterate over faces
    face f = m_comb->firstFace();
    int I = m_comb->maxFaceIndex();
    for (int i = 0; i <= I; i++) {
        if (debug) qDebug() << QString("Face %1:").arg(i);
        adjEntry ae = f->firstAdj();
        int J = f->size();
        QString nodeList = "    ";
        for (int j = 0; j < J; j++) {
            node n = ae->theNode();
            NodeCombStruct *ncs = m_nodeComb.value(n);
            ncs->aesPerFace.insertMulti(f,ae);
            node t = ae->twinNode();
            NodeCombStruct *tcs = m_nodeComb.value(t);
            tcs->aesPerFace.insertMulti(f,ae);
            nodeList += nodeIDString(n) + " ";
            nodeList += nodeIDString(t) + " ";
            // Next ae
            ae = f->nextFaceEdge(ae);
        }
        if (debug) qDebug() << nodeList;
        // Next face.
        f = f->succ();
    }
}

void Planarization::chooseFDTreeFaces(void) {
    // Build ColaFD object.
    node n = NULL;
    int i = 0;
    forall_nodes(n,*m_graph) {
        double w = m_ga->width(n), h = m_ga->height(n);
        double x = m_ga->x(n) - w/2, y = m_ga->y(n) - h/2;
        double X = x + w, Y = y + h;
        rs.push_back( new vpsc::Rectangle(x,X,y,Y) );
        m_nodeIndices.insert(n,i);
        i++;
        //qDebug() << QString("node at: %1, %2").arg(m_ga->x(n)).arg(m_ga->y(n));
        //qDebug() << m_ga->colorNode(n).cstr();
        //qDebug() << QString("node shape: %1").arg(m_ga->shapeNode(n));
    }
    edge e = NULL;
    forall_edges(e,*m_graph) {
        node src = e->source();
        node tgt = e->target();
        int srcIndex = m_nodeIndices.value(src);
        int tgtIndex = m_nodeIndices.value(tgt);
        es.push_back( cola::Edge(srcIndex, tgtIndex) );
    }
    bool op = false;
    double il = m_idealLength;
    cola::ConstrainedFDLayout *fdlayout = new cola::ConstrainedFDLayout(rs,es,il,op);

    // Ask fdlayout for the force vectors.
    QList<node> stubs;
    foreach (node r, m_rootNodes) {
        int u = m_nodeIndices.value(r);
        std::vector<double> ng = fdlayout->negStressGradForNodeOnceRemoved(u);
        // Normalize force vector, then set length equal to twice
        // the dummy node dimension.
        double dx = ng.at(0), dy = ng.at(1);
        double l = sqrt(dx*dx + dy*dy);
        double len = 2*m_avgNodeDim;
        dx = len*dx/l;
        dy = len*dy/l;
        // Add stub node.
        node stub = m_graph->newNode();
        stubs.append(stub);
        m_graph->newEdge(r,stub);
        double x0 = m_ga->x(r), y0 = m_ga->y(r);
        double x1 = x0 + dx, y1 = y0 + dy;
        m_ga->x(stub) = x1;
        m_ga->y(stub) = y1;
        m_ga->width(stub) = m_stubNodeSize.width();
        m_ga->height(stub) = m_stubNodeSize.height();
        //m_ga->labelNode(stub) = ogdf::String('T');
    }
    // Output
    bool writeout = true;
    if (writeout) {
        // Write out the results.
        QString fn = "FD_FC_"+filename+".gml";
        char *c = new char[fn.size()+1];
        for (int i = 0; i < fn.size(); i++) c[i] = fn.at(i).toAscii();
        c[fn.size()] = (char)(0);
        m_ga->writeGML(c);
    }
    // Cleanup
    bool removeStubs = true;
    if (removeStubs) {
        foreach (node s, stubs) {
            m_graph->delNode(s);
        }
    }
}

void Planarization::findHDHVCrossings(void) {
    // Create edge event objects.
    int nH = mH.size(), nD = mD.size(), nV = mV.size();
    int nE = nH + 2*nD + 2*nV;
    EdgeEvent **events = new EdgeEvent*[nE];
    for (int i = 0; i < nH; i++) {
        events[i] = new EdgeEvent(HEDGE,mH.at(i));
    }
    for (int i = 0; i < nD; i++) {
        Edge *d = mD.at(i);
        events[nH+2*i] = new EdgeEvent(DOPENY,d);
        events[nH+2*i+1] = new EdgeEvent(DCLOSEY,d);
    }
    for (int i = 0; i < nV; i++) {
        Edge *v = mV.at(i);
        events[nH+2*nD+2*i] = new EdgeEvent(DOPENY,v);
        events[nH+2*nD+2*i+1] = new EdgeEvent(DCLOSEY,v);
    }
    // Sort them.
    qsort(events,nE,sizeof(EdgeEvent*),cmpEdgeEvent);
    // Find intersections.
    QList<Edge*> openEdges;
    for (int i = 0; i < nE; i++) {
        EdgeEvent *event = events[i];
        switch (event->m_eetype) {
        case DOPENY: {
            Edge *e = event->m_edge;
            e->m_openEdgeIndex = openEdges.size();
            openEdges.append(event->m_edge);
            break;
        }
        case DCLOSEY: {
            Edge *e = event->m_edge;
            int j = e->m_openEdgeIndex;
            openEdges.removeAt(j);
            break;
        }
        case HEDGE: {
            Edge *h = event->m_edge;
            // Check whether h crosses any of the currently open edges.
            foreach (Edge *e, openEdges) {
                if (h->sharesAnEndptWith(*e)) continue; // If they share an endpoint, then they do not cross.
                double y0 = h->constCoord();
                if (!e->coversY(y0)) continue;
                double x0 = e->x(y0);
                if (!h->coversX(x0)) continue;
                // If we reach this point, then h and e cross.
                // Their intersection point is (x0,y0).
                // For now we do not modify any of the existing edges; we merely
                // note the intersection and move on.
                QPointF p(x0,y0);
                h->addIntersection(e,p);
                e->addIntersection(h,p);
            }
        } // end case
        } // end switch
    }
}

void Planarization::Edge::processCrossing(Edge *e) {
    //
}

void Planarization::processCrossings(void) {
    QList<Edge*> newH;
    QList<Edge*> newV;
    QList<Edge*> newD;
    QList<Edge*> allEdges;
    allEdges.append(mH);
    allEdges.append(mV);
    allEdges.append(mD);
    foreach (Edge *e, allEdges) {
        //...
    }
    mH = newH;
    mV = newV;
    mD = newD;
}

void Planarization::findVDCrossings(void) {
    // Create edge event objects.
    int nV = mV.size(), nD = mD.size();
    int nE = nV + 2*nD;
    EdgeEvent **events = new EdgeEvent*[nE];
    for (int i = 0; i < nV; i++) {
        events[i] = new EdgeEvent(VEDGE,mV.at(i));
    }
    for (int i = 0; i < nD; i++) {
        Edge *d = mD.at(i);
        events[nV+2*i] = new EdgeEvent(DOPENY,d);
        events[nV+2*i+1] = new EdgeEvent(DCLOSEY,d);
    }
    // Sort them.
    qsort(events,nE,sizeof(EdgeEvent*),cmpEdgeEvent);
    // Scan once through the sorted list, catching intersections.
    QList<Edge*> openEdges;
    for (int i = 0; i < nE; i++) {
        EdgeEvent event = *events[i];
        switch (event.m_eetype) {
        case DOPENY: {
            Edge *e = event.m_edge;
            e->m_openEdgeIndex = openEdges.size();
            openEdges.append(event.m_edge);
            break;
        }
        case DCLOSEY: {
            int j = event.m_edge->m_openEdgeIndex;
            openEdges.removeAt(j);
            break;
        }
        case VEDGE: {
            Edge *v = event.m_edge;
            // Check whether v crosses any of the currently open diagonal edges.
            foreach (Edge *d, openEdges) {
                if (v->sharesAnEndptWith(*d)) continue; // If they share an endpoint, then they do not cross.
                double x0 = v->constCoord();
                if (!d->coversX(x0)) continue;
                double y0 = d->y(x0);
                if (!v->coversY(y0)) continue;
                // If we reach this point, then v and d cross.
                // Their intersection point is (x0,y0).
                addDummyCross(v,d,QPointF(x0,y0));
                // TODO: Now must push the two new V-edges created onto the
                // front of the event queue.
            }
        }
        }
    }
}

void Planarization::findDDCrossings(void) {
    int nD = mD.size();
    // FIXME: Have to do this with edge queues.
    for (int i = 0; i < nD; i++) {
        Edge *di = mD.at(i);
        for (int j = i+1; j < nD; j++) {
            Edge *dj = mD.at(j);
            QPair<bool,QPointF> X = di->intersectDiagonals(dj);
            if (X.first) {
                QPointF p = X.second;
                addDummyCross(di,dj,p);
                // TODO: add new edges to queues.
            }
        }
    }
}

void Planarization::addCrossing(Edge *e1, Edge *e2, QPointF p) {
    // Add dummy node.
    node dn = m_graph->newNode();
    m_ga->x(dn) = p.x();
    m_ga->y(dn) = p.y();
    m_ga->width(dn) = m_dummyNodeSize.width();
    m_ga->height(dn) = m_dummyNodeSize.height();
    m_dummyNodes.append(dn);
    // Make shape (but do not yet add to canvas).
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    ShapeObj *dummyShape = factory->createShape("org.dunnart.shapes.rect");
    dummyShape->setFillColour(QColor(224,0,0));
    dummyShape->setSize(m_dummyNodeSize);
    dummyShape->setCentrePos(p);
    dummyShape->setLabel(QString::number(m_dummyNodes.size()-1));
    m_dunnartShapes.insert(dn,dummyShape);
    // Tell edges about the crossing.
    e1->addCrossing(dn);
    e2->addCrossing(dn);
}

void Planarization::Edge::addCrossing(node n) {
    int i = 0;
    int N = m_crossings.size();
    double z0 = m_etype==VTYPE ? m_ga->y(n) : m_ga->x(n);
    while (i<N) {
        node m = m_crossings.at(i);
        double z = m_etype==VTYPE ? m_ga->y(m) : m_ga->x(m);
        if (z>z0) break;
        i++;
    }
    if (i<N) {
        m_crossings.insert(i,n);
    } else {
        m_crossings.append(n);
    }
}

void Planarization::Edge::connectCrossings(void) {
    // Add src and tgt of original edge to list of crossing nodes.
    QList<node> nodes(m_crossings);
    double z0 = m_etype==VTYPE ? y0 : x0;
    double z1 = m_etype==VTYPE ? y1 : x1;
    node src = m_ogdfEdge->source(), tgt = m_ogdfEdge->target();
    node first = z0<z1 ? src : tgt;
    node last  = z0<z1 ? tgt : src;
    nodes.prepend(first);
    nodes.append(last);
    // Delete original edge.
    m_plan->delEdge(m_ogdfEdge);
    // Add an edge between each pair of nodes.
    int N = nodes.size();
    for (int i = 0; i + 1 < N; i++) {
        node a = nodes.at(i), b = nodes.at(i+1);
        int af = 0;
        if (m_etype==HTYPE) {
            af = ACAHORIZ;
        } else if (m_etype==VTYPE) {
            af = ACAVERT;
        }
        m_plan->addDummyEdge(a,b,af);
    }
}

QList<Planarization::Edge*> Planarization::addDummyCross(Edge *e1, Edge *e2, QPointF p) {
    edge oe1 = e1->m_ogdfEdge, oe2 = e2->m_ogdfEdge;
    node src1 = oe1->source(), tgt1 = oe1->target(), src2 = oe2->source(), tgt2 = oe2->target();
    // Remove the edges.
    m_graph->delEdge(oe1);
    m_graph->delEdge(oe2);
    // Add dummy node.
    node dn = m_graph->newNode();
    m_ga->x(dn) = p.x();
    m_ga->y(dn) = p.y();
    m_ga->width(dn) = m_dummyNodeSize.width();
    m_ga->height(dn) = m_dummyNodeSize.height();
    // Make shape (but do not yet add to canvas).
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    ShapeObj *dummyShape = factory->createShape("org.dunnart.shapes.rect");
    dummyShape->setFillColour(QColor(224,0,0));
    dummyShape->setSize(m_dummyNodeSize);
    dummyShape->setCentrePos(p);
    m_dunnartShapes.insert(dn,dummyShape);
    // Add dummy edges.
    edge de1s = m_graph->newEdge(src1,dn);
    edge de1t = m_graph->newEdge(tgt1,dn);
    edge de2s = m_graph->newEdge(src2,dn);
    edge de2t = m_graph->newEdge(tgt2,dn);
    // Make records.
    DummyCross *dc = new DummyCross(dn,de1s,de1t,de2s,de2t);
    m_dummyCrosses.append(dc);
    m_dummyNodes.append(dn);
    m_dummyEdges.append(de1s);
    m_dummyEdges.append(de1t);
    m_dummyEdges.append(de2s);
    m_dummyEdges.append(de2t);
    // Split the two Edge structs, and return the halves they reject.
    QList<Edge*> newEdges;
    newEdges.append(e1->rejectHalf(de1s,de1t));
    newEdges.append(e2->rejectHalf(de2s,de2t));
    return newEdges;
}

QPair<bool,QPointF> Planarization::Edge::intersect(Edge *f) {
    // Prepare negative return value.
    QPair<bool,QPointF> noIntersect;
    noIntersect.first = false;
    noIntersect.second = QPointF(0,0);
    // If they share an endpoint, then they do not cross.
    if (sharesAnEndptWith(*f)) return noIntersect;
    // Otherwise they /may/ cross. Prepare positive return value.
    QPair<bool,QPointF> isect;
    isect.first = true;
    // Cases:
    if (m_etype==DTYPE || f->m_etype==DTYPE) {
        // In this case at least one of the edges is a diagonal.
        Edge *d = m_etype==DTYPE ? this : f; // d is DTYPE
        Edge *e = m_etype==DTYPE ? f : this; // e could be any type
        double x,y;
        switch (e->m_etype) {
        case HTYPE: {
            // e is horizontal
            y = e->constCoord();
            if (!d->coversY(y)) return noIntersect;
            x = d->x(y);
            if (!e->coversX(x)) return noIntersect;
            break;
        } // end case
        case VTYPE: {
            // e is vertical
            x = e->constCoord();
            if (!d->coversX(x)) return noIntersect;
            y = d->y(x);
            if (!e->coversY(y)) return noIntersect;
            break;
        } // end case
        case DTYPE: {
            // e and d are both diagonal
            if (e->slope()==d->slope()) return noIntersect; //No intersection if parallel.
            // Else compute point of intersection of lines.
            double m1 = e->slope(), m2 = d->slope(), b1 = e->yInt(), b2 = d->yInt();
            x = (b2-b1)/(m1-m2);
            y = m1*x+b1;
            // Does the point lie on both line /segments/?
            if (x<e->x0 || x>e->x1 || x<d->x0 || x>d->x1) return noIntersect;
            break;
        } // end case
        } // end switch
        isect.second = QPointF(x,y);
        return isect;
    } else {
        // In this case neither edge is diagonal.
        // So if their types are the same then it is HH or VV and they do not intersect.
        if (m_etype==f->m_etype) return noIntersect;
        // In this case they are H and V.
        Edge *h = m_etype==HTYPE ? this : f;
        Edge *v = m_etype==HTYPE ? f : this;
        if (h->coversX(v->constCoord()) && v->coversY(h->constCoord())) {
            isect.second = QPointF(v->constCoord(),h->constCoord());
            return isect;
        } else {
            return noIntersect;
        }
    }
}

/***
  * Return the point of intersection of this edge with another, provided
  * both are diagonal edges.
  */
QPair<bool,QPointF> Planarization::Edge::intersectDiagonals(Edge *d) {
    assert(m_etype==DTYPE && d->m_etype==DTYPE);
    QPair<bool,QPointF> noIntersect;
    noIntersect.first = false;
    noIntersect.second = QPointF(0,0);
    if (sharesAnEndptWith(*d)) return noIntersect; //No intersection if share endpt.
    if (slope()==d->slope()) return noIntersect; //No intersection if parallel.
    // Else compute point of intersection of lines.
    double m1 = slope(), m2 = d->slope(), b1 = yInt(), b2 = d->yInt();
    double x = (b2-b1)/(m1-m2);
    double y = m1*x+b1;
    // Does the point lie on both line /segments/?
    if (x<x0 || x>x1 || x<d->x0 || x>d->x1) return noIntersect;
    QPair<bool,QPointF> intersection;
    intersection.first = true;
    intersection.second = QPointF(x,y);
    return intersection;
}

static int cmpEdgeEvent(const void *p1, const void *p2) {
    QList<ow::Planarization::EdgeEvent> events;
    ow::Planarization::EdgeEvent e1 = ** (ow::Planarization::EdgeEvent **)(p1);
    ow::Planarization::EdgeEvent e2 = ** (ow::Planarization::EdgeEvent **)(p2);
    events.append(e1);
    events.append(e2);
    QList<double> coords;
    foreach (ow::Planarization::EdgeEvent event, events) {
        double z = 0;
        switch (event.m_eetype) {
        case ow::Planarization::HEDGE:
            z = event.m_edge->y0; break;
        case ow::Planarization::VEDGE:
            z = event.m_edge->x0; break;
        case ow::Planarization::DOPENX:
            z = event.m_edge->x0; break;
        case ow::Planarization::DOPENY:
            z = event.m_edge->y0; break;
        case ow::Planarization::DCLOSEX:
            z = event.m_edge->x1; break;
        case ow::Planarization::DCLOSEY:
            z = event.m_edge->y1; break;
        }
        coords.append(z);
    }
    double z1 = coords.at(0), z2 = coords.at(1);
    int c = z1 < z2 ? -1 : z1 > z2 ? 1 : 0;
    return c;
}

}
