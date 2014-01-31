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

namespace dunnart {

Orthowontist::Orthowontist(Canvas *canvas) :
    m_canvas(canvas)
{
}

void Orthowontist::run1(QList<CanvasItem*> items)
{
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

    //...

}

void Orthowontist::buildOGDFGraph(CanvasItemsList items,
        Graph &G, GraphAttributes &GA, shapemap &nodeShapes, connmap &edgeConns)
{
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


}
