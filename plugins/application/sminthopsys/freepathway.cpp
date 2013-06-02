/*
 * Sminthopsys - Dunnart Systems Biology plugin
 *
 * Copyright (C) 2011-2012  Monash University
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
 * Author(s): Steven Kieffer  <http://skieffer.info>
*/

#include <assert.h>

#include <QMap>
#include <QList>
#include <QTime>
#include <QPointF>
#include <QSizeF>

#include "freepathway.h"
#include "dsbreaction.h"
#include "dsbclone.h"

#include "libdunnartcanvas/shape.h"

#include "libogdf/ogdf/basic/Graph_d.h"
#include "libogdf/ogdf/basic/GraphAttributes.h"
#include "libogdf/ogdf/energybased/StressMajorizationSimple.h"

#include "plugins/application/sminthopsys/dsbnode.h"

namespace dunnart{

FreePathway::FreePathway(QList<DSBClone *> clones, QList<DSBReaction *> reacs) :
    m_clones(clones),
    m_reactions(reacs)
{
    //...
}

QSizeF FreePathway::getSize()
{
    return m_size;
}

void FreePathway::setRelPt(QPointF p)
{
    m_relpt = p;
}

void FreePathway::drawRelTo(QPointF q)
{
    QPointF r = m_relpt + q;
    drawAt(r);
}

void FreePathway::redraw()
{
    drawAt(m_basept);
}

QSizeF FreePathway::layout()
{
    foreach (DSBReaction *reac, m_reactions)
    {
        reac->layout();
    }

    nodemap nodeMap;
    ogdf::Graph *G = getOGDFGraph(nodeMap);
    jog(10.0, nodeMap);
    ogdf::GraphAttributes GA(*G);
    extractPosAndSize(nodeMap, GA);
    ogdf::StressMajorization strmaj;
    strmaj.call(GA);
    injectPositions(nodeMap, GA);
    QRectF box = getBbox(nodeMap);
    m_size = box.size();
    return m_size;
}

QRectF FreePathway::getBbox(nodemap &nodeMap)
{
    QList<DSBNode*> nodes = nodeMap.keys();
    assert(nodes.size() > 0);
    DSBNode *f = nodes.first();
    QRectF box = f->getShape()->boundingRect();
    for (int i = 1; i < nodes.size(); i++)
    {
        box = box.united(nodes.at(i)->getShape()->boundingRect());
    }
    return box;
}

void FreePathway::jog(double scale, nodemap &nodeMap)
{
    QTime t = QTime::currentTime();
    int seed = t.msecsTo(QTime(0,0,0,0));
    srand(seed);
    foreach (DSBNode *n, nodeMap.keys())
    {
        ShapeObj *sh = n->getShape();
        QPointF c = sh->centrePos();
        double dx = (rand() / static_cast<double>( RAND_MAX ) - 0.5)*scale;
        double dy = (rand() / static_cast<double>( RAND_MAX ) - 0.5)*scale;
        c += QPointF(dx,dy);
        sh->setCentrePos(c);
    }
}

ogdf::Graph *FreePathway::getOGDFGraph(nodemap &nodeMap)
{
    // You should call layout() before this method, so that the reactions
    // have built their orbits.
    ogdf::Graph *G = new ogdf::Graph;
    // nodes
    foreach (DSBClone *cl, m_clones)
    {
        ogdf::node n = G->newNode();
        nodeMap.insert(cl,n);
    }
    foreach (DSBReaction *re, m_reactions)
    {
        ogdf::node n = G->newNode();
        nodeMap.insert(re,n);
    }
    // edges
    foreach (DSBReaction *re, m_reactions)
    {
        ogdf::node reNode = nodeMap.value(re);
        QList<DSBClone*> sats = re->getAllSatellites();
        foreach (DSBClone *cl, sats)
        {
            ogdf::node clNode = nodeMap.value(cl);
            G->newEdge(reNode,clNode);
        }
    }
    return G;
}

void FreePathway::extractPosAndSize(nodemap &nodeMap, ogdf::GraphAttributes &GA)
{
    foreach (ogdf::node v, nodeMap.values())
    {
        DSBNode *n = nodeMap.key(v);
        ShapeObj *sh = n->getShape();
        GA.width(v) = sh->size().width();
        GA.height(v) = sh->size().height();
        GA.x(v) = sh->centrePos().x();
        GA.y(v) = sh->centrePos().y();
    }
}

void FreePathway::injectPositions(nodemap &nodeMap, ogdf::GraphAttributes &GA)
{
    foreach (ogdf::node v, nodeMap.values())
    {
        DSBNode *n = nodeMap.key(v);
        ShapeObj *sh = n->getShape();
        double cx = GA.x(v) + GA.width(v)/2.0;
        double cy = GA.y(v) + GA.height(v)/2.0;
        sh->setCentrePos(QPointF(cx,cy));
    }
}

QList<CanvasItem*> FreePathway::getAllShapes()
{
    QList<CanvasItem*> items;
    ShapeObj *shape = NULL;
    foreach (DSBClone *cl, m_clones)
    {
        shape = cl->getShape();
        if (shape) { items.append(shape); }
    }
    foreach (DSBReaction *re, m_reactions)
    {
        shape = re->getShape();
        if (shape) { items.append(shape); }
    }
    return items;
}

void FreePathway::acceptCanvasBaseAndRelPts(QPointF parentBasePt)
{
    m_relpt = m_basept - parentBasePt;
    foreach (DSBReaction *reac, m_reactions)
    {
        reac->acceptCanvasBaseAndRelPts(m_basept);
    }
}

void FreePathway::drawAt(QPointF r)
{
    m_basept = r;
    foreach (DSBReaction *reac, m_reactions)
    {
        reac->drawRelTo(m_basept);
    }
}

}
