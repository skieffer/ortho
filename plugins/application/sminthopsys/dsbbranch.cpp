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

#include <QtGui>

#include <stdio.h>
#include <assert.h>

#include "dsbbranch.h"
#include "dsbclone.h"
#include "dsbspecies.h"
#include "dsbreaction.h"
#include "dsbpathway.h"
#include "dsbfork.h"
#include "dsbnode.h"

#include "libdunnartcanvas/connector.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/canvasitem.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/guideline.h"
#include "libdunnartcanvas/separation.h"
#include "libdunnartcanvas/distribution.h"
#include "libdunnartcanvas/relationship.h"

namespace dunnart {

DSBNode *DSBBranch::getPredecessor(DSBNode *node)
{
    DSBNode *pred = NULL;
    int i = nodes.indexOf(node);
    if (i>0)
    {
        pred = nodes.at(i-1);
    }
    return pred;
}

DSBNode *DSBBranch::getSuccessor(DSBNode *node)
{
    DSBNode *succ = NULL;
    int i = nodes.indexOf(node);
    if (i>=0 && i+1 < nodes.size())
    {
        succ = nodes.at(i+1);
    }
    return succ;
}

QList<DSBNode*> DSBBranch::getOwnNodes()
{
    int s = nodes.size();
    int n = cycle ? s - 1 : s;
    return nodes.mid(0,n);
}

QSizeF DSBBranch::layout()
{
    QList<DSBNode*> own = getOwnNodes();
    setMainConnections(own);
    // Put lead node at (0,0).
    // So base point of branch is center of lead node.
    qreal x = 0, y = 0, w = 0, h = 0;
    qreal width = 0, height = 0;
    DSBNode *node = own.first();
    node->setRelPt(QPointF(x,y));
    QSizeF size = node->layout();
    h = size.height();
    y += h/2;
    height += h;
    qreal gap = 20;
    // Now place all other nodes according to their size.
    for (int i = 1; i < own.size(); i++)
    {
        y += gap;
        height += gap;
        node = own.at(i);
        QSizeF size = node->layout();
        h = size.height();
        y += h/2;
        node->setRelPt(QPointF(x,y));
        y += h/2;
        height += h;
        w = size.width();
        width = w > width ? w : width;
    }
    m_size = QSizeF(width, height);
    return m_size;
}

void DSBBranch::acceptCanvasBaseAndRelPts(QPointF parentBasePt)
{
    // TODO
}

QRectF DSBBranch::getBbox()
{
    // layout should have already been called
    QRectF rect = QRectF(0,0,0,0);
    QList<DSBNode*> own = getOwnNodes();
    if (!own.empty())
    {
        DSBNode *n = own.first();
        rect = n->getBbox();
        for (int i = 1; i < own.size(); i++)
        {
            n = own.at(i);
            QRectF r = n->getBbox();
            rect = rect.united(r);
        }
    }
    return rect;
}

void DSBBranch::setMainConnections(QList<DSBNode*> own)
{
    for (int i = 0; i < own.size(); i++)
    {
        DSBNode *n = own.at(i);
        DSBReaction *reac = dynamic_cast<DSBReaction*>(n);
        DSBClone *cl = dynamic_cast<DSBClone*>(n);
        if (reac)
        {
            // predecessor
            n = getPredecessor(reac);
            if (n)
            {
                DSBClone *c = dynamic_cast<DSBClone*>(n);
                assert(c);
                reac->setMainInput(c);
            }
            // successor
            n = getSuccessor(reac);
            if (n)
            {
                DSBClone *c = dynamic_cast<DSBClone*>(n);
                assert(c);
                reac->setMainOutput(c);
            }
        }
        else if (cl)
        {
            DSBFork *fork = cl->getFork();
            if (fork)
            {
                // predecessor
                n = getPredecessor(cl);
                if (n)
                {
                    DSBReaction *r = dynamic_cast<DSBReaction*>(n);
                    assert(r);
                    fork->setMainInput(r);
                }
                // successor
                n = getSuccessor(cl);
                if (n)
                {
                    DSBReaction *r = dynamic_cast<DSBReaction*>(n);
                    assert(r);
                    fork->setMainOutput(r);
                }
            }
        }
        else
        {
            qDebug() << "ERROR: Found Node which was neither Clone nor Reaction: " << n;
        }
    }
}

void DSBBranch::setPathway(DSBPathway *pw)
{
    m_pathway = pw;
    for (int i = 0; i < nodes.size(); i++)
    {
        DSBNode *n = nodes.at(i);
        n->setPathway(pw);
    }
}

void DSBBranch::setCanvas(Canvas *canvas)
{
    m_canvas = canvas;
}

DSBPathway *DSBBranch::getPathway()
{
    return m_pathway;
}

QSizeF DSBBranch::getSize()
{
    return m_size;
}

void DSBBranch::setRelPt(QPointF p)
{
    m_relpt = p;
}

void DSBBranch::drawRelTo(QPointF q)
{
    QPointF r = m_relpt + q;
    drawAt(r);
}

void DSBBranch::redraw()
{
    drawAt(m_basept);
}

void DSBBranch::drawAt(QPointF r)
{
    m_basept = r;
    QList<DSBNode*> own = getOwnNodes();
    for (int i = 0; i < own.size(); i++)
    {
        DSBNode *node = own.at(i);
        node->drawRelTo(r);
    }
}

DSBBranch *DSBBranch::computeChordfreeSubbranch()
{
    if (cycle) {return this;}
    if (nodes.size() < 3) {return this;}
    DSBBranch *cfsb = new DSBBranch;
    cfsb->parent = parent;
    int i = 0;
    DSBNode *node = nodes.at(i);
    cfsb->nodes.append(node);
    while (i+1 < nodes.size())
    {
        int j = nodes.size();
        bool connected = false;
        while (!connected && j > i+1)
        {
            j--;
            connected = node->isConnectedTo(nodes.at(j));
        }
        i = j;
        node = nodes.at(i);
        cfsb->nodes.append(node);
    }
    return cfsb;
}

void DSBBranch::align(bool forward)
{
    if (cycle) { return; } // Shouldn't align cycles.
    if (nodes.size() < 2) { return; } // Can't align fewer than 2 nodes.
    m_guideline = NULL;
    /*
    for (int i = 0; i+1 < nodes.size(); i++)
    {
        ShapeObj *shape1 = nodes.at(i)->getShape();
        ShapeObj *shape2 = nodes.at(i+1)->getShape();
        CanvasItemList items;
        items.append(shape1);
        items.append(shape2);

        // Create separation constraint.
        double sep = 50;
        dtype type = SEP_VERTICAL;
        createSeparation(NULL, type, items, sep);

        // If guideline has not yet been created...
        if (!m_guideline)
        {
            // ...create it.
            m_guideline = createAlignment(ALIGN_CENTER, items);
        }
        else
        {
            // ...otherwise simply align the next node to the guideline.
            items.takeFirst();
            items.push_front(m_guideline);
            createAlignment(ALIGN_CENTER, items);
        }
        // Prompt the layout engine to respond.
        //m_canvas->getActions().clear();
        //m_canvas->interrupt_graph_layout();
    }
    */
    CanvasItemList items;
    if (forward)
    {
        foreach (DSBNode *node, nodes)
        {
            ShapeObj *shape = node->getShape();
            items.append(shape);
        }
    }
    else
    {
        for (int i = nodes.size()-1; i >= 0; i--)
        {
            DSBNode *node = nodes.at(i);
            ShapeObj *shape = node->getShape();
            items.append(shape);
        }
    }
    dtype type = DIST_MIDDLE;
    bool preserveOrder = true;
    //Distribution *dist = createDistribution(NULL, type, items, preserveOrder);
    //dist->setSeparation(100);

    m_guideline = createAlignment(ALIGN_CENTER, items);
}


void DSBBranch::setGuideline()
{
    CanvasItemList shapes;
    QList<DSBNode*> own = getOwnNodes();
    for (int i = 0; i < own.size(); i++)
    {
        DSBNode *node = own.at(i);
        ShapeObj *shape = node->getShape();
        shapes.append(shape);
    }
    m_guideline = createAlignment(ALIGN_CENTER, shapes);
    // Make sure position is set correctly.
    DSBNode *first = own.first();
    QPointF pt = first->getBasePt();
    double x0 = (double) pt.x();
    m_guideline->setPosition(x0);
    // Finally, reinitialize the GraphData object, so that the layout thread knows this
    // guideline is there.
    m_canvas->interrupt_graph_layout();
}

#if 0
void DSBBranch::drawConnectors()
{
    ShapeObj *shp1 = 0, *shp2 = 0;
    if (parent && parent != nodes.first())
    {
        shp1 = parent->getShape();
        shp2 = nodes.first()->getShape();
        connect(shp1,shp2);
    }
    for (int i = 0; i+1 < nodes.size(); i++)
    {
        shp1 = nodes.at(i)->getShape();
        shp2 = nodes.at(i+1)->getShape();
        connect(shp1,shp2);
    }
}
#endif

#if 1
void DSBBranch::drawConnectors()
{
    DSBNode *node1, *node2;
    if (parent && parent != nodes.first())
    {
        node1 = parent;
        node2 = nodes.first();
        connect(node1,node2);
    }
    for (int i = 0; i+1 < nodes.size(); i++)
    {
        node1 = nodes.at(i);
        node2 = nodes.at(i+1);
        connect(node1,node2);
    }
}
#endif

#if 1
void DSBBranch::connect(DSBNode *node1, DSBNode *node2)
{
    DSBReaction *reac1 = dynamic_cast<DSBReaction*>(node1);
    DSBReaction *reac2 = dynamic_cast<DSBReaction*>(node2);
    DSBReaction *reac = NULL;
    DSBClone *cl = NULL;
    if (reac1)
    {
        reac = reac1;
        cl = dynamic_cast<DSBClone*>(node2);
    }
    else
    {
        reac = reac2;
        cl = dynamic_cast<DSBClone*>(node1);
    }
    assert(cl); assert(reac);
    reac->connectTo(cl);
}
#endif

#if 0
void DSBBranch::connect(ShapeObj *shp1, ShapeObj *shp2)
{
    Connector *conn = new Connector();
    conn->initWithConnection(shp1,shp2);
    QUndoCommand *cmd = new CmdCanvasSceneAddItem(m_canvas, conn);
    m_canvas->currentUndoMacro()->addCommand(cmd);
}
#endif

QString DSBBranch::toString()
{
    QString s;
    if (cycle)
    {
        s += "\nCycle:------------------------------------------\n";
    } else {
        s += "\nBranch:-----------------------------------------\n";
    }
    if (parent)
    {
        s += "( ";
        DSBClone *cl = dynamic_cast<DSBClone*>(parent);
        DSBReaction *reac = dynamic_cast<DSBReaction*>(parent);
        if (cl)
        {
            s += cl->getSpecies()->getName();
            s += "-"+cl->getCloneId();
        }
        else if (reac)
        {
            s += reac->getReactionId();
        }
        s += " ) ";
    }
    for (int j = 0; j < nodes.size(); j++)
    {
        DSBClone *cl = dynamic_cast<DSBClone*>(nodes.at(j));
        DSBReaction *reac = dynamic_cast<DSBReaction*>(nodes.at(j));
        if (cl)
        {
            s += cl->getSpecies()->getName();
            s += "-"+cl->getCloneId();
        }
        else if (reac)
        {
            s += " --[ " + reac->getReactionId() + " ]--> ";
        }
    }
    return s;
}

}
