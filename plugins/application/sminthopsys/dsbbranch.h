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

#ifndef DSBBRANCH_H
#define DSBBRANCH_H

#include <QList>
#include <QString>

#include "dsbreclayout.h"

class QRectF;

namespace dunnart {

class DSBNode;
class DSBPathway;
class Canvas;
class ShapeObj;
class Guideline;

class DSBBranch : public DSBRecLayout
{
public:
    // public fields (holdover from implementation as struct)
    bool cycle;
    QList<DSBNode*> nodes;
    DSBNode *parent;
    // Constructor
    DSBBranch() : cycle(false), parent(0), m_pathway(NULL) {}
    // Get and set
    DSBNode *getPredecessor(DSBNode *node);
    DSBNode *getSuccessor(DSBNode *node);
    void setPathway(DSBPathway *pw);
    DSBPathway *getPathway();
    void setCanvas(Canvas *canvas);
    QRectF getBbox();
    QList<DSBNode*> getOwnNodes();
    // RecLayout methods
    QSizeF layout();
    void setRelPt(QPointF p);
    void drawRelTo(QPointF q);
    void drawAt(QPointF r);
    void drawConnectors();
    void setGuideline();
    void redraw();
    QSizeF getSize();
    void acceptCanvasBaseAndRelPts(QPointF parentBasePt);
    // Other methods
    QString toString();
    void align(bool forward);
    DSBBranch *computeChordfreeSubbranch(void);
    // Friends
    friend class DSBFork;

private:
    QPointF m_relpt;
    QPointF m_basept;
    QSizeF m_size;
    Canvas *m_canvas;
    Guideline *m_guideline;

    DSBPathway *m_pathway;

    void setMainConnections(QList<DSBNode*> own);
    //void connect(ShapeObj *shp1, ShapeObj *shp2);
    void connect(DSBNode *node1, DSBNode *node2);
};

}

#endif // DSBBRANCH_H
