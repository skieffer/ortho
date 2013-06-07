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

#ifndef DSBPATHWAY_H
#define DSBPATHWAY_H

#include <QList>
#include <QMap>

#include "dsbreclayout.h"

class QRectF;

namespace dunnart {

class DSBBranch;
class DSBFork;
class DSBNode;
class DSBClone;
class Canvas;
class CanvasItem;

class DSBPathway : public DSBRecLayout
{
public:
    // Constructors
    DSBPathway();
    DSBPathway(DSBNode *head, QList<DSBBranch*> branches);
    // Get and set
    QList<DSBBranch*> getBranches();
    void setCanvas(Canvas *canvas);
    QRectF getBbox();
    DSBClone *getLeadClone();
    DSBBranch *getMainBranch();
    // RecLayout methods
    QSizeF layout();
    void setRelPt(QPointF p);
    void drawRelTo(QPointF q);
    void drawAt(QPointF r);
    void redraw();
    QSizeF getSize();
    void acceptCanvasBaseAndRelPts(QPointF parentBasePt);
    // Other
    QMap<DSBNode*, DSBBranch*> countBranchPoints(QList<DSBBranch*> branches);
    void setFirstBranch(DSBBranch *branch);
    void addBranch(DSBBranch *branch);
    virtual QList<CanvasItem*> getAllShapes();

private:
    QPointF m_relpt;
    QPointF m_basept;
    QSizeF m_size;
    DSBNode *m_headNode;
    Canvas *m_canvas;
    QList<DSBBranch*> m_branches;
    QList<DSBNode*> m_allNodes;

    void addBranchNodes(DSBBranch *branch);
};

}

#endif // DSBPATHWAY_H
