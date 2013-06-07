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

#ifndef DSBNODE_H
#define DSBNODE_H

#include <QList>
#include <QSet>
#include <QString>
#include <QSizeF>

#include "dsbreclayout.h"

class QRectF;

namespace dunnart {

class DSBPathway;
class DSBBranch;
class DSBFork;
class ShapeObj;
class DSBClone;
class DSBReaction;

class DSBNode : public DSBRecLayout
{
public:
    DSBNode() : m_branch(NULL), m_pathway(NULL) {}
    virtual QList<DSBBranch*> findBranchesRec(
            QList<QString>& seen, QList<QString> blacklist,
            bool forward, DSBNode *last = 0) = 0;

    QList<DSBBranch*> findBranches(QList<QString> blacklist, bool forward, bool extended);

    QList<DSBBranch*> mergeSelfWithBranches(
            QList<DSBBranch*> branches, QList<QString> blacklist);

    void setBranchHeadNumber(int n);
    virtual ShapeObj *getShape() = 0;
    virtual void moveShape(qreal dx, qreal dy) = 0;
    virtual QPointF getBasePt() = 0;
    virtual QRectF getBbox() = 0;
    void setBranch(DSBBranch *b);
    DSBBranch *getBranch(void);
    void setPathway(DSBPathway *pw);
    DSBPathway *getPathway(void);
    virtual void connectedComponent(QSet<DSBClone*> &ccClones, QSet<DSBReaction*> &ccReacs) = 0;
    bool isConnectedTo(DSBNode *other);
#if 0
    void addBranch(DSBBranch *branch);
    void addFork(DSBFork *fork);
#endif

    static bool s_followTransporters;

private:
    DSBBranch *findMergeTarget(
            QList<DSBBranch*> branches, QList<QString> blacklist);
    int m_branchHeadNumber;
protected:
    DSBBranch *m_branch;
    DSBPathway *m_pathway;
#if 0
    QList<DSBBranch*> m_branches;
    QList<DSBFork*> m_forks;
#endif
};

}

#endif // DSBNODE_H
