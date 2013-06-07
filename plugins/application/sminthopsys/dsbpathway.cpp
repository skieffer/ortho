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
#include <QRectF>

#include <assert.h>

#include "dsbpathway.h"
#include "dsbnode.h"
#include "dsbbranch.h"
#include "dsbclone.h"
#include "dsbreaction.h"
#include "dsbfork.h"

#include "libdunnartcanvas/canvas.h"

namespace dunnart {

DSBPathway::DSBPathway() :
    m_headNode(NULL),
    m_canvas(NULL)
{}

DSBPathway::DSBPathway(DSBNode *head, QList<DSBBranch *> branches) :
    m_headNode(NULL),
    m_canvas(NULL)
{
    // Find branch starting with selected endpt clone.
    QList<DSBBranch*> otherBranches;
    DSBBranch *mainBranch = 0;
    for (int i = 0; i < branches.size(); i++)
    {
        DSBBranch *b = branches.at(i);
        DSBNode *n = b->nodes.first();
        if (n == head) { mainBranch = b; }
        else { otherBranches.append(b); }
    }
    assert(mainBranch!=0);

    setFirstBranch(mainBranch);

    for (int i = 0; i < otherBranches.size(); i++)
    {
        DSBBranch *b = otherBranches.at(i);
        addBranch(b);
    }

}

QList<CanvasItem*> DSBPathway::getAllShapes()
{
    QList<CanvasItem*> items;
    // TODO
    return items;
}

QList<DSBBranch*> DSBPathway::getBranches()
{
    return m_branches;
}

void DSBPathway::setCanvas(Canvas *canvas)
{
    m_canvas = canvas;
    for (int i = 0; i < m_branches.size(); i++)
    {
        DSBBranch *b = m_branches.at(i);
        b->setCanvas(m_canvas);
    }
}

DSBClone *DSBPathway::getLeadClone()
{
    DSBClone *cl = NULL;
    if (m_headNode)
    {
        cl = dynamic_cast<DSBClone*>(m_headNode);
    }
    return cl;
}

void DSBPathway::setFirstBranch(DSBBranch *branch)
{
    m_branches.clear();
    m_allNodes.clear();
    m_branches.append(branch);
    addBranchNodes(branch);
    branch->setPathway(this);
    m_headNode = branch->nodes.first();
}

DSBBranch *DSBPathway::getMainBranch()
{
    assert(!m_branches.empty());
    return m_branches.first();
}

void DSBPathway::addBranchNodes(DSBBranch *branch)
{
    QList<DSBNode*> nodes = branch->getOwnNodes();
    foreach (DSBNode *node, nodes)
    {
        m_allNodes.append(node);
        node->setBranch(branch);
    }
}

void DSBPathway::addBranch(DSBBranch *branch)
{
    // Find parent node.
    DSBNode *parent = branch->parent;
    if (!m_allNodes.contains(parent))
    {
        qDebug() << "ERROR: could not find node\n" << parent << "\nin pathway\n" << this;
        return;
    }

    // If parent is a reaction, simply add the branch to it.
    DSBReaction *reac = dynamic_cast<DSBReaction*>(parent);
    if (reac)
    {
        DSBNode *first = branch->nodes.first();
        DSBClone *cl = dynamic_cast<DSBClone*>(first);
        assert(cl);
        reac->addOutputBranchHead(cl);
        return;
    }

    // Otherwise parent is a clone.
    DSBClone *cl = dynamic_cast<DSBClone*>(parent);
    assert(cl);
    // Get existing fork, or create one.
    DSBFork *fork = cl->getFork();
    if (!fork)
    {
        fork = new DSBFork(cl);
        fork->setPathway(this);
        cl->setFork(fork);
        // Get branch with parent in it.
        DSBBranch *main = cl->getBranch();
        // Get predecessor and successor of parent in its branch, if any,
        // and set them as upstream resp. downstream reactions.
        DSBNode *n = main->getPredecessor(cl);
        if (n) {
            DSBReaction *reac = dynamic_cast<DSBReaction*>(n);
            assert(reac);
            fork->addUpstream(reac);
        }
        n = main->getSuccessor(cl);
        if (n) {
            DSBReaction *reac = dynamic_cast<DSBReaction*>(n);
            assert(reac);
            fork->addDownstream(reac);
        }
    }
    // Add branch head as downstream member of fork.
    DSBNode *n = branch->nodes.first();
    DSBReaction *downstr = dynamic_cast<DSBReaction*>(n);
    assert(downstr);
    fork->addDownstream(downstr);
    // And set main input to first reaction in branch.
    downstr->setMainInput(cl);

    // Add new branch to list.
    m_branches.append(branch);
    addBranchNodes(branch);
    branch->setPathway(this);
}

/* Let p1, p2, ..., pn be the parent nodes of the passed branches
  b1, b2, ..., bn. We construct the map with pi --> bi for each i.
  Then for any node n, map.count(n) will say how many branches have
  node n as parent, and map.values(n) will return the list of them.
  */
QMap<DSBNode*, DSBBranch*> DSBPathway::countBranchPoints(QList<DSBBranch *> branches)
{
    QMap<DSBNode*, DSBBranch*> map;
    for (int i = 0; i < branches.size(); i++)
    {
        DSBBranch *b = branches.at(i);
        DSBNode *p = b->parent;
        if (p)
        {
            map.insertMulti(p,b);
        }
    }
    return map;
}

QSizeF DSBPathway::layout()
{
    // First layout all branches, so that their sizes are available.
    for (int i = 0; i < m_branches.size(); i++)
    {
        DSBBranch *b = m_branches.at(i);
        b->layout();
    }
    // Set relpt of main branch to (0,0).
    DSBBranch *mainBranch = m_headNode->getBranch();
    mainBranch->setRelPt(QPointF(0,0));
    // Now layout all forks.
    // They will set the relpts of all remaining branches.
    foreach (DSBNode *node, m_allNodes)
    {
        DSBClone *cl = dynamic_cast<DSBClone*>(node);
        if (!cl) { continue; }
        DSBFork *fork = cl->getFork();
        if (!fork) { continue; }
        fork->layout();
    }
    // Get size.
    QRectF rect = getBbox();
    m_size = rect.size();
    return m_size;
}

void DSBPathway::acceptCanvasBaseAndRelPts(QPointF parentBasePt)
{
    // TODO
}

QRectF DSBPathway::getBbox()
{
    // layout methods of all branches should have been called first
    QRectF rect = QRectF(0,0,0,0);
    if (!m_branches.empty())
    {
        DSBBranch *b = m_branches.first();
        rect = b->getBbox();
        for (int i = 1; i < m_branches.size(); i++)
        {
            b = m_branches.at(i);
            QRectF r = b->getBbox();
            rect = rect.united(r);
        }
    }
    return rect;
}

QSizeF DSBPathway::getSize()
{
    return m_size;
}

void DSBPathway::setRelPt(QPointF p)
{
    m_relpt = p;
}

void DSBPathway::drawRelTo(QPointF q)
{
    QPointF r = m_relpt + q;
    drawAt(r);
}

void DSBPathway::redraw()
{
    drawAt(m_basept);
}

void DSBPathway::drawAt(QPointF r)
{
    m_basept = r;
    for (int i = 0; i < m_branches.size(); i++)
    {
        DSBBranch *b = m_branches.at(i);
        b->drawRelTo(m_basept);
    }
    // After all branches have been drawn, can ask them to
    // draw their connectors, and add guidelines.
    for (int i = 0; i < m_branches.size(); i++)
    {
        DSBBranch *b = m_branches.at(i);
        b->drawConnectors();
        //b->setGuideline();
    }

}

}




















