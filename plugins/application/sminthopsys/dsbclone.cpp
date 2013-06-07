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
#include <QSet>
#include <QRectF>

#include "dsbclone.h"
#include "dsbspecies.h"
#include "dsbreaction.h"
#include "dsbbranch.h"
#include "pdepn.h"

#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/pluginshapefactory.h"

namespace dunnart {

DSBClone::DSBClone(DSBSpecies *dsbspec) :
    DSBNode(),
    m_dsbspec(dsbspec),
    m_epn(NULL),
    shapeOnBoard(false),
    m_fork(NULL)
{}

void DSBClone::setCloneNum(int num)
{
    QString snum = QString::number(num);
    QString specId = m_dsbspec->getId();
    m_cloneId = specId + "_clone" + snum;
}

QString DSBClone::getCloneId()
{
    return m_cloneId;
}

void DSBClone::addReactionEntered(DSBReaction *reac)
{
    m_reactionsEntered.append(reac);
}

void DSBClone::addReactionExited(DSBReaction *reac)
{
    m_reactionsExited.append(reac);
}

void DSBClone::addReactionModified(DSBReaction *reac)
{
    m_reactionsModified.append(reac);
}

void DSBClone::setReactionsEntered(QList<DSBReaction *> reacs)
{
    while (!m_reactionsEntered.isEmpty()) { m_reactionsEntered.removeFirst(); }
    for (int i = 0; i < reacs.size(); i++) {
        m_reactionsEntered.append(reacs.at(i));
    }
}

void DSBClone::setReactionsExited(QList<DSBReaction *> reacs)
{
    while (!m_reactionsExited.isEmpty()) { m_reactionsExited.removeFirst(); }
    for (int i = 0; i < reacs.size(); i++) {
        m_reactionsExited.append(reacs.at(i));
    }
}

void DSBClone::setReactionsModified(QList<DSBReaction *> reacs)
{
    while (!m_reactionsModified.isEmpty()) { m_reactionsModified.removeFirst(); }
    for (int i = 0; i < reacs.size(); i++) {
        m_reactionsModified.append(reacs.at(i));
    }
}

void DSBClone::deleteShape()
{
    Canvas *canvas = m_dsbspec->canvas();
    /*
    canvas->deselectAll();
    m_epn->setSelected(true);
    canvas->deleteSelection();
    */
    //m_epn->setSelected(true);
    //canvas->removeItem(m_epn);
    canvas->deleteItem(m_epn);
    m_epn = NULL;
    //m_epn->setSelected(true);
}

QList<Role> DSBClone::getAllRoles()
{
    QList<Role> roles;
    foreach(DSBReaction *r,m_reactionsEntered){roles.append(Role(ENTERING,r));}
    foreach(DSBReaction *r,m_reactionsExited){roles.append(Role(EXITING,r));}
    foreach(DSBReaction *r,m_reactionsModified){roles.append(Role(MODIFYING,r));}
    return roles;
}

void DSBClone::clearRoles()
{
    m_reactionsEntered.clear();
    m_reactionsExited.clear();
    m_reactionsModified.clear();
}

void DSBClone::assign(Role r)
{
    switch(r.type)
    {
    case ENTERING:
        m_reactionsEntered.append(r.reaction);
        break;
    case EXITING:
        m_reactionsExited.append(r.reaction);
        break;
    case MODIFYING:
        m_reactionsModified.append(r.reaction);
        break;
    }
}

DSBSpecies *DSBClone::getSpecies()
{
    return m_dsbspec;
}

void DSBClone::set_is_cloned(bool b)
{
    m_is_cloned = b;
}

QSizeF DSBClone::layout()
{
    // Create shape if don't already have one.
    if (!m_epn)
    {
        PluginShapeFactory *factory = sharedPluginShapeFactory();
        QString type("org.sbgn.pd.00UnspecifiedEPN");
        ShapeObj *shape = factory->createShape(type);
        m_epn = dynamic_cast<PDEPN*>  (shape);
        // Give the epn a pointer to the clone it represents.
        m_epn->setClone(this);
    }
    m_size = m_epn->size();
    return m_size;
}

QSizeF DSBClone::getSize()
{
    return m_size;
}

void DSBClone::setRelPt(QPointF p)
{
    m_relpt = p;
}

void DSBClone::drawRelTo(QPointF q)
{
    QPointF r = m_relpt + q;
    drawAt(r);
}

void DSBClone::redraw()
{
    drawAt(m_basept);
}

void DSBClone::acceptCanvasBaseAndRelPts(QPointF parentBasePt)
{
    if (m_epn)
    {
        m_basept = m_epn->centrePos();
    }
    m_relpt = m_basept - parentBasePt;
}

void DSBClone::drawAt(QPointF r)
{
    m_basept = r;
    // Set shape's properties.
    // Clone marker state
    m_epn->set_is_cloned(m_is_cloned);
    // Size
    m_epn->setSize(m_size);
    // Position
    m_epn->setCentrePos(m_basept);
    // Label
    m_epn->setLabel(m_dsbspec->getName());
    // Add it to the canvas, if necessary.
    if (!shapeOnBoard)
    {
        QUndoCommand *cmd = new CmdCanvasSceneAddItem(m_dsbspec->canvas(), m_epn);
        m_dsbspec->canvas()->currentUndoMacro()->addCommand(cmd);
        shapeOnBoard = true;
    } else {
        // Otherwise, recompute the layout, in case size or position
        // have changed.
        m_dsbspec->canvas()->interrupt_graph_layout();
    }
}

QRectF DSBClone::getBbox()
{
    // layout should have already been called, and relpt set
    QPointF centre = m_relpt;
    QSizeF size = m_size;
    QPointF ulc = QPointF( centre.x()-size.width()/2,
                           centre.y()-size.height()/2);
    return QRectF(ulc,size);
}

ShapeObj *DSBClone::getShape()
{
    return m_epn;
}

void DSBClone::setFork(DSBFork *fork)
{
    m_fork = fork;
}

DSBFork *DSBClone::getFork()
{
    return m_fork;
}

void DSBClone::moveShape(qreal dx, qreal dy)
{
    if (m_epn)
    {
        m_epn->moveBy(dx,dy);
    }
}

QPointF DSBClone::getBasePt()
{
    return m_basept;
}

/* Both reactions entered, and reversible reactions exited, are
  enterable. Compute list of all those.
  */
QList<DSBReaction*> DSBClone::computeEnterableReactions()
{
    QSet<DSBReaction*> enterable;
    // Get reactions entered.
    for (int i = 0; i < m_reactionsEntered.size(); i++)
    {
        DSBReaction *reac = m_reactionsEntered.at(i);
        enterable.insert(reac);
    }
    // Examine reactions exited.
    for (int i = 0; i < m_reactionsExited.size(); i++)
    {
        DSBReaction *reac = m_reactionsExited.at(i);
        // We can only use it if it is reversible.
        if (reac->isReversible()) { enterable.insert(reac); }
    }
    return enterable.toList();
}

/* Both reactions exited, and reversible reactions entered, are
  exitable. Compute list of all those.
  */
QList<DSBReaction*> DSBClone::computeExitableReactions()
{
    QSet<DSBReaction*> exitable;
    // Get reactions exited.
    for (int i = 0; i < m_reactionsExited.size(); i++)
    {
        DSBReaction *reac = m_reactionsExited.at(i);
        exitable.insert(reac);
    }
    // Examine reactions entered.
    for (int i = 0; i < m_reactionsEntered.size(); i++)
    {
        DSBReaction *reac = m_reactionsEntered.at(i);
        // We can only use it if it is reversible.
        if (reac->isReversible()) { exitable.insert(reac); }
    }
    return exitable.toList();
}

QList<DSBBranch*> DSBClone::findBranchesRec(
        QList<QString>& seen, QList<QString> blacklist, bool forward, DSBNode *last)
{
    seen.append(m_cloneId); // Mark self as seen.

    QList<DSBBranch*> branches; // Prepare return value.

    // Now consider all usable reactions.
    QList<DSBReaction*> usable =
            forward ? computeEnterableReactions() : computeExitableReactions();
    for (int i = 0; i < usable.size(); i++)
    {
        DSBReaction *reac = usable.at(i);

        // Do not turn around and go backwards.
        if (reac == last) {continue;}

        // Are we avoiding transporter processes?
        if (!DSBNode::s_followTransporters && reac->isIntercompartmental()) {
            continue;
        }

        // Consider whether this reaction has already been seen or not.
        QString rid = reac->getReactionId();
        if (seen.contains(rid))
        {
            // Reaction has already been seen, so we have found a cycle.
            DSBBranch *b = new DSBBranch;
            b->nodes.append(reac);
            b->cycle = true;
            branches.append(b);
        }
        else
        {
            // No cycle. Recurse.
            QList<DSBBranch*> bb = reac->findBranchesRec(seen, blacklist, forward, this);
            branches.append(bb);
        }
    }
    return mergeSelfWithBranches(branches, blacklist);
}

QList<DSBReaction*> DSBClone::getAllReactions()
{
    QList<DSBReaction*> all;
    all.append(m_reactionsEntered);
    all.append(m_reactionsExited);
    all.append(m_reactionsModified);
    return all;
}

void DSBClone::connectedComponent(QSet<DSBClone *> &ccClones, QSet<DSBReaction *> &ccReacs)
{
    ccClones.insert(this); // add self to component
    // ignore reactions already seen
    QSet<DSBReaction*> rset = getAllReactions().toSet().subtract(ccReacs);
    // ask remaining ones to add to the connected component
    foreach (DSBReaction *reac, rset)
    {
        // skip intercompartmental reactions
        if (reac->isIntercompartmental()) { continue; }
        reac->connectedComponent(ccClones, ccReacs);
    }
}

}





















