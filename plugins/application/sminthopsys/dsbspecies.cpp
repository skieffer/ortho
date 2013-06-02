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

#include "dsbspecies.h"
#include "dsbreaction.h"
#include "dsbclone.h"

#include "libdunnartcanvas/canvas.h"

#include "sbml/SBMLTypes.h"

namespace dunnart {

DSBSpecies::DSBSpecies() : m_nextCloneId(0) {}

DSBSpecies::DSBSpecies(Species *spec) :
    m_nextCloneId(0),
    m_sbmlSpecies(spec)
{
    m_name = QString(spec->getName().c_str());
    m_id = QString(spec->getId().c_str());
    m_compartmentName = QString(spec->getCompartment().c_str());
}

void DSBSpecies::setCanvas(Canvas *canvas)
{
    m_canvas = canvas;
}

Canvas *DSBSpecies::canvas()
{
    return m_canvas;
}

void DSBSpecies::setCompartment(DSBCompartment *comp)
{
    m_compartment = comp;
}

DSBCompartment *DSBSpecies::getCompartment()
{
    return m_compartment;
}

QString DSBSpecies::getCompartmentName()
{
    return m_compartmentName;
}

QString DSBSpecies::getName()
{
    return m_name;
}

QString DSBSpecies::getId()
{
    return m_id;
}

void DSBSpecies::addReactionEntered(DSBReaction *reac)
{
    m_reactionsEntered.append(reac);
}

void DSBSpecies::addReactionExited(DSBReaction *reac)
{
    m_reactionsExited.append(reac);
}

void DSBSpecies::addReactionModified(DSBReaction *reac)
{
    m_reactionsModified.append(reac);
}

QList<DSBClone*> DSBSpecies::getClones()
{
    return m_clones;
}

/* This map must be maintained by all methods that change the
   cloning.

   It associates with each reaction (ID) a list of clones of the present
   species which have been assigned to participate in that reaction, as
   reactant, product, or modifier.

   (SBGN PD diagrams do sometimes feature multiple clones of one
    and the same species playing one and the same role in a reaction,
    be it as reactant or product (or maybe even modifier?).)
  */
DSBCloneAssignment *DSBSpecies::getCloneAssignmentByReactionId(QString rid)
{
    return m_cloneAssignmentsByReactionId.value(rid);
}

void DSBSpecies::deleteClonesAndAssignments()
{
    // Delete shapes from the canvas.
    for (int i = 0; i < m_clones.size(); i++)
    {
        m_clones.at(i)->deleteShape();
    }

    // Delete the clone assignments.
    deleteAssignments();

    // Clear the list of clones, and delete the DSBClone objects themselves.
    while (!m_clones.isEmpty())
    {
        DSBClone *cl = m_clones.takeFirst(); // removes it from the list
        delete cl;
    }
    // Reset nextCloneId to 0.
    m_nextCloneId = 0;
}

void DSBSpecies::deleteAssignments()
{
    // Clear the clone assignment map, and delete the DSBCloneAssignment
    // structs themselves.
    QList<QString> reacIds = m_cloneAssignmentsByReactionId.keys();
    for (int i = 0; i < reacIds.size(); i++)
    {
        QString id = reacIds.at(i);
        DSBCloneAssignment *ca = m_cloneAssignmentsByReactionId.value(id);
        m_cloneAssignmentsByReactionId.remove(id);
        delete ca;
    }
}

/*  Allocate a new DSBClone object, cloning this species,
    give it the next available clone id, add it to the
    list of clones, and return a pointer to it.
  */
DSBClone *DSBSpecies::makeNewClone()
{
    int idNum = m_nextCloneId++;
    DSBClone *cl = new DSBClone(this);
    cl->setCloneNum(idNum);
    m_clones.append(cl);
    return cl;
}

/* The trivial cloning is that in which there is precisely one clone
   of this species, participating in all reactions in which this
   species participates.
  */
void DSBSpecies::setTrivialCloning()
{
    // Delete old clones.
    deleteClonesAndAssignments();

    // Create new clone.
    //DSBClone *cl = new DSBClone(this);
    //m_clones.append(cl);
    DSBClone *cl = makeNewClone();

    // Set it as the sole clone assigned to each reaction.
    for (int i = 0; i < m_reactionsEntered.size(); i++)
    {
        DSBReaction *reac = m_reactionsEntered.at(i);
        cl->addReactionEntered(reac);
        QString rid = reac->getReactionId();
        DSBCloneAssignment *ca =
                m_cloneAssignmentsByReactionId.value( rid, new DSBCloneAssignment() );
        ca->reactants.append(cl);
        m_cloneAssignmentsByReactionId.insert(rid,ca);
    }

    for (int i = 0; i < m_reactionsExited.size(); i++)
    {
        DSBReaction *reac = m_reactionsExited.at(i);
        cl->addReactionExited(reac);
        QString rid = reac->getReactionId();
        DSBCloneAssignment *ca =
                m_cloneAssignmentsByReactionId.value( rid, new DSBCloneAssignment() );
        ca->products.append(cl);
        m_cloneAssignmentsByReactionId.insert(rid,ca);
    }

    for (int i = 0; i < m_reactionsModified.size(); i++)
    {
        DSBReaction *reac = m_reactionsModified.at(i);
        cl->addReactionModified(reac);
        QString rid = reac->getReactionId();
        DSBCloneAssignment *ca =
                m_cloneAssignmentsByReactionId.value( rid, new DSBCloneAssignment() );
        ca->modifiers.append(cl);
        m_cloneAssignmentsByReactionId.insert(rid,ca);
    }

    setCloneMarkers();
}

/* The discrete cloning (named after "the discrete topology") is that
   in which there is one clone for each role played by the species.
  */
void DSBSpecies::setDiscreteCloning()
{
    // Delete old clones.
    deleteClonesAndAssignments();

    // Create a new clone for each role played by this species.
    for (int i = 0; i < m_reactionsEntered.size(); i++)
    {
        DSBClone *cl = makeNewClone();
        DSBReaction *reac = m_reactionsEntered.at(i);
        cl->addReactionEntered(reac);
        QString rid = reac->getReactionId();
        DSBCloneAssignment *ca =
                m_cloneAssignmentsByReactionId.value( rid, new DSBCloneAssignment() );
        ca->reactants.append(cl);
        m_cloneAssignmentsByReactionId.insert(rid,ca);
    }

    for (int i = 0; i < m_reactionsExited.size(); i++)
    {
        DSBClone *cl = makeNewClone();
        DSBReaction *reac = m_reactionsExited.at(i);
        cl->addReactionExited(reac);
        QString rid = reac->getReactionId();
        DSBCloneAssignment *ca =
                m_cloneAssignmentsByReactionId.value( rid, new DSBCloneAssignment() );
        ca->products.append(cl);
        m_cloneAssignmentsByReactionId.insert(rid,ca);
    }

    for (int i = 0; i < m_reactionsModified.size(); i++)
    {
        DSBClone *cl = makeNewClone();
        DSBReaction *reac = m_reactionsModified.at(i);
        cl->addReactionModified(reac);
        QString rid = reac->getReactionId();
        DSBCloneAssignment *ca =
                m_cloneAssignmentsByReactionId.value( rid, new DSBCloneAssignment() );
        ca->modifiers.append(cl);
        m_cloneAssignmentsByReactionId.insert(rid,ca);
    }

    setCloneMarkers();
}

void DSBSpecies::mergeClones(DSBClone *clone, QList<DSBClone *> clones)
{
    // Merge each clone in the list into the singleton.
    foreach (DSBClone *cl, clones)
    {
        QList<Role> roles = cl->getAllRoles();
        foreach (Role role, roles)
        {
            reassign(role,cl,clone);
        }
        m_clones.removeAll(cl);
        cl->deleteShape();
        delete cl;
    }
    setCloneMarkers();
    // Ask the reactions with which this species interacts to clear their connectors,
    // and recompute their orbits.
    QList<DSBReaction*> reacs = getAllReactions();
    foreach (DSBReaction *reac, reacs)
    {
        reac->clearConnectors();
        reac->buildOrbit();
    }
}

void DSBSpecies::fullyClone(DSBClone *cl)
{
    // Save clone's list of reactions now, before it gets changed.
    QList<DSBReaction*> reacs = cl->getAllReactions();
    // Get its list of roles.
    QList<Role> roles = cl->getAllRoles();
    // Clear roles in existing clone, then assign it just the first one.
    cl->clearRoles();
    Role firstRole = roles.first();
    reassign(firstRole,cl,cl);
    // Assign each remaining role to a new clone.
    // Also layout the new clones, and draw them at same place as primary one.
    QPointF bp = cl->getBasePt();
    QPointF rp = cl->m_relpt;
    for (int i = 1; i < roles.size(); i++)
    {
        Role role = roles.at(i);
        DSBClone *c = makeNewClone();
        reassign(role,cl,c);
        c->layout();
        c->m_relpt = rp;
        c->drawAt(bp);
    }
    setCloneMarkers();
    // Ask the reactions with which this clone interacts to clear their connectors,
    // and recompute their orbits.
    foreach (DSBReaction *reac, reacs)
    {
        reac->clearConnectors();
        reac->buildOrbit();
    }
}

void DSBSpecies::setDiscreteCloningUsingExistingClones()
{
    // We assume every reaction in which this species participates already
    // has at least one clone of this species assigned to it.

    // Clear out the clone assignments.
    deleteAssignments();

    // Get copy of clone list, since it will be modified.
    QList<DSBClone*> clones = m_clones;
    foreach (DSBClone *cl, clones)
    {
        // Make list of all roles.

        QList<Role> roles = cl->getAllRoles();
        /*
        QList<Role> roles;
        foreach(DSBReaction *r,cl->m_reactionsEntered){roles.append(Role(ENTERING,r));}
        foreach(DSBReaction *r,cl->m_reactionsExited){roles.append(Role(EXITING,r));}
        foreach(DSBReaction *r,cl->m_reactionsModified){roles.append(Role(MODIFYING,r));}
        */

        // Clear roles in existing clone, then assign it just the first one.
        cl->clearRoles();
        assign(roles.first(),cl);
        // Assign each remaining role to a new clone.
        // Also layout the new clones, and draw them at same place as primary one.
        QPointF bp = cl->getBasePt();
        QPointF rp = cl->m_relpt;
        for (int i = 1; i < roles.size(); i++)
        {
            Role role = roles.at(i);
            DSBClone *c = makeNewClone();
            assign(role,c);
            c->layout();
            c->m_relpt = rp;
            c->drawAt(bp);
        }
    }
    setCloneMarkers();
    // Ask the reactions with which this species interacts to clear their connectors,
    // and recompute their orbits.
    QList<DSBReaction*> reacs = getAllReactions();
    foreach (DSBReaction *reac, reacs)
    {
        reac->clearConnectors();
        reac->buildOrbit();
    }
}

QList<DSBReaction*> DSBSpecies::getAllReactions()
{
    QList<DSBReaction*> all;
    all.append(m_reactionsEntered);
    all.append(m_reactionsExited);
    all.append(m_reactionsModified);
    return all;
}

void DSBSpecies::assign(Role r, DSBClone *cl)
{
    QString rid = r.reaction->getReactionId();
    DSBCloneAssignment *cla = m_cloneAssignmentsByReactionId.value(
                rid, new DSBCloneAssignment() );
    switch(r.type)
    {
    case ENTERING:
        cl->m_reactionsEntered.append(r.reaction);
        cla->reactants.append(cl);
        break;
    case EXITING:
        cl->m_reactionsExited.append(r.reaction);
        cla->products.append(cl);
        break;
    case MODIFYING:
        cl->m_reactionsModified.append(r.reaction);
        cla->modifiers.append(cl);
        break;
    }
    m_cloneAssignmentsByReactionId.insert(rid,cla);
}

void DSBSpecies::reassign(Role r, DSBClone *old, DSBClone *replacement)
{
    QString rid = r.reaction->getReactionId();
    DSBCloneAssignment *cla = m_cloneAssignmentsByReactionId.value(
                rid, new DSBCloneAssignment() );
    QList<DSBClone*> *cloneList;
    switch(r.type)
    {
    case ENTERING:
        replacement->m_reactionsEntered.append(r.reaction);
        cloneList = &cla->reactants;
        break;
    case EXITING:
        replacement->m_reactionsExited.append(r.reaction);
        cloneList = &cla->products;
        break;
    case MODIFYING:
        replacement->m_reactionsModified.append(r.reaction);
        cloneList = &cla->modifiers;
        break;
    }
    int i = (*cloneList).indexOf(old);
    if (i < 0) { (*cloneList).append(replacement); }
    else { (*cloneList).replace(i,replacement); }
    m_cloneAssignmentsByReactionId.insert(rid,cla);
}

void DSBSpecies::setCloneMarkers()
{
    if (m_clones.size() > 1)
    {
        for (int i = 0; i < m_clones.size(); i++)
        {
            m_clones.at(i)->set_is_cloned(true);
        }
    }
    else if (m_clones.size() == 1)
    {
        m_clones.at(0)->set_is_cloned(false);
    }
}

}
