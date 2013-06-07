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

#ifndef DSBSPECIES_H
#define DSBSPECIES_H

#include <QString>
#include <QList>
#include <QMap>

class Species;

namespace dunnart {

class Canvas;
class DSBCompartment;
class DSBReaction;
class DSBClone;

struct DSBCloneAssignment
{
    QList<DSBClone*> reactants;
    QList<DSBClone*> products;
    QList<DSBClone*> modifiers;
};

enum RoleType {ENTERING, EXITING, MODIFYING};

struct Role
{
    Role(RoleType t, DSBReaction *r) : type(t), reaction(r) {}
    RoleType type;
    DSBReaction *reaction;
};

class DSBSpecies
{

public:
    DSBSpecies();
    DSBSpecies(Species *spec);
    QString getName();
    QString getId();
    QString getCompartmentName();
    void addReactionEntered(DSBReaction *reac);
    void addReactionExited(DSBReaction *reac);
    void addReactionModified(DSBReaction *reac);
    void setCanvas(Canvas *canvas);
    Canvas *canvas();
    void setCompartment(DSBCompartment *comp);
    DSBCompartment *getCompartment();
    void setTrivialCloning();
    void setDiscreteCloning();
    void setDiscreteCloningUsingExistingClones();
    void fullyClone(DSBClone *cl);
    QList<DSBClone*> getClones();
    DSBCloneAssignment *getCloneAssignmentByReactionId(QString rid);
    QList<DSBReaction*> getAllReactions(void);
    void mergeClones(DSBClone *clone, QList<DSBClone*> clones);

private:
    int m_nextCloneId;
    Canvas *m_canvas;
    Species *m_sbmlSpecies;
    QString m_name;
    QString m_id;
    QString m_compartmentName;
    DSBCompartment *m_compartment;
    QList<DSBReaction *> m_reactionsEntered;
    QList<DSBReaction *> m_reactionsExited;
    QList<DSBReaction *> m_reactionsModified;
    QList<DSBClone *> m_clones;
    QMap<QString, DSBCloneAssignment*> m_cloneAssignmentsByReactionId;

    void deleteClonesAndAssignments();
    void deleteAssignments();
    void assign(Role r, DSBClone *cl);
    void reassign(Role r, DSBClone *old, DSBClone *replacement);
    void setCloneMarkers();
    DSBClone *makeNewClone();
};

}

#endif // DSBSPECIES_H
