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
#include <QRect>

#include "dsbreaction.h"
#include "dsbspecies.h"
#include "dsbclone.h"
#include "dsbbranch.h"

#include "sbml/SBMLTypes.h"

#include "libdunnartcanvas/pluginshapefactory.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/undo.h"
#include "libdunnartcanvas/connector.h"

namespace dunnart {

DSBReaction::DSBReaction() :
    DSBNode(),
    m_sbmlReaction(NULL),
    m_canvas(NULL),
    m_compartment(NULL),
    m_mainInput(NULL),
    m_mainOutput(NULL),
    m_shape(NULL),
    shapeOnCanvas(false)
{}

DSBReaction::DSBReaction(Reaction *reac) :
    DSBNode(),
    m_sbmlReaction(reac),
    m_canvas(NULL),
    m_compartment(NULL),
    m_mainInput(NULL),
    m_mainOutput(NULL),
    m_shape(NULL),
    shapeOnCanvas(false)
{
    m_name = QString(reac->getName().c_str());
    m_id = QString(reac->getId().c_str());
    m_compartmentName = QString(reac->getCompartment().c_str());
    m_reversible = false; // Default to not reversible.
    if (reac->isSetReversible()) // If reversibility has been stated,
    {
        m_reversible = reac->getReversible(); // then accept the stated value.
    }
}

void DSBReaction::setCompartment(DSBCompartment *comp)
{
    m_compartment = comp;
}

void DSBReaction::setCanvas(Canvas *canvas)
{
    m_canvas = canvas;
}

DSBCompartment *DSBReaction::getCompartment()
{
    return m_compartment;
}

QString DSBReaction::getCompartmentName()
{
    return m_compartmentName;
}

QString DSBReaction::getName()
{
    return m_name;
}

QString DSBReaction::getReactionId()
{
    return m_id;
}

void DSBReaction::addInputBranchHead(DSBClone *head)
{
    m_inputBranchHeads.append(head);
}

void DSBReaction::addOutputBranchHead(DSBClone *head)
{
    m_outputBranchHeads.append(head);
}

bool DSBReaction::isReversible()
{
    return m_reversible;
}

/* Check whether this reaction spans multiple compartments.
  */
bool DSBReaction::isIntercompartmental()
{
    QSet<QString> comps;
    for (int i = 0; i < m_inputs.size(); i++)
    {
        DSBSpecies *spec = m_inputs.at(i);
        QString comp =spec->getCompartmentName();
        comps.insert(comp);
    }
    for (int i = 0; i < m_outputs.size(); i++)
    {
        DSBSpecies *spec = m_outputs.at(i);
        QString comp =spec->getCompartmentName();
        comps.insert(comp);
    }
    for (int i = 0; i < m_modifiers.size(); i++)
    {
        DSBSpecies *spec = m_modifiers.at(i);
        QString comp =spec->getCompartmentName();
        comps.insert(comp);
    }
    return (comps.size() > 1);
}

bool DSBReaction::hasCloneAsInputOrOutput(DSBClone *cl)
{
    // Orbit should have been built before calling this!
    return (m_inputBranchHeads.contains(cl) || m_outputBranchHeads.contains(cl) ||
            m_inSatellites.contains(cl) || m_outSatellites.contains(cl) ||
            m_mainInput == cl || m_mainOutput == cl);
}

/* Give this reaction links to all species involved in it,
   and give those species links to this reaction.
  */
void DSBReaction::doublyLink(QMap<QString, DSBSpecies*> map)
{
    ListOfSpeciesReferences *lsr;
    SimpleSpeciesReference *ssr;
    unsigned int N;
    QString specId;

    // "reactants", or inputs
    lsr = m_sbmlReaction->getListOfReactants();
    N = lsr->size();
    for (unsigned int i = 0; i < N; i++)
    {
        ssr = lsr->get(i);
        specId = QString(ssr->getSpecies().c_str());
        if (!map.contains(specId))
        {
            // TODO: Report error. Reaction is referring to a species that
            // was not declared in the SBML list of species.
            qDebug() << "map does not contain species id " << specId;
        }
        else
        {
            DSBSpecies *dsbspec = map.value(specId);
            m_inputs.append(dsbspec);
            dsbspec->addReactionEntered(this);
        }
    }

    // "products", or outputs
    lsr = m_sbmlReaction->getListOfProducts();
    N = lsr->size();
    for (unsigned int i = 0; i < N; i++)
    {
        ssr = lsr->get(i);
        specId = QString(ssr->getSpecies().c_str());
        if (!map.contains(specId))
        {
            // TODO: Report error. Reaction is referring to a species that
            // was not declared in the SBML list of species.
            qDebug() << "map does not contain species id " << specId;
        }
        else
        {
            DSBSpecies *dsbspec = map.value(specId);
            m_outputs.append(dsbspec);
            dsbspec->addReactionExited(this);
        }
    }

    // modifiers (e.g. catalysts)
    lsr = m_sbmlReaction->getListOfModifiers();
    N = lsr->size();
    for (unsigned int i = 0; i < N; i++)
    {
        ssr = lsr->get(i);
        specId = QString(ssr->getSpecies().c_str());
        if (!map.contains(specId))
        {
            // TODO: Report error. Reaction is referring to a species that
            // was not declared in the SBML list of species.
            qDebug() << "map does not contain species id " << specId;
        }
        else
        {
            DSBSpecies *dsbspec = map.value(specId);
            m_modifiers.append(dsbspec);
            dsbspec->addReactionModified(this);
        }
    }
}

QList<DSBClone*> DSBReaction::getInputClones()
{
    QList<DSBClone*> clones;
    for (int i = 0; i < m_inputs.size(); i++)
    {
        DSBSpecies *spec = m_inputs.at(i);
        DSBCloneAssignment *cla = spec->getCloneAssignmentByReactionId(m_id);
        clones.append(cla->reactants);
    }
    return clones;
}

QList<DSBClone*> DSBReaction::getOutputClones()
{
    QList<DSBClone*> clones;
    for (int i = 0; i < m_outputs.size(); i++)
    {
        DSBSpecies *spec = m_outputs.at(i);
        DSBCloneAssignment *cla = spec->getCloneAssignmentByReactionId(m_id);
        clones.append(cla->products);
    }
    return clones;
}

QList<DSBClone*> DSBReaction::getOpposedClones(DSBClone *clone)
{
    QList<DSBClone*> opp; // Prepare return value.

    // First determine which side, or sides, the clone is on.
    DSBSpecies *spec = clone->getSpecies();
    DSBCloneAssignment *cla = spec->getCloneAssignmentByReactionId(m_id);
    bool isProd = cla->products.contains(clone);
    bool isReac = cla->reactants.contains(clone);

    // If it is a product and not a reactant, then the reactants are opposed.
    if (isProd && !isReac)
    {
        opp = getInputClones();
    }
    // If it is a reactant and not a product, then the products are opposed.
    else if (isReac && !isProd)
    {
        opp = getOutputClones();
    }
    // Otherwise we cannot call anything opposed.

    return opp;
}

QList<DSBBranch*> DSBReaction::findBranchesRec(
        QList<QString>& seen, QList<QString> blacklist, bool forward, DSBNode *last)
{
    seen.append(m_id); // Mark self as seen.

    QList<DSBBranch*> branches; // Prepare return value.

    // Check which side of this reaction the last node 'last' lies on.
    // Then only consider flowing out on opposite side.
    DSBClone *clone = dynamic_cast<DSBClone*>(last);
    QList<DSBClone*> opp = getOpposedClones(clone);

    for (int i = 0; i < opp.size(); i++)
    {
        DSBClone *cl = opp.at(i);

        // Do not turn around and go backwards.
        if (cl == last) {continue;}

        // Consider whether this clone has already been seen or not.
        QString cid = cl->getCloneId();

        if (seen.contains(cid))
        {
            // Clone has already been seen, so we have found a cycle.
            DSBBranch *b = new DSBBranch;
            b->nodes.append(cl);
            b->cycle = true;
            branches.append(b);
        }
        else
        {
            // No cycle. Recurse.
            QList<DSBBranch*> bb = cl->findBranchesRec(seen, blacklist, forward, this);
            branches.append(bb);
        }
    }

    return mergeSelfWithBranches(branches, blacklist);
}

void DSBReaction::buildOrbit()
{
    m_inSatellites.clear();
    m_outSatellites.clear();
    m_modSatellites.clear();
    QList<DSBSpecies*> allSpecies = getAllSpecies();
    foreach (DSBSpecies *spec, allSpecies)
    {
        if (spec->getName()=="melibiose") {
            qDebug()<<"foo";
        }
        DSBCloneAssignment *cla = spec->getCloneAssignmentByReactionId(m_id);
        takeNonBranchHeads(cla->reactants, m_inSatellites);
        takeNonBranchHeads(cla->products, m_outSatellites);
        takeNonBranchHeads(cla->modifiers, m_modSatellites);
    }
}

/* For every clone in src list which is not currently registered as a branch head,
   or as mainInput or mainOutput, append it to the dst list.
   Also, mark all clones appended to the dst list as belonging to the
   same pathway as this reaction.
  */
void DSBReaction::takeNonBranchHeads(QList<DSBClone *> &src, QList<DSBClone *> &dst)
{
    for (int i = 0; i < src.size(); i++)
    {
        DSBClone *cl = src.at(i);
        if (!m_inputBranchHeads.contains(cl) && !m_outputBranchHeads.contains(cl) &&
            m_mainInput != cl && m_mainOutput != cl)
        {
            dst.append(cl);
            cl->setPathway(m_pathway);
        }
    }
}

void DSBReaction::setMainInput(DSBClone *mi)
{
    m_mainInput = mi;
}

void DSBReaction::setMainOutput(DSBClone *mo)
{
    m_mainOutput = mo;
}

QList<DSBSpecies*> DSBReaction::getAllSpecies()
{
    QList<DSBSpecies*> allSpecies;
    allSpecies.append(m_inputs);
    allSpecies.append(m_outputs);
    allSpecies.append(m_modifiers);
    return allSpecies;
}

QList<DSBClone*> DSBReaction::getAllSatellites()
{
    QList<DSBClone*> allSats;
    allSats.append(m_inSatellites);
    allSats.append(m_outSatellites);
    allSats.append(m_modSatellites);
    return allSats;
}

QList<DSBClone*> DSBReaction::getAllClones()
{
    QList<DSBClone*> clones;
    QList<DSBSpecies*> species = getAllSpecies();
    foreach (DSBSpecies *spec, species)
    {
        DSBCloneAssignment *cla = spec->getCloneAssignmentByReactionId(m_id);
        clones.append(cla->reactants);
        clones.append(cla->modifiers);
        clones.append(cla->products);
    }
    return clones;
}

bool DSBReaction::isBranchHead(DSBClone *clone)
{
    return m_inputBranchHeads.contains(clone) || m_outputBranchHeads.contains(clone);
}

QPointF DSBReaction::satPos(int num, int outOf, ReacSide side)
{
    qreal x,y;
    qreal stem = 25;
    qreal conn = 70;
    qreal sqrt2 = 1.414213562;
    qreal sqrt3 = 1.732050808;
    if (outOf == 1 || outOf == 2)
    {
        // Use 45 deg angles
        y = stem + conn/sqrt2;
        x = conn/sqrt2;
        if (num == 1) { x *= -1; }
    }
    else if (outOf == 3 || outOf == 4)
    {
        // Use 30 and 60 deg angles
        qreal A = conn*sqrt3/2;
        qreal a = conn/2;
        if (num == 1 || num == 2) {
            x = a;
            y = stem + A;
        } else {
            x = A;
            y = stem + a;
        }
        if (num == 1 || num == 3) { x *= -1; }
    }
    else
    {
        // This case (more than 4 clones on one side) will be extremely rare.
        // TODO...
    }
    if (side == ABOVE) { y *= -1; }
    return QPointF(x,y);
}

void DSBReaction::acceptCanvasBaseAndRelPts(QPointF parentBasePt)
{
    if (m_shape)
    {
        m_basept = m_shape->centrePos();
    }
    m_relpt = m_basept - parentBasePt;
    QList<DSBClone*> allSats = getAllSatellites();
    foreach (DSBClone *cl, allSats)
    {
        cl->acceptCanvasBaseAndRelPts(m_basept);
    }
}

QSizeF DSBReaction::layout()
{
    clearConnectors();
    buildOrbit();
    // Create a shape, if don't already have one.
    if (!m_shape)
    {
        PluginShapeFactory *factory = sharedPluginShapeFactory();
        QString type("org.sbgn.pd.ProcessNodeVertical");
        ShapeObj *procNode = factory->createShape(type);
        m_shape = procNode;
        procNode->setCentrePos(m_basept);
    }

    // Layout all satellites.
    QList<DSBClone*> allSats = getAllSatellites();
    for (int i = 0; i < allSats.size(); i++)
    {
        DSBClone *cl = allSats.at(i);
        cl->layout();
    }
    // Relpts for satellites above
    int numAbove = m_inSatellites.size();
    for (int i = 0; i < numAbove; i++)
    {
        DSBClone *sat = m_inSatellites.at(i);
        QPointF p = satPos(i+1,numAbove,ABOVE);
        QString foo = m_id;
        sat->setRelPt(p);
    }
    // Relpts for satellites below
    int numBelow = m_outSatellites.size();
    for (int i = 0; i < numBelow; i++)
    {
        DSBClone *sat = m_outSatellites.at(i);
        QPointF p = satPos(i+1,numBelow,BELOW);
        sat->setRelPt(p);
    }

    // For now, process nodes have the fixed size of:
    //   16x16 box
    //   stems of length 17
    // Hence, 50x16 (horiz.), or 16x50 (vert.).
    qreal height = 50;
    if (numAbove > 0) { height += 90; }
    if (numAbove > 2) { height += 10; }
    if (numBelow > 0) { height += 90; }
    if (numBelow > 2) { height += 10; }
    m_size = QSizeF(216,height);
    return m_size;
}

QRectF DSBReaction::getBbox()
{
    // layout should have already been called, and relpt set
    QPointF centre = m_relpt;
    QSizeF minsize = QSizeF(16,50);
    QPointF ulc = QPointF( centre.x()-minsize.width()/2,
                           centre.y()-minsize.height()/2);
    QRectF rect = QRectF(ulc,minsize);
    QList<DSBClone*> allSats = getAllSatellites();
    for (int i = 0; i < allSats.size(); i++)
    {
        DSBClone *cl = allSats.at(i);
        QRectF clRect = cl->getBbox();
        rect = rect.united(clRect);
    }
    return rect;
}

QSizeF DSBReaction::getSize()
{
    return m_size;
}

QPointF DSBReaction::getBasePt()
{
    return m_basept;
}

void DSBReaction::setRelPt(QPointF p)
{
    m_relpt = p;
}

void DSBReaction::drawRelTo(QPointF q)
{
    QPointF r = m_relpt + q;
    drawAt(r);
}

void DSBReaction::redraw()
{
    drawAt(m_basept);
}

void DSBReaction::drawAt(QPointF r)
{
    m_basept = r;
    // Add shape to the canvas if not already there.
    if (!shapeOnCanvas)
    {
        QUndoCommand *cmd = new CmdCanvasSceneAddItem(m_canvas, m_shape);
        m_canvas->currentUndoMacro()->addCommand(cmd);
        shapeOnCanvas = true;
    }
    // Draw all satellites.
    QList<DSBClone*> allSats = getAllSatellites();
    for (int i = 0; i < allSats.size(); i++)
    {
        DSBClone *cl = allSats.at(i);
        cl->drawRelTo(r);
        // connector
        /*
        Connector *conn = new Connector();
        conn->initWithConnection(procNode,cl->getShape());
        QUndoCommand *cmd = new CmdCanvasSceneAddItem(m_canvas, conn);
        m_canvas->currentUndoMacro()->addCommand(cmd);
        */
        connectTo(cl);
    }
}

void DSBReaction::connectTo(DSBClone *cl)
{
    // Do we already have a connection to this clone?
    if (m_connectors.contains(cl)) { return; }
    // If not, then create one.
    Connector *conn = new Connector();
    conn->initWithConnection(m_shape,cl->getShape());
    if (cl==m_mainInput || m_inputBranchHeads.contains(cl)
            || m_inSatellites.contains(cl))
    {
        if (m_reversible) {
            conn->setDirected(true);
        }
    }
    else if (cl==m_mainOutput || m_outputBranchHeads.contains(cl)
             || m_outSatellites.contains(cl))
    {
        conn->setDirected(true);
    }
    m_connectors.insert(cl,conn);
    QUndoCommand *cmd = new CmdCanvasSceneAddItem(m_canvas, conn);
    m_canvas->currentUndoMacro()->addCommand(cmd);
}

void DSBReaction::clearConnectors()
{
    m_canvas->stop_graph_layout();
    QList<DSBClone*> clones = m_connectors.keys();
    foreach (DSBClone *cl, clones)
    {
        Connector *conn = m_connectors.value(cl);
        m_connectors.remove(cl);
        //m_canvas->deleteItem(conn);
        m_canvas->clearSelection();
        conn->setSelected(true);
        m_canvas->deleteSelection();
        //
        delete conn;
    }
    m_canvas->restart_graph_layout();
}

ShapeObj *DSBReaction::getShape()
{
    return m_shape;
}

void DSBReaction::moveShape(qreal dx, qreal dy)
{
    if (m_shape)
    {
        m_shape->moveBy(dx,dy);
    }
}

void DSBReaction::connectedComponent(QSet<DSBClone *> &ccClones, QSet<DSBReaction *> &ccReacs)
{
    ccReacs.insert(this); // add self to component
    // ignore clones already seen
    QSet<DSBClone*> cset = getAllClones().toSet().subtract(ccClones);
    // ask remaining ones to add to the connected component
    foreach (DSBClone *cl, cset)
    {
        cl->connectedComponent(ccClones, ccReacs);
    }
}

}

























