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

//diag:
#include <iostream>
#include "libdunnartcanvas/templates.h"
#include "libdunnartcanvas/graphlayout.h"
#include "libdunnartcanvas/pluginshapefactory.h"
//

#include "findBranchesDialog.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/canvasitem.h"
#include "libdunnartcanvas/shape.h"

#include "pdepn.h"
#include "dsbclone.h"
#include "dsbspecies.h"
#include "dsbcompartment.h"
#include "dsbbranch.h"
#include "dsbpathway.h"

namespace dunnart {

FindBranchesDialog::FindBranchesDialog(Canvas *canvas, QWidget *parent)
    : QDialog(parent)
{
    m_canvas = canvas;

    // Build the dialog.
    QLabel *endpointLabel = new QLabel(this);
    endpointLabel->setText(tr("Endpoint species:"));

    m_endpointEdit = new QLineEdit(this);

    m_endpointCBox = new QComboBox(this);
    m_endpointCBox->addItem(tr("first"));
    m_endpointCBox->addItem(tr("last"));

    QLabel *layoutLabel = new QLabel(this);
    layoutLabel->setText(tr("Layout:"));

    m_layoutCBox = new QComboBox(this);
    m_layoutCBox->addItem(tr("Vertical"));
    m_layoutCBox->addItem(tr("Horizontal"));
    m_layoutCBox->addItem(tr("Clockwise"));
    m_layoutCBox->addItem(tr("Counterclockwise"));

    QPushButton *findButton = new QPushButton(tr("Find"),this);
    findButton->setDefault(true);
    QPushButton *cancelButton = new QPushButton(tr("Cancel"),this);
    QPushButton *testButton = new QPushButton(tr("Test"),this);

    connect(findButton, SIGNAL(clicked()), this, SLOT(findBranches()));
    connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
    connect(testButton, SIGNAL(clicked()), this, SLOT(test()));

    QGridLayout *layout = new QGridLayout;
    layout->addWidget(endpointLabel,  0, 0);
    layout->addWidget(m_endpointEdit, 0, 1);
    layout->addWidget(m_endpointCBox,   0, 2);
    layout->addWidget(layoutLabel,    1, 0);
    layout->addWidget(m_layoutCBox,     1, 1);
    layout->addWidget(cancelButton,   2, 0);
    layout->addWidget(findButton,     2, 1);
    layout->addWidget(testButton,     2, 2);
    setLayout(layout);

    setWindowTitle(tr("Find Branches"));

    // Populate endpoint box.
    getSelectedClone();

    // Watch change of selection.
    connect(m_canvas, SIGNAL(selectionChanged()),
            this, SLOT(canvasSelectionChanged()));
}

void FindBranchesDialog::test(){
    PluginShapeFactory *factory = sharedPluginShapeFactory();
    ShapeObj *r1 = factory->createShape("org.dunnart.shapes.rect");
    ShapeObj *r2 = factory->createShape("org.dunnart.shapes.rect");
    r1->setSize(QSizeF(100,200));
    r2->setSize(QSizeF(25,25));
    //r1->addContainedShape(r2);
    m_canvas->stop_graph_layout();
    m_canvas->addItem(r1);
    m_canvas->addItem(r2);
    r1->addContainedShape(r2);
    m_canvas->interrupt_graph_layout();
    /*
    if (m_endpointClone)
    {
        //DSBCompartment *comp = m_endpointClone->getSpecies()->getCompartment();
        //comp->dumpPathwayNodePositions();
        //comp->jogPathways();
        //LinearTemplate *lintemp = new LinearTemplate(0,0,TEMPLATE_LINEAR_VERT,m_canvas);
        //m_canvas->m_graphlayout->initThread();
        //m_canvas->m_graphlayout = new GraphLayout(m_canvas);
        //m_canvas->m_graphlayout->freeShiftFromDunnart = false;
    }
    */
}

/* Respond to a change in the canvas selection.
 */
void FindBranchesDialog::canvasSelectionChanged()
{
    getSelectedClone();
}

/* Consult the canvas selection, and check whether exactly one
   shape node is selected. If so, populate edit box with its name.
   Otherwise, put a message requesting selection of just one shape.
 */
void FindBranchesDialog::getSelectedClone()
{
    CanvasItemList selection = m_canvas->selectedItems();
    int n = selection.size();
    if (n != 1)
    {
        m_endpointEdit->setText(tr("not exactly one node selected"));
    }
    else // Exactly one object was selected.
    {
        CanvasItem *item = selection.first();
        //ShapeObj *shape = qobject_cast<ShapeObj *> (item);
        PDEPN *epn = dynamic_cast<PDEPN *>(item);
        if (epn)
        {
            m_endpointEdit->setText(epn->getLabel());
            m_endpointIDString = epn->idString();
            m_endpointClone = epn->getClone();
        }
        else
        {
            m_endpointEdit->setText(tr("selection is not an EPN"));
        }

    }

}

/* Carry out the find-branches action, as specified in the dialog box.
  */
void FindBranchesDialog::findBranches()
{
    // If an endpoint has been selected, then go ahead.
    if (m_endpointClone)
    {
        // Get the compartment in which the endpoint clone lives.
        DSBCompartment *comp = m_endpointClone->getSpecies()->getCompartment();

        // Determine whether we're searching forward, or backward.
        int whichEnd = m_endpointCBox->currentIndex();
        bool forward = (whichEnd == 0);

        // Ask the compartment to clone all blacklisted species except for
        // that of the selected clone.
        QSet<QString> names = comp->m_default_blacklist.toSet();
        QString name = m_endpointClone->getSpecies()->getName();
        names.remove(name);
        comp->setDiscreteCloningsByName(names.toList());

        // Find branches.
        QList<DSBBranch*> branches = comp->findBranches(m_endpointClone, forward);

        // Build pathway.
        if (branches.size()>0)
        {
            DSBPathway *pathway = new DSBPathway(m_endpointClone, branches);
            pathway->setCanvas(m_canvas);
            comp->addPathway(pathway);
        }
        // Redisplay compartment.
        comp->redisplay();
    }

    // Finally, "accept" the click on the dialog's OK button.
    accept();
}

}
