/*
 * Dunnart - Constraint-based Diagram Editor
 *
 * Copyright (C) 2003-2007  Michael Wybrow  <mjwybrow@users.sourceforge.net>
 * Copyright (C) 2006-2008  Monash University
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
*/

#include <QSignalMapper>

#include "layoutproperties.h"
#include "canvas.h"
#include "canvasview.h"


namespace dunnart {


LayoutPropertiesDialog::LayoutPropertiesDialog(Canvas *canvas, QWidget *parent)
    : QDockWidget(parent),
      m_canvas(NULL),
      m_mode_signal_mapper(NULL)
{
    setupUi(this);

    m_mode_signal_mapper = new QSignalMapper(this);
    m_mode_signal_mapper->setMapping(
            organicStructureButton, Canvas::OrganicLayout);
    m_mode_signal_mapper->setMapping(
            flowStructureButton, Canvas::FlowLayout);
    m_mode_signal_mapper->setMapping(
            layeredStuctureButton, Canvas::LayeredLayout);

    connect(m_mode_signal_mapper, SIGNAL(mapped(int)),
            this, SIGNAL(setOptStructuralLayoutMode(int)));

    connect(organicStructureButton, SIGNAL(clicked()),
            m_mode_signal_mapper, SLOT(map()));
    connect(flowStructureButton, SIGNAL(clicked()),
            m_mode_signal_mapper, SLOT(map()));
    connect(layeredStuctureButton, SIGNAL(clicked()),
            m_mode_signal_mapper, SLOT(map()));

    connect(layoutButton, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptAutomaticGraphLayout(bool)));

    connect(preventOverlapsCheckBox, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptPreventOverlaps(bool)));
    connect(this, SIGNAL(optChangedPreventOverlaps(bool) ),
            preventOverlapsCheckBox, SLOT(setChecked(bool)));

    connect(preventOverlapsCheckBox2, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptPreventOverlaps(bool)));
    connect(this, SIGNAL(optChangedPreventOverlaps(bool) ),
            preventOverlapsCheckBox2, SLOT(setChecked(bool)));

    connect(snapToCheckBox, SIGNAL(clicked(bool)),
            this, SIGNAL(setOptSnapTo(bool)));
    connect(this, SIGNAL(optChangedSnapTo(bool)),
            snapToCheckBox, SLOT(setChecked(bool)));

    connect(gridSnapCheckBox, SIGNAL(clicked(bool)),
            this, SIGNAL(setOptGridSnap(bool)));
    connect(this, SIGNAL(optChangedGridSnap(bool)),
            gridSnapCheckBox, SLOT(setChecked(bool)));

    connect(relaxCheckBox, SIGNAL(clicked(bool)),
            this, SIGNAL(setOptRelax(bool)));
    connect(this, SIGNAL(optChangedRelax(bool)),
            relaxCheckBox, SLOT(setChecked(bool)));

    connect(preserveTopologyCheckBox, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptPreserveTopology(bool)));
    connect(this, SIGNAL(optChangedPreserveTopology(bool) ),
            preserveTopologyCheckBox, SLOT(setChecked(bool)));

    connect(rubberBandRoutingCheckBox, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptRubberBandRouting(bool)));
    connect(this, SIGNAL(optChangedRubberBandRouting(bool) ),
            rubberBandRoutingCheckBox, SLOT(setChecked(bool)));

    connect(rubberBandRoutingCheckBox2, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptRubberBandRouting(bool)));
    connect(this, SIGNAL(optChangedRubberBandRouting(bool) ),
            rubberBandRoutingCheckBox2, SLOT(setChecked(bool)));

    connect(pageBoundaryCheckBox, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptFitWithinPage(bool)));
    connect(this, SIGNAL(optChangedFitWithinPage(bool) ),
            pageBoundaryCheckBox, SLOT(setChecked(bool)));

    connect(pageBoundaryCheckBox2, SIGNAL(clicked(bool) ),
            this, SIGNAL(setOptFitWithinPage(bool)));
    connect(this, SIGNAL(optChangedFitWithinPage(bool) ),
            pageBoundaryCheckBox2, SLOT(setChecked(bool)));

    connect(spacingSlider, SIGNAL(sliderMoved(int)),
            this, SIGNAL(setOptShapeNonoverlapPadding(int)));
    connect(spacingSlider2, SIGNAL(sliderMoved(int)),
            this, SIGNAL(setOptShapeNonoverlapPadding(int)));
    connect(this, SIGNAL(optChangedShapeNonoverlapPadding(int)),
            spacingSlider, SLOT(setValue(int)));
    connect(this, SIGNAL(optChangedShapeNonoverlapPadding(int)),
            spacingSlider2, SLOT(setValue(int)));

    // Preserve topology control should only be enabled when overlaps are.
    connect(this, SIGNAL(optChangedPreventOverlaps(bool) ),
            preserveTopologyCheckBox, SLOT(setEnabled(bool)));

    changeCanvas(canvas);
}


void LayoutPropertiesDialog::changeCanvas(Canvas *canvas)
{
    if (m_canvas)
    {
        disconnect(m_canvas, 0, this, 0);
        disconnect(this, 0, m_canvas, 0);
        disconnect(idealLengthSlider, 0, m_canvas, 0);
        disconnect(snapDistanceSlider, 0, m_canvas, 0);
        disconnect(snapStrengthSlider, 0, m_canvas, 0);
        disconnect(gridWidthSlider, 0, m_canvas, 0);
        disconnect(gridHeightSlider, 0, m_canvas, 0);
        disconnect(relaxThresholdSlider, 0, m_canvas, 0);
        disconnect(flowSeparationSlider, 0, m_canvas, 0);
        disconnect(flowDirectionDial, 0, m_canvas, 0);
    }
    m_canvas = canvas;

    connect(this, SIGNAL(setOptAutomaticGraphLayout(bool)),
            m_canvas, SLOT(setOptAutomaticGraphLayout(bool)));
    connect(m_canvas, SIGNAL(optChangedAutomaticLayout(bool) ),
            this, SLOT(changeAutomaticLayoutMode(bool)));

    connect(this, SIGNAL(setOptStructuralLayoutMode(int)),
            m_canvas, SLOT(setOptLayoutModeFromInt(int)));
    connect(m_canvas, SIGNAL(optChangedLayoutMode(int)),
            this, SLOT(changeStructuralLayoutMode(int)));

    connect(this, SIGNAL(setOptPreventOverlaps(bool)),
            m_canvas, SLOT(setOptPreventOverlaps(bool)));
    connect(m_canvas, SIGNAL(optChangedPreventOverlaps(bool) ),
            this, SIGNAL(optChangedPreventOverlaps(bool)));

    connect(this, SIGNAL(setOptSnapTo(bool)),
            m_canvas, SLOT(setOptSnapTo(bool)));
    connect(m_canvas, SIGNAL(optChangedSnapTo(bool)),
            this, SIGNAL(optChangedSnapTo(bool)));

    connect(this, SIGNAL(setOptGridSnap(bool)),
            m_canvas, SLOT(setOptGridSnap(bool)));
    connect(m_canvas, SIGNAL(optChangedGridSnap(bool)),
            this, SIGNAL(optChangedGridSnap(bool)));

    connect(this, SIGNAL(setOptRelax(bool)),
            m_canvas, SLOT(setOptRelax(bool)));
    connect(m_canvas, SIGNAL(optChangedRelax(bool)),
            this, SIGNAL(optChangedRelax(bool)));

    connect(this, SIGNAL(setOptPreserveTopology(bool) ),
            m_canvas, SLOT(setOptPreserveTopology(bool) ));
    connect(m_canvas, SIGNAL(optChangedPreserveTopology(bool) ),
            this, SIGNAL(optChangedPreserveTopology(bool) ));

    connect(this, SIGNAL(setOptRubberBandRouting(bool)),
            m_canvas, SLOT(setOptRubberBandRouting(bool)));
    connect(m_canvas, SIGNAL(optChangedRubberBandRouting(bool) ),
            this, SIGNAL(optChangedRubberBandRouting(bool) ));

    connect(this, SIGNAL(setOptFitWithinPage(bool)),
            m_canvas, SLOT(setOptFitWithinPage(bool)));
    connect(m_canvas, SIGNAL(optChangedFitWithinPage(bool) ),
            this, SIGNAL(optChangedFitWithinPage(bool) ));

    connect(idealLengthSlider, SIGNAL(sliderMoved(int)),
            m_canvas, SLOT(setOptIdealEdgeLengthModifierFromSlider(int)));
    connect(m_canvas, SIGNAL(optChangedIdealEdgeLengthModifier(double)),
            this, SLOT(changeIdealEdgeLength(double)));

    connect(snapDistanceSlider, SIGNAL(sliderMoved(int)),
            m_canvas, SLOT(setOptSnapDistanceModifierFromSlider(int)));
    connect(m_canvas, SIGNAL(optChangedSnapDistanceModifier(double)),
            this, SLOT(changeSnapDistance(double)));

    connect(snapStrengthSlider, SIGNAL(sliderMoved(int)),
            m_canvas, SLOT(setOptSnapStrengthModifierFromSlider(int)));
    connect(m_canvas, SIGNAL(optChangedSnapStrengthModifier(double)),
            this, SLOT(changeSnapStrength(double)));

    connect(gridWidthSlider, SIGNAL(sliderMoved(int)),
            m_canvas, SLOT(setOptGridWidthModifierFromSlider(int)));
    connect(m_canvas, SIGNAL(optChangedGridWidthModifier(double)),
            this, SLOT(changeGridWidth(double)));

    connect(gridHeightSlider, SIGNAL(sliderMoved(int)),
            m_canvas, SLOT(setOptGridHeightModifierFromSlider(int)));
    connect(m_canvas, SIGNAL(optChangedGridHeightModifier(double)),
            this, SLOT(changeGridHeight(double)));

    connect(relaxThresholdSlider, SIGNAL(sliderMoved(int)),
            m_canvas, SLOT(setOptRelaxThresholdModifierFromSlider(int)));
    connect(m_canvas, SIGNAL(optChangedRelaxThresholdModifier(double)),
            this, SLOT(changeRelaxThreshold(double)));

    connect(flowSeparationSlider, SIGNAL(sliderMoved(int)),
            m_canvas, SLOT(setOptFlowSeparationModifierFromSlider(int)));
    connect(m_canvas, SIGNAL(optChangedDirectedEdgeSeparationModifier(double)),
            this, SLOT(changeDirectedEdgeSeparationModifier(double)));

    connect(flowDirectionDial, SIGNAL(valueChanged(int)),
            m_canvas, SLOT(setOptFlowDirectionFromDial(int)));
    connect(m_canvas, SIGNAL(optChangedFlowDirection(int)),
            flowDirectionDial, SLOT(setValue(int)));

    connect(this, SIGNAL(setOptShapeNonoverlapPadding(int)),
            m_canvas, SLOT(setOptShapeNonoverlapPadding(int)));
    connect(m_canvas, SIGNAL(optChangedShapeNonoverlapPadding(int)),
            this, SIGNAL(optChangedShapeNonoverlapPadding(int)));

    // Set initial control values.
    idealLengthSlider->setSliderPosition(
            m_canvas->optIdealEdgeLengthModifier() * 100);

    double snap_distance = m_canvas->optSnapDistanceModifier();
    snapDistanceSlider->setSliderPosition(snap_distance);
    snapDistanceLabel->setText(QString("Distance: %1").arg(snap_distance));

    double snap_strength = m_canvas->optSnapStrengthModifier();
    snapStrengthSlider->setSliderPosition(snap_strength);
    snapStrengthLabel->setText(QString("Strength: %1").arg(snap_strength));

    double grid_width = m_canvas->optGridWidthModifier();
    gridWidthSlider->setSliderPosition(grid_width);
    gridWidthSlider->setText(QString("Width: %1").arg(grid_width));

    double grid_height = m_canvas->optGridHeightModifier();
    gridHeightSlider->setSliderPosition(grid_height);
    gridHeightSlider->setText(QString("Height: %1").arg(grid_height));

    double relax_threshold = m_canvas->optRelaxThresholdModifier();
    relaxThresholdSlider->setSliderPosition(relax_threshold * 100);
    relaxStressLabel->setText(QString("Stress Threshold: %1").arg(relax_threshold));

    flowSeparationSlider->setSliderPosition(
            m_canvas->optFlowSeparationModifier() * 100);

    bool value = m_canvas->optPreventOverlaps();
    preventOverlapsCheckBox->setChecked(value);
    preventOverlapsCheckBox2->setChecked(value);
    preserveTopologyCheckBox->setEnabled(value);

    snapToCheckBox->setChecked(m_canvas->optSnapTo());
    gridSnapCheckBox->setChecked(m_canvas->optGridSnap());
    relaxCheckBox->setChecked(m_canvas->optRelax());

    flowDirectionDial->setSliderPosition(m_canvas->optFlowDirection());

    spacingSlider->setSliderPosition(m_canvas->optShapeNonoverlapPadding());
    spacingSlider2->setSliderPosition(m_canvas->optShapeNonoverlapPadding());

    value = m_canvas->optPreserveTopology();
    preserveTopologyCheckBox->setChecked(value);

    value = m_canvas->optRubberBandRouting();
    rubberBandRoutingCheckBox->setChecked(value);
    rubberBandRoutingCheckBox2->setChecked(value);

    value = m_canvas->optFitWithinPage();
    pageBoundaryCheckBox->setChecked(value);
    pageBoundaryCheckBox2->setChecked(value);

    changeAutomaticLayoutMode(m_canvas->optAutomaticGraphLayout());
    changeStructuralLayoutMode(m_canvas->optLayoutMode());
}

void LayoutPropertiesDialog::changeStructuralLayoutMode(int mode)
{
    switch(mode)
    {
        case Canvas::OrganicLayout:
            organicStructureButton->setChecked(true);
            flowStructureButton->setChecked(false);
            layeredStuctureButton->setChecked(false);
            flowOptionsLabel->setEnabled(false);
            flowSeparationSlider->setEnabled(false);
            flowDirectionDial->setEnabled(false);
            break;
        case Canvas::FlowLayout:
            organicStructureButton->setChecked(false);
            flowStructureButton->setChecked(true);
            layeredStuctureButton->setChecked(false);
            flowOptionsLabel->setEnabled(true);
            flowSeparationSlider->setEnabled(true);
            flowDirectionDial->setEnabled(true);
            break;
        case Canvas::LayeredLayout:
            organicStructureButton->setChecked(false);
            flowStructureButton->setChecked(false);
            layeredStuctureButton->setChecked(true);
            flowOptionsLabel->setEnabled(true);
            flowSeparationSlider->setEnabled(true);
            flowDirectionDial->setEnabled(true);
            break;
        default:
            break;
    }
}

void LayoutPropertiesDialog::changeAutomaticLayoutMode(bool auto_layout)
{
    layoutButton->setChecked(auto_layout);
    if (auto_layout)
    {
        layoutModeLabel->setText(tr("Graph layout mode:\nAutomatic layout"));
        layoutOptionPages->setCurrentIndex(1);
    }
    else
    {
        layoutModeLabel->setText(tr("Graph layout mode:\nManual layout"));
        layoutOptionPages->setCurrentIndex(0);
    }
}

void LayoutPropertiesDialog::changeDirectedEdgeSeparationModifier(double value)
{
    flowSeparationSlider->setSliderPosition(value * 100);
}

void LayoutPropertiesDialog::changeIdealEdgeLength(double value)
{
    idealLengthSlider->setSliderPosition(value * 100);
}

void LayoutPropertiesDialog::changeSnapDistance(double value)
{
    snapDistanceSlider->setSliderPosition(value);
    snapDistanceLabel->setText(QString("Distance: %1").arg(value));
}

void LayoutPropertiesDialog::changeSnapStrength(double value)
{
    snapStrengthSlider->setSliderPosition(value);
    snapStrengthLabel->setText(QString("Strength: %1").arg(value));
}

void LayoutPropertiesDialog::changeGridWidth(double value)
{
    gridWidthSlider->setSliderPosition(value);
    gridWidthSlider->setText(QString("Width: %1").arg(value));
}

void LayoutPropertiesDialog::changeGridHeight(double value)
{
    gridHeightSlider->setSliderPosition(value);
    gridHeightSlider->setText(QString("Height: %1").arg(value));
}

void LayoutPropertiesDialog::changeRelaxThreshold(double value)
{
    relaxThresholdSlider->setSliderPosition(value * 100);
    relaxStressLabel->setText(QString("Stress Threshold: %1").arg(value));
}

}
// vim: filetype=cpp ts=4 sw=4 et tw=0 wm=0 cindent
