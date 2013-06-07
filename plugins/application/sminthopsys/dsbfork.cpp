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

#include "dsbfork.h"
#include "dsbclone.h"
#include "dsbpathway.h"
#include "dsbbranch.h"
#include "dsbreaction.h"

namespace dunnart {

DSBFork::DSBFork(DSBClone *c) :
    m_centre(c),
    m_pathway(NULL),
    m_mainInput(NULL),
    m_mainOutput(NULL)
{}

void DSBFork::addUpstream(DSBReaction *reac)
{
    m_upstreamReacs.append(reac);
}

void DSBFork::addDownstream(DSBReaction *reac)
{
    m_downstreamReacs.append(reac);
}

void DSBFork::setMainInput(DSBReaction *mi)
{
    m_mainInput = mi;
}

void DSBFork::setMainOutput(DSBReaction *mo)
{
    m_mainOutput = mo;
}

void DSBFork::setPathway(DSBPathway *pw)
{
    m_pathway = pw;
}

DSBPathway *DSBFork::getPathway()
{
    return m_pathway;
}

QSizeF DSBFork::layout()
{
    // TODO: Write proper method.
    // For now, we just set relpt of branches headed by off-main downstreamReacs.
    // We assume our centre clone has already been given a relpt within its branch,
    // and that its branch has as well, relative to its pathway.
    QPointF crel = m_centre->m_relpt;
    DSBBranch *b = m_centre->getBranch();
    QPointF brel = b->m_relpt;
    QPointF frel = brel + crel; // relpt of this fork within its pathway,
                                // where basept of fork is centre of centre clone
    qreal x = frel.x(); qreal y = frel.y();
    qreal jump = 250;
    qreal disp = -jump;
    x += disp; y += jump/2;
    for (int i = 0; i < m_downstreamReacs.size(); i++)
    {
        DSBReaction *reac = m_downstreamReacs.at(i);
        if (reac == m_mainOutput) { continue; }
        DSBBranch *b = reac->getBranch();
        b->setRelPt(QPointF(x,y));
        // Prepare next x value.
        if (disp<0) { disp *= -1; }
        else { disp = -disp-jump; }
    }
    m_size = QSizeF(0,0);
    return m_size;
}

void DSBFork::acceptCanvasBaseAndRelPts(QPointF parentBasePt)
{
    // TODO
}

QSizeF DSBFork::getSize()
{
    return m_size;
}

void DSBFork::setRelPt(QPointF p)
{
    m_relpt = p;
}

void DSBFork::drawRelTo(QPointF q)
{
    QPointF r = m_relpt + q;
    drawAt(r);
}

void DSBFork::redraw()
{
    drawAt(m_basept);
}

void DSBFork::drawAt(QPointF r)
{
    m_basept = r;
}


}
