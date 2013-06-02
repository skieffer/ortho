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

#include "sbgnArc.h"

namespace dunnart {

SBGNArc::SBGNArc(QDomNode node)
{
    QDomNamedNodeMap attrs = node.attributes();
    m_srcId = attrs.namedItem("source").toAttr().value();
    m_tgtId = attrs.namedItem("target").toAttr().value();
    m_class = attrs.namedItem("class").toAttr().value();

    QDomNode start = node.namedItem("start");
    QDomNamedNodeMap startAttrs = start.attributes();
    float x = startAttrs.namedItem("x").toAttr().value().toFloat();
    float y = startAttrs.namedItem("y").toAttr().value().toFloat();
    m_startPt = QPointF(x,y);

    QDomNode end = node.namedItem("end");
    QDomNamedNodeMap endAttrs = end.attributes();
    x = endAttrs.namedItem("x").toAttr().value().toFloat();
    y = endAttrs.namedItem("y").toAttr().value().toFloat();
    m_endPt = QPointF(x,y);

}

QString SBGNArc::srcId()
{
    return m_srcId;
}

QString SBGNArc::tgtId()
{
    return m_tgtId;
}

QString SBGNArc::toString()
{
    QString s = "";
    s += "SBGN Arc:\n";
    s += QString("  src: %1, tgt: %2, class: %3\n").arg(m_srcId).arg(m_tgtId).arg(m_class);
    //
    return s;
}

}









