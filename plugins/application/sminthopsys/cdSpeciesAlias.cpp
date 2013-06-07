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

#include "cdSpeciesAlias.h"

#include <QtGui>
#include <QString>
#include <QDomNode>
#include <QDomAttr>
#include <QDomNamedNodeMap>

CDSpeciesAlias::CDSpeciesAlias(QDomNode node)
{
    QDomNamedNodeMap attrs = node.attributes();
    QDomNode idNode = attrs.namedItem("id");
    QDomAttr idAttr = idNode.toAttr();
    m_id = idAttr.value();

    m_speciesName = attrs.namedItem("species").toAttr().value();
    m_compartmentAlias = attrs.namedItem("compartmentAlias").toAttr().value();

    QDomNode bounds = node.namedItem("bounds");
    QDomNamedNodeMap battrs = bounds.attributes();
    m_x = battrs.namedItem("x").toAttr().value().toDouble();
    m_y = battrs.namedItem("y").toAttr().value().toDouble();
    m_width = battrs.namedItem("w").toAttr().value().toDouble();
    m_height = battrs.namedItem("h").toAttr().value().toDouble();

    //m_x = battrs.namedItem("x").toAttr().value();
    //m_y = battrs.namedItem("y").toAttr().value();
    //m_width = battrs.namedItem("w").toAttr().value();
    //m_height = battrs.namedItem("h").toAttr().value();

    qDebug() << toString();
    //qDebug() << "found species alias with id " << idAttr.value();
}

QString CDSpeciesAlias::toString()
{
    QString s = "";
    s += m_id + "," + m_speciesName + "," + m_compartmentAlias;
    s += "," + QString::number(m_x);
    s += "," + QString::number(m_y);
    s += "," + QString::number(m_width);
    s += "," + QString::number(m_height);

    //s += "," + m_x;
    //s += "," + m_y;
    //s += "," + m_width;
    //s += "," + m_height;

    return s;
}
