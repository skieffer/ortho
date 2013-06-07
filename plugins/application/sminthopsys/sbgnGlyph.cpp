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

#include <QMap>

#include "sbgnGlyph.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/pluginshapefactory.h"

namespace dunnart {

SBGNPort::SBGNPort(QDomNode node)
{
    QDomNamedNodeMap attrs = node.attributes();

    float x = attrs.namedItem("x").toAttr().value().toFloat();
    float y = attrs.namedItem("y").toAttr().value().toFloat();
    m_pos = QPointF(x,y);

    m_id = attrs.namedItem("id").toAttr().value();
}


SBGNGlyph::SBGNGlyph(QDomNode node) :
    m_shape(NULL), m_orient(HORIZ)
{
    QDomNamedNodeMap attrs = node.attributes();
    m_id = attrs.namedItem("id").toAttr().value();
    m_class = attrs.namedItem("class").toAttr().value();

    // Get orientation, if defined.
    QDomNode orient = attrs.namedItem("orientation");
    if (!orient.isNull())
    {
        QString ostr = orient.toAttr().value();
        if (ostr == "vertical")
        {
            m_orient = VERT;
        }
        else if (ostr == "horizontal")
        {
            m_orient = HORIZ;
        }
    }

    // Get any ports that may be defined.
    QDomNode port = node.namedItem("port");
    while (!port.isNull()) {
        addPort(port);
        port = port.nextSiblingElement("port");
    }

    // Label
    QDomNode label = node.namedItem("label");
    if (!label.isNull()) {
        QDomNamedNodeMap labelAttrs = label.attributes();
        m_labelText = labelAttrs.namedItem("text").toAttr().value();
        // Get bbox of label, if it is specified.
        if (label.hasChildNodes())
        {
            QDomNode labelBboxNode = label.namedItem("bbox");
            if (!labelBboxNode.isNull())
            {
                QDomNamedNodeMap lbboxAttrs = labelBboxNode.attributes();
                float x = lbboxAttrs.namedItem("x").toAttr().value().toFloat();
                float y = lbboxAttrs.namedItem("y").toAttr().value().toFloat();
                float w = lbboxAttrs.namedItem("w").toAttr().value().toFloat();
                float h = lbboxAttrs.namedItem("h").toAttr().value().toFloat();
                m_labelBbox = QRectF(x,y,w,h);
            }
        }
    }

    // Is it marked as a clone?
    QDomNode cloneMarker = node.namedItem("clone");
    m_isCloned = !cloneMarker.isNull();

    // Bbox of node
    QDomNode bbox = node.namedItem("bbox");
    QDomNamedNodeMap bboxAttrs = bbox.attributes();
    float x = bboxAttrs.namedItem("x").toAttr().value().toFloat();
    float y = bboxAttrs.namedItem("y").toAttr().value().toFloat();
    float w = bboxAttrs.namedItem("w").toAttr().value().toFloat();
    float h = bboxAttrs.namedItem("h").toAttr().value().toFloat();
    m_Bbox = QRectF(x,y,w,h);

    // Make the shape.
    makeShape();
}

void SBGNGlyph::addNeighbour(SBGNGlyph *nbr)
{
    m_nbrs.append(nbr);
}

/* Build the set C of SBGNGlyphs which (a) are connected to this one, and
   (b) are in the passed set R.
  */
void SBGNGlyph::getRestrConnComp(QSet<SBGNGlyph*> &R, QSet<SBGNGlyph*> &C)
{
    // Add self to component.
    C.insert(this);
    // Compute set of eligible neighbours.
    QSet<SBGNGlyph*> N = m_nbrs.toSet().intersect(R);
    // Recurse on those neighbours which are not already in the set C.
    foreach (SBGNGlyph *n, N) {
        if (!C.contains(n)) {
            n->getRestrConnComp(R,C);
        }
    }
}

QString SBGNGlyph::id()
{
    return m_id;
}

void SBGNGlyph::makeShape()
{
    if (m_class == "compartment")
    {
        // TODO
        return;
    }

    QString type = "org.sbgn.pd.00UnspecifiedEPN";
    if (m_class == "process")
    {
        switch (m_orient)
        {
        case VERT:
            type = "org.sbgn.pd.ProcessNodeVertical";
            break;
        case HORIZ:
            type = "org.sbgn.pd.ProcessNodeHorizontal";
            break;
        }
    }
    else
    {
        QMap<QString,QString> map;
        map.insert("simple chemical","org.sbgn.pd.01SimpleChemEPN");
        map.insert("macromolecule","org.sbgn.pd.02MacromolEPN");
        //map.insert("","");
        type = map.value(m_class,"org.sbgn.pd.00UnspecifiedEPN");
    }

    PluginShapeFactory *factory = sharedPluginShapeFactory();
    m_shape = factory->createShape(type);

    //qDebug() << "centre: " << m_Bbox.center().x() << "," << m_Bbox.center().y();
    m_shape->setCentrePos(m_Bbox.center());
    m_shape->setSize(m_Bbox.size());

    if (m_labelText.length() > 0) {
        m_shape->setLabel(m_labelText);
    }

}

ShapeObj *SBGNGlyph::shape()
{
    return m_shape;
}

void SBGNGlyph::addPort(QDomNode port)
{
    SBGNPort *p = new SBGNPort(port);
    m_ports.append(p);
}

qreal SBGNGlyph::cx()
{
    QPointF c = m_Bbox.center();
    return c.x();
}

qreal SBGNGlyph::cy()
{
    QPointF c = m_Bbox.center();
    return c.y();
}

QString SBGNGlyph::toString()
{
    QString s = "";
    s += "SBGN Glyph:\n";
    s += QString("  id: %1, class: %2\n").arg(m_id).arg(m_class);
    s += (m_isCloned ? "  clone\n" : "");
    s += QString("  label: %1\n").arg(m_labelText);
    //
    return s;
}

}









