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

#ifndef SBGNGLYPH_H
#define SBGNGLYPH_H

#include <QList>
#include <QSet>
#include <QDomNode>
#include <QRectF>
#include <QPointF>

namespace dunnart {

class ShapeObj;

struct SBGNPort
{
    SBGNPort(QDomNode node);
    QString m_id;
    QPointF m_pos;
};

enum ProcOrient { VERT, HORIZ };

class SBGNGlyph
{
public:
    SBGNGlyph(QDomNode node);
    QString toString(void);
    ShapeObj *shape(void);
    QString id(void);
    qreal cx(void);
    qreal cy(void);
    void addNeighbour(SBGNGlyph* nbr);
    void getRestrConnComp(QSet<SBGNGlyph *> &R, QSet<SBGNGlyph *> &C);

private:
    void makeShape(void);
    void addPort(QDomNode port);

    QString m_id;
    QString m_class;
    bool m_isCloned;
    QString m_labelText;
    QRectF m_labelBbox;
    QRectF m_Bbox;
    ShapeObj *m_shape;
    ProcOrient m_orient;
    QList<SBGNPort*> m_ports;
    QList<SBGNGlyph*> m_nbrs;
};

}

#endif // SBGNGLYPH_H
