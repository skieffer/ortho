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

#ifndef SBGNARC_H
#define SBGNARC_H

#include <QDomNode>
#include <QPointF>

namespace dunnart {

class SBGNArc
{
public:
    SBGNArc(QDomNode node);
    QString toString(void);
    QString srcId(void);
    QString tgtId(void);

private:
    QString m_srcId;
    QString m_tgtId;
    QString m_class;
    QPointF m_startPt;
    QPointF m_endPt;
};

}

#endif // SBGNARC_H
