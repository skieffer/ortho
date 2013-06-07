/*
 * Dunnart - Constraint-based Diagram Editor
 *
 * Copyright (C) 2010-2011  Monash University
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
 * Author(s): Sarah Boyd  <Sarah.Boyd@monash.edu>
*/


#include "pdsourcesink.h"

using namespace dunnart;

QPainterPath SourceOrSink::buildPainterPath(void)
{
    QPainterPath p;

    // make sure this glyph stays cirular
    qreal draw_width = width();
    if ( height() < width()) draw_width = height();

    QRectF r(-draw_width/2, -draw_width/2, draw_width, draw_width);

    p.addEllipse(r);
    p.moveTo(r.bottomLeft());
    p.lineTo(r.topRight());

    return p;
}
// vim: filetype=cpp ts=4 sw=4 et tw=0 wm=0 cindent
