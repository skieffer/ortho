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

#ifndef SMINTHOPSYSFREEPATHWAY_H
#define SMINTHOPSYSFREEPATHWAY_H

#include <QList>
#include <QMap>
#include <QRectF>

#include "dsbpathway.h"
#include "libogdf/ogdf/basic/Graph_d.h"

namespace ogdf {
class GraphAttributes;
}
namespace dunnart {

class DSBNode;
class DSBClone;
class DSBReaction;
class CanvasItem;

typedef QMap<DSBNode*,ogdf::node> nodemap;

class FreePathway : public DSBPathway
{
public:
    // Constructor
    FreePathway(QList<DSBClone*> clones, QList<DSBReaction*> reacs);
    // RecLayout methods
    QSizeF layout();
    void setRelPt(QPointF p);
    void drawRelTo(QPointF q);
    void drawAt(QPointF r);
    void redraw();
    QSizeF getSize();
    void acceptCanvasBaseAndRelPts(QPointF parentBasePt);
    QList<CanvasItem*> getAllShapes();
    ogdf::Graph *getOGDFGraph(nodemap& nodeMap);
    void extractPosAndSize(nodemap& nodeMap, ogdf::GraphAttributes& GA);
    void injectPositions(nodemap& nodeMap, ogdf::GraphAttributes& GA);
    static void jog(double scale, nodemap& nodeMap);
    static QRectF getBbox(nodemap& nodeMap);

private:
    QList<DSBClone*> m_clones;
    QList<DSBReaction*> m_reactions;

    QPointF m_relpt;
    QPointF m_basept;
    QSizeF m_size;
};

}

#endif // SMINTHOPSYSFREEPATHWAY_H
