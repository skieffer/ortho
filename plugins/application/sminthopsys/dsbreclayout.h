/*
 * Sminthopsys - Dunnart Systems Biology plugin
 *
 * Copyright (C) 2011  Monash University
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
 * Author(s): Steven Kieffer <http://skieffer.info>
*/

#ifndef DSBRECLAYOUT_H
#define DSBRECLAYOUT_H

#include <QSizeF>
#include <QPointF>

namespace dunnart {

//! Classes that implement this interface conform to a recursive layout
//! paradigm for objects that contain other objects, as detailed in the
//! documentation for the methods below.
class DSBRecLayout
{
public:

    //! The first step of the layout method should be to recurse, i.e.
    //! to ask all contained objects to lay themselves out. Second,
    //! based on the sizes that those contained objects return, each of
    //! them should be assigned a location within the present object,
    //! using their setRelPt method. Finally, it should now be possible
    //! to return the size of the present object.
    virtual QSizeF layout() = 0;

    //! Each reclayout object must have an established base point,
    //! or origin for its local coordinate system. For a box, this might
    //! be the upper-left corner; for a disc it might be the center point.
    //! The "relpt" of the present object should be set by its parent
    //! container, to give its location relative to the origin of that
    //! parent container.
    virtual void setRelPt(QPointF p) = 0;

    //! After the layout phase, the drawing phase can begin. Typically,
    //! the top-level container will be asked to draw itself relative to
    //! the point (0,0), so that with its relpt set to p, it will draw
    //! itself with its base point at p. Thereafter, each
    //! parent container drawn rel-to q, and having its own relpt set to p,
    //! computes the absolute location of its own base point as q + p,
    //! draws itself with its own base point at q + p, and then
    //! and asks each of its children to draw themselves
    //! rel-to the point q + p.
    virtual void drawRelTo(QPointF q) = 0;

    //! Ignore own rel-pt, and simply draw with base point at r.
    virtual void drawAt(QPointF r) = 0;

    //! Draw again, with base point at same place as last time.
    virtual void redraw() = 0;

    //! Convenience method for drawing top-level container, relative to
    //! the origin of the overall coordinate system.
    void draw()
    {
        drawRelTo(QPointF(0,0));
    }

    //! Retrieve the size, after layout has been called.
    virtual QSizeF getSize() = 0;

    virtual void acceptCanvasBaseAndRelPts(QPointF parentBasePt) = 0;

};

}

#endif // DSBRECLAYOUT_H
