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

#ifndef DSBCOMPARTMENT_H
#define DSBCOMPARTMENT_H

#include <QString>
#include <QList>
#include <QMap>
#include <QObject>
#include <QGraphicsItem>

#include "dsbreclayout.h"
#include "libdunnartcanvas/shape.h"

#define CONTAINEDSHAPES

class QRectF;
class QPainter;
class QPainterPath;
class QStyleOptionGraphicsItem;
class QWidget;
class QGraphicsSceneMouseEvent;
class QAction;
class QMenu;

namespace dunnart {

class DSBSpecies;
class DSBReaction;
class DSBClone;
class DSBNode;
class DSBBranch;
class DSBPathway;
class Canvas;
class Cluster;
class CanvasItem;

class DSBCompartment;

class CompShape2 : public ShapeObj
{
public:
    CompShape2(DSBCompartment *comp, qreal x, qreal y, qreal w, qreal h);
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
protected:
    QAction *buildAndExecContextMenu(QGraphicsSceneMouseEvent *event, QMenu& menu);
private:
    DSBCompartment *m_compartment;
    qreal m_penWidth;
    qreal m_cornerRadius;
};

#ifdef CONTAINEDSHAPES
class CompartmentShape : public ShapeObj
{
public:
    CompartmentShape(DSBCompartment *comp, qreal x, qreal y, qreal w, qreal h);
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    void resize(qreal w, qreal h);
protected:
    //void mousePressEvent(QGraphicsSceneMouseEvent *event);
    QAction *buildAndExecContextMenu(QGraphicsSceneMouseEvent *event, QMenu& menu);
private:
    DSBCompartment *m_compartment;
    qreal m_width;
    qreal m_height;
    qreal m_penWidth;
    qreal m_cornerRadius;
};
#else
class CompartmentShape : public QGraphicsItem
{
public:
    CompartmentShape(DSBCompartment *comp, qreal x, qreal y, qreal w, qreal h);
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    void resize(qreal w, qreal h);
protected:
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    QAction *buildAndExecContextMenu(QGraphicsSceneMouseEvent *event, QMenu& menu);
private:
    DSBCompartment *m_compartment;
    qreal m_width;
    qreal m_height;
    qreal m_penWidth;
    qreal m_cornerRadius;
};
#endif

class DSBCompartment : public DSBRecLayout, public QObject
{
public:
    // Constructors
    DSBCompartment(QString compartmentName);
    // Building methods
    void addSpecies(DSBSpecies *spec);
    void addReaction(DSBReaction *reac);
    void addCompartment(DSBCompartment *comp);
    void addCompartments(QList<DSBCompartment*> comps);
    void addPathway(DSBPathway *pw);
    QList<DSBBranch*> findBranches(DSBClone *endpt, bool forward);
    QList<DSBBranch*> findBranches(DSBClone *endpt, bool forward, QList<QString> blacklist);
    void buildConnectedPathways(void);
    void setTrivialCloning(void);
    void setDiscreteCloningsByName(QList<QString> names);
    void cloneCurrencyMolecules(void);
    // Various layout methods
    QSizeF rowLayout(void);
    QSizeF layoutSquareCloneArray(QList<DSBClone*> clones, int ulx, int uly);
    // RecLayout methods
    QSizeF layout();
    void setRelPt(QPointF p);
    void drawRelTo(QPointF q);
    void drawAt(QPointF r);
    void redraw();
    QSizeF getSize();
    void acceptCanvasBaseAndRelPts(QPointF parentBasePt);
    void redisplay(bool reLayout = true);
    // Misc get and set
    QString getName(void);
    void setParent(DSBCompartment *comp);
    void setCanvas(Canvas *canvas);
    Canvas *getCanvas(void);
    void setBoundaryVisible(bool b);
    void adjustSize(void);

    void dumpPathwayNodePositions(void);
    void dumpAllClonePositions(void);

    QList<QString> m_default_blacklist;

public slots:
    void jogPathways(void);

private:
    QString m_compartmentName;
    QPointF m_relpt;
    QPointF m_basept;
    QSizeF m_size;
    Canvas *m_canvas;
    DSBCompartment *m_parentCompartment;
    QList<DSBSpecies *> m_species;
    QList<DSBReaction *> m_reactions;
    QList<DSBCompartment*> m_compartments;
    QList<DSBPathway*> m_pathways;
    bool m_show_reactions;
    bool m_boundaryVisible;
    CompShape2 *m_boundaryShape;
    Cluster *m_cluster;

    QList<DSBClone*> getAllClones(void);
    QList<DSBClone*> getLooseClones(void);
    QList<CanvasItem*> getAllShapes(void);


};

}

#endif // DSBCOMPARTMENT_H
