/*
 * Dunnart - Constraint-based Diagram Editor
 *
 * Copyright (C) 2003-2007  Michael Wybrow
 * Copyright (C) 2006-2008  Monash University
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
 * Author(s): Michael Wybrow  <http://michael.wybrow.info/>
*/
#include <QtGui>
#include <QtSvg>
#include <QParallelAnimationGroup>

#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/shape.h"
#include "libdunnartcanvas/canvasitem.h"
#include "libdunnartcanvas/visibility.h"
#include "libdunnartcanvas/canvas.h"
#include "libdunnartcanvas/separation.h"
#include "libdunnartcanvas/distribution.h"
#include "libdunnartcanvas/guideline.h"
#include "libdunnartcanvas/connector.h"
#include "libdunnartcanvas/canvasitem.h"
#include "libdunnartcanvas/gmlgraph.h"
#include "libdunnartcanvas/connectionpininfo.h"
#include "libdunnartcanvas/pluginfileiofactory.h"

#include "libdunnartcanvas/graphlayout.h"
#include "libdunnartcanvas/oldcanvas.h"
#include "libdunnartcanvas/undo.h"
#include "libdunnartcanvas/templates.h"
#include "libdunnartcanvas/ui/createseparation.h"
#include "libdunnartcanvas/ui/createdistribution.h"
#include "libdunnartcanvas/ui/createtemplate.h"

#include "libdunnartcanvas/BCLayout.h"

#include "libavoid/libavoid.h"
#include "libtopology/orthogonal_topology.h"
#include "libdunnartcanvas/graphdata.h"
//#include "libcola/compound_constraints.h"
#include "libcola/shortest_paths.h"


namespace dunnart {

enum {
    PASTE_UPDATEOBJIDS = 1,
    PASTE_UPDATEIDPROPS,
    PASTE_FINDBADDISTROS,
    PASTE_REMOVEBADDISTROS,
    PASTE_SELECTSHAPES
};

Actions::Actions()
{
    clear();
}


void Actions::clear(void)
{
    flags = ACTION_NONE;
    moveList.clear();
    resizeList.clear();
}

bool Actions::empty(void) const
{
    return moveList.empty() && resizeList.empty() &&
           !(flags & (ACTION_ADDITIONS | ACTION_MODIFICATIONS | ACTION_DELETIONS));
}

// Resize Handle Types
enum {
    HAND_TOP_LEFT = 0,
    HAND_TOP_CENTRE,
    HAND_TOP_RIGHT,
    HAND_RIGHT_CENTRE,
    HAND_BOTTOM_RIGHT,
    HAND_BOTTOM_CENTRE,
    HAND_BOTTOM_LEFT,
    HAND_LEFT_CENTRE
};

class SelectionResizeHandle : public Handle {
    public:
        SelectionResizeHandle(int index, double xpos,
                double ypos) :
            Handle(NULL, index, 0),
            m_conn(NULL),
            m_pos(xpos, ypos)
        {
            // Position it in front of other objects.
            this->setZValue(1000000);

            switch (index)
            {
            case HAND_TOP_LEFT:
            case HAND_BOTTOM_RIGHT:
                setCursor(Qt::SizeFDiagCursor);
                break;
            case HAND_TOP_CENTRE:
            case HAND_BOTTOM_CENTRE:
                setCursor(Qt::SizeVerCursor);
                break;
            case HAND_TOP_RIGHT:
            case HAND_BOTTOM_LEFT:
                setCursor(Qt::SizeBDiagCursor);
                break;
            case HAND_RIGHT_CENTRE:
            case HAND_LEFT_CENTRE:
                setCursor(Qt::SizeHorCursor);
                break;
            default:
                break;
            }

            setHoverMessage("Selection Resize Handle - Click "
                    "and drag to resize the selected shapes. .");
        }
    protected:
        Connector *m_conn;
        QPointF m_pos;
        virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event)
        {
            int index = this->handleFlags();
            QPointF scenePosition = event->scenePos();

            Canvas *canvas = dynamic_cast<Canvas *> (scene());

            canvas->moveSelectionResizeHandle(index, scenePosition);
        }
        void mousePressEvent(QGraphicsSceneMouseEvent *event)
        {
            Canvas *canvas = dynamic_cast<Canvas *> (scene());
            canvas->beginUndoMacro(tr("Resize"));
            canvas->storeSelectionResizeInfo();
            event->accept();
            Handle::mousePressEvent(event);
        }
};

Canvas::Canvas()
    : QGraphicsScene(),
      gd2013_batch_mode(false),
      gd2013_step(0),
      gd2013_type(0),
      m_action_finished_timer(NULL),
      m_visual_page_buffer(3.0),
      m_layout_update_timer(NULL),
      m_layout_finish_timer(NULL),
      m_processing_layout_updates(false),
      m_graphlayout(NULL),
      m_router(NULL),
      m_svg_renderer(NULL),
      m_status_bar(NULL),
      m_max_string_id(0),
      m_max_internal_id(0),
      m_gml_graph(NULL),
      m_use_gml_clusters(true),
      m_most_recent_stress(0),
      m_most_recent_ortho_obj_func(0),
      m_most_recent_crossing_count(0),
      m_most_recent_coincidence_count(0),
      m_stress_bar_maximum(5000),
      m_obliquity_bar_maximum(5000),
      m_connector_nudge_distance(0),
      m_why_is_it_triggered_twice(true),
      m_trying_alignments(false),
      m_max_align_tries(180),
      m_num_align_tries(0),
      m_align_pairs_tried(NULL),
      m_apply_alignments_epsilon(-DBL_MAX),
      m_opt_draw_separation_indicators(true),
      m_opt_edge_node_repulsion(false),
      m_opt_ideal_edge_length_modifier(1.0),
      m_opt_snap_distance_modifier(50),
      m_opt_snap_strength_modifier(20),
      m_opt_snap_grid_width(100.0),
      m_opt_snap_grid_height(100.0),
      m_opt_relax_threshold_modifier(0.1),
      m_opt_shape_nonoverlap_padding(0.0),
      m_dragged_item(NULL),
      m_dragged_with_force(false),
      m_lone_selected_item(NULL),
      m_undo_stack(NULL),
      m_current_undo_macro(NULL),
      m_hide_selection_handles(false),
      m_overlay_router_raw_routes(false),
      m_overlay_router_display_routes(false),
      m_overlay_router_obstacles(false),
      m_overlay_router_visgraph(false),
      m_overlay_router_orthogonal_visgraph(false),
      m_rendering_for_printing(false),
      m_edit_mode(ModeSelection),
      m_routing_event_posted(false),
      m_canvas_font(NULL),
      m_canvas_font_size(DEFAULT_CANVAS_FONT_SIZE),
      m_animation_group(NULL),
      m_bclayout(NULL),
      //m_why_is_it_triggered_twice(true),
      //m_trying_alignments(false),
      //m_max_align_tries(180),
      //m_num_align_tries(0),
      //m_align_pairs_tried(NULL),
      /*
      m_most_recent_stress(0),
      m_most_recent_ortho_obj_func(0),
      m_most_recent_crossing_count(0),
      m_most_recent_coincidence_count(0),
      m_stress_bar_maximum(5000),
      m_obliquity_bar_maximum(5000),
      m_apply_alignments_epsilon(-DBL_MAX),
      m_opt_draw_separation_indicators(true),
      m_opt_edge_node_repulsion(false),
            */
      m_tentative_guideline_timestamp(0)
{
    m_ideal_connector_length = 100;
    m_flow_separation_modifier = 0.5;
    m_sticky_nodes = false;
    m_downward_edges = false;
    m_avoid_connector_crossings = false;
    m_nudge_orthogonal_routes = true;
    m_avoid_cluster_crossings = false;
    m_rectangle_constraint_test = false;
    m_batch_diagram_layout = false;
    m_simple_paths_during_layout = true;
    m_force_orthogonal_connectors = false;
    m_infer_tentative_alignments = true;

    m_undo_stack = new QUndoStack(this);

    // This is faster for dynamic scenes.
    setItemIndexMethod(QGraphicsScene::NoIndex);

    // Options for controlling behaviour of Constraint-based layout:
    m_opt_automatic_graph_layout = false;
    m_opt_prevent_overlaps       = false;
    m_opt_snap_to                = false;
    m_opt_grid_snap              = false;
    m_opt_relax                  = false;
    m_opt_preserve_topology      = false;
    m_opt_rubber_band_routing    = false;
    m_opt_fit_within_page        = false;
    m_opt_colour_interfering_connectors = false;
    m_opt_connector_rounding_distance = 5;
    m_opt_stuctural_editing_disabled = false;
    m_opt_flow_direction = FlowDown;
    m_opt_layered_alignment_position = ShapeMiddle;

    // Set default canvas background colour.
    m_opt_canvas_background_colour = QColor(189, 189, 223);

    // Default list of connector colors. Use only dark colors (and full
    // opacity) since connectors are drawn as thin lines so light
    // colors will be hard to see.
    m_default_connector_colours <<
            QColor(0, 0, 0) <<             // default (black)
            QColor(0, 0, 255) <<           // blue
            QColor(139, 0, 0) <<           // darkred
            QColor(143, 188, 143) <<       // dark sea green
            QColor(85, 26, 139) <<         // purple4
            QColor(139, 0, 139) <<         // magenta4
            QColor(139, 35, 35) <<         // brown4
            QColor(25, 25, 112) <<         // midnight blue
            QColor(148, 0, 211);           // dark violet

    m_graphlayout = new GraphLayout(this);

    // Avoid::PolyLineRouting
    m_router = new Avoid::Router(Avoid::OrthogonalRouting |
            Avoid::PolyLineRouting);
    m_router->SimpleRouting = true;
    m_router->setRoutingParameter(Avoid::shapeBufferDistance, 4.0);
    m_router->setRoutingOption(
            Avoid::nudgeOrthogonalSegmentsConnectedToShapes, true);

    m_router->setRoutingParameter(Avoid::segmentPenalty, 50);
    m_router->setRoutingParameter(Avoid::clusterCrossingPenalty, 0);
    //m_router->setRoutingParameter(Avoid::fixedSharedPathPenalty);

    m_animation_group = new QParallelAnimationGroup();

    m_selection_resize_handles = QVector<SelectionResizeHandle *>(8);
    for (int i = 0; i < 8; ++i)
    {
        m_selection_resize_handles[i] = new SelectionResizeHandle(i, 0, 0);
        addItem(m_selection_resize_handles[i]);
        m_selection_resize_handles[i]->setVisible(false);
    }

    connect(this, SIGNAL(selectionChanged()), this,
            SLOT(selectionChangeTriggers()));

    m_bclayout = new BCLayout(this);

    computeNeighbourhoods();

#ifdef FPSTIMER
    m_convergence_update_count = 0;
    m_convergence_timer_running = false;
#endif
}


Canvas::~Canvas()
{
    delete m_graphlayout;
    delete m_router;
    delete m_animation_group;

    if (m_svg_renderer)
    {
        delete m_svg_renderer;
    }

    // Free selection resize handles if they are not currently displayed
    // and thus owned by the scene.
    for (int i = 0; i < m_selection_resize_handles.size(); ++i)
    {
        if (m_selection_resize_handles[i]->scene() == NULL)
        {
            // Not part of the scene, so we need to free.
            delete m_selection_resize_handles[i];
            m_selection_resize_handles[i] = NULL;
        }
    }
}


bool Canvas::loadGmlDiagram(const QFileInfo& fileInfo)
{
    setOptFitWithinPage(true);
    setOptAutomaticGraphLayout(true);
    setOptShapeNonoverlapPadding(10);
    int cxoff, cyoff;
    m_gml_graph = new gml::Graph(this, fileInfo.absolutePath().toStdString(),
            gml::Page(this), gml::COff(cxoff, cyoff));
    return true;
}


bool Canvas::loadDiagram(const QString& filename)
{
    if (filename.isEmpty())
    {
        return false;
    }

    QString errorMessage;
    QFileInfo fileInfo(filename);
    PluginFileIOFactory *fileIOFactory = sharedPluginFileIOFactory();
    bool successful = fileIOFactory->loadDiagramFromFile(this, fileInfo,
            errorMessage);

    if (successful)
    {
        this->setFilename(filename);
    }
    else
    {
        // We weren't successful loading, so show an error message.
        QString warning = QString(
                QObject::tr("<p><b>The document \"%1\" could not be loaded.</b></p>"
                "<p>%2</p>")).arg(fileInfo.fileName()).arg(errorMessage);

        QWidget *window = views().first()->window();
        QMessageBox message(QMessageBox::Warning, "Error Loading File",
                            warning, QMessageBox::Ok, window);
        message.setWindowModality(Qt::WindowModal);
        message.exec();
    }
    return successful;
}

void Canvas::setSvgRendererForFile(const QString& filename)
{
    m_svg_renderer = new QSvgRenderer(filename);
}

void Canvas::postDiagramLoad(void)
{
    setOptAutomaticGraphLayout(m_opt_automatic_graph_layout);

    // Update on-screen representation of dependant indicators.
    foreach (CanvasItem *item, items())
    {
        Guideline *gobj =    dynamic_cast<Guideline *> (item);
        Distribution *dobj = dynamic_cast<Distribution *> (item);
        Separation *sobj =   dynamic_cast<Separation *> (item);

        if (dobj)
        {
            dobj->updateFromLayout(dobj->getSeparation());
        }
        else if (sobj)
        {
            sobj->updateFromLayout(sobj->getSeparation());
        }
        else  if (gobj)
        {
            gobj->updateFromLayout(gobj->position(), true);
        }
    }
    if (!m_opt_preserve_topology)
    {
        bool lastSimpleRouting = m_router->SimpleRouting;
        m_router->SimpleRouting = false;
        if (!m_batch_diagram_layout)
        {
            //reroute_connectors(this);
        }
        m_router->SimpleRouting = lastSimpleRouting;
    }

    // QT clear_undo_stack();
}


CanvasItem *Canvas::getItemByID(QString ID) const
{
    assert(this != NULL);
    foreach (CanvasItem *item, items())
    {
        if (ID == item->idString())
        {
            return item;
        }
    }
    return NULL;
}

CanvasItem *Canvas::getItemByInternalId(uint internalId) const
{
    assert(this != NULL);
    foreach (CanvasItem *item, items())
    {
        if (internalId == item->internalId())
        {
            return item;
        }
    }
    return NULL;
}

QFont &Canvas::canvasFont(void)
{
    if (m_canvas_font == NULL)
    {
        QFontDatabase database;
        database.addApplicationFont(":/resources/DejaVuSans.ttf");

        m_canvas_font = new QFont("DejaVu Sans", m_canvas_font_size);
    }
    return *m_canvas_font;
}

QUndoStack *Canvas::undoStack(void) const
{
    return m_undo_stack;
}

UndoMacro *Canvas::currentUndoMacro(void)
{
    if (m_current_undo_macro == NULL)
    {
        m_current_undo_macro = new UndoMacro(this);
    }
    return m_current_undo_macro;
}

UndoMacro *Canvas::beginUndoMacro(const QString& text)
{
    m_current_undo_macro = new UndoMacro(this);
    m_current_undo_macro->setText(text);
    m_undo_stack->push(m_current_undo_macro);

    return m_current_undo_macro;
}

void Canvas::endUndoMacro(void)
{
    m_current_undo_macro = NULL;
}


Actions& Canvas::getActions(void)
{
    return m_actions;
}

uint Canvas::assignInternalId(void)
{
    return ++m_max_internal_id;
}


QString Canvas::assignStringId(QString id)
{
    if (id.isEmpty())
    {
        id.setNum(++m_max_string_id);
    }
    else
    {
        if (idIsUnique(id))
        {
            // If id is a uint then count it in m_max_id, since we use integer
            // ids when we get clashes and reassign ids.
            bool isUInt = false;
            uint idUInt = id.toUInt (&isUInt);
            if (isUInt)
            {
                m_max_string_id = qMax(idUInt, m_max_string_id);
            }
        }
        else
        {
            ++m_max_string_id;
            qWarning("Clashing id \"%s\", reassigned \"%d\"", qPrintable(id),
                     m_max_string_id);
            id.setNum(m_max_string_id);
        }

    }
    return id;
}

bool Canvas::idIsUnique(QString id) const
{
    unsigned int count = 0;
    CanvasItem *item;
    foreach (item, items())
    {
        if (item->idString() == id)
        {
            ++count;
        }
    }
    return (count == 1);
}


void Canvas::setExpandedPage(const QRectF newExpandedPage)
{
    if (newExpandedPage != this->m_expanded_page)
    {
        QRectF updateArea = this->m_expanded_page;
        this->m_expanded_page = newExpandedPage;
        updateArea |= newExpandedPage;

        // Update each of the four page-ballon rectangles, triggering redraw.
        // Empty rectanges, if the drawing is contained on that side, do nothing.
        QRectF updateRect;

        updateRect = updateArea;
        updateRect.setTop(this->m_page.bottom());
        update(updateRect);

        updateRect = updateArea;
        updateRect.setBottom(this->m_page.top());
        update(updateRect);

        updateRect = updateArea;
        updateRect.setLeft(this->m_page.right());
        update(updateRect);

        updateRect = updateArea;
        updateRect.setRight(this->m_page.left());
        update(updateRect);
    }
}

void Canvas::drawBackground(QPainter *painter, const QRectF& rect)
{
    if ( m_rendering_for_printing )
    {
        if (m_opt_grid_snap) drawGridLines(painter, rect);
        return;
    }

    if ( m_expanded_page.isNull() )
    {
        // No expanded page: effectively show just the normal page.
        m_expanded_page = m_page;
    }

    // Draws purple background and the white page (if it is set).
    painter->fillRect(rect, m_opt_canvas_background_colour);
    painter->fillRect(m_expanded_page, QColor(200, 200, 200));
    painter->fillRect(m_page, QColor(255, 255, 255));
    painter->setPen(QColor(110, 110, 110));
    painter->drawRect(m_expanded_page);

    // Draw snap grid lines
    if (m_opt_grid_snap) drawGridLines(painter, rect);
}

void Canvas::drawGridLines(QPainter *painter, const QRectF &rect)
{
    QPen pen(QColor(128,128,128,128));
    pen.setWidth(1);
    pen.setCosmetic(true);
    painter->setPen(pen);

    double W = m_opt_snap_grid_width;
    double H = m_opt_snap_grid_height;

    int n0 = ceil(rect.left()/W);
    int n1 = floor(rect.right()/W);
    int n = n1 - n0 + 1;
    double x = n0*W;

    for (int j = 0; j < n; j++) {
        painter->drawLine(x,rect.top(),x,rect.bottom());
        x += W;
    }

    int m0 = ceil(rect.top()/H);
    int m1 = floor(rect.bottom()/H);
    int m = m1 - m0 + 1;
    double y = m0*H;

    for (int i = 0; i < m; i++) {
        painter->drawLine(rect.left(),y,rect.right(),y);
        y += H;
    }
}


void Canvas::drawForeground(QPainter *painter, const QRectF& rect)
{
    Q_UNUSED(rect);

    if ( m_rendering_for_printing )
    {
        // Don't draw any foreground at all.
        return;
    }

    if (m_overlay_router_obstacles)
    {
        QPen pen(Qt::red);
        pen.setCosmetic(true);
        painter->setPen(pen);
        foreach (CanvasItem *item, items())
        {
            ShapeObj *shape = dynamic_cast<ShapeObj *> (item);

            if (shape && shape->avoidRef)
            {
                // Draw the rectangular box used for orthogonal routing.
                Avoid::Box bBox = shape->avoidRef->routingBox();
                Avoid::Point topLeft = bBox.min;
                Avoid::Point bottomRight = bBox.max;

                QRectF rect(QPointF(topLeft.x, topLeft.y),
                        QPointF(bottomRight.x, bottomRight.y));
                painter->drawRect(rect);

                // Draw the polygon used for polyline routing.
                Avoid::Polygon poly = shape->avoidRef->routingPolygon();
                QPolygonF polygon;

                for (size_t i = 0; i < poly.size(); ++i)
                {
                    const Avoid::Point& point = poly.at(i);
                    polygon << QPointF(point.x, point.y);
                }
                painter->drawPolygon(polygon);
            }
        }
    }

    if (m_overlay_router_visgraph)
    {
        Avoid::EdgeList& visList = router()->visGraph;

        Avoid::EdgeInf *finish = visList.end();
        for (Avoid::EdgeInf *t = visList.begin(); t != finish; t = t->lstNext)
        {
            std::pair<Avoid::Point, Avoid::Point> ptpair = t->points();

            QPointF pt1(ptpair.first.x, ptpair.first.y);
            QPointF pt2(ptpair.second.x, ptpair.second.y);

            std::pair<Avoid::VertID, Avoid::VertID> ids = t->ids();

            if (ids.first.isConnPt() || ids.second.isConnPt())
            {
                // Endpt
                QColor colour(Qt::blue);
                colour.setAlpha(50);
                QPen pen(colour);
                pen.setCosmetic(true);
                painter->setPen(pen);
                painter->drawLine(pt1, pt2);
            }
            else
            {
                // Shape
                QColor colour(Qt::red);
                colour.setAlpha(50);
                QPen pen(colour);
                pen.setCosmetic(true);
                painter->setPen(pen);
                painter->drawLine(pt1, pt2);
            }
            QPen pen(Qt::black);
            pen.setCosmetic(true);
            painter->setPen(pen);
            painter->drawPoint(pt1);
            painter->drawPoint(pt2);
        }
    }

    if (m_overlay_router_orthogonal_visgraph)
    {
        Avoid::EdgeList& visList = router()->visOrthogGraph;

        Avoid::EdgeInf *finish = visList.end();
        for (Avoid::EdgeInf *t = visList.begin(); t != finish; t = t->lstNext)
        {
            std::pair<Avoid::Point, Avoid::Point> ptpair = t->points();

            QPointF pt1(ptpair.first.x, ptpair.first.y);
            QPointF pt2(ptpair.second.x, ptpair.second.y);

            std::pair<Avoid::VertID, Avoid::VertID> ids = t->ids();

            if (ids.first.isConnPt() || ids.second.isConnPt())
            {
                // Endpt
                QColor colour(Qt::blue);
                colour.setAlpha(50);
                QPen pen(colour);
                pen.setCosmetic(true);
                painter->setPen(pen);
                painter->drawLine(pt1, pt2);
            }
            else
            {
                // Shape
                QColor colour(Qt::red);
                colour.setAlpha(50);
                QPen pen(colour);
                pen.setCosmetic(true);
                painter->setPen(pen);
                painter->drawLine(pt1, pt2);
            }
            QPen pen(Qt::black);
            pen.setCosmetic(true);
            painter->setPen(pen);
            painter->drawPoint(pt1);
            painter->drawPoint(pt2);
        }
    }

    const int routeWidth = 3;
    if (m_overlay_router_raw_routes)
    {
        QColor colour(Qt::blue);
        colour.setAlpha(50);
        QPen pen(colour);
        pen.setWidth(routeWidth);
        pen.setCosmetic(true);
        painter->setPen(pen);
        for (Avoid::ConnRefList::const_iterator ci = router()->connRefs.begin();
             ci != router()->connRefs.end(); ++ci)
        {
            Avoid::Polygon path = (*ci)->route();
            QPolygonF points;
            for (unsigned int i = 0; i < path.size(); ++i)
            {
                points << QPointF(path.ps[i].x, path.ps[i].y);
            }
            painter->drawPolyline(points);
        }
    }

    if (m_overlay_router_display_routes)
    {
        QColor colour(Qt::blue);
        colour.setAlpha(50);
        QPen pen(colour);
        pen.setWidth(routeWidth);
        pen.setCosmetic(true);
        painter->setPen(pen);
        for (Avoid::ConnRefList::const_iterator ci = router()->connRefs.begin();
                ci != router()->connRefs.end(); ++ci)
        {
            Avoid::Polygon path = (*ci)->displayRoute();
            QPolygonF points;
            for (unsigned int i = 0; i < path.size(); ++i)
            {
                points << QPointF(path.ps[i].x, path.ps[i].y);
            }
            painter->drawPolyline(points);
        }
    }

    /*
        if (! router->clusterRefs.empty() )
        {
            // There are clusters so do cluster routing.
            for (Avoid::ClusterRefList::const_iterator cl =
                    router->clusterRefs.begin();
                    cl != router->clusterRefs.end(); ++cl)
            {
                Avoid::ReferencingPolygon& cBoundary = (*cl)->polygon();
                printf("PRINT: Points: %lu\n",  cBoundary.size());
                for (size_t j = 0; j < cBoundary.size(); ++j)
                {
                    Point p1 = cBoundary.at(j);
                    Point p2 = cBoundary.at((j + 1) % cBoundary.size());

                    unoffsetPoint(p1);
                    unoffsetPoint(p2);
                    aalineRGBA(screen, (int) p1.x, (int) p1.y, (int) p2.x,
                            (int) p2.y, 0, 0, 255, 255);
                }
            }
        }
    */
}


void Canvas::processSelectionDropEvent(QGraphicsSceneMouseEvent *event)
{
    Q_UNUSED (event)

    if (!isLayoutSuspended())
    {
        // Don't do this processing if the layout is running.
        return;
    }

    foreach (CanvasItem *selectedItem, selectedItems())
    {
        ShapeObj *selectedShape =
                dynamic_cast<ShapeObj *> (selectedItem);

        if (!selectedShape)
        {
            continue;
        }

        QRectF selectedShapeRect = selectedShape->boundingRect().translated(
                selectedShape->scenePos());
        
        foreach (CanvasItem *item, items())
        {
            ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
            
            if (shape)
            {
                shape->removeContainedShape(selectedShape);
            }
        }

        foreach (CanvasItem *item, items())
        {
            ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
            if (!shape || (shape == selectedShape))
            {
                continue;
            }

            QRectF shapeRect = 
                    shape->boundingRect().translated(shape->scenePos());

            if (shape && shapeRect.contains(selectedShapeRect))
            {
                shape->addContainedShape(selectedShape);
            }
        }
    }
}

void Canvas::setDraggedItem(CanvasItem *item, bool withForce)
{
    if (item == NULL)
    {
        // Finished dragging.
        QApplication::restoreOverrideCursor();
        if (m_dragged_item)
        {
            glueObjectsToIndicators();
            clearIndicatorHighlights(true);
            m_hide_selection_handles = false;
            repositionAndShowSelectionResizeHandles(true);
        }        
        m_dragged_item = NULL;
        if (m_dragged_with_force == false)
        {
            // If the dragging wasn't forceful, then we call
            // processResponseTasks() so it notices there is no longer
            // a dragged item.  If the dragging was forceful, this action
            // is performed once the layout converges.
            processResponseTasks();
        }
    }

    if ((m_dragged_item == NULL) && item)
    {
        // Started dragging.
        m_dragged_item = item;
        clearIndicatorHighlights(true);
        createIndicatorHighlightCache();
        m_hide_selection_handles = true;
        hideSelectionResizeHandles();
        beginUndoMacro(tr("Move"));
    }
    else
    {
        clearIndicatorHighlights();
        assert(item == m_dragged_item);
    }

    m_dragged_with_force |= withForce;
}


bool Canvas::layoutRunningAndNotProcessingUpdates(void) const
{
    return m_graphlayout->isRunning() && !m_processing_layout_updates;
}

void Canvas::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    // For some reason if a shape is dragged off the canvas under other 
    // widget and back then it doesn't recieve a mouse release event, so 
    // effectively pass on the event here.
    if (event->button() == Qt::LeftButton)
    {
        if (m_dragged_item)
        {
            m_dragged_item->dragReleaseEvent(event);
        }
    }

    QGraphicsScene::mouseReleaseEvent(event);
}


void Canvas::customEvent(QEvent *event)
{
    if (dynamic_cast<LayoutUpdateEvent *> (event))
    {
        //qDebug() << "Update " << (long long) this;
        this->startLayoutUpdateTimer();
    }
    else if (dynamic_cast<LayoutFinishedEvent *> (event))
    {
        //qDebug() << "Finish " << (long long) this;
        //qDebug() << "Layout finished.";
        this->startLayoutFinishTimer();
    }
    else if (dynamic_cast<RoutingRequiredEvent *> (event))
    {
        // Call libavoid's processTransaction and reroute connectors.
        m_routing_event_posted = false;
        reroute_connectors(this);
    }
    else if (dynamic_cast<ConstraintRejectedEvent *> (event))
    {
        qDebug() << "Received ConstraintRejectedEvent.";
        ConstraintRejectedEvent *cre = dynamic_cast<ConstraintRejectedEvent *> (event);
        if (cre->m_unsat) {
            appliedAlignmentWasUnsat(cre);
        } else if (m_rejecting_alignments) {
            rejectAlignmentsCallback(cre);
        } else {
            Guideline *gdln = cre->m_guideline;
            if (gdln->canvas())
            {
                stop_graph_layout();
                UndoMacro *undoMacro = beginUndoMacro(tr("Delete"));
                CanvasItemSet dummySet;
                gdln->deactivateAll(dummySet);
                QUndoCommand *cmd = new CmdCanvasSceneRemoveItem(this, gdln);
                undoMacro->addCommand(cmd);
            }
        }
    }
    else if (dynamic_cast<TrialAlignmentEvent*>(event))
    {
        TrialAlignmentEvent *tae = dynamic_cast<TrialAlignmentEvent*>(event);
        int flag = tae->m_flag;
        switch(flag) {
        case 0:
            applyAlignments();
            break;
        case 1:
            initRejectAlignments();
            break;
        case 2:
            rejectAlignments();
            break;
        }
    }
    else
    {
        QGraphicsScene::customEvent(event);
    }
}


void Canvas::toggleSelectedShapePinning(void)
{
    foreach (CanvasItem *item, selectedItems())
    {
        if (ShapeObj *shape = isShapeForLayout(item))
        {
            // Toggle pinned position setting.
            if (shape->isPinned())
            {
                shape->setPinned(false);
            }
            else
            {
                shape->setPinned(true);
            }
        }
    }
}


int Canvas::editMode(void) const
{
    return m_edit_mode;
}

void Canvas::setEditMode(int mode)
{
    m_edit_mode = mode;
    emit editModeChanged(mode);

    // Update selection.
    selectionChangeTriggers();

    // Repaint selected objects (for changed selection cues).
    foreach (QGraphicsItem *item, QGraphicsScene::selectedItems())
    {
        item->update();
    }
}


void Canvas::setStatusBar(QStatusBar *statusBar)
{
    m_status_bar = statusBar;

    if (m_status_bar && !m_status_messages.empty())
    {
        m_status_bar->showMessage(m_status_messages.top());
    }
}


void Canvas::pushStatusMessage(const QString& message)
{
    m_status_messages.push(message);

    if (m_status_bar)
    {
        m_status_bar->showMessage(m_status_messages.top());
    }
}


void Canvas::popStatusMessage(void)
{
    m_status_messages.pop();

    if (m_status_bar)
    {
        if (m_status_messages.empty())
        {
            m_status_bar->clearMessage();
        }
        else
        {
            m_status_bar->showMessage(m_status_messages.top());
        }
    }
}

void Canvas::postRoutingRequiredEvent(void)
{
    if (!m_routing_event_posted)
    {
        QCoreApplication::postEvent(this, new RoutingRequiredEvent(),
                Qt::NormalEventPriority);
        m_routing_event_posted = true;
    }
}

void Canvas::selectAll(void)
{
    QPainterPath selectionArea;
    selectionArea.addPolygon(sceneRect());
    setSelectionArea(selectionArea);
}


void Canvas::templateFromSelection(int type)
{
    CanvasItemList selected_items = selectedItems();

    Indicator *indicator = NULL;
//    Guideline *guide = NULL;
    switch (type)
    {
        case TEMPLATE_LINEAR_HORI:
        {
            qDebug("creating horizontal linear template");

            double xpos = 150;
            double ypos = 150;
//            guide = new Guideline(GUIDE_TYPE_HORI, 150);
//            guide = createAlignment(ALIGN_CENTER, selected_items);

            indicator = new LinearTemplate(xpos, ypos, type, this);
//            Canvas::addItem(guide);
            Canvas::addItem(indicator);
            break;
        }
        case TEMPLATE_LINEAR_VERT:
        {
            double xpos = 50;
            double ypos = 50;
            indicator = new LinearTemplate(xpos, ypos, type, this);
            break;
        }
        case TEMPLATE_BRANCHED:
        {
            double xpos = 50;
            double ypos = 150;
            indicator = new BranchedTemplate(xpos, ypos);
            break;
        }
        default:
            break;
    }

    // Delselect shapes so they can be moved by layout solver.
    deselectAll();
    // Clear previous moves.
    Actions& actions = getActions();
    actions.clear();
//    add_undo_record(DELTA_ADD, indicator);
//add_undo_record(DELTA_MOVE, guide);
    // Relayout.
    interrupt_graph_layout();
}

void Canvas::alignSelection(int type)
{
    CanvasItemList selected_items = selectedItems();
    beginUndoMacro(tr("Create Alignment"));
    Guideline *guide = createAlignment((atypes) type, selected_items);

    // Delselect shapes so they can be moved by layout solver.
    deselectAll();
    // Clear previous moves.
    Actions& actions = getActions();
    actions.clear();
    // UNDO add_undo_record(DELTA_MOVE, guide);
    // Relayout.
    interrupt_graph_layout();
}


void Canvas::distributeSelection(int type)
{
    CreateDistributionDialog *sepDialog =
            dynamic_cast<CreateDistributionDialog *> (sender()->parent());
    CanvasItemList selected_items = selectedItems();
    QWidget *window = (sepDialog) ? sepDialog->window() : NULL;
    beginUndoMacro(tr("Create Distribution"));
    createDistribution(window, (dtype) type, selected_items);

    // Delselect shapes so they can be moved by layout solver.
    deselectAll();
    // Clear previous moves.
    Actions& actions = getActions();
    actions.clear();
    // Relayout.
    interrupt_graph_layout();
}

void Canvas::separateSelection(int type)
{
    CreateSeparationDialog *sepDialog =
            dynamic_cast<CreateSeparationDialog *> (sender()->parent());
    double minSeparationDist = (sepDialog) ? sepDialog->separationDistance() : 50.0;
    CanvasItemList selected_items = selectedItems();

    QWidget *window = (sepDialog) ? sepDialog->window() : NULL;
    beginUndoMacro(tr("Create Separation"));
    createSeparation(window, (dtype) type, selected_items, minSeparationDist);

    // Deselect shapes so they can be moved by layout solver.
    deselectAll();
    // Clear previous moves.
    Actions& actions = getActions();
    actions.clear();
    // Relayout.
    fully_restart_graph_layout();
}


QRectF Canvas::pageRect(void) const
{
    return m_page;
}

double Canvas::visualPageBuffer(void) const
{
    return m_visual_page_buffer;
}

GraphLayout *Canvas::layout(void) const
{
    return m_graphlayout;
}


Avoid::Router *Canvas::router(void) const
{
    return m_router;
}


void Canvas::setPageRect(const QRectF &rect)
{
    if (!rect.isEmpty())
    {
        QRectF updateArea = m_expanded_page;
        updateArea |= m_page;

        m_page = rect;
        m_expanded_page = QRectF();

        updateArea |= m_page;

        // Schedule repaint.
        update(updateArea);
    }
}

QString Canvas::saveConstraintInfoToString(void) const
{
    QDomDocument doc("svg");

    QDomElement svg = doc.createElement("svg");
    doc.appendChild(svg);

    QDomElement options = writeLayoutOptionsToDomElement(doc);
    doc.appendChild(options);

    // Put things into a multimap before outputting them,
    // so that they will be sorted in the correct Z-order.
    foreach (CanvasItem *item, items())
    {
        Indicator *indicator = dynamic_cast<Indicator *> (item);

        if (indicator)
        {
            QDomElement node = indicator->to_QDomElement(XMLSS_ALL, doc);
            QDomNode parent = node.parentNode();
            if (!parent.isNull())
            {
                // The constructed node will have a parent if it was placed
                // in a group.  In this case, add the group node to the tree.
                svg.appendChild(node.parentNode());
            }
            else
            {
                svg.appendChild(node);
            }
        }
    }

    qDebug() << doc.toString();
    return doc.toString();
}


void Canvas::loadConstraintInfoFromString(const QString& constraintInfo)
{
    QDomDocument doc("svg");

    bool parseNamespaces = true;
    doc.setContent(constraintInfo, parseNamespaces);

    QDomElement root = doc.documentElement();

    // Actually do the pasting, in correct order.
    for (int pass = 0; pass < PASS_LAST; ++pass)
    {
        this->recursiveReadSVG(root, x_dunnartNs, pass);
    }

    this->fully_restart_graph_layout();
}

void Canvas::setIdealConnectorLength(const double length)
{
    m_ideal_connector_length = length;
}

double Canvas::idealConnectorLength(void) const
{
    return m_ideal_connector_length;
}

bool Canvas::avoidConnectorCrossings(void) const
{
    return m_avoid_connector_crossings;
}

bool Canvas::avoidClusterCrossings(void) const
{
    return m_avoid_cluster_crossings;
}


bool Canvas::optAutomaticGraphLayout(void) const
{
    return m_opt_automatic_graph_layout;
}

Canvas::LayoutMode Canvas::optLayoutMode(void) const
{
    return (LayoutMode) m_graphlayout->mode;
}

Canvas::FlowDirection Canvas::optFlowDirection(void) const
{
    return (FlowDirection) m_opt_flow_direction;
}

bool Canvas::optPreventOverlaps(void) const
{
    return m_opt_prevent_overlaps;
}

bool Canvas::optSnapTo(void) const
{
    return m_opt_snap_to;
}

bool Canvas::optGridSnap(void) const
{
    return m_opt_grid_snap;
}

bool Canvas::optEdgeNodeRepulsion(void) const
{
    return m_opt_edge_node_repulsion;
}

bool Canvas::optDrawSeparationIndicators(void) const
{
    return m_opt_draw_separation_indicators;
}

bool Canvas::optRelax(void) const
{
    return m_opt_relax;
}

bool Canvas::optPreserveTopology(void) const
{
    return m_opt_preserve_topology;
}

bool Canvas::optRubberBandRouting(void) const
{
    return m_opt_rubber_band_routing;
}

bool Canvas::optFitWithinPage(void) const
{
    return m_opt_fit_within_page;
}

bool Canvas::optColourInterferingConnectors(void) const
{
    return m_opt_colour_interfering_connectors;
}

void Canvas::setDebugCOLAOutput(const bool value)
{
    m_graphlayout->setOutputDebugFiles(value);
}

QColor Canvas::optCanvasBackgroundColour(void) const
{
    return m_opt_canvas_background_colour;
}


int Canvas::optRoutingShapePadding(void) const
{
    return (int) m_router->routingParameter(Avoid::shapeBufferDistance);
}

int Canvas::optShapeNonoverlapPadding(void) const
{
    return (int) m_opt_shape_nonoverlap_padding;
}

void Canvas::setOptAutomaticGraphLayout(const bool value)
{
    // Remember previous value.
    bool had_auto_layout = m_opt_automatic_graph_layout;

    // Set new value.
    m_opt_automatic_graph_layout = value;
    
    if (m_opt_automatic_graph_layout)
    {
        m_router->SimpleRouting = m_simple_paths_during_layout;
    }
    else
    {
        m_router->SimpleRouting = true;
    }
    emit optChangedAutomaticLayout(m_opt_automatic_graph_layout);
    
    fully_restart_graph_layout();

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    if (!m_opt_automatic_graph_layout && had_auto_layout)
    {
        // If autolayout was set previously, then update all the positions
        // of obstacles with libavoid and reroute connectors.
        QList<CanvasItem *> canvas_items = items();
        for (int i = 0; i < canvas_items.size(); ++i)
        {
            ShapeObj *shape = dynamic_cast<ShapeObj *> (canvas_items.at(i));
            if (shape)
            {
                router()->moveShape(shape->avoidRef, 0, 0);
            }
        }
        reroute_connectors(this, true);
    }
    QApplication::restoreOverrideCursor();
}

void Canvas::setOptLayoutMode(const LayoutMode mode)
{
    m_graphlayout->setLayoutMode((GraphLayout::Mode) mode);
    emit optChangedLayoutMode(mode);
    fully_restart_graph_layout();
}

void Canvas::setOptLayoutModeFromInt(const int mode)
{
    setOptLayoutMode((LayoutMode) mode);
}

bool Canvas::optStructuralEditingDisabled(void) const
{
    return m_opt_stuctural_editing_disabled;
}


void Canvas::setOptStructuralEditingDisabled(const bool value)
{
    m_opt_stuctural_editing_disabled = value;

    clearSelection();

    emit optChangedStructuralEditingDisabled(m_opt_stuctural_editing_disabled);
}

void Canvas::setSelection(const QList<CanvasItem *>& newSelection)
{
    clearSelection();
    foreach (CanvasItem *item, newSelection)
    {
        item->setSelected(true);
    }
}

void Canvas::setOptPreserveTopology(const bool value)
{
    m_opt_preserve_topology = value;
    emit optChangedPreserveTopology(m_opt_preserve_topology);
    fully_restart_graph_layout();
}


void Canvas::setOptRubberBandRouting(const bool value)
{
    m_opt_rubber_band_routing = value;

    m_router->RubberBandRouting = m_opt_rubber_band_routing;
    emit optChangedRubberBandRouting(m_opt_rubber_band_routing);
    restart_graph_layout();
}


void Canvas::setOptFitWithinPage(const bool value)
{
    m_opt_fit_within_page = value;
    if (!m_opt_fit_within_page)
    {
        setExpandedPage(QRectF());
    }
    emit optChangedFitWithinPage(m_opt_fit_within_page);
    interrupt_graph_layout();
}


void Canvas::setOptPreventOverlaps(const bool value)
{
    m_opt_prevent_overlaps = value;
    emit optChangedPreventOverlaps(m_opt_prevent_overlaps);
    fully_restart_graph_layout();
}


void Canvas::setOptSnapTo(const bool value)
{
    m_opt_snap_to = value;
    emit optChangedSnapTo(m_opt_snap_to);
    fully_restart_graph_layout();
}


void Canvas::setOptGridSnap(const bool value)
{
    m_opt_grid_snap = value;
    setOptIdealEdgeLengthModifier(m_opt_snap_grid_width/100.0);
    emit optChangedGridSnap(m_opt_grid_snap);
    fully_restart_graph_layout();
    this->update(combinedViewsRect());
}

void Canvas::setOptEdgeNodeRepulsion(const bool value)
{
    m_opt_edge_node_repulsion = value;
    emit optChangedEdgeNodeRepulsion(m_opt_edge_node_repulsion);
    fully_restart_graph_layout();
}

void Canvas::setOptDrawSeparationIndicators(const bool value)
{
    m_opt_draw_separation_indicators = value;
    this->update(combinedViewsRect());
}

QRectF Canvas::combinedViewsRect(void) const
{
    QRectF viewsRect;

    // Determine bounds by taking the union of all view bounds.
    foreach (QGraphicsView *view, views())
    {
        viewsRect = viewsRect.united(QRectF(view->mapToScene(0,0),
                view->mapToScene(view->width(), view->height())));
    }
    return viewsRect;
}


void Canvas::setOptRelax(const bool value)
{
    m_opt_relax = value;
    emit optChangedRelax(m_opt_relax);
    fully_restart_graph_layout();
}


void Canvas::setOptRoutingPenaltySegment(const int value)
{
    m_router->setRoutingParameter(Avoid::segmentPenalty, (double) value);

    reroute_all_connectors(this);
}


void Canvas::setOptConnRoundingDist(const int value)
{
    m_opt_connector_rounding_distance = value;

    redraw_connectors(this);
}


void Canvas::setOptRoutingShapePadding(const int value)
{
    m_router->setRoutingParameter(Avoid::shapeBufferDistance, (double) value);
    emit optChangedRoutingShapePadding(value);

    m_router->markAllObstaclesAsMoved();
    reroute_all_connectors(this);
    update();
}

void Canvas::setOptShapeNonoverlapPadding(const int value)
{
    m_opt_shape_nonoverlap_padding = value;
    emit optChangedShapeNonoverlapPadding(value);
    fully_restart_graph_layout();
}

void Canvas::setOptIdealEdgeLengthModifierFromSlider(int int_modifier)
{
    double double_modifier = int_modifier / 100.0;
    setOptIdealEdgeLengthModifier(double_modifier);
}


void Canvas::setOptIdealEdgeLengthModifier(double modifier)
{
    m_opt_ideal_edge_length_modifier = modifier;
    emit optChangedIdealEdgeLengthModifier(modifier);
    interrupt_graph_layout();
}


void Canvas::setOptSnapDistanceModifierFromSlider(int int_modifier)
{
    //qDebug() << "Snap Distance int modifier: " << int_modifier;
    double double_modifier = int_modifier / 1.0;
    setOptSnapDistanceModifier(double_modifier);
}

void Canvas::setOptSnapDistanceModifier(double modifier)
{
    m_opt_snap_distance_modifier = modifier;
    emit optChangedSnapDistanceModifier(modifier);
    fully_restart_graph_layout();
}


void Canvas::setOptSnapStrengthModifierFromSlider(int int_modifier)
{
    double double_modifier = int_modifier / 1.0;
    setOptSnapStrengthModifier(double_modifier);
}

void Canvas::setOptSnapStrengthModifier(double modifier)
{
    m_opt_snap_strength_modifier = modifier;
    emit optChangedSnapStrengthModifier(modifier);
    fully_restart_graph_layout();
}

void Canvas::setOptGridSizeFromSlider(int intValue)
{
    double doubleValue = intValue / 1.0;

    m_opt_snap_grid_width = doubleValue;
    m_opt_snap_grid_height = doubleValue;
    emit optChangedGridWidth(doubleValue);
    emit optChangedGridHeight(doubleValue);
    emit optChangedGridSize(doubleValue);
    fully_restart_graph_layout();
    this->update(combinedViewsRect());
}

void Canvas::setOptGridWidth(double value)
{
    m_opt_snap_grid_width =value;
    emit optChangedGridWidth(value);
    fully_restart_graph_layout();
    this->update(combinedViewsRect());
}

void Canvas::setOptGridHeight(double value)
{
    m_opt_snap_grid_height = value;
    emit optChangedGridHeight(value);
    fully_restart_graph_layout();
    this->update(combinedViewsRect());
}

void Canvas::setOptRelaxThresholdModifierFromSlider(int int_modifier)
{
    double double_modifier = int_modifier / 100.0;
    setOptRelaxThresholdModifier(double_modifier);
}

void Canvas::setOptRelaxThresholdModifier(double modifier)
{
    m_opt_relax_threshold_modifier = modifier;
    emit optChangedRelaxThresholdModifier(modifier);
    fully_restart_graph_layout();
}


void Canvas::setOptLayeredAlignmentPosition(const LayeredAlignment pos)
{
    m_opt_layered_alignment_position = pos;
    emit optChangedLayeredAlignmentPosition(pos);
    interrupt_graph_layout();
}

void Canvas::setOptFlowSeparationModifier(const double value)
{
    m_flow_separation_modifier = value;
    emit optChangedDirectedEdgeSeparationModifier(value);
    interrupt_graph_layout();
}

void Canvas::setOptFlowSeparationModifierFromSlider(const int intValue)
{
    double doubleValue = intValue / 100.0;
    setOptFlowSeparationModifier(doubleValue);
}

void Canvas::setOptFlowDirection(const FlowDirection value)
{
    m_opt_flow_direction = value;
    emit optChangedFlowDirection(value);
    interrupt_graph_layout();
}

void Canvas::setOptFlowDirectionFromDial(const int value)
{
    setOptFlowDirection((FlowDirection) value);
}

void Canvas::setOptCanvasBackgroundColour(const QColor colour)
{
   m_opt_canvas_background_colour = colour;
   emit canvasDrawingChanged();
}

bool Canvas::hasVisibleOverlays(void) const
{
    return m_overlay_router_raw_routes|| m_overlay_router_display_routes ||
            m_overlay_router_obstacles || m_overlay_router_visgraph ||
            m_overlay_router_orthogonal_visgraph;
}

void Canvas::setRenderingForPrinting(const bool printingMode)
{
    m_rendering_for_printing = printingMode;
}

bool Canvas::isRenderingForPrinting(void) const
{
    return m_rendering_for_printing;
}
void Canvas::setOverlayRouterObstacles(const bool value)
{
    m_overlay_router_obstacles = value;
    emit debugOverlayEnabled(hasVisibleOverlays());
    this->update();
}


bool Canvas::overlayRouterObstacles(void) const
{
    return m_overlay_router_obstacles;
}

void Canvas::setOverlayRouterRawRoutes(const bool value)
{
    m_overlay_router_raw_routes = value;
    emit debugOverlayEnabled(hasVisibleOverlays());
    this->update();
}


bool Canvas::overlayRouterRawRoutes(void) const
{
    return m_overlay_router_raw_routes;
}

void Canvas::setOverlayRouterDisplayRoutes(const bool value)
{
    m_overlay_router_display_routes = value;
    emit debugOverlayEnabled(hasVisibleOverlays());
    this->update();
}


bool Canvas::overlayRouterDisplayRoutes(void) const
{
    return m_overlay_router_display_routes;
}

void Canvas::setOverlayRouterVisGraph(const bool value)
{
    m_overlay_router_visgraph = value;
    emit debugOverlayEnabled(hasVisibleOverlays());
    this->update();
}


bool Canvas::overlayRouterVisGraph(void) const
{
    return m_overlay_router_visgraph;
}


void Canvas::setOverlayRouterOrthogonalVisGraph(const bool value)
{
    m_overlay_router_orthogonal_visgraph = value;
    emit debugOverlayEnabled(hasVisibleOverlays());
    this->update();
}


bool Canvas::overlayRouterOrthogonalVisGraph(void) const
{
    return m_overlay_router_orthogonal_visgraph;
}


double Canvas::optIdealEdgeLengthModifier(void) const
{
    return m_opt_ideal_edge_length_modifier;
}


double Canvas::optSnapDistanceModifier(void) const
{
    return m_opt_snap_distance_modifier;
}

double Canvas::optSnapStrengthModifier(void) const
{
    return m_opt_snap_strength_modifier;
}

double Canvas::optGridWidth(void) const
{
    return m_opt_snap_grid_width;
}

double Canvas::optGridHeight(void) const
{
    return m_opt_snap_grid_height;
}

double Canvas::optRelaxThresholdModifier(void) const
{
    return m_opt_relax_threshold_modifier;
}

int Canvas::optConnectorRoundingDistance(void) const
{
    return m_opt_connector_rounding_distance;
}

double Canvas::optStressBarMaximum(void) const
{
    return m_stress_bar_maximum;
}

double Canvas::optObliquityBarMaximum(void) const
{
    return m_obliquity_bar_maximum;
}


int Canvas::optRoutingPenaltySegment(void) const
{
    return (int) m_router->routingParameter(Avoid::segmentPenalty);
}

double Canvas::optFlowSeparationModifier(void) const
{
    return m_flow_separation_modifier;
}

Canvas::LayeredAlignment Canvas::optLayeredAlignmentPosition(void) const
{
    return (LayeredAlignment) m_opt_layered_alignment_position;
}

void Canvas::bringToFront(void)
{
    if (selectedItems().isEmpty())
    {
        return;
    }

    foreach (CanvasItem *item, selectedItems())
    {
        item->bringToFront();
    }
}


void Canvas::sendToBack(void)
{
    if (selectedItems().isEmpty())
    {
        return;
    }

    foreach (CanvasItem *item, selectedItems())
    {
        item->sendToBack();
    }
}


void Canvas::deselectAll(void)
{
    foreach (CanvasItem *item, selectedItems())
    {
        item->setSelected(false);
    }
}


void Canvas::cutSelection(void)
{
    if (selectedItems().empty())
    {
        return;
    }
    copySelection();
    deleteSelection();
}


void Canvas::copySelection(void)
{
    if (selectedItems().empty())
    {
        return;
    }
    QDomDocument clipboard;

    QDomElement svg = clipboard.createElement("svg");
    clipboard.appendChild(svg);
    newProp(svg, "xmlns", "http://www.w3.org/2000/svg");
    newProp(svg, "xmlns:dunnart", x_dunnartURI);

    //QT selection_box_props(&clip_x, &clip_y, &clip_w, &clip_h);

    foreach (CanvasItem *item, selectedItems())
    {
        QDomElement elem = item->to_QDomElement(XMLSS_ALL, clipboard);

        svg.appendChild(elem);
    }
    m_clipboard = clipboard.toString();
    emit clipboardContentsChanged();
}


void Canvas::pasteSelection(void)
{
    QDomDocument clipboard;
    bool parseNamespaces = true;
    clipboard.setContent(m_clipboard, parseNamespaces);

    if (!clipboard.hasChildNodes())
    {
        // No children, so clipboard is empty.
        return;
    }

    beginUndoMacro(tr("Paste"));

    this->deselectAll();

    QString dunnartNs = x_dunnartNs;

    m_paste_id_map.clear();
    m_paste_bad_constraint_ids.clear();

    // Assign new clipboard IDs.
    recursiveMapIDs(clipboard, dunnartNs, PASTE_UPDATEOBJIDS);

    // Update IDs for connectors and relationships.
    recursiveMapIDs(clipboard, dunnartNs, PASTE_UPDATEIDPROPS);

    // Find bad distributions and separations.
    recursiveMapIDs(clipboard, dunnartNs, PASTE_FINDBADDISTROS);

    // Remove bad distributions and separations.
    recursiveMapIDs(clipboard, dunnartNs, PASTE_REMOVEBADDISTROS);

    qDebug() << clipboard.toString();
    // Actually do the pasting, in correct order.
    for (int pass = 0; pass < PASS_LAST; ++pass)
    {
        this->recursiveReadSVG(clipboard, dunnartNs, pass);
    }

    // Select all new shapes.
    recursiveMapIDs(clipboard, dunnartNs, PASTE_SELECTSHAPES);

    // Find the centre of pasted items, so we know how much to move them.
    QPointF oldCentrePos = diagramBoundingRect(selectedItems()).center();

    // If the cursor is inside the canvas, then paste the objects centred
    // under the cursor, otherwise paste to the centre of the visible canvas.
    QGraphicsView *view = views().first();
    QList<CanvasItem *> selected_items = selectedItems();
    QPointF pastePosition;
    if (view->underMouse())
    {
        pastePosition = view->mapToScene(
                view->mapFromGlobal(QCursor::pos()));
    }
    else
    {
        pastePosition = QRectF(view->mapToScene(0,0),
                view->mapToScene(view->width(), view->height())).center();
    }

    // Move the new selection to paste position.
    QPointF difference = pastePosition - oldCentrePos;
    foreach (CanvasItem *item, selected_items)
    {
        item->moveBy(difference.x(), difference.y());
    }

    // Put the distribution indicators at their default positions:
    foreach (CanvasItem *item, selected_items)
    {
        Distribution *distro = dynamic_cast<Distribution *> (item);
        Separation *sep = dynamic_cast<Separation *> (item);
        bool store_undo = false;

        if (distro)
        {
            distro->moveToDefaultPos(store_undo);
        }
        else if (sep)
        {
            sep->moveToDefaultPos(store_undo);
        }
    }
    interrupt_graph_layout();
    restart_graph_layout();
}

void Canvas::deleteSelection(void)
{
    QList<CanvasItem *> sel_copy = this->selectedItems();
    deleteItems(sel_copy);
    restart_graph_layout();
}

void Canvas::deleteItem(CanvasItem *item)
{
    qDebug() << "deleting item: " << item;
    QList<CanvasItem*> items;
    items.append(item);
    deleteItems(items);
}

void Canvas::deleteItems(QList<CanvasItem*> items)
{
    if (items.empty())
    {
        return;
    }

    // Finish dragging if we were in the middle of that.
    if (m_dragged_item)
    {
        m_dragged_item->dragReleaseEvent(NULL);
    }

    Actions& actions = getActions();
    actions.clear();

    UndoMacro *undoMacro = beginUndoMacro(tr("Delete"));

    // Disconnect all shape--connector connections where one is inside the
    // set and one is outside the set.
    CanvasItemSet selSet;
    QList<ShapeObj *> sel_shapes;

    stop_graph_layout();
    if (m_opt_stuctural_editing_disabled)
    {
        // If structural editing is disabled then we should only allow
        // deletion of constraint indicators.
        for (QList<CanvasItem *>::iterator sh = items.begin();
                sh != items.end(); )
        {
            Indicator *indicator = dynamic_cast<Indicator *> (*sh);
            if (indicator)
            {
                // Indicator, so leave in list.
                sh++;
            }
            else
            {
                // Not an idicator, so remove from list.
                sh = items.erase(sh);
            }
        }
    }

    // Do we need to deselectAll? This is a vestige from when this code
    // belonged to the deleteSelection method.
    this->deselectAll();
    //

    // Do distro's first, in case they have guidelines and distros selected.
    for (QList<CanvasItem *>::iterator sh = items.begin();
            sh != items.end(); ++sh)
    {
        Distribution *distro = dynamic_cast<Distribution *> (*sh);
        Separation *sep = dynamic_cast<Separation *> (*sh);

        if (distro || sep)
        {
            // Deactivate constraints:
            (*sh)->deactivateAll(selSet);
        }
    }

    for (QList<CanvasItem *>::iterator sh = items.begin();
            sh != items.end(); ++sh)
    {
        Connector  *conn  = dynamic_cast<Connector *>  (*sh);
        ShapeObj *shape = dynamic_cast<ShapeObj *> (*sh);
        Guideline *guide = dynamic_cast<Guideline *> (*sh);

        if (conn || guide || shape)
        {
            (*sh)->deactivateAll(selSet);
        }
        if (shape)
        {
            sel_shapes.push_back(shape);
        }
    }

    // Remove containment relationships.
    QList<CanvasItem *> canvas_items = this->items();
    for (int i = 0; i < canvas_items.size(); ++i)
    {
        ShapeObj *shape = dynamic_cast<ShapeObj *> (canvas_items.at(i));
        if (shape)
        {
            shape->removeContainedShapes(sel_shapes);
        }
    }

    for (QList<CanvasItem *>::iterator sh = items.begin();
            sh != items.end(); ++sh)
    {
        QUndoCommand *cmd = new CmdCanvasSceneRemoveItem(this, *sh);
        undoMacro->addCommand(cmd);
    }

    //restart_graph_layout();
}



/*
void Canvas::deleteSelection(void)
{
    if (selectedItems().empty())
    {
        return;
    }

    // Finish dragging if we were in the middle of that.
    if (m_dragged_item)
    {
        m_dragged_item->dragReleaseEvent(NULL);
    }

    Actions& actions = getActions();
    actions.clear();

    UndoMacro *undoMacro = beginUndoMacro(tr("Delete"));

    // Disconnect all shape--connector connections where one is inside the
    // selection and one is outside the selection.
    CanvasItemSet selSet;
    // Use a copy since the deselect action will change the selection.
    QList<CanvasItem *> sel_copy = this->selectedItems();
    QList<ShapeObj *> sel_shapes;

    stop_graph_layout();
    if (m_opt_stuctural_editing_disabled)
    {
        // If structural editing is disabled then we should only allow
        // deletion of constraint indicators.
        for (QList<CanvasItem *>::iterator sh = sel_copy.begin();
                sh != sel_copy.end(); )
        {
            Indicator *indicator = dynamic_cast<Indicator *> (*sh);
            if (indicator)
            {
                // Indicator, so leave in list.
                sh++;
            }
            else
            {
                // Not an idicator, so remove from list.
                sh = sel_copy.erase(sh);
            }
        }
    }

    this->deselectAll();

    // Do distro's first, incase they have guidelines and distros selected.
    foreach (CanvasItem *item, sel_copy)
    {
        Distribution *distro = dynamic_cast<Distribution *> (item);
        Separation *sep = dynamic_cast<Separation *> (item);

        if (distro || sep)
        {
            // Deactivate constraints:
            item->deactivateAll(selSet);
        }
    }
    foreach (CanvasItem *item, sel_copy)
    {
        Connector  *conn  = dynamic_cast<Connector *>  (item);
        ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
        Guideline *guide = dynamic_cast<Guideline *> (item);
        
        if (conn || guide || shape)
        {
            item->deactivateAll(selSet);
        }
        if (shape)
        {
            sel_shapes.push_back(shape); 
        }
    }
    
    // Remove containment relationships.
    foreach (CanvasItem *item, items())
    {
        ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
        if (shape)
        {
            shape->removeContainedShapes(sel_shapes);
        }
    }

    foreach (CanvasItem *item, sel_copy)
    {
        QUndoCommand *cmd = new CmdCanvasSceneRemoveItem(this, item);
        undoMacro->addCommand(cmd);
    }
    restart_graph_layout();
}
*/

// Delay, in milliseconds, to give the event loop time to respond to normal
// events like mouse movements.  A delay of zero will flood the event queue.
static const unsigned int updateEventDelay = 15;

void Canvas::startLayoutUpdateTimer(void)
{
    if (!m_layout_update_timer)
    {
        m_layout_update_timer = new QTimer(this);
        connect(m_layout_update_timer, SIGNAL(timeout()), this,
                SLOT(processLayoutUpdateEvent()));
    }
    if (!m_layout_update_timer->isActive())
    {
        m_layout_update_timer->start(updateEventDelay);
    }
}


void Canvas::startLayoutFinishTimer(void)
{
    if (!m_layout_finish_timer)
    {
        m_layout_finish_timer = new QTimer(this);
        connect(m_layout_finish_timer, SIGNAL(timeout()), this,
                SLOT(processLayoutFinishedEvent()));
    }
    m_layout_finish_timer->start(updateEventDelay);
}


void Canvas::setFilename(QString filename)
{
    m_filename = filename;
    QFileInfo info(m_filename);
    emit diagramFilenameChanged(info);
}


QString Canvas::filename(void)
{
    return m_filename;
}


QList<CanvasItem *> Canvas::items(void) const
{
    QList<CanvasItem *> canvasItems;

    foreach (QGraphicsItem *item, QGraphicsScene::items())
    {
        CanvasItem *canvasItem = dynamic_cast<CanvasItem *> (item);
        if (canvasItem)
        {
            canvasItems.push_back(canvasItem);
        }
    }

    return canvasItems;
}


QList<CanvasItem *> Canvas::selectedItems(void) const
{
    QList<CanvasItem *> filteredSelection;
    
    // Filter and return just the CanvasItem-based objects.
    foreach (QGraphicsItem *item, QGraphicsScene::selectedItems())
    {
        CanvasItem *canvasItem = dynamic_cast<CanvasItem *> (item);
        if (canvasItem)
        {
            filteredSelection.push_back(canvasItem);
        }
    }

    if (m_edit_mode == ModeConnection)
    {
        // In connection mode we want to allow selection of lone connectors
        // for editing purposes, but not selection of other objects.
        if ((filteredSelection.size() == 1) &&
                dynamic_cast<Connector *> (filteredSelection.first()))
        {
            // Return the single selected connector
            return filteredSelection;
        }
        else
        {
            // Otherwise, return no selection.
            QList<CanvasItem *> emptySelection;
            return emptySelection;
        }
    }

    return filteredSelection;
}

bool Canvas::useGmlClusters(void) const
{
    return m_use_gml_clusters;
}

void Canvas::setNudgeDistance(const double dist)
{
    m_connector_nudge_distance = dist;
}

void Canvas::updateConnectorsForLayout(void)
{
    if ((!m_opt_preserve_topology || (m_graphlayout->runLevel != 1)) &&
            !m_gml_graph && !m_batch_diagram_layout)
    {
        // Don't reroute connectors in the case of topology preserving layout.
        reroute_connectors(this);
    }
}

void Canvas::processLayoutUpdateEvent(void)
{
    if (m_action_finished_timer)
    {
        m_action_finished_timer->start();
    }

    //qDebug("LayoutUpdateEvent");
#ifdef FPSTIMER
    if (!m_convergence_timer_running)
    {
        m_convergence_timer.start();
        m_convergence_timer_running = true;
    }
    m_convergence_update_count++;
#endif
    foreach (CanvasItem *item, items())
    {
        if (item && item->constraintConflict())
        {
            item->setConstraintConflict(false);
        }
    }

    bool success = m_graphlayout->processReturnPositions();
    if (success == false)
    {
        // finished processing for the moment.
        m_layout_update_timer->stop();
    }
    updateConnectorsForLayout();

    //qDebug("processLayoutUpdateEvent %7d", ++layoutUpdates);
}


static QRectF expandRect(const QRectF& origRect, double amount)
{
    return (origRect.isEmpty()) ?
            origRect : origRect.adjusted(-amount, -amount, amount, amount);
}

#define ACTION_TIMER_DELAY 500

void Canvas::processLayoutFinishedEvent(void)
{
    m_layout_finish_timer->stop();

    //qDebug("LayoutFinishedEvent");
#ifdef FPSTIMER
    if (m_convergence_timer_running)
    {
        double elapsedSecs = m_convergence_timer.elapsed() / 1000.0;

        qDebug("************** Avg Framerate: %g\n", m_convergence_update_count / elapsedSecs);
        qDebug("************** Time to converge: %g\n", elapsedSecs);
        m_convergence_update_count = 0;
        m_convergence_timer_running = false;
    }
#endif
    /*
    if (!straighten_bends &&
            automatic_graph_layout && simple_paths_during_layout)
    {
        // Only bother doing this if automatic graph layout is on
        // and we would have been drawing simple paths up until this
        // point.
        bool lastSimpleRouting = router->SimpleRouting;
        router->SimpleRouting = false;
        reroute_connectors(this);
        router->SimpleRouting = lastSimpleRouting;
    }
    repaint_canvas();
    */
    GraphLayout* gl=m_graphlayout;
    if (!m_opt_automatic_graph_layout && !m_batch_diagram_layout)
    {
        // Do connector post-processing.
        reroute_connectors(this, false, true);
    }
    bool changes = !getActions().empty();
    if (gl->runLevel == 0)
    {
        gl->runLevel=1;
        //qDebug("runLevel=1");
        interrupt_graph_layout();
        changes = true;
    }

    if ((m_dragged_item == NULL) &&  m_dragged_with_force)
    {
        m_dragged_with_force = false;

        // Let the layout notice there is no longer a dragged item
        // and let it converge again.
        processResponseTasks();
        changes = true;
    }

    emit layoutHasConverged();

    if (m_batch_diagram_layout && !changes)
    {
        // Reroute connectors.
        reroute_connectors(this, true, true);
        // Nudge connectors if requested (-w option)
        if (m_connector_nudge_distance > 0)
        {
            nudgeConnectors(this, m_connector_nudge_distance, true);
        }
        // redo the connector interference coloring after nudging
        if (m_opt_colour_interfering_connectors)
        {
            colourInterferingConnectors(this);
        }
        // Fit the page size to the entire diagram.
        //QT getPageSize(NULL, NULL, NULL, NULL, BUT_FITPAGETODIAGRAM);
        // Save the SVG and exit.
        //QT saveDiagramAsSVG(this, filename());
        exit(EXIT_SUCCESS);
    }
//#define showObjectiveBars
#ifdef  showObjectiveBars
    if (!changes) {
        computeOrthoObjective();
    }
#endif

    if (m_trying_alignments) {
#define newApplyAlignments
#ifdef newApplyAlignments
        applyAlignmentsCallback();
#else
        tryAlignments();
#endif
    }
    else if (gd2013_batch_mode)
    {
        if (!m_action_finished_timer)
        {
            m_action_finished_timer = new QTimer(this);
            connect(m_action_finished_timer, SIGNAL(timeout()), this,
                    SLOT(actionFinished()));
        }
        m_action_finished_timer->start(ACTION_TIMER_DELAY);
    }
}

void Canvas::actionFinished(void)
{
    m_action_finished_timer->stop();
    ++gd2013_step;
    qDebug() << "Type: " << gd2013_type << "  Step: " << gd2013_step;
    if (gd2013_type == 1)
    {
        m_opt_snap_strength_modifier = 150;

        if (gd2013_step == 1)
        {
            // Clear the output file.
            QFileInfo fileBase(filename());
            QString outputFile = QString("/Users/mjwybrow/Desktop/Orthogonal-paper/%1.txt").arg(fileBase.completeBaseName());

            QFile file(outputFile);
            if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
                return;
            }
            m_action_timer.start();
            setOptAutomaticGraphLayout(true);
        }
        else if (gd2013_step == 2)
        {
            setOptPreventOverlaps(true);
        }
        else if (gd2013_step == 3)
        {
            QFileInfo file(filename());
            paperExport(file.completeBaseName(), "fd", 2);
            setOptGridSnap(true);
            m_action_timer.start();
        }
        else if (gd2013_step == 4)
        {
            QFileInfo file(filename());
            paperExport(file.completeBaseName(), "grid");
            exit(0);
        }
    }
    if (gd2013_type == 2)
    {
        if (gd2013_step == 1)
        {
            setOptAutomaticGraphLayout(true);
        }
        else if (gd2013_step == 2)
        {
            setOptPreventOverlaps(true);
        }
        else if (gd2013_step == 3)
        {
            setOptSnapTo(true);
            m_action_timer.start();
        }
        else if (gd2013_step == 4)
        {
            QFileInfo file(filename());
            paperExport(file.completeBaseName(), "node");
            setOptGridSnap(true);
            m_action_timer.start();
        }
        else if (gd2013_step == 5)
        {
            QFileInfo file(filename());
            paperExport(file.completeBaseName(), "node+grid");
            exit(0);
        }
    }
    if (gd2013_type == 3)
    {
        if (gd2013_step == 1)
        {
            setOptAutomaticGraphLayout(true);
        }
        else if (gd2013_step == 2)
        {
            initTryAlignments();
            m_action_timer.start();
        }
        else if (gd2013_step == 3)
        {
            QFileInfo file(filename());
            paperExport(file.completeBaseName(), "aca");
            setOptGridSnap(true);
            m_action_timer.start();
        }
        else if (gd2013_step == 4)
        {
            setOptPreventOverlaps(true);
        }
        else if (gd2013_step == 5)
        {
            QFileInfo file(filename());
            paperExport(file.completeBaseName(), "aca+grid");
            exit(0);
        }
    }
}

void Canvas::paperExport(QString diagram, QString type, int actions)
{
    double elapsedSecs = (m_action_timer.elapsed() -
            (ACTION_TIMER_DELAY * actions)) / 1000.0;

    QString fileBase = QString("%1-%2").arg(diagram).arg(type);
    setPageRect(
            expandRect(diagramBoundingRect(
                    items()), 8));
    QString outputFile = QString("/Users/mjwybrow/Desktop/Orthogonal-paper/%1.pdf").arg(fileBase);
    qDebug() << outputFile;
//            saveDiagram("");
    exportDiagramToFile(outputFile);
    outputFile = QString("/Users/mjwybrow/Desktop/Orthogonal-paper/%1.svg").arg(fileBase);
    saveDiagram(outputFile);

    outputFile = QString("/Users/mjwybrow/Desktop/Orthogonal-paper/%1.txt").arg(diagram);
    QFile file(outputFile);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Text))
    {
        return;
    }

    uint nodes = 0;
    uint edges = 0;

    foreach (CanvasItem *item, items())
    {
        ShapeObj *node = dynamic_cast<ShapeObj *> (item);
        Connector *edge = dynamic_cast<Connector *> (item);

        if (node)
        {
            ++nodes;
        }
        else if (edge)
        {
            ++edges;
        }
    }

    computeOrthoObjective();

    QTextStream out(&file);
#if 1
    out << fileBase << ","
        << nodes << ","
        << edges << ","
        << elapsedSecs << ","
        << m_most_recent_stress << ","
        << m_most_recent_coincidence_count << ","
        << m_most_recent_coincidence_group_count << ","
        << m_most_recent_crossing_count << ","
        << m_most_recent_angle_res_score << ","
        << m_most_recent_ang_res_by_degree.at(0) << ","
        << m_most_recent_ang_res_by_degree.at(1) << ","
        << m_most_recent_ang_res_by_degree.at(2) << ","
        << m_most_recent_average_obliqueness << ","
        << m_most_recent_avg_grid_distance << ","
        << m_most_recent_edge_node_overlap_count << ","
        << m_most_recent_ortho_obj_func << "\n";
#else
    out << "\n" << fileBase << "\n\n";
    out << "Nodes:\t\t" << nodes << "\n";
    out << "Edges:\t\t" << edges << "\n";
    out << "Time:\t\t" << elapsedSecs << "\n";
    out << "Stress:\t\t" << m_most_recent_stress << "\n";
    out << "Coincidences:\t" << m_most_recent_coincidence_count << "\n";
    out << "Coincidence groups:\t" << m_most_recent_coincidence_group_count << "\n";
    out << "Crossings:\t" << m_most_recent_crossing_count << "\n";
    out << "Average Angular Resolution:\t" << m_most_recent_angle_res_score << "\n";
    out << "Average Ang. Res. for deg-2 nodes:\t" << m_most_recent_ang_res_by_degree.at(0) << "\n";
    out << "Average Ang. Res. for deg-3 nodes:\t" << m_most_recent_ang_res_by_degree.at(1) << "\n";
    out << "Average Ang. Res. for deg-4 nodes:\t" << m_most_recent_ang_res_by_degree.at(2) << "\n";
    out << "Average Obliqueness:\t" << m_most_recent_average_obliqueness << "\n";
    out << "Average Grid Distance:\t" << m_most_recent_avg_grid_distance << "\n";
    out << "Edge-Node Overlaps:\t" << m_most_recent_edge_node_overlap_count << "\n";
    out << "Score:\t\t" << m_most_recent_ortho_obj_func << "\n";
#endif
}

void Canvas::exportDiagramToFile(QString filename)
{
    QFileInfo file(filename);
    if (file.suffix() == "svg")
    {
        QSvgGenerator generator;
        generator.setFileName(filename);

        generator.setSize(pageRect().size().toSize());
        QRectF targetRect(QPointF(0, 0), QSizeF(sceneRect().size()));

        QRectF viewbox(pageRect().topLeft() - sceneRect().topLeft(),
                pageRect().size());
        generator.setViewBox(viewbox);

        generator.setTitle(QFileInfo(filename).fileName());
        generator.setDescription(tr("This file was exported from Dunnart.  "
                                    "http://www.dunnart.org/"));

        QPainter painter;
        if (painter.begin(&generator))
        {
            painter.setRenderHint(QPainter::Antialiasing);
            setRenderingForPrinting(true);
            render(&painter, targetRect, sceneRect(), Qt::IgnoreAspectRatio);
            setRenderingForPrinting(false);

            painter.end();
        }
        else
        {
            qDebug("Export SVG painter failed to begin.");
        }
    }
    else
    {
        // Use QPrinter for PDF and PS.
        QPrinter printer;
        printer.setOutputFileName(filename);
        printer.setPaperSize(pageRect().size(), QPrinter::Millimeter);
        QPainter painter;
        if (painter.begin(&printer))
        {
            painter.setRenderHint(QPainter::Antialiasing);
            setRenderingForPrinting(true);
            render(&painter, QRectF(), pageRect().adjusted(+3, +3, -3, -3),
                    Qt::IgnoreAspectRatio);
            setRenderingForPrinting(false);
        }
        else
        {
            qDebug("Export PDF/PS painter failed to begin.");
        }
    }
}

QSvgRenderer *Canvas::svgRenderer(void) const
{
    return m_svg_renderer;
}


void Canvas::setLayoutSuspended(bool suspend)
{
    if (suspend)
    {
        // Stop the automatic layout.
        m_graphlayout->setLayoutSuspended(true);
    }
    else
    {
        // Restart the automatic layout.
        // Make it look like things just moved.
        Actions& actions = getActions();
        foreach (CanvasItem *item, selectedItems())
        {
            actions.moveList.push_back(item);
        }

        m_graphlayout->setLayoutSuspended(false);
        this->interrupt_graph_layout();
        this->restart_graph_layout();
    }
}

bool Canvas::isLayoutSuspended(void) const
{
    return m_graphlayout->isFreeShiftFromDunnart();
}

void Canvas::createIndicatorHighlightCache(void)
{
    if (!isLayoutSuspended())
    {
        // User is not holding ALT, so don't automatically attach to guidelines.
        return;
    }

    foreach (CanvasItem *item, items())
    {
        Guideline *g = dynamic_cast<Guideline *> (item);
        if (g && (g->isSelected() == false))
        {
            // We don't want guides that are being moved, because they
            // are attached via multi-way constraints to selected shapes.
            bool invalid = false;
            for (RelsList::iterator r = g->relationships.begin(); r != g->relationships.end();
                    r++)
            {
                if ((*r)->shape && (*r)->shape->isSelected())
                {
                    invalid = true;
                    break;
                }
            }
            if (invalid)
            {
                continue;
            }

            if (g->get_dir() == GUIDE_TYPE_VERT)
            {
                double pos = g->x();
                m_vguides[pos - 1] = g;
                m_vguides[pos] = g;
                m_vguides[pos + 1] = g;
            }
            else
            {
                double pos = g->y();
                m_hguides[pos - 1] = g;
                m_hguides[pos] = g;
                m_hguides[pos + 1] = g;
            }
        }
    }
}


void Canvas::highlightIndicatorsForItemMove(CanvasItem *item)
{
    bool vfound = false, hfound = false;
    ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
    if (shape && (shape->canvasItemFlags() & CanvasItem::ItemIsAlignable))
    {
        for (int i = 0; i < 6; i++)
        {
            if (!(shape->rels[i]))
            {
                double pos = shape->attachedGuidelinePosition((atypes) i);

                if (!hfound && (i < 3))
                {
                    if (m_hguides.find(pos) != m_hguides.end())
                    {
                        m_hguides[pos]->setHighlighted(true);
                        m_hguides[pos]->update();
                    }
                }
                else if (!vfound && (i >= 3))
                {
                    if (m_vguides.find(pos) != m_vguides.end())
                    {
                        m_vguides[pos]->setHighlighted(true);
                        m_vguides[pos]->update();
                    }
                }
            }
        }
    }
}

bool Canvas::inSelectionMode(void) const
{
    return (m_edit_mode == ModeSelection);
}

void Canvas::clearIndicatorHighlights(const bool clearCache)
{
    foreach (CanvasItem *item, items())
    {
        Guideline *g = dynamic_cast<Guideline *> (item);
        if (g)
        {
            if (g->isHighlighted())
            {
                g->setHighlighted(false);
                g->update();
            }
        }
    }

    if (clearCache)
    {
        m_hguides.clear();
        m_vguides.clear();
    }
}


void Canvas::glueObjectsToIndicators(void)
{
    if (!isLayoutSuspended())
    {
        // User is not holding ALT, so don't automatically attach to guidelines.
        return;
    }

    foreach (CanvasItem *item, selectedItems())
    {
        bool vfound = false, hfound = false;
        ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
        if (shape && (shape->canvasItemFlags() & CanvasItem::ItemIsAlignable))
        {
            for (int i = 0; i < 6; i++)
            {
                if (!(shape->rels[i]))
                {
                    double pos = shape->attachedGuidelinePosition((atypes) i);

                    if (!hfound && (i < 3))
                    {
                        if (m_hguides.find(pos) != m_hguides.end())
                        {
                            new Relationship(m_hguides[pos], shape, (atypes) i);
                        }
                    }
                    else if (!vfound && (i >= 3))
                    {
                        if (m_vguides.find(pos) != m_vguides.end())
                        {
                            new Relationship(m_vguides[pos], shape, (atypes) i);
                        }
                    }
                }
            }
        }
    }
    m_hguides.clear();
    m_vguides.clear();
}

void Canvas::processResponseTasks(void)
{
    GraphLayout* gl = m_graphlayout;
    gl->setRestartFromDunnart();
    gl->apply(!m_opt_automatic_graph_layout);

    getActions().clear();
}

void Canvas::processUndoResponseTasks(void)
{
    getActions().clear();
    repositionAndShowSelectionResizeHandles(true);
    stop_graph_layout();
    reroute_connectors(this);
}

void Canvas::fully_restart_graph_layout(void)
{
    GraphLayout* gl = m_graphlayout;
    gl->setInterruptFromDunnart();
    gl->setRestartFromDunnart();
    gl->runLevel=0;
    if (gl->mode == GraphLayout::LAYERED) 
    {
        //gl->runLevel=1;
    }
    gl->apply(!m_opt_automatic_graph_layout);
}

void Canvas::restart_graph_layout(void)
{
    GraphLayout* gl = m_graphlayout;
    gl->setRestartFromDunnart();
    gl->apply(!m_opt_automatic_graph_layout);
}

void Canvas::stop_graph_layout(void)
{
    GraphLayout* gl = m_graphlayout;
    gl->setInterruptFromDunnart();
}


void Canvas::interrupt_graph_layout(void)
{
    GraphLayout* gl = m_graphlayout;
    gl->setInterruptFromDunnart();
    gl->apply(!m_opt_automatic_graph_layout);
}


void Canvas::selectionChangeTriggers(void)
{
    int shapeCount = 0;
    int indicatorCount = 0;

    this->hideSelectionResizeHandles();

    // Update the area covered by the selection and the resize handles.
    this->update(m_selection_shapes_bounding_rect.adjusted(
            -BOUNDINGRECTPADDING, -BOUNDINGRECTPADDING,
            +BOUNDINGRECTPADDING, BOUNDINGRECTPADDING));
    // Then reset the selection shape's bounding Rectangle.
    m_selection_shapes_bounding_rect = QRectF();

    QList<CanvasItem *> selected_items = selectedItems();
    foreach (CanvasItem *item, selected_items)
    {
        ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
        if (shape)
        {
            if (!shape->sizeLocked())
            {
                // Build the bounding rectangle from the union of all
                // shape's boundingRects in the selection.
                m_selection_shapes_bounding_rect =
                        m_selection_shapes_bounding_rect.united(
                        shape->boundingRect().translated(shape->scenePos()));
            }
            shapeCount++;
        }
        else if (dynamic_cast<Indicator *> (item))
        {
            indicatorCount++;
        }
    }

    // Remove the boundingRect padding.
    m_selection_shapes_bounding_rect = m_selection_shapes_bounding_rect.adjusted(
            +BOUNDINGRECTPADDING, +BOUNDINGRECTPADDING,
            -BOUNDINGRECTPADDING, -BOUNDINGRECTPADDING);

    if (shapeCount > 0)
    {
        this->repositionAndShowSelectionResizeHandles();
    }

    CanvasItem *new_lone_selected_item = NULL;
    if (selected_items.size() == 1)
    {
        new_lone_selected_item = selected_items.front();
    }
    if (m_lone_selected_item != new_lone_selected_item)
    {
        if (m_lone_selected_item)
        {
            m_lone_selected_item->loneSelectedChange(false);
        }
        m_lone_selected_item = new_lone_selected_item;
        if (m_lone_selected_item)
        {
            m_lone_selected_item->loneSelectedChange(true);
        }
    }

#if 0
    if (queryMode)
    {
        // Nicely handle selection changes while we are in pair_query mode.
        if (selection.size() == 1)
        {
            // We have one object selected. Set queryObj if it is a shape.
            queryObj = dynamic_cast<CanvasItem *> (selection.front());
            if (queryObj)
            {
                if (dynamic_cast<CanvasItem *> (active_obj))
                {
                    // Query the shape under the mouse if there is one.
                    pair_query(active_obj);
                }
                else if (active_obj)
                {
                    // This case copes with us being over a handle.
                    // Query the parent of the handle, i.e, the shape.
                    if (dynamic_cast<CanvasItem *> (active_obj->get_parent()))
                    {
                        pair_query(active_obj->get_parent());
                    }
                }
            }
            else
            {
                // Not a queryObj, so reset illumination.
                resetQueryModeIllumination(false);
            }
        }
        else
        {
            // Not a single object in the selection, so reset illumination.
            bool resetQueryObj = true;
            resetQueryModeIllumination(resetQueryObj);
        }
    }
#endif
}


void Canvas::storeSelectionResizeInfo(void)
{
    // Store info about selected shape position and dimensions in
    // relation to the boundingRect of the selection.  This way, when the
    // selection is resized, we can use this information to resize each
    // of the selected shapes.
    QPointF selectionTopLeft = m_selection_shapes_bounding_rect.topLeft();
    QPointF selectionDimensions =
            m_selection_shapes_bounding_rect.bottomRight() - selectionTopLeft;
    assert(selectionDimensions.x() >= 0);
    assert(selectionDimensions.y() >= 0);
    QList<CanvasItem *> selected_items = selectedItems();
    m_selection_resize_info = QVector<QRectF>(selected_items.size());
    size_t index = 0;
    foreach (CanvasItem *item, selected_items)
    {
        ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
        if (shape && !shape->sizeLocked())
        {
            QRectF shapeBR = shape->boundingRect().adjusted(
                    +BOUNDINGRECTPADDING, +BOUNDINGRECTPADDING,
                    -BOUNDINGRECTPADDING, -BOUNDINGRECTPADDING);
            shapeBR = shapeBR.translated(shape->scenePos());
            QPointF topLeft = (shapeBR.topLeft() - selectionTopLeft);
            topLeft = QPointF(topLeft.x() / selectionDimensions.x(),
                              topLeft.y() / selectionDimensions.y());
            QPointF bottomRight =
                    (shapeBR.bottomRight() - selectionTopLeft);
            bottomRight = QPointF(bottomRight.x() / selectionDimensions.x(),
                                  bottomRight.y() / selectionDimensions.y());
            m_selection_resize_info[index] = QRectF(topLeft, bottomRight);
        }
        ++index;
    }
}

void Canvas::moveSelectionResizeHandle(const int index, const QPointF pos)
{

    bool aroundCentre =
            (QApplication::keyboardModifiers() & Qt::MetaModifier);
    QPointF oppositePos = m_selection_shapes_bounding_rect.center() - (pos -
            m_selection_shapes_bounding_rect.center());
    switch (index)
    {
    case HAND_TOP_LEFT:
        m_selection_shapes_bounding_rect.setTopLeft(pos);
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setBottomRight(oppositePos);
        }
        break;
    case HAND_TOP_CENTRE:
        m_selection_shapes_bounding_rect.setTop(pos.y());
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setBottom(oppositePos.y());
        }
        break;
    case HAND_TOP_RIGHT:
        m_selection_shapes_bounding_rect.setTopRight(pos);
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setBottomLeft(oppositePos);
        }
        break;
    case HAND_RIGHT_CENTRE:
        m_selection_shapes_bounding_rect.setRight(pos.x());
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setLeft(oppositePos.x());
        }
        break;
    case HAND_BOTTOM_RIGHT:
        m_selection_shapes_bounding_rect.setBottomRight(pos);
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setTopLeft(oppositePos);
        }
        break;
    case HAND_BOTTOM_CENTRE:
        m_selection_shapes_bounding_rect.setBottom(pos.y());
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setTop(oppositePos.y());
        }
        break;
    case HAND_BOTTOM_LEFT:
        m_selection_shapes_bounding_rect.setBottomLeft(pos);
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setTopRight(oppositePos);
        }
        break;
    case HAND_LEFT_CENTRE:
        m_selection_shapes_bounding_rect.setLeft(pos.x());
        if (aroundCentre)
        {
            m_selection_shapes_bounding_rect.setRight(oppositePos.x());
        }
        break;
    default:
        break;
    }

    m_selection_shapes_bounding_rect = m_selection_shapes_bounding_rect.normalized();

    // Update the resize handle positions.
    this->repositionAndShowSelectionResizeHandles();

    // Calculate and apply new positions and dimensions for all shapes in
    // the selection.
    QPointF selectionTopLeft = m_selection_shapes_bounding_rect.topLeft();
    QPointF selectionDimensions =
            m_selection_shapes_bounding_rect.bottomRight() - selectionTopLeft;
    size_t ind = 0;
    foreach (CanvasItem *item, selectedItems())
    {
        ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
        if (shape && !shape->sizeLocked())
        {
            QPointF topLeft = m_selection_resize_info[ind].topLeft();
            QPointF bottomRight = m_selection_resize_info[ind].bottomRight();

            topLeft = QPointF(topLeft.x() * selectionDimensions.x(),
                              topLeft.y() * selectionDimensions.y());
            bottomRight = QPointF(bottomRight.x() * selectionDimensions.x(),
                                  bottomRight.y() * selectionDimensions.y());

            topLeft += selectionTopLeft;
            bottomRight += selectionTopLeft;

            QRectF newSize(topLeft, bottomRight);
            shape->setPosAndSize(newSize.center(),
                    QSizeF(qMax(newSize.width(), MIN_SHAPE_SIZE),
                           qMax(newSize.height(), MIN_SHAPE_SIZE)));
        }
        ++ind;
    }

    // Cause the layout engine to notice changes to shape sizes.
    this->interrupt_graph_layout();
    this->restart_graph_layout();
}


void Canvas::hideSelectionResizeHandles(void)
{
    foreach (SelectionResizeHandle *handle, m_selection_resize_handles)
    {
        handle->setVisible(false);
    }
}

void Canvas::repositionAndShowSelectionResizeHandles(bool calculatePosition)
{
    if (m_opt_stuctural_editing_disabled)
    {
        // Don't show resize handles if structural editing is disabled.
        return;
    }

    if (m_hide_selection_handles)
    {
        return;
    }

    if (calculatePosition)
    {
        m_selection_shapes_bounding_rect = QRectF();

        foreach (CanvasItem *item, selectedItems())
        {
            ShapeObj *shape = dynamic_cast<ShapeObj *> (item);
            if (shape && !shape->sizeLocked())
            {
                // Build the bounding rectangle from the union of all shape's
                // boundingRects in the selection.
                m_selection_shapes_bounding_rect = m_selection_shapes_bounding_rect.united(
                        shape->boundingRect().translated(shape->scenePos()));
            }
        }
        // Remove the boundingRect padding.
        m_selection_shapes_bounding_rect = m_selection_shapes_bounding_rect.adjusted(
                +BOUNDINGRECTPADDING, +BOUNDINGRECTPADDING,
                -BOUNDINGRECTPADDING, -BOUNDINGRECTPADDING);
    }

    if (m_selection_shapes_bounding_rect.isEmpty())
    {
        return;
    }

    // Reposition resize handles.
    m_selection_resize_handles[HAND_TOP_LEFT]->setPos(
            m_selection_shapes_bounding_rect.topLeft());
    m_selection_resize_handles[HAND_TOP_CENTRE]->setPos(
            m_selection_shapes_bounding_rect.center().x(),
            m_selection_shapes_bounding_rect.top());
    m_selection_resize_handles[HAND_TOP_RIGHT]->setPos(
            m_selection_shapes_bounding_rect.topRight());
    m_selection_resize_handles[HAND_RIGHT_CENTRE]->setPos(
            m_selection_shapes_bounding_rect.right(),
            m_selection_shapes_bounding_rect.center().y());
    m_selection_resize_handles[HAND_BOTTOM_RIGHT]->setPos(
            m_selection_shapes_bounding_rect.bottomRight());
    m_selection_resize_handles[HAND_BOTTOM_CENTRE]->setPos(
            m_selection_shapes_bounding_rect.center().x(),
            m_selection_shapes_bounding_rect.bottom());
    m_selection_resize_handles[HAND_BOTTOM_LEFT]->setPos(
            m_selection_shapes_bounding_rect.bottomLeft());
    m_selection_resize_handles[HAND_LEFT_CENTRE]->setPos(
            m_selection_shapes_bounding_rect.left(),
            m_selection_shapes_bounding_rect.center().y());

    // Show resize handles.
    foreach (SelectionResizeHandle *handle, m_selection_resize_handles)
    {
        handle->setVisible(true);
    }
}


bool Canvas::singlePropUpdateID(QDomElement& node, const QString& prop,
            const QString ns)
{
    bool wasSuccessful = false;
    QString oldId = nodeAttribute(node, ns, prop);
    if (!oldId.isNull())
    {
        QString propertyName = (ns.isEmpty()) ? prop : ns + ":" + prop;
        if (m_paste_id_map.find(oldId) != m_paste_id_map.end())
        {
            // The object this ID refers to is in the selection.
            node.setAttribute(propertyName, m_paste_id_map[oldId]);
            wasSuccessful = true;
        }
        else
        {
            node.removeAttribute(prop);
        }
    }
    return wasSuccessful;
}



void Canvas::recursiveMapIDs(QDomNode start, const QString& ns, int pass)
{

    for (QDomNode curr = start; !curr.isNull(); curr = curr.nextSibling())
    {
        if (curr.isElement())
        {
            QDomElement element = curr.toElement();
            QString idVal = element.attribute(x_id);

            if (!idVal.isEmpty())
            {
                if (pass == PASTE_UPDATEOBJIDS)
                {
                    QString newId = QString::number(++m_max_string_id);
                    m_paste_id_map[idVal] = newId;
                    element.setAttribute(x_id, newId);
                }
                else if (pass == PASTE_SELECTSHAPES)
                {
                    CanvasItem *obj = this->getItemByID(idVal);
                    if (obj)
                    {
                        obj->setSelected(true);
                    }
                }
                else if (pass == PASTE_REMOVEBADDISTROS)
                {
                    if (m_paste_bad_constraint_ids.contains(idVal))
                    {
                        // This is a bad distribution or separation, so
                        // remove its dunnart type so that it is ignored.
                        element.removeAttribute(x_type);
                    }
                }
            }

            if (pass == PASTE_UPDATEIDPROPS)
            {
                // Update single properties that refer to IDs.
                singlePropUpdateID(element, x_srcID);
                singlePropUpdateID(element, x_dstID);
                singlePropUpdateID(element, x_constraintID);
                singlePropUpdateID(element, x_objOneID);
                singlePropUpdateID(element, x_objTwoID);
            }
            else if (pass == PASTE_FINDBADDISTROS)
            {
                if (nodeAttribute(element, ns, x_type) == x_constraint)
                {
                    QString relType = nodeAttribute(element, ns, x_relType);
                    if ((relType ==  x_distribution) || (relType == x_separation))
                    {
                        QString indicatorID = essentialProp<QString>(element, x_constraintID);
                        QString objID1 = element.attribute(x_objOneID);
                        QString objID2 = element.attribute(x_objTwoID);

                        if (objID1.isEmpty() || objID2.isEmpty())
                        {
                            // This is a distribution or separation
                            // relationship for a guideline or guidelines
                            // that are not in the selection.

                            // Add distribution or separation to bad list.
                            m_paste_bad_constraint_ids.push_back(indicatorID);
                            // Remove the type attribute so that this node
                            // is effectively ignored in later processing.
                            element.removeAttribute(x_type);
                        }
                    }
                    if (relType == x_alignment)
                    {
                        QString objID1 = element.attribute(x_objOneID);

                        if (objID1.isEmpty())
                        {
                            // This is an alignment relationship for a
                            // shape that is not in the selection.

                            // Remove the type attribute so that this node
                            // is effectively ignored in later processing.
                            element.removeAttribute(x_type);
                        }
                    }
                }
            }
        }
        recursiveMapIDs(curr.firstChild(), ns, pass);
    }
}


// Return true for namespaces that are not used by Dunnart (since we
// will be copy elements/properties in such namespaces straight through).
static bool is_external_ns(const QString& ns)
{
    if (ns.isEmpty() || (ns == x_dunnartNs) || (ns == "inkscape") ||
            (ns ==  "sodipodi") )
    {
        return false;
    }
    else
    {
        return true;
    }
}

void Canvas::loadSVGRootNodeAttributes(const QDomElement& svgRoot)
{
    if (svgRoot.hasAttribute("viewBox"))
    {
        QString value = svgRoot.attribute("viewBox");
        if (!value.isEmpty())
        {
            double pageX, pageY, pageW, pageH;
            sscanf(value.toLatin1().data(), "%lf %lf %lf %lf",
                    &pageX, &pageY, &pageW, &pageH);
            this->setPageRect(QRectF(pageX, pageY, pageW, pageH));
        }
    }
}

void Canvas::recursiveReadSVG(const QDomNode& start, const QString& dunnartNS,
        int pass)
{
    if (pass == PASS_CLUSTERS)
    {
        // Cause shapes to be added before clusters try and reference them.
        router()->processTransaction();
    }

    for (QDomNode curr = start; !curr.isNull(); curr = curr.nextSibling())
    {
        if (!curr.prefix().isEmpty())
        {
            if (is_external_ns(curr.prefix()))
            {
                m_extra_namespaces_map[curr.prefix()] = curr.namespaceURI();
            }
        }

        if (curr.isElement())
        {
            const QDomElement element = curr.toElement();

            if (pass == PASS_SHAPES)
            {
                if ((element.localName() == "options") &&
                    (element.prefix() == x_dunnartNs))
                {
                    this->loadLayoutOptionsFromDomElement(element);
                }
                else if ((element.localName() == "identification") &&
                        (element.prefix() == "proorigami"))
                {
                    // For Pro-origami diagrams, use orthogonal connectors.
                    this->m_force_orthogonal_connectors = true;
                    // Don't allow the user to change diagram structure.
                    this->setOptStructuralEditingDisabled(true);
                    // Prevent overlaps.
                    this->m_opt_prevent_overlaps = true;
                }
                else if ((element.localName() == "svg") &&
                         element.prefix().isEmpty())
                {
                    this->loadSVGRootNodeAttributes(element);
                }
                else if (is_external_ns(element.prefix()))
                {
                    // Save nodes for external namespaces to output
                    // unchanged on saving.
                    // [ADS] FIXME: the tree structure of these external
                    // nodes will be lost, we are just storing them in
                    // a list regardless of depth.
                    QDomNode nodecopy = curr.cloneNode();
                    m_external_node_list.push_back(nodecopy);
                }
            }

            // Read other entities.
            if (nodeHasAttribute(element, dunnartNS, x_type))
            {
                // We have found a non-Dunnart node with a "dunnart:type"
                // attribute, thus we look for other attributes on this node
                // that are in the Dunnart namespace.
                CanvasItem::create(this, element, dunnartNS, pass);
            }
            if ((element.localName() == "node") &&
                (element.prefix() == x_dunnartNs))
            {
                // We have found a standard dunnart:node node.  We can read
                // attributes from this without any namespace.
                CanvasItem::create(this, element, "", pass);
            }
        }
        this->recursiveReadSVG(curr.firstChild(), dunnartNS, pass);
    }
}


bool Canvas::forceOrthogonalConnectors(void) const
{
    return m_force_orthogonal_connectors;
}

// Atributes:
static const char *x_EXPERIMENTAL_rect = "EXP-rectConstraint";
static const char *x_layoutMethod = "layoutMethod";
static const char *x_layoutMode = "layoutMode";
static const char *x_layoutBeautification = "layoutBeautify";
static const char *x_preventOverlaps = "preventOverlaps";
static const char *x_automaticGraphLayout =
        "automaticGraphLayout";
static const char *x_avoidBuffer = "avoidBuffer";
static const char *x_routingBuffer = "routingBuffer";
static const char *x_flowSeparation = "flowSeparation";
static const char *x_flowDirection = "flowDirection";
static const char *x_layeredAlignment = "layeredAlignment";
static const char *x_defaultIdealConnectorLength =
        "defaultIdealConnectorLength";
static const char *x_pageBoundaryConstraints =
        "pageBoundaryConstraints";
static const char *x_penaliseCrossings = "penaliseCrossings";
static const char *x_segmentPenalty = "segmentPenalty";
static const char *x_colourInterferingConnectors =
        "colourInterferingConnectors";
static const char *x_rubberBandRouting =
        "rubberBandRouting";
static const char *x_interferingConnectorColours =
        "interferingConnectorColours";

void Canvas::saveDiagram(const QString& outputFilename)
{
    QFileInfo fileInfo(outputFilename);
    PluginFileIOFactory *fileIOFactory = sharedPluginFileIOFactory();
    QString errorMessage;
    bool successful = fileIOFactory->saveDiagramToFile(this, fileInfo,
            errorMessage);

    if (successful)
    {
        undoStack()->setClean();
    }
    else
    {
        // We weren't successful saving, so show an error message.
        QString warning = QString(
                QObject::tr("<p><b>The document \"%1\" could not be saved.</b></p>"
                "<p>%2</p>")).arg(fileInfo.fileName()).arg(errorMessage);

        QWidget *window = views().first()->window();
        QMessageBox message(QMessageBox::Warning, "Error Saving File",
                            warning, QMessageBox::Ok, window);
        message.setWindowModality(Qt::WindowModal);
        message.exec();
    }
}


void Canvas::loadLayoutOptionsFromDomElement(const QDomElement& options)
{
    GraphLayout *gl = this->layout();
    assert(gl);

    int method = gl->optimizationMethod;
    if(optionalProp(options,x_layoutMethod,method)) {
        gl->setOptimizationMethod((GraphLayout::OptimizationMethod)method);
    }
    int mode = gl->mode;
    if (optionalProp(options,x_layoutMode,mode))
    {
        setOptLayoutMode((Canvas::LayoutMode) mode);
    }

    bool booleanVal = false;
    if (optionalProp(options,x_automaticGraphLayout,booleanVal))
    {
        setOptAutomaticGraphLayout(booleanVal);
    }

    if (optionalProp(options,x_layoutBeautification,booleanVal))
    {
        setOptPreserveTopology(booleanVal);
    }
    gl->runLevel = optPreserveTopology();

    if (optionalProp(options,x_preventOverlaps,booleanVal))
    {
        setOptPreventOverlaps(booleanVal);
    }

    if (optionalProp(options,x_pageBoundaryConstraints,booleanVal))
    {
        setOptFitWithinPage(booleanVal);
    }

    if (optionalProp(options,x_rubberBandRouting,booleanVal))
    {
        setOptRubberBandRouting(booleanVal);
    }

    optionalProp(options,x_EXPERIMENTAL_rect,m_rectangle_constraint_test);
    optionalProp(options,x_avoidBuffer,m_opt_shape_nonoverlap_padding);

    double routingBuffer;
    if (optionalProp(options,x_routingBuffer,routingBuffer))
    {
        setOptRoutingShapePadding(routingBuffer);
    }
    else
    {
        setOptRoutingShapePadding(m_opt_shape_nonoverlap_padding);
    }
    optionalProp(options,x_flowSeparation,m_flow_separation_modifier);
    optionalProp(options,x_flowDirection,m_opt_flow_direction);
    optionalProp(options,x_layeredAlignment,m_opt_layered_alignment_position);

    double ideal_connector_length_modifier;
    if (optionalProp(options,x_defaultIdealConnectorLength,
           ideal_connector_length_modifier))
    {
        setOptIdealEdgeLengthModifier(ideal_connector_length_modifier);
    }

    optionalProp(options,x_penaliseCrossings,m_avoid_connector_crossings);
    double segment_penalty;
    if (optionalProp(options,x_segmentPenalty, segment_penalty))
    {
        router()->setRoutingParameter(Avoid::segmentPenalty, segment_penalty);
    }
    optionalProp(options,x_colourInterferingConnectors,
            m_opt_colour_interfering_connectors);

    if (options.hasAttribute(x_interferingConnectorColours))
    {
        setInterferingConnectorColours(
                options.attribute(x_interferingConnectorColours));
    }

    unsigned int fileFontSize;
    if (optionalProp(options,x_fontSize,fileFontSize) && fileFontSize>0)
    {
        m_canvas_font_size = fileFontSize;
    }
}



QDomElement Canvas::writeLayoutOptionsToDomElement(QDomDocument& doc) const
{
    // Store Dunnart options.
    QDomElement dunOpts = doc.createElement("dunnart:options");

    // layout properties
    GraphLayout* gl = this->layout();
    newProp(dunOpts, x_automaticGraphLayout, optAutomaticGraphLayout());
    newProp(dunOpts, x_layoutMethod, (int)gl->optimizationMethod);
    newProp(dunOpts, x_layoutMode, (int)gl->mode);
    newProp(dunOpts, x_layoutBeautification, optPreserveTopology());
    newProp(dunOpts, x_preventOverlaps, optPreventOverlaps());
    newProp(dunOpts, x_avoidBuffer, m_opt_shape_nonoverlap_padding);
    newProp(dunOpts, x_routingBuffer, optRoutingShapePadding());
    newProp(dunOpts, x_flowSeparation, m_flow_separation_modifier);
    if (m_opt_flow_direction != FlowDown)
    {
        newProp(dunOpts, x_flowDirection, m_opt_flow_direction);
    }
    if (m_opt_layered_alignment_position != ShapeMiddle)
    {
        newProp(dunOpts, x_layeredAlignment, m_opt_layered_alignment_position);
    }
    newProp(dunOpts, x_pageBoundaryConstraints, optFitWithinPage());
    newProp(dunOpts, x_defaultIdealConnectorLength,
            optIdealEdgeLengthModifier());
    newProp(dunOpts, x_penaliseCrossings, m_avoid_connector_crossings);
    newProp(dunOpts, x_segmentPenalty,
            router()->routingParameter(Avoid::segmentPenalty));
    newProp(dunOpts, x_colourInterferingConnectors,
            m_opt_colour_interfering_connectors);
    newProp(dunOpts, x_rubberBandRouting, optRubberBandRouting());
    if (!m_interfering_connector_colours.isEmpty() &&
            (m_interfering_connector_colours != m_default_connector_colours))
    {
        QStringList colourStrings;
        QColor colour;
        foreach (colour, m_interfering_connector_colours)
        {
            // name is #RRGGBB, but just want RRGGBB
            colourStrings.append(colour.name().right(6));
        }
        // Build string of colours separated by commmas.
        newProp(dunOpts, x_interferingConnectorColours,
                colourStrings.join(",").toLatin1().data());
    }

    if (m_canvas_font_size != DEFAULT_CANVAS_FONT_SIZE)
    {
        newProp(dunOpts,x_fontSize,m_canvas_font_size);
    }

    return dunOpts;
}


//
// setInterferingConnectorColors() - set list of colors to color interfering
//                                   connectors
//
// Parameters:
//  color_list_str - comma-delmited list of hex RRGGBB color values (no '#' prefix)
//
void Canvas::setInterferingConnectorColours(const QString colourListString)
{
    QStringList colourList = colourListString.split(QChar(','));

    m_interfering_connector_colours.clear();

    QString colourStr;
    foreach (colourStr, colourList)
    {
        // Put in the form "#RRGGBB"
        colourStr.prepend("#");
        QColor colour(colourStr);
        if (colour.isValid())
        {
            m_interfering_connector_colours.append(colour);
        }
    }
}

//
// interferingConnectorColors() - get list of colors to color intefering
//                                   connectors
//
// Returns a list of QColors previously set up with the
// setInterferiongConenctorColours method, or if unset then
// the default list of colors.
//
const QList<QColor> Canvas::interferingConnectorColours(void) const
{
    if (m_interfering_connector_colours.isEmpty())
    {
        return m_default_connector_colours;
    }
    else
    {
        return m_interfering_connector_colours;
    }
}

QRectF diagramBoundingRect(const QList<CanvasItem *>& list)
{
    QRectF rect;

    foreach (CanvasItem *item, list)
    {
        if (!dynamic_cast<Indicator *> (item))
        {
            rect |= item->sceneBoundingRect();
        }
    }

    return rect;
}

void Canvas::improveOrthogonalTopology()
{
    router()->processActions();

    using namespace Avoid;
    // Create a router and set its parameters.
    Router *router2 = new Router(OrthogonalRouting);
    router2->setRoutingParameter((RoutingParameter)0, 50);
    router2->setRoutingParameter((RoutingParameter)1, 0);
    router2->setRoutingParameter((RoutingParameter)2, 0);
    router2->setRoutingParameter((RoutingParameter)3, 4000);
    router2->setRoutingParameter((RoutingParameter)4, 0);
    router2->setRoutingParameter((RoutingParameter)5, 100);
    router2->setRoutingParameter((RoutingParameter)6, 0);
    router2->setRoutingParameter((RoutingParameter)7, 4);
    router2->setRoutingOption((RoutingOption)0, true);
    router2->setRoutingOption((RoutingOption)1, true);
    router2->setRoutingOption((RoutingOption)2, false);
    router2->setRoutingOption((RoutingOption)3, false);

    // Keep track of router references to canvas shapes and connectors.
    cola::VariableIDMap idMap;
    QMap<ShapeObj*,ShapeRef*> d2rShapes;
    QMap<Connector*,ConnRef*> d2rConns;

    // Create router references.
    foreach (CanvasItem *item, items())
    {
        if (ShapeObj *shape = isShapeForLayout(item))
        {
            Polygon poly = shape->polygon();
            int id = shape->internalId();
            idMap.addMappingForVariable(id,id);
            ShapeRef *sr = new ShapeRef(router2, poly, id);
            d2rShapes.insert(shape,sr);
            //qDebug() << shape->internalId();
        }
        else if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            conn->applyNewRoute(conn->avoidRef->displayRoute(), true);

            int id = conn->internalId();
            idMap.addMappingForVariable(id,id);
            ConnRef *cref = new ConnRef(router2, id);
            d2rConns.insert(conn,cref);
            //qDebug() << conn->internalId();
            QPair<CPoint, CPoint> connpts = conn->get_connpts();
            ConnEnd srcPt(Point(connpts.first.x, connpts.first.y),ConnDirAll);
            ConnEnd dstPt(Point(connpts.second.x, connpts.second.y),ConnDirAll);
            cref->setEndpoints(srcPt, dstPt);
            cref->set_route(conn->avoidRef->displayRoute());
        }
    }

    // Build a GraphData object to get all the rectangles and constraints on
    // the canvas, plus the cluster hierarchy.
    GraphData *graph = new GraphData(this, true, m_graphlayout->mode, false, 10000);

    // Create the topology addon.
    topology::AvoidTopologyAddon topologyAddon(
                graph->rs, graph->ccs, &(graph->clusterHierarchy), idMap, DBL_MAX);

    //router2->outputInstanceToSVG("test-jog2-before");

    // Do the routing
    router2->processActions();
    foreach (Connector *conn, d2rConns.keys())
    {
        ConnRef *cr = d2rConns.value(conn);
        cr->set_route(conn->avoidRef->displayRoute().simplify());
    }
    router2->setTopologyAddon(&topologyAddon);
    router2->improveOrthogonalTopology();

    //router2->processTransaction();

    //router2->outputInstanceToSVG("test-jog2-after");

    // Apply the changes to the canvas.
    foreach (ShapeObj *shape, d2rShapes.keys())
    {
        ShapeRef *sr = d2rShapes.value(shape);
        Point pt = sr->position();
        shape->setCentrePos(QPointF(pt.x,pt.y));
    }
    foreach (Connector *conn, d2rConns.keys())
    {
        ConnRef *cr = d2rConns.value(conn);
        conn->applyNewRoute(cr->displayRoute(),true);
    }

    // Get the GraphLayout object to recognize the changes.
    stop_graph_layout();

    // Clean up.
    delete graph;
    delete router2;
}

/* Build the set C of ShapeObjs which (a) are connected to sh (according to nbrs), and
   (b) are in the passed set R.
  */
void Canvas::getRestrConnComp(QMap<ShapeObj *, ShapeObj *> nbrs, ShapeObj *sh,
                              QSet<ShapeObj *> &R, QSet<ShapeObj *> &C)
{
    // Add sh to component.
    C.insert(sh);
    // Compute set of eligible neighbours.
    QSet<ShapeObj*> N = nbrs.values(sh).toSet().intersect(R);
    // Recurse on those neighbours which are not already in the set C.
    foreach (ShapeObj *n, N) {
        if (!C.contains(n)) {
            getRestrConnComp(nbrs,n,R,C);
        }
    }
}

struct AlignDesc {
    AlignDesc(atypes at, CanvasItemList cil) {
        alignType = at;
        items = cil;
    }
    AlignDesc(Connector *conn, Canvas::Dimension d) {
        connector = conn;
        dim = d;
    }
    atypes alignType;
    CanvasItemList items;
    Connector *connector;
    Canvas::Dimension dim;
    double goalDelta;
};

QList<AlignDesc> Canvas::inferAlignOneDim(QList<ShapeObj *> shapes, QMap<ShapeObj *, ShapeObj *> nbrs,
                                          Dimension dim, qreal tolerance)
{
    // Sort by x-coord for vertical alignment, y-coord for horizontal.
    QMap<qreal,ShapeObj*> coordmap;
    foreach(ShapeObj *sh, shapes)
    {
        qreal coord = (dim == VERT ? sh->centrePos().x() : sh->centrePos().y());
        coordmap.insertMulti(coord,sh);
    }
    // Partition into equivalence classes by coord
    qreal eps = tolerance; // tolerance
    QList< QList<ShapeObj*> > classes;
    QList<qreal> keys = coordmap.uniqueKeys();
    int N = keys.length();
    QList<ShapeObj*> curr;
    qreal k = keys.at(0);
    curr.append(coordmap.values(k));
    qreal lastKey = k;
    for (int i = 1; i < N; i++)
    {
        qreal k = keys.at(i);
        if (k - lastKey >= eps) {
            // The current key is not within tolerance of the previous one.
            // So record the current class, and then start a new one.
            classes.append(curr);
            curr.clear();
        }
        curr.append(coordmap.values(k));
        if (i == N - 1) {
            classes.append(curr);
        }
        lastKey = k;
    }
    // Prepare return value.
    QList<AlignDesc> aligns;
    // Now partition each list into connected components.
    foreach (QList<ShapeObj*> list, classes)
    {
        QList< QSet<ShapeObj*> > ccs;
        QSet<ShapeObj*> S = list.toSet();
        QSet<ShapeObj*> R = list.toSet();
        while (!S.isEmpty()) {
            ShapeObj *sh = S.toList().first();
            QSet<ShapeObj*> C;
            getRestrConnComp(nbrs, sh, R, C);
            ccs.append(C);
            S = S.subtract(C);
        }
        // Create an alignment constraint for CCs of 2 or more elements.
        foreach (QSet<ShapeObj*> cc, ccs)
        {
            if (cc.size() < 2) { continue; }
            CanvasItemList items;
            foreach (ShapeObj *sh, cc) {
                items.append(sh);
            }
            atypes alignType = (dim == VERT ? ALIGN_CENTER : ALIGN_MIDDLE);
            //createAlignment(alignType, items);
            aligns.append(AlignDesc(alignType,items));
        }
    }
    return aligns;
}

QList<AlignDesc> Canvas::inferAlignments(qreal tolerance)
{
    // Build list of shapes, and map from each shape to list of all of its neighbours.
    QList<ShapeObj*> shapes;
    QMap<ShapeObj*,ShapeObj*> nbrs;
    foreach (CanvasItem *item, items())
    {
        if (ShapeObj *shape = dynamic_cast<ShapeObj*>(item))
        {
            shapes.append(shape);
        }
        else if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            QPair<ShapeObj*, ShapeObj*> endpts = conn->getAttachedShapes();
            nbrs.insertMulti(endpts.first,endpts.second);
            nbrs.insertMulti(endpts.second,endpts.first);
        }
    }

    QList<AlignDesc> aligns;

    // Vertical alignments
    aligns.append(inferAlignOneDim(shapes, nbrs, VERT, tolerance));
    // Horizontal alignments
    aligns.append(inferAlignOneDim(shapes, nbrs, HORIZ, tolerance));

    return aligns;
}

void Canvas::inferAndApplyAlignments()
{
    qreal tolerance = m_opt_snap_distance_modifier;
//#define reach
#ifdef reach
    tolerance = m_opt_snap_grid_width * 1.2;
#endif
    QList<AlignDesc> aligns = inferAlignments(tolerance);
    stop_graph_layout();
    foreach (AlignDesc d, aligns) {
        Guideline *gdln = createAlignment(d.alignType, d.items);
        if (m_infer_tentative_alignments) {
            gdln->setTentative(true);
        }
    }
    restart_graph_layout();
}

void Canvas::computeNeighbourhoods()
{
    m_align_nbrs.clear();
    m_shapes_by_id.clear();
    QList<Connector*> conns;
    foreach (CanvasItem *item, items())
    {
        if (ShapeObj *shape = dynamic_cast<ShapeObj*>(item))
        {
            uint id = shape->internalId();
            m_shapes_by_id.insert(id,shape);
        }
        else if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            conns.append(conn);
            QPair<ShapeObj*, ShapeObj*> endpts = conn->getAttachedShapes();
            m_align_nbrs.insertMulti(endpts.first,endpts.second);
            m_align_nbrs.insertMulti(endpts.second,endpts.first);
        }
    }
}

void Canvas::initTryAlignments()
{
    // Determine the range of IDs of shapes.
    // Also build neighbour sets.
    uint maxID = 0;
    m_align_nbrs.clear();
    m_shapes_by_id.clear();
    setOptRelax(false);
    QList<Connector*> conns;
    QList<Guideline*> guidelines;
    foreach (CanvasItem *item, items())
    {
        if (ShapeObj *shape = dynamic_cast<ShapeObj*>(item))
        {
            uint id = shape->internalId();
            maxID = qMax(id, maxID);
            m_shapes_by_id.insert(id,shape);
#if 0
            // Debugging, change shape labels to be IDs.
            QString label;
            label.setNum(id);
            shape->setLabel(label);
#endif
        }
        else if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            conns.append(conn);
            QPair<ShapeObj*, ShapeObj*> endpts = conn->getAttachedShapes();
            m_align_nbrs.insertMulti(endpts.first,endpts.second);
            m_align_nbrs.insertMulti(endpts.second,endpts.first);
        }
        else if (Guideline *g = dynamic_cast<Guideline*>(item))
        {
            guidelines.append(g);
        }
    }
    maxID++; // increment, so IDs themselves can be used as indices into array
    m_max_shape_id = maxID;
    uint size = maxID*maxID;

    // Array that just says whether we've tried aligning a pair of nodes or not:
    if (m_align_pairs_tried) {
        delete[] m_align_pairs_tried;
    }
    m_align_pairs_tried = new bool[size];
    for (uint i = 0; i < size; i++) {
        m_align_pairs_tried[i] = false;
    }

    // Array that indicates the "alignment state" of each pair of nodes, using
    // Canvas::AlignmentFlags

    m_alignment_state = Matrix2d<int>(maxID,maxID);
    for (uint i = 0; i < maxID; i++) {
        for (uint j = 0; j < maxID; j++) {
            m_alignment_state(i,j) = 0;
        }
    }
    foreach (Connector *conn, conns) {
        ShapeObj *s = conn->getAttachedShapes().first, *t = conn->getAttachedShapes().second;
        int sID = s->internalId(), tID = t->internalId();
        m_alignment_state(sID,tID) = Connected;
        m_alignment_state(tID,sID) = Connected;
    }
    foreach (Guideline *gl, guidelines) {
        AlignmentFlags af = gl->get_dir()==GUIDE_TYPE_VERT ? Vertical : Horizontal;
        CanvasItemSet objSet;
        gl->addAttachedShapesToSet(objSet);
        QList<CanvasItem*> items = objSet.toList();
        int n = items.size();
        for (int i = 0; i+1 < n; i++) {
            ShapeObj *s = dynamic_cast<ShapeObj*>(items.at(i));
            if (!s) continue;
            int sID = s->internalId();
            for (int j = i+1; j<n; j++) {
                ShapeObj *t = dynamic_cast<ShapeObj*>(items.at(j));
                if (!t) continue;
                int tID = t->internalId();
                m_alignment_state(sID,tID) |= af;
                m_alignment_state(tID,sID) |= af;
            }
        }
    }
    m_num_align_tries = 0;
    m_max_actual_dG = -DBL_MAX;
#ifdef newApplyAlignments
        applyAlignments();
#else
        tryAlignments();
#endif
}

/*
* Update alignment state matrix for the alignment of shapes s and t in the
* dimension indicated by the alignment flag.
*/
void Canvas::updateAlignmentStates(ShapeObj *s, ShapeObj *t, AlignmentFlags a)
{
    // Get sets of indices of shapes already aligned with each of s and t.
    int sID = s->internalId(), tID = t->internalId();
    QList<int> sA, tA;
    for (int j = 0; j < m_max_shape_id; j++) {
        if (m_alignment_state(sID, j) & a) sA.append(j);
        if (m_alignment_state(tID, j) & a) tA.append(j);
    }
    // s and t are aligned deliberately
    m_alignment_state(sID,tID) |= a | Deliberate;
    m_alignment_state(tID,sID) |= a | Deliberate;
    // s is now aligned with each element of tA, and vice versa
    foreach(int k, tA) {
        m_alignment_state(sID,k) |= a;
        m_alignment_state(k,sID) |= a;
    }
    foreach(int k, sA) {
        m_alignment_state(tID,k) |= a;
        m_alignment_state(k,tID) |= a;
    }
}

bool Canvas::predictCoincidence(Connector *conn, Dimension dim)
{
    ShapeObj *s = conn->getAttachedShapes().first, *t = conn->getAttachedShapes().second;
    double sz = dim==HORIZ ? s->centrePos().x() : s->centrePos().y();
    double tz = dim==HORIZ ? t->centrePos().x() : t->centrePos().y();
    double lowCoord = sz<tz ? sz : tz;
    double highCoord = sz<tz ? tz : sz;
    int lowID = sz<tz ? s->internalId() : t->internalId();
    int highID = sz<tz ? t->internalId() : s->internalId();
    // Let L and H be the low and high shapes respectively.
    // We consider each node U which is already aligned with either L or H.
    // Any such node must have lower coord than L if it is connected to L, and
    // higher coord than H if it is connected to H. If either of those conditions
    // fails, then we predict coincidence.
    bool coincidence = false;
    AlignmentFlags a = dim==HORIZ ? Horizontal : Vertical;
    for (int j = 0; j < m_max_shape_id; j++) {
        if (j==lowID || j==highID) continue;
        if (!m_shapes_by_id.contains(j)) continue;
        ShapeObj *r = m_shapes_by_id.value(j);
        int lj = m_alignment_state(lowID, j);
        int hj = m_alignment_state(highID, j);
        if (lj&a || hj&a) {
            double z = dim==HORIZ ? r->centrePos().x() : r->centrePos().y();
            // low shape
            if ( lj&Connected && lowCoord < z ) {
                coincidence = true; break;
            }
            // high shape
            if ( hj&Connected && z < highCoord ) {
                coincidence = true; break;
            }
        }
    }
    return coincidence;
}

/// Say whether named side of shape is "clear", meaning that there is not
/// currently a neighbour of this shape which is on that side and roughly
/// aligned, up to the named tolerance.
///
/// If a shape object is passed as the optional fourth argument (default
/// value is NULL), then this shape is considered an exception; i.e. if this
/// turns out to be the only neighbour on the named side, true is still returned.
///
/// Should make 'side' into an enum. For now:
/// 0 = right, 1 = top, 2 = left, 3 = bottom.
bool Canvas::sideIsClear(ShapeObj *s, int side, double tolerance, ShapeObj* except) {
    bool clear = true;
    QList<ShapeObj*> nbrs = m_align_nbrs.values(s);
    double sx = s->centrePos().x(), sy = s->centrePos().y();
    foreach (ShapeObj *nbr, nbrs) {
        if (nbr == except) continue;
        double nx = nbr->centrePos().x(), ny = nbr->centrePos().y();
        if (side==0) {
            // right side
            if (nx > sx && fabs(ny-sy) < tolerance) { clear = false; break; }
        } else if (side==1) {
            // top side
            if (ny < sy && fabs(nx-sx) < tolerance) { clear = false; break; }
        } else if (side==2) {
            // left side
            if (nx < sx && fabs(ny-sy) < tolerance) { clear = false; break; }
        } else { // side == 3
            // bottom side
            if (ny > sy && fabs(nx-sx) < tolerance) { clear = false; break; }
        }
    }
    return clear;
}

void Canvas::applyAlignments()
{
    m_trying_alignments = false;
    //m_previous_ortho_goal_value = computeOrthoObjective();
    // Build list of potential alignments.
    QList<AlignDesc*> ads;
    foreach (CanvasItem *item, items())
    {
        if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            ShapeObj *s = conn->getAttachedShapes().first;
            ShapeObj *t = conn->getAttachedShapes().second;
            int sID = s->internalId(), tID = t->internalId();
            // If not already aligned, then try both alignments.
            if ( !( m_alignment_state(sID,tID) & (Horizontal|Vertical) ) ) {
                ads.append( new AlignDesc(conn,HORIZ) );
                ads.append( new AlignDesc(conn,VERT) );
            }
        }
    }
    // Compute predicted goal function changes.
    predictOrthoObjectiveChange(ads);

#define chooseMostNegativeChange
#ifdef  chooseMostNegativeChange
    // Find most negative change.
    double dG = DBL_MAX; AlignDesc *a = NULL;
    foreach (AlignDesc *b, ads) {
        if (b->goalDelta < dG) {
            dG = b->goalDelta;
            a = b;
        }
    }
    //bool proceed = dG<0;
    bool proceed = dG<1000000;
#else
    // Find largest estimated increase in goal function.
    double dG = DBL_MIN; AlignDesc *a = NULL;
    foreach (AlignDesc *b, ads) {
        if (b->goalDelta > dG) {
            dG = b->goalDelta;
            a = b;
        }
    }
    bool proceed = dG>-1000000;
#endif

    // Try another alignment only if it looks like one will
    // result in a decrease of the goal function.
    if (proceed) {
        // Get the shapes involved.
        ShapeObj *s = a->connector->getAttachedShapes().first;
        ShapeObj *t = a->connector->getAttachedShapes().second;
        //qDebug() << "%%%aligning:" << s->internalId() << t->internalId() << a->dim << s->centrePos().x() << s->centrePos().y() << t->centrePos().x() << t->centrePos().y();
        // Update various records.
        m_trying_alignments = true; // (enables callback trigger in Canvas::processLayoutFinishedEvent)
        m_num_align_tries++;
        AlignmentFlags af = a->dim==VERT ? Vertical : Horizontal;
        updateAlignmentStates(s,t,af);
        int sID = s->internalId(), tID = t->internalId();
        m_align_pairs_tried[sID*m_max_shape_id+tID] = true;
        m_align_pairs_tried[tID*m_max_shape_id+sID] = true;
        // Create new guideline and restart graph layout.
        CanvasItemList items;
        items.append(s); items.append(t);
#define alliedSeps
#ifdef  alliedSeps
        // First create a separation constraint to preserve their orthogonal ordering.
        dtype dt = a->dim==VERT ? SEP_VERTICAL : SEP_HORIZONTAL;
        double minDist = (s->width()+t->width())/2.0;
        bool sort = true;
        createSeparation(NULL,dt,items,minDist,sort);
#endif
        // Now create alignment.
        atypes at = a->dim==VERT ? ALIGN_CENTER : ALIGN_MIDDLE;
        Guideline *gdln = createAlignment(at,items);
        gdln->setTentative(true);
        long timestamp = ++m_tentative_guideline_timestamp;
        //qDebug() << "created guideline number" << timestamp;
        gdln->setTimestamp(timestamp);
        //qDebug() << "Trying alignment...";
        fully_restart_graph_layout();
    } else {
        qDebug() << "Finished applying alignments.";
        //qDebug() << m_alignment_state.toString();
        TrialAlignmentEvent *tae = new TrialAlignmentEvent(1);
        QCoreApplication::postEvent(this, tae, Qt::LowEventPriority);
    }
}

void Canvas::appliedAlignmentWasUnsat(ConstraintRejectedEvent *cre)
{
    // TODO: cre now has m_shape field pointing to the specific
    // ShapeObj on the Guideline whose alignment constraint was rejected,
    // so we should try just removing this shape from the guideline.
    /* Refer to ShapeObj::buildAndExecContextMenu on this.
      Looks like we need to figure out which Relationship object in
      the ShapeObj's rels array corresponds to this Guideline, and then
      call that Relationship's deactivate() method.
      */
    ShapeObj *shape = cre->m_shape;
    Guideline *reject = cre->m_guideline;
    qDebug() << "Found UNSAT tentative constraint.";
    qDebug() << "Guideline:" << reject->idString() << "Shape:" << shape->idString();
    // Delete guideline.
    stop_graph_layout();
    UndoMacro *undoMacro = beginUndoMacro(tr("Delete"));
    CanvasItemSet dummySet;
    reject->deactivateAll(dummySet);
    QUndoCommand *cmd = new CmdCanvasSceneRemoveItem(this, reject);
    undoMacro->addCommand(cmd);
    fully_restart_graph_layout();
}

void Canvas::applyAlignmentsCallback()
{
    TrialAlignmentEvent *tae = new TrialAlignmentEvent(0);
    QCoreApplication::postEvent(this, tae, Qt::LowEventPriority);
    return;

    m_trying_alignments = false;
    double G  = m_previous_ortho_goal_value;
    double Gp = computeOrthoObjective();
    double dG = Gp - G;
    //qDebug() << "Assessing result of alignment. dG =" << dG;
    m_max_actual_dG = std::max(m_max_actual_dG, dG);
    //qDebug() << "Max actual dG:" << m_max_actual_dG;
    if (dG < -m_apply_alignments_epsilon) {
        // Will continue applying alignments.
        TrialAlignmentEvent *tae = new TrialAlignmentEvent(0);
        QCoreApplication::postEvent(this, tae, Qt::LowEventPriority);

    } else {
        qDebug() << "Finished applying alignments.";
        //qDebug() << m_alignment_state.toString();
        TrialAlignmentEvent *tae = new TrialAlignmentEvent(1);
        QCoreApplication::postEvent(this, tae, Qt::LowEventPriority);
    }
}

void Canvas::initRejectAlignments()
{
//#define doRejectionPhase
#ifdef  doRejectionPhase
    qDebug() << "Rejecting alignments...";
    setOptRelaxThresholdModifier(0);
    setOptRelax(true);
    m_reject_phase_previous_goal_value = 0;
    m_reject_phase_previous_stress_value = 0;
    m_reject_phase_previous_align = NULL;
    rejectAlignments();
#endif
}

void Canvas::rejectAlignments()
{
    computeOrthoObjective();
    m_rejecting_alignments = true;
    fully_restart_graph_layout();
}

void Canvas::rejectAlignmentsCallback(ConstraintRejectedEvent *cre)
{
    m_rejecting_alignments = false;
    bool proceed = true;
    double currentGoal = computeOrthoObjective();
    double currentStress = computeStress();
    // Was there a previous run?
    if (m_reject_phase_previous_align!=NULL) {
        // Check whether it went well.
        double threshold = 20;
        double dG = currentGoal - m_reject_phase_previous_goal_value;
        double dS = currentStress - m_reject_phase_previous_stress_value;
        qDebug() << "Rejection dG" << dG << "dS" << dS;
        if (dG > threshold) {
            // It did not go well enough.
            // Restore the last alignment, and quit.
            proceed = false;
            setOptRelax(false);
            AlignDesc *ad = m_reject_phase_previous_align;
            Guideline *gdln = createAlignment(ad->alignType, ad->items);
            gdln->setTentative(true);
            fully_restart_graph_layout();
        }
    }

    // Try again?
    if (proceed) {
        m_reject_phase_previous_goal_value = currentGoal;
        m_reject_phase_previous_stress_value = currentStress;
        Guideline *gdln = cre->m_guideline;
        if (gdln->canvas())
        {
            // Save a record of the alignment.
            CanvasItemSet itemSet;
            gdln->addAttachedShapesToSet(itemSet);
            atypes at = gdln->get_dir()==GUIDE_TYPE_VERT ? ALIGN_CENTER : ALIGN_MIDDLE;
            m_reject_phase_previous_align = new AlignDesc(at,itemSet.toList());
            // Delete guideline.
            stop_graph_layout();
            UndoMacro *undoMacro = beginUndoMacro(tr("Delete"));
            CanvasItemSet dummySet;
            gdln->deactivateAll(dummySet);
            QUndoCommand *cmd = new CmdCanvasSceneRemoveItem(this, gdln);
            undoMacro->addCommand(cmd);
        }
        // Run it again.
        TrialAlignmentEvent *tae = new TrialAlignmentEvent(2);
        QCoreApplication::postEvent(this, tae, Qt::LowEventPriority);
    }
}

void Canvas::tryAlignments()
{
    m_trying_alignments = false;
    //qDebug() << "Num align tries:" << m_num_align_tries;
    //if (m_num_align_tries >= m_max_align_tries) return;
    double score = computeOrthoObjective();
    //qDebug() << "Objective function:" << score;
    double eps = 10;
    double sig = m_opt_snap_distance_modifier;
    //qDebug() << "snap distance:" << sig;
    // For now, try simply aligning neighbours which are not already aligned.
    foreach (CanvasItem *item, items())
    {
        if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            ShapeObj *s = conn->getAttachedShapes().first;
            ShapeObj *t = conn->getAttachedShapes().second;

            // Have we tried this pair before?
            int sid = s->idString().toInt(), tid = t->idString().toInt();
            int i = sid <= tid ? sid : tid;
            int j = sid <= tid ? tid : sid;
            int pairIndex = i*m_max_shape_id + j;
            if (m_align_pairs_tried[pairIndex]) continue; // do not try a second time

            double sx=s->centrePos().x(), sy=s->centrePos().y();
            double tx=t->centrePos().x(), ty=t->centrePos().y();
            double ady=fabs(ty-sy), adx=fabs(tx-sx);

            //if (adx < eps || ady < eps) continue; // already aligned

            int plan = 0; // 0 = plan no alignment; 1 = horizontal; 2 = vertical

            // Would a horizontal alignment be suitable?
            if (ady <= sig) {
                // Check whether the appropriate sides of the shapes are open.
                ShapeObj *left  = sx < tx ? s : t;
                ShapeObj *right = sx < tx ? t : s;
                if ( sideIsClear(left,0,eps,right) && sideIsClear(right,2,eps,left) ) plan = 1;
            }
            // Would a vertical alignment be suitable?
            if (adx <= sig) {
                ShapeObj *above = sy < ty ? s : t;
                ShapeObj *below = sy < ty ? t : s;
                if ( sideIsClear(above,3,eps,below) && sideIsClear(below,1,eps,above) ) {
                    plan = plan==0 ? 2 : ( adx < ady ? 2 : 1 );
                }
            }

            bool doOverlapPrevention = true;
            if (!doOverlapPrevention) {
                plan = 0;
                if (adx <= sig || ady <= sig) {
                    plan = adx < ady ? 2 : 1;
                }
            }

            if ((i==35&&j==71) || (i==71&&j==35)) {
                qDebug() << "i,j" << i << j << "plan:" << plan;
            }

            if (plan > 0) {
                // Will try an alignment.
                CanvasItemList items;
                items.append(s); items.append(t);
                atypes a = plan==2 ? ALIGN_CENTER : ALIGN_MIDDLE;
                //qDebug() << "Trying alignment" << a << "adx=" << adx << "ady=" << ady << "s:" << s->idString() << "t:" << t->idString();
                Guideline *gdln = createAlignment(a,items);
                //qDebug() << "guideline is" << (gdln==NULL?"NULL":"not null");
                gdln->setTentative(true);
                AlignmentFlags af = plan==2 ? Vertical : Horizontal;
                updateAlignmentStates(s,t,af);
                m_trying_alignments = true;
                m_num_align_tries++;
                m_align_pairs_tried[pairIndex] = true; // Mark this pair of shapes as "tried".
                break;
            }
        }
    }
    if (m_trying_alignments) {
        fully_restart_graph_layout();
    } else {
        if (sig < 250) {
            m_opt_snap_distance_modifier *= 2; // experimental!
            qDebug() << "Doubling snap distance, and trying alignments again.";
            tryAlignments();
        }
    }
}

LineSegment::LineSegment(Connector *conn) : connector(conn)
{
    ShapeObj *s = conn->getAttachedShapes().first;
    ShapeObj *t = conn->getAttachedShapes().second;
    double sx=s->centrePos().x(), sy=s->centrePos().y();
    double tx=t->centrePos().x(), ty=t->centrePos().y();
    computeParameters(sx,sy,tx,ty);
}

LineSegment::LineSegment(double sx, double sy, double tx, double ty) : connector(NULL)
{
    computeParameters(sx,sy,tx,ty);
}

void LineSegment::computeParameters(double sx, double sy, double tx, double ty)
{
    p1 = Avoid::Point(sx, sy);
    p2 = Avoid::Point(tx, ty);

    // If the points are coincident, set angle to -1 and quit.
    if (sx==tx && sy==ty) { angle = -1; return; }
    double a;
    if (ty>sy) {
        a = atan2(ty-sy,tx-sx);
    } else if (ty<sy) {
        a = atan2(sy-ty,sx-tx);
    } else {
        a = 0;
    }
    angle = (int)(round(a*180/3.1415927)) % 180;
    if (angle==0) {
        intercept = (sy+ty)/2.0;
        t0 = sx < tx ? sx : tx;
        t1 = sx < tx ? tx : sx;
    } else if (angle==90) {
        intercept = (sx+tx)/2.0;
        t0 = sy < ty ? sy : ty;
        t1 = sy < ty ? ty : sy;
    } else {
        intercept = sx - sy*(tx-sx)/(ty-sy);
        a = intercept;
        double st = sqrt((sx-a)*(sx-a)+sy*sy) * (sy < 0 ? -1 : 1);
        double tt = sqrt((tx-a)*(tx-a)+ty*ty) * (ty < 0 ? -1 : 1);
        t0 = st < tt ? st : tt;
        t1 = st < tt ? tt : st;
    }
}

bool LineSegment::intersects(LineSegment *other, double tolerance)
{
    Q_UNUSED (tolerance)

    bool ans = Avoid::segmentIntersect(p1, p2, other->p1, other->p2);
    if (ans) {
        connector->addIntersector(other->connector);
        other->connector->addIntersector(connector);
    } else {
        connector->removeIntersector(other->connector);
        other->connector->removeIntersector(connector);
    }
    return ans;
}

bool LineSegment::coincidesWith(LineSegment *other, double angleTolerance, double interceptTolerance)
{
    bool ans = false;
    if ( fabs(angle-other->angle) < angleTolerance && fabs(intercept-other->intercept) < interceptTolerance ) {
        // In this case, we consider the /lines/ to be the same, so we must say
        // whether the line /segments/ overlap.
        double u0 = other->t0, u1 = other->t1;
        // The question is whether the intervals (t0,t1) and (u0,u1) intersect.
        // We should also ensure some fair bit of overlap.
        //return (t1>u0 && u1>t0);
        ans = u0 + 3*interceptTolerance < t1 && t0 + 3*interceptTolerance < u1;
    }
    if (ans) {
        connector->addCoincider(other->connector);
        other->connector->addCoincider(connector);
    } else {
        connector->removeCoincider(other->connector);
        other->connector->removeCoincider(connector);
    }
    return ans;
}

double LineSegment::obliquityScore()
{
    // Put angle in range from 1 to 89.
    int a = angle > 90 ? angle - 90 : angle;
    // If close enough to orthogonal, then zero obliquity.
    int t = 3;
    if (a<=t || a>=90-t) return 0;
    // Else base score on distance from 45 deg.
    double d = fabs(a-45); // 0 <= d <= 45-t
    double deg45_cost = 1;
    return deg45_cost + d;

}

void Canvas::updateStress(double stress) {
    m_most_recent_stress = stress;
    int stressBarValue = (int)(round(stress));
    emit newStressBarValue(stressBarValue);
}

/*
* Predicts what the change in the ortho objective function would be if you were
* to align the endpoints of the connector in the named dimension.
*
* Returns  DBL_MAX if this alignment is predicted to result in creation
* of an edge coincidence.
*
* Leaves out crossing detection -- just uses the current value.
*
*/
double Canvas::predictOrthoObjectiveChange(Connector *conn, Dimension dim)
{
    // Would it create a coincidence?
    if (predictCoincidence(conn,dim)) return DBL_MAX;

    // Obliquity change
    LineSegment seg(conn);
    double dOb = -seg.obliquityScore();

    // Stress change
    GraphData *graph = new GraphData(this, false, m_graphlayout->mode, false, 10000);
    ShapeObj *s1 = conn->getAttachedShapes().first;
    ShapeObj *s2 = conn->getAttachedShapes().second;
    unsigned id1 = graph->getNodeID(s1);
    unsigned id2 = graph->getNodeID(s2);
    QPointF p1 = s1->centrePos(), p2 = s2->centrePos();
    vpsc::Rectangle *r1 = graph->rs.at(id1);
    vpsc::Rectangle *r2 = graph->rs.at(id2);
    // Now translate rectangles according to proposed alignment.
    if (dim == HORIZ) {
        double yav = ( p1.y() + p2.y() ) / 2.0;
        r1->moveCentreY(yav);
        r2->moveCentreY(yav);
    } else { // dim == VERT
        double xav = ( p1.x() + p2.x() ) / 2.0;
        r1->moveCentreX(xav);
        r2->moveCentreX(xav);
    }
    double stress = m_graphlayout->computeStressOfGraph(graph);
    double dStress = stress-m_most_recent_stress;

    // Sum up
    OrthoWeights o = m_ortho_weights;
    double predictedDelta = o.wob*dOb + o.wst*dStress;
    return predictedDelta;
}

double Canvas::computeStress()
{
    GraphData *graph = new GraphData(this, false, m_graphlayout->mode, false, 10000);

    vpsc::Rectangles rs = graph->rs;
    std::vector<cola::Edge> es = graph->edges;
    double idealLength = 1.0;
    std::valarray<double> eLengths;
    graph->getEdgeLengths(eLengths);
    unsigned n = rs.size();

    std::valarray<double> X = std::valarray<double>(n);
    std::valarray<double> Y = std::valarray<double>(n);
    unsigned i=0;
    for(vpsc::Rectangles::const_iterator ri=rs.begin();ri!=rs.end();++ri,++i) {
        X[i]=(*ri)->getCentreX();
        Y[i]=(*ri)->getCentreY();
    }

    // Compute D and G matrices.
    double** D=new double*[n];
    unsigned short** G=new unsigned short*[n];
    for(unsigned i=0;i<n;i++) {
        D[i]=new double[n];
        G[i]=new unsigned short[n];
    }

    //shortest_paths::johnsons(n,D,es,&eLengths); // Getting weird crash on this one.
    shortest_paths::johnsons(n,D,es,&eLengths);

    for(unsigned i=0;i<n;i++) {
        for(unsigned j=0;j<n;j++) {
            if(i==j) continue;
            double& d=D[i][j];
            unsigned short& p=G[i][j];
            p=2;
            if(d==DBL_MAX) {
                // i and j are in disconnected subgraphs
                p=0;
            } else if(!&eLengths) {
                D[i][j]*=idealLength;
            }
        }
    }
    for(std::vector<cola::Edge>::const_iterator e=es.begin();e!=es.end();++e) {
        unsigned u=e->first, v=e->second;
        G[u][v]=G[v][u]=1;
    }

    // Compute stress
    double stress=0;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            if(p==0) continue;
            double rx=X[u]-X[v], ry=Y[u]-Y[v];
            double l=sqrt(rx*rx+ry*ry);
            double d=D[u][v];
            if(l>d && p>1) continue;
            double d2=d*d;
            double rl=d-l;
            double s=rl*rl/d2;
            stress+=s;
        }
    }
    // Note: we are missing any stress that would be added by a topologyAddon.

    return stress;
}

void Canvas::predictStressChange(QList<AlignDesc*> ads)
{
    GraphData *graph = new GraphData(this, false, m_graphlayout->mode, false, 10000);

    vpsc::Rectangles rs = graph->rs;
    std::vector<cola::Edge> es = graph->edges;
    double idealLength = 1.0;
    std::valarray<double> eLengths;
    graph->getEdgeLengths(eLengths);
    unsigned n = rs.size();

    std::valarray<double> X = std::valarray<double>(n);
    std::valarray<double> Y = std::valarray<double>(n);
    unsigned i=0;
    for(vpsc::Rectangles::const_iterator ri=rs.begin();ri!=rs.end();++ri,++i) {
        X[i]=(*ri)->getCentreX();
        Y[i]=(*ri)->getCentreY();
    }

    // Compute D and G matrices.
    double** D=new double*[n];
    unsigned short** G=new unsigned short*[n];
    for(unsigned i=0;i<n;i++) {
        D[i]=new double[n];
        G[i]=new unsigned short[n];
    }

    //shortest_paths::johnsons(n,D,es,&eLengths); // Getting weird crash on this one.
    shortest_paths::johnsons(n,D,es,&eLengths);

    for(unsigned i=0;i<n;i++) {
        for(unsigned j=0;j<n;j++) {
            if(i==j) continue;
            double& d=D[i][j];
            unsigned short& p=G[i][j];
            p=2;
            if(d==DBL_MAX) {
                // i and j are in disconnected subgraphs
                p=0;
            } else if(!&eLengths) {
                D[i][j]*=idealLength;
            }
        }
    }
    for(std::vector<cola::Edge>::const_iterator e=es.begin();e!=es.end();++e) {
        unsigned u=e->first, v=e->second;
        G[u][v]=G[v][u]=1;
    }

    // Compute stress
    double stress0=0;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            if(p==0) continue;
            double rx=X[u]-X[v], ry=Y[u]-Y[v];
            double l=sqrt(rx*rx+ry*ry);
            double d=D[u][v];
            if(l>d && p>1) continue;
            double d2=d*d;
            double rl=d-l;
            double s=rl*rl/d2;
            stress0+=s;
        }
    }
    // Note: we are missing any stress that would be added by a topologyAddon.

    double minDelta=DBL_MAX, maxDelta=DBL_MIN;
    foreach (AlignDesc *ad, ads) {
        // Now compute stress if the alignment were to be put in effect.
        Connector *conn = ad->connector;
        Dimension dim   = ad->dim;
        ShapeObj *s1 = conn->getAttachedShapes().first;
        ShapeObj *s2 = conn->getAttachedShapes().second;
        unsigned id1 = graph->getNodeID(s1);
        unsigned id2 = graph->getNodeID(s2);
        QPointF p1 = s1->centrePos(), p2 = s2->centrePos();
        // Translate coords according to proposed alignment.
        double c1=0,c2=0;
        if (dim == HORIZ) {
            double yav = ( p1.y() + p2.y() ) / 2.0;
            c1 = Y[id1]; c2 = Y[id2];
            Y[id1] = yav; Y[id2] = yav;
        } else { // dim == VERT
            double xav = ( p1.x() + p2.x() ) / 2.0;
            c1 = X[id1]; c2 = X[id2];
            X[id1] = xav; X[id2] = xav;
        }
        // Recompute stress.
        double astress = 0;
        for(unsigned u=0;(u + 1)<n;u++) {
            for(unsigned v=u+1;v<n;v++) {
                unsigned short p=G[u][v];
                if(p==0) continue;
                double rx=X[u]-X[v], ry=Y[u]-Y[v];
                double l=sqrt(rx*rx+ry*ry);
                double d=D[u][v];
                if(l>d && p>1) continue;
                double d2=d*d;
                double rl=d-l;
                double s=rl*rl/d2;
                astress+=s;
            }
        }
        // Restore coords.
        if (dim == HORIZ) {
            Y[id1] = c1; Y[id2] = c2;
        } else { // dim == VERT
            X[id1] = c1; X[id2] = c2;
        }

        double dstress = astress-stress0;

        //qDebug() << "dstress" << s1->internalId() << s2->internalId() << dstress;
        minDelta = dstress<minDelta?dstress:minDelta;
        maxDelta = dstress>maxDelta?dstress:maxDelta;

        OrthoWeights o = m_ortho_weights;
        ad->goalDelta = o.wst*dstress;
    }
    //qDebug() << "max and min estimated stress change:" << maxDelta << minDelta;

}

/*
* This version handles a whole list of Connector-Dimension pairs at once.
*/
void Canvas::predictOrthoObjectiveChange(QList<AlignDesc *> &ads)
{
    // First compute predicted change in stress.
    predictStressChange(ads);
    // Now consider other predictions.
    foreach (AlignDesc *ad, ads) {
        Connector *conn = ad->connector;
        Dimension dim   = ad->dim;
        ShapeObj *s1    = conn->getAttachedShapes().first;
        ShapeObj *s2    = conn->getAttachedShapes().second;

        // Would it create a coincidence?
        if (predictCoincidence(conn,dim)) {
            ad->goalDelta = DBL_MAX;
            continue;
        }

        // Obliquity change
        LineSegment seg(conn);
        double dOb = -seg.obliquityScore();

        // Incidence score
        int inc1 = m_align_nbrs.values(s1).size();
        int inc2 = m_align_nbrs.values(s2).size();
        double inc = inc1+inc2;

        // Stress change
        //double dStress = 0;
        /*
        double stress = computeStress();
        double dStress = stress-m_most_recent_stress;
        qDebug() << "computed stress" << stress;
        qDebug() << "difference:" << dStress;
        */
        //dStress = 0;

        // Sum up
        OrthoWeights o = m_ortho_weights;
#define considerObliquity
#ifdef  considerObliquity
        ad->goalDelta += o.wob*dOb;
#endif

//#define justIncidenceScore
#ifdef  justIncidenceScore
        ad->goalDelta = -inc;
#endif

#define bendPointPenalty
//#define SBGN
#ifdef  SBGN
        double score = 0;
        int npd1 = nonPendantDegree(conn->getAttachedShapes().first);
        int npd2 = nonPendantDegree(conn->getAttachedShapes().second);
        //if (npd1!=2 || npd2!=2) score = 100000;

        score = 1000;
        if (npd1>=2 && npd2>=2 && (npd1==2 || npd2==2)) score=0;
        if (npd1==1 || npd2==1) score=DBL_MAX;
        //if (npd1==1 && npd2==1) score=DBL_MAX;

        ad->goalDelta += score;

        // Would this alignment make either endpoint into a degree-2 bend point?
        if (npd1==2) {
            if (dim==HORIZ && numVAligns(s1)>0 || dim==VERT && numHAligns(s1)>0) ad->goalDelta += 1000;
        }
        if (npd2==2) {
            if (dim==HORIZ && numVAligns(s2)>0 || dim==VERT && numHAligns(s2)>0) ad->goalDelta += 1000;
        }
#else
        // Would this alignment make either endpoint into a degree-2 bend point?
        if (inc1==2) {
            if (dim==HORIZ && numVAligns(s1)>0 || dim==VERT && numHAligns(s1)>0) ad->goalDelta += 1000;
        }
        if (inc2==2) {
            if (dim==HORIZ && numVAligns(s2)>0 || dim==VERT && numHAligns(s2)>0) ad->goalDelta += 1000;
        }
#endif

    }
}

int Canvas::numHAligns(ShapeObj *s) {
    int sid = s->internalId();
    int n = 0;
    for (int j = 0; j < m_max_shape_id; j++) {
        if (m_alignment_state(sid,j) & Horizontal) n++;
    }
    return n;
}

int Canvas::numVAligns(ShapeObj *s) {
    int sid = s->internalId();
    int n = 0;
    for (int j = 0; j < m_max_shape_id; j++) {
        if (m_alignment_state(sid,j) & Vertical) n++;
    }
    return n;
}

/// How many shapes is shape s connected to which are not pendants?
int Canvas::nonPendantDegree(ShapeObj *s)
{
    QList<ShapeObj*> nbrs = m_align_nbrs.values(s);
    int d = 0;
    foreach (ShapeObj *t, nbrs) {
        int n = m_align_nbrs.values(t).size();
        if (n>1) d++;
    }
    return d;
}

/*
* Predicts what the ortho objective function would be if you were
* to align the endpoints of the connector in the named dimension.
*
* Returns -1 if this alignment is predicted to result in creation
* of an edge coincidence.
*
* Leaves out crossing detection -- just uses the current value.
*
* Only uses a heuristic to check if a coincidence will be created.
*
*/
double Canvas::predictOrthoObjective(Connector *conn, Dimension dim)
{
    // Shouldn't be called according to Steve.
    abort();

    // Obliquity
    QList<LineSegment*> segs;
    double obliquity = 0;
    foreach (CanvasItem *item, items())
    {
        if (Connector *c = dynamic_cast<Connector*>(item))
        {
            if (c==conn) continue; // (If conn is aligned, it will have 0 obliquity.)
            LineSegment *s = new LineSegment(c);
            segs.append(s);
            obliquity += s->obliquityScore();
        }
    }

    // Stress
    GraphData *graph = new GraphData(this, false, m_graphlayout->mode, false, 10000);
    ShapeObj *s1 = conn->getAttachedShapes().first;
    ShapeObj *s2 = conn->getAttachedShapes().second;
    unsigned id1 = graph->getNodeID(s1);
    unsigned id2 = graph->getNodeID(s2);
    QPointF p1 = s1->centrePos(), p2 = s2->centrePos();

    vpsc::Rectangle *r1 = graph->rs.at(id1);
    vpsc::Rectangle *r2 = graph->rs.at(id2);
    // Now translate rectangles according to proposed alignment.
    if (dim == HORIZ) {
        double yav = ( p1.y() + p2.y() ) / 2.0;
        r1->moveCentreY(yav);
        r2->moveCentreY(yav);
    } else { // dim == VERT
        double xav = ( p1.x() + p2.x() ) / 2.0;
        r1->moveCentreX(xav);
        r2->moveCentreX(xav);
    }

    double stress = m_graphlayout->computeStressOfGraph(graph);

    // Crossings
    int crossings = m_most_recent_crossing_count;

    // Coincidences
    int coincidences = m_most_recent_coincidence_count;

    // Sum up
    OrthoWeights o = m_ortho_weights;
    double prediction = o.wcr*crossings + o.wco*coincidences + o.wob*obliquity + o.wst*stress;
    return prediction;
}

class NeighbourAngleLessThan {
public:
    NeighbourAngleLessThan(ShapeObj *shape)
        : m_centre_shape_for_angles(shape)
    {

    }

    bool operator()(const ShapeObj *t, const ShapeObj *u)
    {
        ShapeObj *s = m_centre_shape_for_angles;
        double sx = s->centrePos().x(), sy = s->centrePos().y();
        double tx = t->centrePos().x(), ty = t->centrePos().y();
        double ux = u->centrePos().x(), uy = u->centrePos().y();
        double tAngle = atan2(ty-sy,tx-sx), uAngle = atan2(uy-sy,ux-sx);
        return tAngle < uAngle;
    }
private:
    ShapeObj *m_centre_shape_for_angles;
};


double Canvas::computeOrthoObjective()
{
    // Build a LineSegment for each Connector.
    // Sum obliquity now.
    QList<LineSegment*> segs;
    double obliquity = 0;
    foreach (CanvasItem *item, items())
    {
        if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            LineSegment *s = new LineSegment(conn);
            segs.append(s);
            obliquity += s->obliquityScore();
        }
    }
    emit newObliquityBarValue( (int)(round(obliquity)) );
    if (segs.size()>0) {
        m_most_recent_average_obliqueness = obliquity/segs.size();
    }

    // Compare each pair of edges to see whether they cross and/or are coincident.
    // For now we do the naive quadratic run-time approach. Improve later if needed.
    int m = segs.size();
    int crossings = 0, coincidences = 0;

    double angleTolerance = 3;
    double interceptTolerance = 3;

    QMap<LineSegment*,QSet<LineSegment*> > coincidenceGroups;

    for (int i = 0; i+1 < m; i++) {
        LineSegment *s1 = segs.at(i);
        for (int j = i+1; j < m; j++) {
            LineSegment *s2 = segs.at(j);
            if (s1->intersects(s2,interceptTolerance)) {
                crossings++;
            }
            if (s1->coincidesWith(s2,angleTolerance,interceptTolerance)) {
                coincidences++;
                if (!coincidenceGroups.contains(s1) && !coincidenceGroups.contains(s2)) {
                    QSet<LineSegment*> set;
                    set.insert(s1);
                    set.insert(s2);
                    coincidenceGroups.insert(s1,set);
                    coincidenceGroups.insert(s2,set);
                } else if (coincidenceGroups.contains(s1) && !coincidenceGroups.contains(s2)) {
                    QSet<LineSegment*> set = coincidenceGroups.value(s1);
                    set.insert(s2);
                    coincidenceGroups.insert(s2,set);
                } else if (!coincidenceGroups.contains(s1) && coincidenceGroups.contains(s2)) {
                    QSet<LineSegment*> set = coincidenceGroups.value(s2);
                    set.insert(s1);
                    coincidenceGroups.insert(s1,set);
                } else { // Both s1 and s2 already have a set. We must merge them.
                    QSet<LineSegment*> set1 = coincidenceGroups.value(s1);
                    QSet<LineSegment*> set2 = coincidenceGroups.value(s2);
                    set1.unite(set2);
                    foreach (LineSegment *ls, set1) {
                        coincidenceGroups.insert(ls,set1);
                    }
                }
            }
        }
    }
    int coinGpCount = 0;
    QSet<LineSegment*> seen;
    foreach (LineSegment *ls, coincidenceGroups.keys()) {
        if (seen.contains(ls)) continue;
        coinGpCount++;
        seen.unite(coincidenceGroups.value(ls));
    }
    m_most_recent_coincidence_group_count = coinGpCount;
    qDebug() << "coincidences:" << coincidences;
    qDebug() << "coincidence groups:" << coinGpCount;

    m_most_recent_crossing_count = crossings;
    m_most_recent_coincidence_count = coincidences;
    emit newCrossingCount(crossings);
    emit newCoincidenceCount(coincidences);

    // Angular resolution score
    computeNeighbourhoods();
    double angres = 0; int nodeCount = 0;
    double angres2 = 0; int node2Count = 0;
    double angres3 = 0; int node3Count = 0;
    double angres4 = 0; int node4Count = 0;
    foreach (CanvasItem *item, items())
    {
        if (ShapeObj *s = dynamic_cast<ShapeObj*>(item))
        {
            nodeCount++;
            QList<ShapeObj*> nbrs = m_align_nbrs.values(s);
            int n = nbrs.size();
            //qDebug() << "Shape" << s->internalId() << "has" << n << "nbrs.";
            if (n<2) continue;
            // Sort nbrs by angle
            NeighbourAngleLessThan compare(s);
            qSort(nbrs.begin(), nbrs.end(), compare);

            double ideal = 2*3.1415926/n;
            for (int i = 0; i < n; i++) {
                ShapeObj *t = nbrs.at(i);
                ShapeObj *u = nbrs.at((i+1)%n);
                double sx = s->centrePos().x(), sy = s->centrePos().y();
                double tx = t->centrePos().x(), ty = t->centrePos().y();
                double ux = u->centrePos().x(), uy = u->centrePos().y();
                double vx = tx-sx, vy = ty-sy;
                double wx = ux-sx, wy = uy-sy;
                double L = sqrt(vx*vx+vy*vy)*sqrt(wx*wx+wy*wy);
                if (L==0) continue;
                double c = (vx*wx+vy*wy)/L;
                if (c <= -1) c=-0.9999999;
                if (c >= 1)  c=0.9999999;
                //if (fabs(c) > 1) qDebug() << "acos arg outside [-1,1]!" << c;
                double theta = acos(c);
                double a = fabs(theta-ideal);
                angres += a;
                if (n==2) {
                    node2Count++;
                    angres2 += a;
                } else if (n==3) {
                    node3Count++;
                    angres3 += a;
                } else if (n==4) {
                    node4Count++;
                    angres4 += a;
                }
            }
        }
    }
    angres/=nodeCount;
    m_most_recent_angle_res_score = angres;
    qDebug() << "Angular resolution score:" << angres;

    if (node2Count>0) angres2/=node2Count;
    if (node3Count>0) angres3/=node3Count;
    if (node4Count>0) angres4/=node4Count;

    qDebug() << "Ang res for deg 2, 3, 4:" << angres2 << angres3 << angres4;

    m_most_recent_ang_res_by_degree.clear();
    m_most_recent_ang_res_by_degree.append(angres2);
    m_most_recent_ang_res_by_degree.append(angres3);
    m_most_recent_ang_res_by_degree.append(angres4);


    // Grid distance score
    double avggriddist = 0; nodeCount = 0;
    foreach (CanvasItem *item, items())
    {
        if (ShapeObj *s = dynamic_cast<ShapeObj*>(item))
        {
            nodeCount++;
            double W = m_opt_snap_grid_width, H = m_opt_snap_grid_height;
            double sx = s->centrePos().x(), sy = s->centrePos().y();
            double qx,rx,qy,ry;
            rx = fabs(modf(sx/W,&qx));
            ry = fabs(modf(sy/H,&qy));
            double dx=rx<=0.5?rx*W:(1-rx)*W;
            double dy=ry<=0.5?ry*H:(1-ry)*H;
            avggriddist += sqrt(dx*dx+dy*dy);
        }
    }
    avggriddist/=nodeCount;
    m_most_recent_avg_grid_distance = avggriddist;
    qDebug() << "Average grid distance:" << avggriddist;

    // Edge-node intersection
    int edgeNodeCrossings = 0;
    foreach (CanvasItem *item, items())
    {
        if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            ShapeObj *s = conn->getAttachedShapes().first, *t = conn->getAttachedShapes().second;
            double sx=s->centrePos().x(), sy=s->centrePos().y();
            double tx=t->centrePos().x(), ty=t->centrePos().y();
            Avoid::Point sp(sx,sy);
            Avoid::Point tp(tx,ty);
            //qDebug() << "sp,tp" << sx << sy << tx << ty;
            foreach (CanvasItem *item2, items())
            {
                if (ShapeObj *u = dynamic_cast<ShapeObj*>(item2))
                {
                    if (u==s || u==t) continue;
                    QRectF r = u->boundingRect();
                    r.moveCenter(u->centrePos());
                    double ax=r.left(), ay=r.top();
                    double bx=r.right(), by=r.bottom();
                    //qDebug() << "tl,br" << ax << ay << bx << by;
                    Avoid::Point tl(ax,ay);
                    Avoid::Point tr(bx,ay);
                    Avoid::Point bl(ax,by);
                    Avoid::Point br(bx,by);
                    if (Avoid::segmentIntersect(sp,tp,tl,tr) ||
                        Avoid::segmentIntersect(sp,tp,tr,br) ||
                        Avoid::segmentIntersect(sp,tp,br,bl) ||
                        Avoid::segmentIntersect(sp,tp,bl,tl)
                            ) {
                        edgeNodeCrossings++;
                    }
                }
            }
        }
    }
    m_most_recent_edge_node_overlap_count = edgeNodeCrossings;
    qDebug() << "Edge-node overlaps:" << edgeNodeCrossings;

    // Clean up
    foreach (LineSegment *s, segs) delete s;
    segs.clear();

    // Compute the score.
    OrthoWeights o = m_ortho_weights;
    //qDebug() << "Stress:" << m_most_recent_stress;
    double score = o.wcr*crossings + o.wco*coincidences + o.wob*obliquity + o.wst*m_most_recent_stress;
    emit newOrthoGoalBarValue( (int)(round(score)) );
    m_most_recent_ortho_obj_func = score;
    return score;
}

void Canvas::arrangePendants()
{
    // Why does this method get called twice when we press the button?
    if (!m_why_is_it_triggered_twice) {
        m_why_is_it_triggered_twice = true;
        return;
    }
    m_why_is_it_triggered_twice = false;

    // Build list of shapes, and map from each shape to list of all of its neighbours.
    // Keep track of max and min x and y coords.
    // Also record all connectors.
    QList<ShapeObj*> shapes;
    QMap<ShapeObj*,ShapeObj*> nbrs;
    QList<Connector*> conns;
    double xmin = DBL_MAX , xmax = DBL_MIN, ymin = DBL_MAX, ymax = DBL_MIN;
    foreach (CanvasItem *item, items())
    {
        if (ShapeObj *shape = dynamic_cast<ShapeObj*>(item))
        {
            shapes.append(shape);
            QPointF p = shape->centrePos();
            xmin = p.x() < xmin ? p.x() : xmin;
            xmax = p.x() > xmax ? p.x() : xmax;
            ymin = p.y() < ymin ? p.y() : ymin;
            ymax = p.y() > ymax ? p.y() : ymax;
        }
        else if (Connector *conn = dynamic_cast<Connector*>(item))
        {
            conns.append(conn);
            QPair<ShapeObj*, ShapeObj*> endpts = conn->getAttachedShapes();
            nbrs.insertMulti(endpts.first,endpts.second);
            nbrs.insertMulti(endpts.second,endpts.first);
        }
    }
    //qDebug() << "xmin" << xmin << "xmax" << xmax << "ymin" << ymin << "ymax" << ymax;
    // Now work out which grid locations are occupied.
    double W = m_opt_snap_grid_width, H = m_opt_snap_grid_height;
    int j0 = (int)(floor(xmin/W))-1, j1 = (int)(floor(xmax/W))+1;
    int i0 = (int)(floor(ymin/H))-1, i1 = (int)(floor(ymax/H))+1;
    int cols = j1 - j0 + 1, rows = i1 - i0 + 1;
    bool occupied[rows][cols];
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            occupied[i][j] = false;
        }
    }
    double tolerance = 0.25;
    foreach (ShapeObj *shape, shapes) {
        QPointF p = shape->centrePos();
        double j = p.x()/W, i = p.y()/H;
        double qj, rj, qi, ri;
        rj = modf(j,&qj);
        ri = modf(i,&qi);
        if (rj < -0.5) {
            qj -= 1; rj += 1;
        } else if (rj > 0.5) {
            qj += 1; rj -= 1;
        }
        if (ri < -0.5) {
            qi -= 1; ri += 1;
        } else if (ri > 0.5) {
            qi += 1; ri -= 1;
        }
        //qDebug() << "rj" << rj << "ri" << ri;
        if (fabs(rj) < tolerance && fabs(ri) < tolerance) {
            // E.g. if tolerance = 0.25 then we get here iff shape lies within
            // a quarter of a grid length in each dimension.
            int jz = (int)(qj), iz = (int)(qi);
            int ja = jz-j0, ia = iz-i0;
            occupied[ia][ja] = true;
        }
    }

    // check
    QString a = "Just nodes:\n";
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            a += occupied[i][j] ? "1 " : "0 ";
        }
        a += "\n";
    }
    qDebug() << a;

    // Also mark locations as occupied if an edge runs through.
    foreach (Connector *conn, conns) {
        ShapeObj *s = conn->getAttachedShapes().first;
        ShapeObj *t = conn->getAttachedShapes().second;
        double sx=s->centrePos().x(), sy=s->centrePos().y();
        double tx=t->centrePos().x(), ty=t->centrePos().y();
        double dy=sy-ty, dx=tx-sx;
        int a = (int)(ceil(atan2(dy,dx)*180/3.1415927));
        a = a < 0 ? a + 360 : a; // direction from s to t, from 0 to 360
        int eps = 5;
        if (eps < a%45 && a%45 < 45-eps) continue; // not within eps of a multiple of 45 deg

        // Get nearest grid point [si][sj] to s.
        double qsx, rsx, qsy, rsy;
        int sj, si;
        rsx = modf(sx/W,&qsx); rsy = modf(sy/H,&qsy);
        sj = fabs(rsx) < 0.5 ? (int)(qsx) : ( rsx > 0 ? (int)(qsx)+1 : (int)(qsx)-1 );
        sj -= j0;
        si = fabs(rsy) < 0.5 ? (int)(qsy) : ( rsy > 0 ? (int)(qsy)+1 : (int)(qsy)-1 );
        si -= i0;
        // Get nearest grid point [ti][tj] to t.
        double qtx, rtx, qty, rty;
        int tj, ti;
        rtx = modf(tx/W,&qtx); rty = modf(ty/H,&qty);
        tj = fabs(rtx) < 0.5 ? (int)(qtx) : ( rtx > 0 ? (int)(qtx)+1 : (int)(qtx)-1 );
        tj -= j0;
        ti = fabs(rty) < 0.5 ? (int)(qty) : ( rty > 0 ? (int)(qty)+1 : (int)(qty)-1 );
        ti -= i0;

        if (a < eps || a > 360-eps) {
            //
            //   s t
            //
            occupied[si  ][sj+1] = true;
            occupied[ti  ][tj-1] = true;
        } else if (a < 45+eps) {
            //     t
            //   s
            //
            occupied[si-1][sj+1] = true;
            occupied[ti+1][tj-1] = true;
        } else if (a < 90+eps) {
            //   t
            //   s
            //
            occupied[si-1][sj  ] = true;
            occupied[ti+1][tj  ] = true;
        } else if (a < 135+eps) {
            // t
            //   s
            //
            occupied[si-1][sj-1] = true;
            occupied[ti+1][tj+1] = true;
        } else if (a < 180+eps) {
            //
            // t s
            //
            occupied[si  ][sj-1] = true;
            occupied[ti  ][tj+1] = true;
        } else if (a < 225+eps) {
            //
            //   s
            // t
            occupied[si+1][sj-1] = true;
            occupied[ti-1][tj+1] = true;
        } else if (a < 270+eps) {
            //
            //   s
            //   t
            occupied[si+1][sj  ] = true;
            occupied[ti-1][tj  ] = true;
        } else if (a < 315+eps) {
            //
            //   s
            //     t
            occupied[si+1][sj+1] = true;
            occupied[ti-1][tj-1] = true;
        }
    }

    // check again
    a = "Wtih edges:\n";
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            a += occupied[i][j] ? "1 " : "0 ";
        }
        a += "\n";
    }
    qDebug() << a;

    // Now try to reposition each pendant node.
    QList<ShapePosInfo*> posinfos;
    foreach (ShapeObj *shape, shapes) {
        // Is it a pendant?
        QList<ShapeObj*> nbhd = nbrs.values(shape);
        if (nbhd.size() > 1 || nbhd.size() < 1) continue;
        // Get the owning shape.
        ShapeObj *owner = nbhd.at(0);
        // Get actual coords of p (pendant) and o (owner).
        double px = shape->centrePos().x(), py = shape->centrePos().y();
        double ox = owner->centrePos().x(), oy = owner->centrePos().y();
        // Get nearest grid point [pi][pj] to p.
        double qpx, rpx, qpy, rpy;
        int pj, pi;
        rpx = modf(px/W,&qpx); rpy = modf(py/H,&qpy);
        pj = fabs(rpx) < 0.5 ? (int)(qpx) : ( rpx > 0 ? (int)(qpx)+1 : (int)(qpx)-1 );
        pj -= j0;
        pi = fabs(rpy) < 0.5 ? (int)(qpy) : ( rpy > 0 ? (int)(qpy)+1 : (int)(qpy)-1 );
        pi -= i0;
        // Get nearest grid point [oi][oj] to o.
        double qox, rox, qoy, roy;
        int oj, oi;
        rox = modf(ox/W,&qox); roy = modf(oy/H,&qoy);
        oj = fabs(rox) < 0.5 ? (int)(qox) : ( rox > 0 ? (int)(qox)+1 : (int)(qox)-1 );
        oj -= j0;
        oi = fabs(roy) < 0.5 ? (int)(qoy) : ( roy > 0 ? (int)(qoy)+1 : (int)(qoy)-1 );
        oi -= i0;
        // If already ortho adjacent, skip.
        if ( ( pi==oi && (pj==oj-1 || pj==oj+1) ) ||
             ( pj==oj && (pi==oi-1 || pi==oi+1) )
                ) continue;
        // First consider whether any orthogonally adjacent position is available.
        QList<int> open;
        if (!occupied[oi  ][oj+1]) open.append(0);
        if (!occupied[oi-1][oj  ]) open.append(1);
        if (!occupied[oi  ][oj-1]) open.append(2);
        if (!occupied[oi+1][oj  ]) open.append(3);
        if (open.size()==0) {
            // There were no ortho. adj. positions open. Try diagonal.
            open.clear();
            if (!occupied[oi-1][oj+1]) open.append(4);
            if (!occupied[oi-1][oj-1]) open.append(5);
            if (!occupied[oi+1][oj-1]) open.append(6);
            if (!occupied[oi+1][oj+1]) open.append(7);
        }
        if (open.size()>0) {
            // At least one position is open. We will take one that is closest to (px,py).
            double bestX = px, bestY = py;
            int bestI = pi, bestJ = pj;
            double minDistSq = DBL_MAX;
            foreach (int pos, open) {
                double x=0, y=0;
                int i=0, j=0;
                switch (pos) {
                case 0:
                    y = oy  ; x = ox+W;
                    i = oi  ; j = oj+1;
                    break;
                case 1:
                    y = oy-H; x = ox  ;
                    i = oi-1; j = oj  ;
                    break;
                case 2:
                    y = oy  ; x = ox-W;
                    i = oi  ; j = oj-1;
                    break;
                case 3:
                    y = oy+H; x = ox  ;
                    i = oi+1; j = oj  ;
                    break;
                case 4:
                    y = oy-H; x = ox+W;
                    i = oi-1; j = oj+1;
                    break;
                case 5:
                    y = oy-H; x = ox-W;
                    i = oi-1; j = oj-1;
                    break;
                case 6:
                    y = oy+H; x = ox-W;
                    i = oi+1; j = oj-1;
                    break;
                case 7:
                    y = oy+H; x = ox+W;
                    i = oi+1; j = oj+1;
                    break;
                }
                double distSq = (x-px)*(x-px) + (y-py)*(y-py);
                if (distSq < minDistSq) {
                    bestX = x; bestY = y;
                    bestI = i; bestJ = j;
                    minDistSq = distSq;
                }
            }
            if (bestI != pi || bestJ != pj) {
                shape->setCentrePos(QPointF(bestX,bestY));
                //posinfos.append(m_graphlayout->makeShapePosInfo(shape,bestX,bestY));
                occupied[pi][pj] = 0;
                occupied[bestI][bestJ] = 1;
            }
        }
    }
    //m_graphlayout->processShapePosInfos(posinfos);
    fully_restart_graph_layout();

    //redraw_connectors(this);
    //reroute_all_connectors(this);
    //QCoreApplication::postEvent(this, new LayoutUpdateEvent(), Qt::LowEventPriority);


    //reroute_connectors(this);

}

void Canvas::applyFM3()
{
    m_bclayout->applyFM3();
}

void Canvas::layoutBCTrees()
{
    m_bclayout->layoutBCTrees();
}

void Canvas::BCWithFM3()
{
    m_bclayout->orthoLayout(1);
}

void Canvas::BCWithSpringEmbedder()
{
    m_bclayout->orthoLayout(2);
}

void Canvas::BCWithKamadaKawai()
{
    m_bclayout->orthoLayout(3);
}

void Canvas::BCWithPlanarization()
{
    m_bclayout->orthoLayout(4);
}

void Canvas::BCWithPlanarizationGrid()
{
    m_bclayout->orthoLayout(5);
}

void Canvas::BCOther()
{
    m_bclayout->orthoLayout(6);
}

}
// vim: filetype=cpp ts=4 sw=4 et tw=0 wm=0 cindent

