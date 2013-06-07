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

//! @file
//! Canvas class.  This is DunnartCanvas' equivalent of QGraphicsScene.

#ifndef CANVAS_H_
#define CANVAS_H_

#include <QGraphicsScene>
#include <QAction>
#include <QEvent>
#include <QList>
#include <QStack>
#include <QString>
#include <QDomDocument>
#include <QUndoCommand>
#include <QColor>
#include <QElapsedTimer>

#include "libavoid/geometry.h"

class QToolBar;
class QStatusBar;
class QSvgRenderer;
class QUndoStack;
class QUndoCommand;
class QFileInfo;
class QParallelAnimationGroup;

class BuiltinLayoutFileIOPlugin;
class BuiltinSVGFileIOPlugin;


namespace Avoid {
class Router;
}

namespace dunnart {

namespace gml {
class Graph;
}

class CanvasItem;
class Guideline;
class Connector;
class GraphLayout;
class SelectionResizeHandle;
class UndoMacro;

class BCLayout;
class ShapeObj;

class ConstraintRejectedEvent;

typedef QList<CanvasItem *> CanvasItemsList;

class Actions {
    public:
        unsigned int flags;
        CanvasItemsList moveList;
        CanvasItemsList resizeList;

        Actions();
        void clear(void);
        bool empty(void) const;
};

struct AlignDesc;

struct OrthoWeights {
    OrthoWeights() :
        wcr(22.0),
        wco(45.0),
        wob(1.0),
        wst(1.0) {}
    double wcr;
    double wco;
    double wob;
    double wst;
};

template<typename T>
struct Matrix2d
{
    int rows, cols;
    std::vector<T> data;
    Matrix2d() : rows(0), cols(0) {}
    Matrix2d(int rows, int cols) : rows(rows), cols(cols), data(rows*cols)
    { }

    T operator()(int i, int j) const
    {
        Q_ASSERT(i < rows);
        Q_ASSERT(j < cols);
        return data[i*cols+j];
    }
    T& operator()(int i, int j)
    {
        Q_ASSERT(i < rows);
        Q_ASSERT(j < cols);
        return data[i*cols+j];
    }

    QString toString() {
        QString s = "";
        s += "\n  ";
        for (int j=0; j<cols; j++) s += QString(" %1").arg(j,2);
        for (int i=0; i<rows; i++) {
            s += "\n";
            s += QString("%1").arg(i,2);
            for (int j=0; j<cols; j++) {
                s += QString(" %1").arg(data[i*cols+j],2);
            }
        }
        return s;
    }

};

static const unsigned int DEFAULT_CANVAS_FONT_SIZE = 11;

static const uint MESSAGEBOX_PIXMAP_SIZE = 70;

static const unsigned int ACTION_NONE           = 0;
static const unsigned int ACTION_ADDITIONS      = 1;
static const unsigned int ACTION_DELETIONS      = 2;
static const unsigned int ACTION_MODIFICATIONS  = 4;

static const int ModeSelection   = 1;
static const int ModeConnection  = 2;

enum loadPass
{
    PASS_SHAPES,
    PASS_CLUSTERS,
    PASS_CONNECTORS,
    PASS_RELATIONSHIPS,
    PASS_LAST
};

class Canvas : public QGraphicsScene
{
    Q_OBJECT

    Q_PROPERTY (bool automaticGraphLayout READ optAutomaticGraphLayout WRITE setOptAutomaticGraphLayout)
    Q_PROPERTY (LayoutMode layoutMode READ optLayoutMode WRITE setOptLayoutMode)
    Q_PROPERTY (FlowDirection flowDirection READ optFlowDirection WRITE setOptFlowDirection)
    Q_PROPERTY (double flowSeparationModifier READ optFlowSeparationModifier WRITE setOptFlowSeparationModifier)
    Q_PROPERTY (LayeredAlignment layeredAlignmentPosition READ optLayeredAlignmentPosition WRITE setOptLayeredAlignmentPosition)
    Q_PROPERTY (bool preventOverlaps READ optPreventOverlaps WRITE setOptPreventOverlaps)
    Q_PROPERTY (int shapeNonOverlapPadding READ optShapeNonoverlapPadding WRITE setOptShapeNonoverlapPadding)
    Q_PROPERTY (bool preserveTopology READ optPreserveTopology WRITE setOptPreserveTopology)
    Q_PROPERTY (bool rubberBandRouting READ optRubberBandRouting WRITE setOptRubberBandRouting)
    Q_PROPERTY (bool fitDiagramWithinPage READ optFitWithinPage WRITE setOptFitWithinPage)
    //Q_PROPERTY (bool colourInterferingConnectors READ optColourInterferingConnectors)
    Q_PROPERTY (double ideaEdgeLengthModifier READ optIdealEdgeLengthModifier WRITE setOptIdealEdgeLengthModifier)
    Q_PROPERTY (int routingShapePadding READ optRoutingShapePadding WRITE setOptRoutingShapePadding)
    Q_PROPERTY (int connectorRoundingDistance READ optConnectorRoundingDistance WRITE setOptConnRoundingDist)
    Q_PROPERTY (int routingSegmentPenalty READ optRoutingPenaltySegment WRITE setOptRoutingPenaltySegment)
    Q_PROPERTY (bool structuralEditingDisabled READ optStructuralEditingDisabled WRITE setOptStructuralEditingDisabled)
    Q_PROPERTY (QColor canvasBackgroundColour READ optCanvasBackgroundColour WRITE setOptCanvasBackgroundColour)
    Q_ENUMS (FlowDirection)
    Q_ENUMS (LayoutMode)
    Q_ENUMS (LayeredAlignment)

    public:

        Canvas();
        virtual ~Canvas();

        enum FlowDirection
        {
            FlowDown  = 0,
            FlowLeft  = 1,
            FlowUp    = 2,
            FlowRight = 3
        };

        enum LayeredAlignment
        {
            ShapeMiddle = 0,
            ShapeStart = 1,
            ShapeEnd = 2
        };

        enum LayoutMode
        {
            OrganicLayout = 0,
            FlowLayout    = 1,
            LayeredLayout = 2
        };

        enum Dimension
        {
            VERT,
            HORIZ
        };

        enum AlignmentFlags
        {
            Horizontal = 1,
            Vertical   = 2,
            Deliberate = 4,
            Connected  = 8
        };

        bool loadGmlDiagram(const QFileInfo& fileInfo);
        void loadSVGRootNodeAttributes(const QDomElement& svgRoot);
        void startLayoutUpdateTimer(void);
        void startLayoutFinishTimer(void);
        void setFilename(QString filename);
        QString filename(void);
        QList<CanvasItem *> items(void) const;
        QList<CanvasItem *> selectedItems(void) const;
        void setSelection(const QList<CanvasItem *>& newSelection);
        void postDiagramLoad(void);
        void setExpandedPage(const QRectF newExpandedPage);
        QRectF pageRect(void) const;
        void setPageRect(const QRectF &rect);
        GraphLayout *layout(void) const;
        Avoid::Router *router(void) const;
        void setStatusBar(QStatusBar *statusBar);
        void pushStatusMessage(const QString& message);
        void popStatusMessage(void);
        void highlightIndicatorsForItemMove(CanvasItem *item);
        void moveSelectionResizeHandle(const int index, const QPointF pos);
        void storeSelectionResizeInfo(void);
        QFont& canvasFont(void);
        QRectF combinedViewsRect(void) const;

        QString saveConstraintInfoToString(void) const;
        void loadConstraintInfoFromString(const QString& constraintInfo);

        // XXX These need to be reworked and given meaningful names:
        void interrupt_graph_layout(void);
        void stop_graph_layout(void);
        void restart_graph_layout(void);
        void fully_restart_graph_layout(void);

        CanvasItem *getItemByID(QString ID) const;
        CanvasItem *getItemByInternalId(uint internalId) const;
        void processSelectionDropEvent(QGraphicsSceneMouseEvent *event);
        
        bool optAutomaticGraphLayout(void) const;
        bool optPreventOverlaps(void) const;
        bool optSnapTo(void) const;
        bool optGridSnap(void) const;
        bool optRelax(void) const;
        bool optPreserveTopology(void) const;
        bool optRubberBandRouting(void) const;
        bool optFitWithinPage(void) const;
        bool optColourInterferingConnectors(void) const;
        bool optStructuralEditingDisabled(void) const;
        double optIdealEdgeLengthModifier(void) const;
        double optSnapDistanceModifier(void) const;
        double optSnapStrengthModifier(void) const;
        double optGridWidth(void) const;
        double optGridHeight(void) const;
        double optRelaxThresholdModifier(void) const;
        int optConnectorRoundingDistance(void) const;
        int optRoutingPenaltySegment(void) const;
        int optRoutingShapePadding(void) const;
        LayoutMode optLayoutMode(void) const;
        FlowDirection optFlowDirection(void) const;
        double optFlowSeparationModifier(void) const;
        int optShapeNonoverlapPadding(void) const;
        LayeredAlignment optLayeredAlignmentPosition(void) const;
        double optStressBarMaximum(void) const;
        double optObliquityBarMaximum(void) const;
        QColor optCanvasBackgroundColour(void) const;
        bool optDrawSeparationIndicators(void) const;
        bool optEdgeNodeRepulsion(void) const;

        bool overlayRouterRawRoutes(void) const;
        bool overlayRouterDisplayRoutes(void) const;
        bool overlayRouterObstacles(void) const;
        bool overlayRouterVisGraph(void) const;
        bool overlayRouterOrthogonalVisGraph(void) const;
        Actions& getActions(void);
        QString assignStringId(QString id);
        uint assignInternalId(void);
        int editMode(void) const;
        void setEditMode(int mode);

        void setLayoutSuspended(bool suspend);
        bool isLayoutSuspended(void) const;

        void setDraggedItem(CanvasItem *item, bool withForce = false);
        bool layoutRunningAndNotProcessingUpdates(void) const;

        QSvgRenderer *svgRenderer(void) const;
        QUndoStack *undoStack(void) const;
        UndoMacro *currentUndoMacro(void);
        UndoMacro *beginUndoMacro(const QString& text);
        void endUndoMacro(void);

        void saveDiagram(const QString& outputFilename);
        const QList<QColor> interferingConnectorColours(void) const;
        double visualPageBuffer(void) const;
        bool useGmlClusters(void) const;
        void setNudgeDistance(const double dist);
        void setIdealConnectorLength(const double length);
        double idealConnectorLength(void) const;
        bool avoidConnectorCrossings(void) const;
        bool avoidClusterCrossings(void) const;
        bool forceOrthogonalConnectors(void) const;
        void repositionAndShowSelectionResizeHandles(
                const bool calculatePosition = false);
        void exportDiagramToFile(QString filename);
        void paperExport(QString diagram, QString type, int actions = 1);
        bool tryingAlignments(void) const
        {
            return m_trying_alignments;
        }

        bool gd2013_batch_mode;
        uint gd2013_step;
        uint gd2013_type;
        QTimer *m_action_finished_timer;
        QElapsedTimer m_action_timer;

    public slots:
        void actionFinished(void);
        void bringToFront(void);
        void sendToBack(void);
        void deselectAll(void);
        void cutSelection(void);
        void copySelection(void);
        void pasteSelection(void);
        void deleteSelection(void);
        void deleteItem(CanvasItem *item);
        void deleteItems(QList<CanvasItem*> items);

        void toggleSelectedShapePinning(void);
        void selectAll(void);
        void templateFromSelection(int type);
        void alignSelection(int type);
        void distributeSelection(int type);
        void separateSelection(int type);

        void improveOrthogonalTopology(void);
        void inferAndApplyAlignments(void);
        void initTryAlignments(void);
        void tryAlignments(void);
        void applyAlignments(void);
        void rejectAlignments(void);
        void arrangePendants(void);
        void applyFM3(void);
        void layoutBCTrees(void);
        void BCWithFM3(void);
        void BCWithSpringEmbedder(void);
        void BCWithKamadaKawai(void);
        void BCWithPlanarization(void);
        void BCWithPlanarizationGrid(void);
        void BCOther(void);

        void customEvent(QEvent *event);
        void setOptIdealEdgeLengthModifierFromSlider(int int_modifier);
        void setOptIdealEdgeLengthModifier(double modifier);
        void setOptSnapDistanceModifierFromSlider(int int_modifier);
        void setOptSnapDistanceModifier(double modifier);
        void setOptSnapStrengthModifierFromSlider(int int_modifier);
        void setOptSnapStrengthModifier(double modifier);
        void setOptGridSizeFromSlider(int intValue);
        void setOptGridWidth(double value);
        void setOptGridHeight(double value);
        void setOptRelaxThresholdModifierFromSlider(int int_modifier);
        void setOptRelaxThresholdModifier(double modifier);
        void updateStress(double stress);

        void setDebugCOLAOutput(const bool value);
        void setOptAutomaticGraphLayout(const bool value);
        void setOptPreventOverlaps(const bool value);
        void setOptSnapTo(const bool value);
        void setOptGridSnap(const bool value);
        void setOptRelax(const bool value);
        void setOptPreserveTopology(const bool value);
        void setOptRubberBandRouting(const bool value);
        void setOptFitWithinPage(const bool value);
        void setOptRoutingPenaltySegment(const int value);
        void setOptRoutingShapePadding(const int value);
        void setOptConnRoundingDist(const int value);
        void setOptStructuralEditingDisabled(const bool value);
        void setOptLayoutMode(const LayoutMode mode);
        void setOptLayoutModeFromInt(const int mode);
        void setOptFlowSeparationModifier(const double value);
        void setOptFlowSeparationModifierFromSlider(const int intValue);
        void setOptFlowDirection(const FlowDirection value);
        void setOptFlowDirectionFromDial(const int value);
        void setOptShapeNonoverlapPadding(const int value);
        void setOptLayeredAlignmentPosition(const LayeredAlignment pos);
        void setOptCanvasBackgroundColour(const QColor colour);
        void setOptDrawSeparationIndicators(const bool value);
        void setOptEdgeNodeRepulsion(const bool value);

        void processResponseTasks(void);
        void processUndoResponseTasks(void);

        void setOverlayRouterRawRoutes(const bool value);
        void setOverlayRouterDisplayRoutes(const bool value);
        void setOverlayRouterObstacles(const bool value);
        void setOverlayRouterVisGraph(const bool value);
        void setOverlayRouterOrthogonalVisGraph(const bool value);

        // When rendering for printing, constraint indicators, selection
        // cues and other decorations are not painted.  This mode is used
        // for printing documents as well as exporting SVG, PDF and PS files.
        bool isRenderingForPrinting(void) const;
        void setRenderingForPrinting(const bool printingMode);
        bool inSelectionMode(void) const;
        void postRoutingRequiredEvent(void);

        double computeOrthoObjective(void);
        double predictOrthoObjective(Connector *conn, Dimension dim);
        double predictOrthoObjectiveChange(Connector *conn, Dimension dim);
        void predictOrthoObjectiveChange(QList<AlignDesc*> &ads);
        bool predictCoincidence(Connector *conn, Dimension dim);

    signals:
        void diagramFilenameChanged(const QFileInfo& title);
        void canvasDrawingChanged(void);
        void debugOverlayEnabled(bool enabled);
        void clipboardContentsChanged(void);
        void editModeChanged(const int mode);
        void layoutHasConverged(void);

        void optChangedAutomaticLayout(bool checked);
        void optChangedPreserveTopology(bool checked);
        void optChangedPreventOverlaps(bool checked);
        void optChangedSnapTo(bool checked);
        void optChangedGridSnap(bool checked);
        void optChangedRelax(bool checked);
        void optChangedRubberBandRouting(bool checked);
        void optChangedFitWithinPage(bool checked);
        void optChangedStructuralEditingDisabled(bool checked);
        void optChangedIdealEdgeLengthModifier(double value);
        void optChangedSnapDistanceModifier(double value);
        void optChangedSnapStrengthModifier(double value);
        void optChangedGridWidth(double value);
        void optChangedGridHeight(double value);
        void optChangedGridSize(double value);
        void optChangedRelaxThresholdModifier(double value);
        void optChangedLayoutMode(int mode);
        void optChangedDirectedEdgeSeparationModifier(double modifier);
        void optChangedFlowDirection(int direction);
        void optChangedRoutingShapePadding(int padding);
        void optChangedShapeNonoverlapPadding(int padding);
        void optChangedLayeredAlignmentPosition(LayeredAlignment pos);
        void optChangedDrawSeparationIndicators(bool checked);
        void optChangedEdgeNodeRepulsion(bool checked);
        void newStressBarValue(int value);
        void newCrossingCount(int value);
        void newCoincidenceCount(int value);
        void newObliquityBarValue(int value);
        void newOrthoGoalBarValue(int value);

    private slots:
        void processLayoutUpdateEvent(void);
        void processLayoutFinishedEvent(void);
        void selectionChangeTriggers(void);
    protected:
        virtual void drawBackground(QPainter *painter, const QRectF& rect);
        virtual void drawForeground(QPainter *painter, const QRectF& rect);
        virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

    private:
        void drawGridLines(QPainter *painter, const QRectF& rect);
        bool loadDiagram(const QString& filename);
        bool idIsUnique(QString id) const;
        void recursiveMapIDs(QDomNode start, const QString& ns, int pass);
        bool singlePropUpdateID(QDomElement& node, const QString& prop,
                const QString ns = QString());
        QDomElement writeLayoutOptionsToDomElement(QDomDocument& doc) const;
        void loadLayoutOptionsFromDomElement(const QDomElement& options);
        void setSvgRendererForFile(const QString& filename);
        void recursiveReadSVG(const QDomNode& start, const QString& dunnartNS,
                int pass);
        void setInterferingConnectorColours(const QString colourListString);
        void hideSelectionResizeHandles(void);
        void createIndicatorHighlightCache(void);
        void clearIndicatorHighlights(const bool clearCache = false);
        void glueObjectsToIndicators(void);
        bool hasVisibleOverlays(void) const;
        void updateConnectorsForLayout(void);

        QList<AlignDesc> inferAlignments(qreal tolerance=1.0);
        QList<AlignDesc> inferAlignOneDim(QList<ShapeObj*> shapes,
                              QMap<ShapeObj*,ShapeObj*> nbrs,
                              Dimension dim, qreal tolerance=1.0);
        void getRestrConnComp(QMap<ShapeObj*,ShapeObj*> nbrs, ShapeObj *sh,
                              QSet<ShapeObj*> &R, QSet<ShapeObj*> &C);


        double m_visual_page_buffer;
        QString m_filename;
        QTimer *m_layout_update_timer;
        QTimer *m_layout_finish_timer;
        bool m_processing_layout_updates;
        QString m_clipboard;
        QRectF m_page;
        QRectF m_expanded_page;
        GraphLayout* m_graphlayout;
        Avoid::Router *m_router;
        QSvgRenderer *m_svg_renderer;
        QStatusBar *m_status_bar;
        QStack<QString> m_status_messages;
        uint m_max_string_id;
        uint m_max_internal_id;
        gml::Graph *m_gml_graph;
        bool m_use_gml_clusters;

        double m_most_recent_stress;
        double m_most_recent_ortho_obj_func;
        int m_most_recent_crossing_count;
        int m_most_recent_coincidence_count;
        double m_stress_bar_maximum;
        double m_obliquity_bar_maximum;
        OrthoWeights m_ortho_weights;

        double m_connector_nudge_distance;
        double m_ideal_connector_length;
        double m_flow_separation_modifier;
        bool m_rectangle_constraint_test;
        bool m_sticky_nodes;
        bool m_downward_edges;
        bool m_avoid_connector_crossings;
        bool m_avoid_cluster_crossings;
        bool m_nudge_orthogonal_routes;
        bool m_simple_paths_during_layout;
        bool m_batch_diagram_layout;
        bool m_force_orthogonal_connectors;
        bool m_infer_tentative_alignments;
        bool m_why_is_it_triggered_twice;
        bool m_trying_alignments;
        bool m_rejecting_alignments;
        int m_max_align_tries;
        int m_num_align_tries;
        int m_max_shape_id;
        bool * m_align_pairs_tried;
        QMap<ShapeObj*,ShapeObj*> m_align_nbrs;
        QMap<uint,ShapeObj*> m_shapes_by_id;
        bool sideIsClear(ShapeObj* s, int side, double tolerance, ShapeObj* except=NULL);
        //int **m_alignment_state;
        Matrix2d<int> m_alignment_state;
        void updateAlignmentStates(ShapeObj *s, ShapeObj *t, AlignmentFlags a);
        void appliedAlignmentWasUnsat(ConstraintRejectedEvent *cre);
        void applyAlignmentsCallback(void);
        void initRejectAlignments(void);
        void rejectAlignmentsCallback(ConstraintRejectedEvent *cre);
        double m_previous_ortho_goal_value;
        double m_apply_alignments_epsilon;
        double m_max_actual_dG;
        bool m_opt_draw_separation_indicators;
        bool m_opt_edge_node_repulsion;
        double computeStress(void);
        void predictStressChange(QList<AlignDesc*> ads);
        double m_reject_phase_previous_goal_value;
        double m_reject_phase_previous_stress_value;
        AlignDesc *m_reject_phase_previous_align;
        int nonPendantDegree(ShapeObj *s);
        int numHAligns(ShapeObj *s);
        int numVAligns(ShapeObj *s);
        long m_tentative_guideline_timestamp;
        //bool neighbourAngleLessThan(ShapeObj *t, ShapeObj *u);
        ShapeObj *m_centre_shape_for_angles;
        void computeNeighbourhoods(void);
        double m_most_recent_angle_res_score;
        double m_most_recent_avg_grid_distance;
        int m_most_recent_edge_node_overlap_count;
        double m_most_recent_average_obliqueness;

        double m_opt_ideal_edge_length_modifier;
        double m_opt_snap_distance_modifier;
        double m_opt_snap_strength_modifier;
        double m_opt_snap_grid_width;
        double m_opt_snap_grid_height;
        double m_opt_relax_threshold_modifier;
        double m_opt_shape_nonoverlap_padding;
        int  m_opt_connector_rounding_distance;
        bool m_opt_automatic_graph_layout;
        bool m_opt_prevent_overlaps;
        bool m_opt_snap_to;
        bool m_opt_grid_snap;
        bool m_opt_relax;
        bool m_opt_preserve_topology;
        bool m_opt_rubber_band_routing;
        bool m_opt_fit_within_page;
        bool m_opt_colour_interfering_connectors;
        bool m_opt_stuctural_editing_disabled;
        int  m_opt_flow_direction;
        int m_opt_layered_alignment_position;
        QColor m_opt_canvas_background_colour;
        Actions m_actions;

        QMap<QString, QString> m_paste_id_map;
        QList<QString> m_paste_bad_constraint_ids;
        // list of nodes in other namespaces
        QList<QDomNode> m_external_node_list;
        QMap<QString, QString> m_extra_namespaces_map;

        // Default list of connector colors.
        QList<QColor> m_default_connector_colours;

        // List of connector colors.
        // This list is used for automatic coloring of connectors that cross or
        // have shared paths, by coloring interference graph.
        // It will be either as sepcified through the interferingConnectorColours
        // dunnart option in XML, set with the setInterferingConnectorColours()
        // fuction, or defaults to m_interfering_connector_colours.
        // Access via interferingConnectorColours().
        QList<QColor> m_interfering_connector_colours;

        QMap<int, Guideline *> m_vguides, m_hguides;
        CanvasItem *m_dragged_item;
        bool m_dragged_with_force;
        CanvasItem *m_lone_selected_item;
        QUndoStack *m_undo_stack;
        UndoMacro *m_current_undo_macro;

        QRectF m_selection_shapes_bounding_rect;
        QVector<SelectionResizeHandle *> m_selection_resize_handles;
        QVector<QRectF> m_selection_resize_info;
        bool m_hide_selection_handles;

        bool m_overlay_router_raw_routes;
        bool m_overlay_router_display_routes;
        bool m_overlay_router_obstacles;
        bool m_overlay_router_visgraph;
        bool m_overlay_router_orthogonal_visgraph;
        bool m_rendering_for_printing;
        int m_edit_mode;
        bool m_routing_event_posted;
        QFont *m_canvas_font;
        unsigned int m_canvas_font_size;
        QParallelAnimationGroup *m_animation_group;

        BCLayout *m_bclayout;

#ifdef FPSTIMER
        QElapsedTimer m_convergence_timer;
        unsigned int m_convergence_update_count;
        bool m_convergence_timer_running;
#endif

        friend class GraphLayout;
        friend class GraphData;
        friend class UndoMacro;
        friend class MainWindow;
        friend struct ShapePosInfo;
        friend class ObjectsRepositionedAnimation;
        friend class CanvasTabWidget;

        friend class ::BuiltinLayoutFileIOPlugin;
        friend class ::BuiltinSVGFileIOPlugin;
};


class LayoutUpdateEvent : public QEvent
{
    public:
        LayoutUpdateEvent() :
            QEvent((QEvent::Type) (QEvent::User + 1))
        {
        }
};
        
class LayoutFinishedEvent : public QEvent
{
    public:
        LayoutFinishedEvent() :
            QEvent((QEvent::Type) (QEvent::User + 2))
        {
        }
};

class RoutingRequiredEvent : public QEvent
{
    public:
        RoutingRequiredEvent() :
            QEvent((QEvent::Type) (QEvent::User + 3))
        {
        }
};

class ConstraintRejectedEvent : public QEvent
{
    public:
        ConstraintRejectedEvent() :
            QEvent((QEvent::Type) (QEvent::User + 4)),
            m_guideline(NULL),
            m_unsat(false)
        {
        }
        Guideline *m_guideline;
        ShapeObj *m_shape;
        bool m_unsat;
};

class TrialAlignmentEvent : public QEvent
{
public:
    TrialAlignmentEvent(int flag) :
        QEvent((QEvent::Type) (QEvent::User + 5)),
        m_flag(flag)
    {
    }
    int m_flag;
};

class Connector;

struct LineSegment {
    LineSegment() : connector(NULL) {}
    LineSegment(Connector *conn);
    LineSegment(double sx, double sy, double tx, double ty);
    void computeParameters(double sx, double sy, double tx, double ty);
    bool intersects(LineSegment *other, double tolerance);
    bool coincidesWith(LineSegment *other, double angleTolerance, double interceptTolerance);
    double obliquityScore(void);
    Connector *connector;
    // Angles are integers from 0 to 179 inclusive.
    int angle;
    double m_cos;
    double m_sin;
    // Intercept is an x-intercept unless angle == 0,
    // in which case it is a y-intercept.
    double intercept;
    // t0 and t1 are the parameters giving the
    // start and end points of the line segment,
    // if the line is parameterized by motion at
    // unit speed in angle direction from intercept point.
    double t0;
    double t1;

    // libavoid segment endpoints for segmentIntersect testing.
    Avoid::Point p1;
    Avoid::Point p2;
};

extern QRectF diagramBoundingRect(const QList<CanvasItem *>& list);


}

Q_DECLARE_METATYPE (dunnart::Canvas::FlowDirection)
Q_DECLARE_METATYPE (dunnart::Canvas::LayeredAlignment)
Q_DECLARE_METATYPE (dunnart::Canvas::LayoutMode)

#endif // CANVAS_H_
// vim: filetype=cpp ts=4 sw=4 et tw=0 wm=0 cindent

