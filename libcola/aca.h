/*
 * vim: ts=4 sw=4 et tw=0 wm=0
 *
 * libcola - A library providing force-directed network layout using the 
 *           stress-majorization method subject to separation constraints.
 *
 * Copyright (C) 2006-2014  Monash University
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * See the file LICENSE.LGPL distributed with the library.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author(s):  Tim Dwyer
 *             Michael Wybrow
 *             Steve Kieffer
*/

#ifndef ACA_H
#define ACA_H

#include <vector>

#include "libcola/cola.h"

namespace cola {

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
        assert(i < rows);
        assert(j < cols);
        return data[i*cols+j];
    }
    T& operator()(int i, int j)
    {
        assert(i < rows);
        assert(j < cols);
        return data[i*cols+j];
    }

    std::string toString() {
        std::string s = "";
        s += "\n  ";
        char buffer [10];
        for (int j=0; j<cols; j++) {
            sprintf(buffer," %2d",j);
            s += std::string(buffer);
        }
        for (int i=0; i<rows; i++) {
            s += "\n";
            sprintf(buffer,"%2d",i);
            s += std::string(buffer);
            for (int j=0; j<cols; j++) {
                sprintf(buffer," %2d",data[i*cols+j]);
                s += std::string(buffer);
            }
        }
        return s;
    }

};

enum ACAFlags {
    ACAHORIZ = 1,
    ACAVERT  = 2,
    ACADELIB = 4,
    ACACONN  = 8
};

enum ACASepFlags {
    ACANOSEP     =  0,
    ACANORTH     =  1,
    ACAEAST      =  2,
    ACASOUTH     =  4,
    ACAWEST      =  8,
    ACANORTHEAST =  3,
    ACASOUTHEAST =  6,
    ACANORTHWEST =  9,
    ACASOUTHWEST = 12
};

struct OrderedAlignment {
    cola::SeparationConstraint* separation;
    cola::AlignmentConstraint* alignment;
    int rect1;
    int rect2;
    ACAFlags af;
    ACASepFlags sf;
};


/**
 * @brief Implements the Adaptive Constrained Alignment (ACA) algorithm.
 *
 * See
 * Kieffer, Steve, Tim Dwyer, Kim Marriott, and Michael Wybrow.
 * "Incremental grid-like layout using soft and hard constraints." In Graph
 * Drawing 2013, pp. 448-459. Springer International Publishing, 2013.
 */
class ACALayout {
public:
    /**
     * @brief Constructs an adaptive constrained alignment layout instance.
     *
     * Parameters are the same as for the ConstrainedFDLayout constructor,
     * with the addition of a vector of constraints passed by reference.
     *
     * If the vector of constraints is non-empty, these constraints will be
     * applied throughout the ACA process, and the new constraints created
     * by ACA will not conflict with any of these.
     *
     * Conversely, if no constraints are passed, the empty vector is still
     * used to store the new constraints created by ACA.
     *
     * @param[in] rs  Bounding boxes of nodes at their initial positions.
     * @param[in] es  Simple pair edges, giving indices of the start and end 
     *                nodes in rs.
     * @param[in] ccs  Vector of any pre-existing constraints, and holder where
     *                 new constraints created by ACA will be recorded.
     * @param[in] idealLength  A scalar modifier of ideal edge lengths in 
     *                         eLengths or 1 if no ideal lengths are specified..
     * @param preventOverlaps  Causes non-overlap constraints to be generated 
     *                          for all rectangles, if it is set to true.
     * @param[in] eLengths  Individual ideal lengths for edges.
     *                      The actual ideal length used for the ith edge is 
     *                      idealLength*eLengths[i], or if eLengths is NULL a
     *                      then just idealLength is used (i.e., eLengths[i] 
     *                      is assumed to be 1).
     * @param[in] done  A test of convergence operation called at the end of 
     *                  each iteration (optional).  If not given, uses a
     *                  default TestConvergence object.
     * @param[in] preIteration  An operation called before each iteration
     *                          (optional).
     */
    ACALayout(
        const vpsc::Rectangles& rs,
        const std::vector<cola::Edge>& es,
        CompoundConstraints& ccs,
        const double idealLength,
        const bool preventOverlaps,
        const EdgeLengths& eLengths = StandardEdgeLengths,
        TestConvergence* doneTest = NULL,
        PreIteration* preIteration=NULL);
    ~ACALayout();
    /**
     * @brief Creates alignments.
     *
     * This is the main functionality of ACA. This function should be called
     * on an existing layout in order to greedily align edges until any further
     * alignments would create edge overlaps.
     *
     * If the graph does not have an initial layout already, then the 'layout'
     * function may be called instead.
     */
    void createAlignments(void);
    /**
     * @brief Do an initial layout, and then create alignments.
     *
     * This is a convenience function which first does a constrained force-directed
     * layout of the given graph, and then calls the 'createAlignments' function.
     */
    void layout(void);

    // For debugging:
    std::string writeAlignmentTable(void);
    std::string writeSeparationTable(void);
    std::string aStateBeforeChop;
    std::string sStateBeforeChop;

private:
    /**
     * Used by the constructor to compute the degrees of nodes in the graph, and
     * determine which are the leaves, and which are the nodes of degree 2, both
     * in the sense of ordinary degree and non-leaf degree.
     */
    void computeDegrees(void);
    /**
     * Used by the constructor to convert the incoming cola constraints into vpsc
     * constraints, and the generate any additional variables that are needed by
     * those constraints.
     */
    void generateVPSCConstraints(void);
    /**
     * Used only during initialisation of the state tables, to handle any additional
     * variables generated by the incoming cola constraints.
     */
    int adjustVarNumForExtraVars(vpsc::Dim dim, int k);
    /**
     * Used by the constructor to initialise the state tables, computing the transitive
     * closure of all incoming alignments and separations, based on the vpsc constraints
     * to which the incoming cola constraints compile.
     */
    void initStateTables(void);
    /**
     * Record the specified alignment between rectangles i and j.
     * Also record all additional alignments arising from the transitive
     * closure.
     *
     * The numCols variable is only there to manage the state table
     * initialisation process, where we have to deal with the extra variables
     * generated by the incoming constraints.
     */
    void recordAlignmentWithClosure(int i, int j, ACAFlags af, int numCols = 0);
    /**
     * Record the specified separation between rectangles i and j.
     * The ACASepFlag sf names a compass direction D, and the understanding
     * is that the direction from rectangle i to rectangle j is D.
     * For example, to record that the x-coord of j is >= that of i,
     * the appropriate ACASepFlag is ACAEAST.
     *
     * If there is an existing separation in the alternate dimension
     * to the one passed, the two will be combined.
     *
     * Also record all additional separations arising from the transitive
     * closure.
     *
     * The numCols variable is only there to manage the state table
     * initialisation process, where we have to deal with the extra variables
     * generated by the incoming constraints.
     */
    void recordSeparationWithClosure(int i, int j, ACASepFlags sf, int numCols = 0);


    // TODO ? ---------------------------------------------
    void initialPositions(void);
    void moveCoincidentNodes(void);
    void initialLayout(void);
    void acaLoopOneByOne(void);
    void acaLoopAllAtOnce(void);
    void finalLayout(void);
    void updateAlignmentState(OrderedAlignment *oa);
    OrderedAlignment *chooseOA(void);
    bool createsCoincidence(int src, int tgt, ACAFlags af);
    double deflection(int src, int tgt, ACAFlags af);
    double bendPointPenalty(int src, int tgt, ACAFlags af);
    double leafPenalty(int src, int tgt);
    // -----------------------------------------------------


    int m_n; // number of nodes
    int m_m; // number of edges
    int m_numExtraXVars;
    int m_numExtraYVars;
    vpsc::Rectangles m_rs;
    std::vector<cola::Edge> m_es;
    cola::CompoundConstraints m_ccs;

    vpsc::Variables m_xvs;
    vpsc::Variables m_yvs;
    vpsc::Constraints m_xEqCs;
    vpsc::Constraints m_xIneqCs;
    vpsc::Constraints m_yEqCs;
    vpsc::Constraints m_yIneqCs;

    double m_idealLength;
    bool m_preventOverlaps;
    EdgeLengths m_edgeLengths;
    TestConvergence *m_doneTest;
    PreIteration *m_preIteration;
    double m_nodePadding;

    // Configuration:
    bool m_addBendPointPenalty;
    bool m_postponeLeaves;
    bool m_useNonLeafDegree;
    bool m_allAtOnce;

    std::multimap<int,int> m_nbrs; // neighbours
    std::multimap<int,int> m_nlnbrs; // non-leaf neighbours
    std::set<int> m_leaves; // indices of rectangles that are leaves
    std::set<int> m_deg2Nodes; // degree-2 nodes w.r.t. ordinary nbrs
    std::set<int> m_nldeg2Nodes; // degree-2 nodes w.r.t. only non-leaf neighbours

    Matrix2d<int> *m_alignmentState;
    Matrix2d<int> *m_separationState;
    std::vector<OrderedAlignment*> m_ordAligns;
};

} // namespace cola


#endif // ACA_H
