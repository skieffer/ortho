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

};

enum ACAFlags {
    ACAHORIZ = 1,
    ACAVERT  = 2,
    ACADELIB = 4,
    ACACONN  = 8
};

enum ACASepFlags {
    NONE      =  0,
    NORTH     =  1,
    NORTHEAST =  2,
    EAST      =  3,
    SOUTHEAST =  4,
    SOUTH     = -1,
    SOUTHWEST = -2,
    WEST      = -3,
    NORTHWEST = -4
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
 * Drawing, pp. 448-459. Springer International Publishing, 2013.
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

private:
    void computeDegrees(void);
    void initStateTables(void);
    void generateVPSCConstraints(void);


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

    unsigned m_n; // number of nodes
    unsigned m_m; // number of edges
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
    double *m_edgeLengths;
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
    std::set<int> m_nldeg2Nodes; // same but only non-leaf neighbours

    Matrix2d<int> m_alignmentState;
    Matrix2d<int> m_separationState;
    std::vector<OrderedAlignment*> m_ordAligns;
};

} // namespace cola


#endif // ACA_H
