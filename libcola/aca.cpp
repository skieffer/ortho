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

#include "libcola/aca.h"

using namespace std;
using namespace vpsc;

namespace cola {

ACALayout::ACALayout(
        const vpsc::Rectangles& rs,
        const std::vector& es,
        CompoundConstraints& ccs,
        const double idealLength,
        const bool preventOverlaps,
        const EdgeLengths& eLengths,
        TestConvergence *doneTest,
        PreIteration *preIteration)
    : m_n(rs.size()),
      m_m(es.size()),
      m_rs(rs),
      m_es(es),
      m_ccs(ccs),
      m_idealLength(idealLength),
      m_preventOverlaps(preventOverlaps),
      m_edgeLengths(eLengths),
      m_doneTest(doneTest),
      m_preIteration(preIteration),
      m_nodePadding(0),
      m_addBendPointPenalty(true),
      m_postponeLeaves(true),
      m_useNonLeafDegree(true),
      m_allAtOnce(false)
{
    computeDegrees();
    generateVPSCConstraints();
    initStateTables();
}

void ACALayout::computeDegrees(void)
{
    // Map node indices to indices of their neighbours.
    for (int j = 0; j < m_m; j++) {
        cola::Edge e = m_es.at(j);
        m_nbrs.insert(std::pair<int,int>(e.first,e.second));
        m_nbrs.insert(std::pair<int,int>(e.second,e.first));
    }
    // Find the leaves and degree-2 nodes.
    for (int i = 0; i < m_n; i++) {
        int c = m_nbrs.count(i);
        if (c == 1) m_leaves.insert(i);
        if (c == 2) m_deg2Nodes.insert(i);
    }
    // Compute non-leaf neighbours.
    for (int i = 0; i < m_n; i++) {
        // Iterate over the neighbours of node i, retaining only the non-leaves.
        std::pair< std::multimap<int,int>::iterator, std::multimap<int,int>::iterator > range;
        range = m_nbrs.equal_range(i);
        for (std::multimap<int,int>::iterator it=range.first; it!=range.second; ++it) {
            int j = it->second;
            if (m_leaves.count(j)==0) {
                m_nlnbrs.insert(std::pair<int,int>(i,j));
            }
        }
    }
    // Find the non-leaf degree-2 nodes.
    for (int i = 0; i < m_n; i++) {
        if (m_nlnbrs.count(i) == 2) m_nldeg2Nodes.insert(i);
    }
}

void ACALayout::generateVPSCConstraints(void)
{
    // First generate x-variables and y-variables.
    for (int i = 0; i < m_n; i++) {
        m_xvs.push_back(new vpsc::Variable(i));
        m_yvs.push_back(new vpsc::Variable(i));
    }
    // Generate all VPSC constraints.
    vpsc::Constraints xcs, ycs;
    for (int k = 0; k < m_ccs.size(); k++) {
        cola::CompoundConstraint *cc = m_ccs.at(k);
        cc->generateSeparationConstraints(vpsc::XDIM, m_xvs, xcs, m_rs);
        cc->generateSeparationConstraints(vpsc::YDIM, m_yvs, ycs, m_rs);
    }
    // Store by dimension and by equality vs. inequality constraints.
    for (int k = 0; k < xcs.size(); k++) {
        vpsc::Constraint *c = xcs.at(k);
        if (c->equality) {
            m_xEqCs.push_back(c);
        } else {
            m_xIneqCs.push_back(c);
        }
    }
    for (int k = 0; k < ycs.size(); k++) {
        vpsc::Constraint *c = ycs.at(k);
        if (c->equality) {
            m_yEqCs.push_back(c);
        } else {
            m_yIneqCs.push_back(c);
        }
    }
}

void ACALayout::initStateTables(void)
{
    m_alignmentState = Matrix2d<int>(m_n,m_n);
    m_separationState = Matrix2d<int>(m_n,m_n);
    // Initialise with zeroes.
    for (int i = 0; i < m_n; i++) {
        for (int j = 0; j < m_n; j++) {
            m_alignmentState(i,j) = 0;
            m_separationState(i,j) = 0;
        }
    }
    // Note connections in alignment state table.
    for (int j = 0; j < m_m; j++) {
        cola::Edge e = m_es.at(j);
        int src = e.first, tgt = e.second;
        m_alignmentState(src,tgt) = ACACONN;
        m_alignmentState(tgt,src) = ACACONN;
    }
    // Consider equality constraints in the x-dimension.
    for (int k = 0; k < m_xEqCs.size(); k++) {
        vpsc::Constraint *c = m_xEqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        if (gap==0) {
            // It is an alignment.
            recordAlignmentWithClosure(l,r,ACAVERT);
        } else {
            // It is a disalignment, or separation.
            ACASepFlags sf = gap > 0 ? ACAEAST : ACAWEST;
            recordSeparationWithClosure(l,r,sf);
        }
    }
    // Consider equality constraints in the y-dimension.
    for (int k = 0; k < m_yEqCs.size(); k++) {
        vpsc::Constraint *c = m_yEqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        if (gap==0) {
            // It is an alignment.
            recordAlignmentWithClosure(l,r,ACAHORIZ);
        } else {
            // It is a disalignment, or separation.
            ACASepFlags sf = gap > 0 ? ACASOUTH : ACANORTH;
            recordSeparationWithClosure(l,r,sf);
        }
    }
    // Consider inequality constraints in the x-dimension.
    for (int k = 0; k < m_xIneqCs.size(); k++) {
        vpsc::Constraint *c = m_xIneqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        if (gap < 0) continue; // does not constrain r to be on one side of l
        recordSeparationWithClosure(l,r,ACAEAST);
    }
    // Consider inequality constraints in the y-dimension.
    for (int k = 0; k < m_yIneqCs.size(); k++) {
        vpsc::Constraint *c = m_yIneqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        if (gap < 0) continue; // does not constrain r to be on one side of l
        recordSeparationWithClosure(l,r,ACASOUTH);
    }
}

void ACALayout::recordAlignmentWithClosure(int i, int j, ACAFlags af)
{
    // TODO
}

void ACALayout::recordSeparationWithClosure(int i, int j, ACASepFlags sf)
{
    // TODO
}


} // namespace cola
































