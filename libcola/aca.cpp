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
        const std::vector<cola::Edge>& es,
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

ACALayout::~ACALayout(void)
{
    delete m_alignmentState;
    delete m_separationState;
}

void ACALayout::createAlignments(void)
{
    if (m_allAtOnce) {
        acaLoopAllAtOnce();
    } else {
        acaLoopOneByOne();
    }
}

void ACALayout::layout(void)
{
    layoutWithCurrentConstraints();
    createAlignments();
}

void ACALayout::addBendPointPenalty(bool b)
{
    m_addBendPointPenalty = b;
}

void ACALayout::postponeLeaves(bool b)
{
    m_postponeLeaves = b;
}

void ACALayout::useNonLeafDegree(bool b)
{
    m_useNonLeafDegree = b;
}

void ACALayout::allAtOnce(bool b)
{
    m_allAtOnce = b;
}

std::string ACALayout::writeAlignmentTable(void)
{
    std::string s = m_alignmentState->toString();
    return s;
}

std::string ACALayout::writeSeparationTable(void)
{
    std::string s = m_separationState->toString();
    return s;
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
    for (unsigned k = 0; k < m_ccs.size(); k++) {
        cola::CompoundConstraint *cc = m_ccs.at(k);
        cc->generateVariables(vpsc::XDIM, m_xvs);
        cc->generateVariables(vpsc::YDIM, m_yvs);
        cc->generateSeparationConstraints(vpsc::XDIM, m_xvs, xcs, m_rs);
        cc->generateSeparationConstraints(vpsc::YDIM, m_yvs, ycs, m_rs);
    }
    // Some extra variables may have been generated in each dimension, on top of
    // the one we have for each rectangle.
    m_numExtraXVars = m_xvs.size() - m_n;
    m_numExtraYVars = m_yvs.size() - m_n;
    // Store constraints by dimension and by equality vs. inequality constraints.
    for (unsigned k = 0; k < xcs.size(); k++) {
        vpsc::Constraint *c = xcs.at(k);
        if (c->equality) {
            m_xEqCs.push_back(c);
        } else {
            m_xIneqCs.push_back(c);
        }
    }
    for (unsigned k = 0; k < ycs.size(); k++) {
        vpsc::Constraint *c = ycs.at(k);
        if (c->equality) {
            m_yEqCs.push_back(c);
        } else {
            m_yIneqCs.push_back(c);
        }
    }
}

/**
 * This is for use only during initialisation of the state tables, when we have
 * to deal with the extra variables generated by the incoming constraints.
 * For that time we arrange the variables like this:
 *
 *     ... rectangles ... | ... extra x-vars ... | ... extra y-vars ...
 *
 * So for example if there were 5 rectangles, 2 extra x-vars, and 3 extra y-vars,
 * then the columns (and rows) in the state table corresponding to the 7 x- and 8 y-vars
 * would be:
 *     x: 0 1 2 3 4 | 5 6
 *     y: 0 1 2 3 4 | 7 8 9
 */
int ACALayout::adjustVarNumForExtraVars(vpsc::Dim dim, int k)
{
    if (dim==vpsc::YDIM && k >= m_n) k += m_numExtraXVars;
    return k;
}

void ACALayout::initStateTables(void)
{
    // Start by building tables large enough to handle each rectangle, as
    // well as each extra X-var, and each extra Y-var.
    int N = m_n + m_numExtraXVars + m_numExtraYVars;
    m_alignmentState  = new Matrix2d<int>(N,N);
    m_separationState = new Matrix2d<int>(N,N);
    // Initialise with zeroes.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            (*m_alignmentState)(i,j) = 0;
            (*m_separationState)(i,j) = 0;
        }
    }
    // Note connections in alignment state table.
    for (int j = 0; j < m_m; j++) {
        cola::Edge e = m_es.at(j);
        int src = e.first, tgt = e.second;
        (*m_alignmentState)(src,tgt) = ACACONN;
        (*m_alignmentState)(tgt,src) = ACACONN;
    }
    // Consider equality constraints in the x-dimension.
    for (unsigned k = 0; k < m_xEqCs.size(); k++) {
        vpsc::Constraint *c = m_xEqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        // Adjust for extra variables.
        l = adjustVarNumForExtraVars(vpsc::XDIM,l);
        r = adjustVarNumForExtraVars(vpsc::XDIM,r);
        if (gap==0) {
            // It is an alignment.
            recordAlignmentWithClosure(l,r,ACAVERT,N);
        } else {
            // It is a disalignment, or separation.
            ACASepFlags sf = gap > 0 ? ACAEAST : ACAWEST;
            recordSeparationWithClosure(l,r,sf,N);
        }
    }
    // Consider equality constraints in the y-dimension.
    for (unsigned k = 0; k < m_yEqCs.size(); k++) {
        vpsc::Constraint *c = m_yEqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        // Adjust for extra variables.
        l = adjustVarNumForExtraVars(vpsc::YDIM,l);
        r = adjustVarNumForExtraVars(vpsc::YDIM,r);
        if (gap==0) {
            // It is an alignment.
            recordAlignmentWithClosure(l,r,ACAHORIZ,N);
        } else {
            // It is a disalignment, or separation.
            ACASepFlags sf = gap > 0 ? ACASOUTH : ACANORTH;
            recordSeparationWithClosure(l,r,sf,N);
        }
    }
    // Consider inequality constraints in the x-dimension.
    for (unsigned k = 0; k < m_xIneqCs.size(); k++) {
        vpsc::Constraint *c = m_xIneqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        if (gap < 0) continue; // does not constrain r to be on one side of l
        // Adjust for extra variables.
        l = adjustVarNumForExtraVars(vpsc::XDIM,l);
        r = adjustVarNumForExtraVars(vpsc::XDIM,r);
        recordSeparationWithClosure(l,r,ACAEAST,N);
    }
    // Consider inequality constraints in the y-dimension.
    for (unsigned k = 0; k < m_yIneqCs.size(); k++) {
        vpsc::Constraint *c = m_yIneqCs.at(k);
        int l = c->left->id, r = c->right->id;
        double gap = c->gap;
        if (gap < 0) continue; // does not constrain r to be on one side of l
        // Adjust for extra variables.
        l = adjustVarNumForExtraVars(vpsc::YDIM,l);
        r = adjustVarNumForExtraVars(vpsc::YDIM,r);
        recordSeparationWithClosure(l,r,ACASOUTH,N);
    }
    // Record snapshot for debugging purposes.
    aStateBeforeChop = m_alignmentState->toString();
    sStateBeforeChop = m_separationState->toString();
    // Now that we have computed the transitive closure of all the passed
    // constraints, we will no longer need the rows and columns for the extra
    // variables, so we chop those off now.
    Matrix2d<int> *aState = new Matrix2d<int>(m_n,m_n);
    Matrix2d<int> *sState = new Matrix2d<int>(m_n,m_n);
    for (int i = 0; i < m_n; i++) {
        for (int j = 0; j < m_n; j++) {
            (*aState)(i,j) = (*m_alignmentState)(i,j);
            (*sState)(i,j) = (*m_separationState)(i,j);
        }
    }
    delete m_alignmentState;
    delete m_separationState;
    m_alignmentState  = aState;
    m_separationState = sState;
}

void ACALayout::recordAlignmentWithClosure(int i, int j, ACAFlags af, int numCols)
{
    if (numCols == 0) numCols = m_n;
    // Get the set of all indices already aligned with i, including i itself.
    // Do likewise for j.
    std::set<int> Ai, Aj;
    Ai.insert(i);
    Aj.insert(j);
    for (int k = 0; k < numCols; k++) {
        if ((*m_alignmentState)(i,k) & af) Ai.insert(k);
        if ((*m_alignmentState)(j,k) & af) Aj.insert(k);
    }
    // Now record that everything in Ai is aligned with everything in Aj.
    for (std::set<int>::iterator it=Ai.begin(); it!=Ai.end(); ++it) {
        for (std::set<int>::iterator jt=Aj.begin(); jt!=Aj.end(); ++jt) {
            (*m_alignmentState)(*it,*jt) |= af;
            (*m_alignmentState)(*jt,*it) |= af;
        }
    }
}

ACASepFlags negateSepFlag(ACASepFlags sf) {
    unsigned short c = (unsigned short) sf;
    c += 16*c;
    c &= 60; // 00111100
    unsigned short b = c >> 2;
    ACASepFlags nf = (ACASepFlags) b;
    return nf;
}

void ACALayout::recordSeparationWithClosure(int i, int j, ACASepFlags sf, int numCols)
{
    if (numCols == 0) numCols = m_n;
    // Reduce to the case where sf is either ACAEAST or ACASOUTH.
    switch (sf) {
    case ACANOSEP:
        // The separation table is monotonic; you cannot remove separations
        // that already exist.
        return;
    case ACANORTH:
        recordSeparationWithClosure(j,i,ACASOUTH);
        return;
    case ACAWEST:
        recordSeparationWithClosure(j,i,ACAEAST);
        return;
    case ACANORTHEAST:
        recordSeparationWithClosure(i,j,ACANORTH);
        recordSeparationWithClosure(i,j,ACAEAST);
        return;
    case ACASOUTHEAST:
        recordSeparationWithClosure(i,j,ACASOUTH);
        recordSeparationWithClosure(i,j,ACAEAST);
        return;
    case ACASOUTHWEST:
        recordSeparationWithClosure(i,j,ACASOUTH);
        recordSeparationWithClosure(i,j,ACAWEST);
        return;
    case ACANORTHWEST:
        recordSeparationWithClosure(i,j,ACANORTH);
        recordSeparationWithClosure(i,j,ACAWEST);
        return;
    case ACAEAST:
    case ACASOUTH:
        // We drop down to this point only if sf is ACAEAST or ACASOUTH.
        // The code is the same for both cases.
        // For simplicity we express the comments for the case ACAEAST only.
        // Let L be the set of all indices west or equal to i;
        // let U be the set of all indices east or equal to j.
        std::set<int> L, U;
        L.insert(i);
        U.insert(j);
        for (int k = 0; k < numCols; k++) {
            if ((*m_separationState)(k,i) & sf) L.insert(k);
            if ((*m_separationState)(j,k) & sf) U.insert(k);
        }
        // Now record that everything in L is west of everything in U.
        ACASepFlags nf = negateSepFlag(sf);
        for (std::set<int>::iterator it=L.begin(); it!=L.end(); ++it) {
            for (std::set<int>::iterator jt=U.begin(); jt!=U.end(); ++jt) {
                (*m_separationState)(*it,*jt) |= sf;
                (*m_separationState)(*jt,*it) |= nf;
            }
        }
    }
}

void ACALayout::layoutWithCurrentConstraints(void)
{
    cola::ConstrainedFDLayout *fdlayout =
#define OLDCOLA
#ifdef OLDCOLA
            new cola::ConstrainedFDLayout(m_rs,m_es,m_idealLength,
                                          m_preventOverlaps);
#else
            new cola::ConstrainedFDLayout(m_rs,m_es,m_idealLength,
                                          m_preventOverlaps,m_edgeLengths,
                                          m_doneTest,m_preIteration);
#endif
    fdlayout->setConstraints(m_ccs);
    fdlayout->run(true,true);
    delete fdlayout;
}

void ACALayout::acaLoopOneByOne(void)
{
    // Choose a first alignment.
    OrderedAlignment *oa = chooseOA();
    while (oa) {
        // Add the new separated alignment constraints.
        m_ccs.push_back(oa->separation);
        m_ccs.push_back(oa->alignment);
        // Redo the layout, with the new constraints.
        layoutWithCurrentConstraints();
        // Update state tables.
        updateStateTables(oa);
        // Choose next ordered alignment.
        oa - chooseOA();
    }
}

void ACALayout::acaLoopAllAtOnce(void)
{
    // Choose a first alignment.
    OrderedAlignment *oa = chooseOA();
    while (oa) {
        // Add the new separated alignment constraints.
        m_ccs.push_back(oa->separation);
        m_ccs.push_back(oa->alignment);
        // Update state tables.
        updateStateTables(oa);
        // Choose next ordered alignment.
        oa - chooseOA();
    }
    // Redo the layout, with the new constraints.
    layoutWithCurrentConstraints();
}

void ACALayout::updateStateTables(OrderedAlignment *oa)
{
    recordAlignmentWithClosure(oa->rect1,oa->rect2,oa->af);
    recordSeparationWithClosure(oa->rect1,oa->rect2,oa->sf);
}

OrderedAlignment *ACALayout::chooseOA(void)
{
    // TODO
}


} // namespace cola
































