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

ACASepFlags negateSepFlag(ACASepFlags sf)
{
    unsigned short c = (unsigned short) sf;
    c += 16*c;
    c &= 60; // 00111100
    unsigned short b = c >> 2;
    ACASepFlags nf = (ACASepFlags) b;
    return nf;
}

ACAFlags sepToAlignFlag(ACASepFlags sf)
{
    return sf==ACANORTH || sf==ACASOUTH ? ACAVERT : ACAHORIZ;
}

ACAFlags perpAlignFlag(ACAFlags af)
{
    return af==ACAHORIZ ? ACAVERT : ACAHORIZ;
}

ACASepFlags vectorToSepFlag(double dx, double dy)
{
    int f = 0;
    f |= dx > 0 ? ACAEAST : dx < 0 ? ACAWEST : 0;
    f |= dy > 0 ? ACASOUTH : dy < 0 ? ACANORTH : 0;
    return (ACASepFlags) f;
}

bool propsedSepConflictsWithExistingPosition(ACASepFlags pro, ACASepFlags ex)
{
    // Proposed separation pro conflicts with existing position ex
    // if the union of their bits contains both north and south, or
    // both east and west.
    int u = pro | ex;
    return ( ( (u&5) == 5 ) || ( (u&10) == 10 ) );
}

bool propsedSepConflictsWithExistingSepCo(ACASepFlags pro, ACASepFlags ex)
{
    // Proposed separation pro conflicts with existing separation ex
    // if ex has any of the complementary bits of pro.
    int proComp = 15 - pro;
    return (proComp & ex);
}

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
      m_allAtOnce(false),
      m_aggressiveOrdering(false)
{
    computeDegrees();
    generateVPSCConstraints();
    initStateTables();
}

ACALayout::~ACALayout(void)
{
    delete m_alignmentState;
    delete m_separationState;
    delete m_fdlayout;
}

void ACALayout::createAlignments(void)
{
    if (m_allAtOnce) {
        acaLoopAllAtOnce();
    } else {
        acaLoopOneByOne();
    }
}

bool ACALayout::createOneAlignment(void)
{
    return acaLoopOnce();
}

void ACALayout::layout(void)
{
    layoutWithCurrentConstraints();
    createAlignments();
}

cola::ConstrainedFDLayout *ACALayout::getFDLayout(void)
{
    return m_fdlayout;
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

void ACALayout::aggressiveOrdering(bool b)
{
    m_aggressiveOrdering = b;
}

void ACALayout::setAlignmentOffsetsForCompassDirection(ACASepFlags sf, EdgeOffsets offsets)
{
    assert(offsets.size()==(size_t)m_m); // There should be one offset for each edge.
    m_edgeOffsets.insert( std::pair<ACASepFlags,EdgeOffsets>(sf,offsets) );
}

void ACALayout::setAllowedSeparations(std::vector<ACASepFlags> seps)
{
    assert(seps.size()==(size_t)m_m); // There should be one flag for each edge.
    m_allowedSeps = seps;
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
    // Get the set Ai of all indices already aligned with i, including i itself.
    // Do likewise for j.
    // Also if af is ACAVERT then let L be the set of all indices west of i and
    // U the set of all indices east of j. Do similarly if af is ACAHORIZ.
    std::set<int> Ai, Aj;
    Ai.insert(i);
    Aj.insert(j);
    std::set<int> U, L;
    ACASepFlags sf = af == ACAVERT ? ACAEAST : ACASOUTH;
    for (int k = 0; k < numCols; k++) {
        if ( (*m_alignmentState)(i,k) & af ) Ai.insert(k);
        if ( (*m_alignmentState)(j,k) & af ) Aj.insert(k);
        if ( (*m_separationState)(k,i) & sf ) L.insert(k);
        if ( (*m_separationState)(j,k) & sf ) U.insert(k);
    }
    // Record that everything in Ai is aligned with everything in Aj.
    // This is the transitive closure of the new alignment for the alignments table.
    for (std::set<int>::iterator it=Ai.begin(); it!=Ai.end(); ++it) {
        for (std::set<int>::iterator jt=Aj.begin(); jt!=Aj.end(); ++jt) {
            (*m_alignmentState)(*it,*jt) |= af;
            (*m_alignmentState)(*jt,*it) |= af;
        }
    }
    // Record that everything in L goes in direction sf to everything in U.
    // This is the transitive closure of the new alignment for the separation table.
    ACASepFlags nf = negateSepFlag(sf);
    for (std::set<int>::iterator it=L.begin(); it!=L.end(); ++it) {
        for (std::set<int>::iterator jt=U.begin(); jt!=U.end(); ++jt) {
            (*m_separationState)(*it,*jt) |= sf;
            (*m_separationState)(*jt,*it) |= nf;
        }
    }
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
    case ACAEAST: // Drop to next case.
    case ACASOUTH: {
        // We drop down to this point only if sf is ACAEAST or ACASOUTH.
        // The code is the same for both cases.
        // For simplicity we express the comments for the case ACAEAST only.
        // Let L be the set of all indices west, aligned, or equal to i;
        // let U be the set of all indices east, aligned, or equal to j.
        std::set<int> L, U;
        L.insert(i);
        U.insert(j);
        ACAFlags af = perpAlignFlag(sepToAlignFlag(sf));
        for (int k = 0; k < numCols; k++) {
            if ( ( (*m_separationState)(k,i) & sf ) || ( (*m_alignmentState)(k,i) & af ) ) L.insert(k);
            if ( ( (*m_separationState)(j,k) & sf ) || ( (*m_alignmentState)(j,k) & af ) ) U.insert(k);
        }
        // Now record that everything in L is west of everything in U.
        // This is the transitive closure of the new separation.
        ACASepFlags nf = negateSepFlag(sf);
        for (std::set<int>::iterator it=L.begin(); it!=L.end(); ++it) {
            for (std::set<int>::iterator jt=U.begin(); jt!=U.end(); ++jt) {
                (*m_separationState)(*it,*jt) |= sf;
                (*m_separationState)(*jt,*it) |= nf;
            }
        }
        return;
    }
    default:
        return;
    }
}

void ACALayout::layoutWithCurrentConstraints(void)
{
    if (m_fdlayout) delete m_fdlayout;
    m_fdlayout =
#define OLDCOLA
#ifdef OLDCOLA
            new cola::ConstrainedFDLayout(m_rs,m_es,m_idealLength,
                                          m_preventOverlaps);
#else
            new cola::ConstrainedFDLayout(m_rs,m_es,m_idealLength,
                                          m_preventOverlaps,m_edgeLengths,
                                          m_doneTest,m_preIteration);
#endif
    m_fdlayout->setConstraints(m_ccs);
    m_fdlayout->run(true,true);
}

bool ACALayout::acaLoopOnce(void)
{
    OrderedAlignment *oa = chooseOA();
    if (oa) {
        // Add the new separated alignment constraints.
        m_ccs.push_back(oa->separation);
        m_ccs.push_back(oa->alignment);
        // Redo the layout, with the new constraints.
        layoutWithCurrentConstraints();
        // Update state tables.
        updateStateTables(oa);
        return true;
    } else {
        return false;
    }
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
        oa = chooseOA();
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
        oa = chooseOA();
    }
    // Redo the layout, with the new constraints.
    layoutWithCurrentConstraints();
}

void ACALayout::updateStateTables(OrderedAlignment *oa)
{
    ACASepFlags sf = oa->af==ACAHORIZ ? ACAEAST : ACASOUTH;
    recordAlignmentWithClosure(oa->left,oa->right,oa->af);
    recordSeparationWithClosure(oa->left,oa->right,sf);
}

OrderedAlignment *ACALayout::chooseOA(void)
{
    OrderedAlignment *oa = NULL;
    // Initialise minPenalty to a value exceeding the maximum penalty
    // any edge can actually be assigned (so effectively infinity).
    double minPenalty = PENALTY_BOUND;
    // Consider each edge for potential alignment.
    for (int j = 0; j < m_m; j++) {
        cola::Edge e = m_es.at(j);
        int src = e.first, tgt = e.second;
        // If already aligned, then skip this edge.
        int astate = (*m_alignmentState)(src,tgt);
        if (astate & (ACAHORIZ|ACAVERT)) continue;
        // Otherwise consider placing tgt in each of the compass directions from src.
        ACASepFlags sf;
        int sn = 1;
        while (sn < 16) {
            sf = (ACASepFlags) sn;
            // shift left to next cardinal direction
            // (do it now in case of 'continue's below)
            sn = sn << 1;
            // Check feasibility.
            if (badSeparation(j,sf)) continue;
            if (createsOverlap(src,tgt,sf)) continue;
            // Passed feasibility tests.
            // Now compute penalty.
            double p = 0;
            p += deflection(src,tgt,sf);
            if (m_addBendPointPenalty) p += bendPointPenalty(src,tgt,sf);
            if (m_postponeLeaves) p += leafPenalty(src,tgt);
            // Replace current ordered alignment if penalty is lower.
            if (p < minPenalty) {
                minPenalty = p;
                if (!oa) oa = new OrderedAlignment;
                oa->af = sepToAlignFlag(sf);
                oa->left  = sf==ACANORTH || sf==ACAWEST ? tgt : src;
                oa->right = oa->left == tgt ? src : tgt;
                EdgeOffset offset = getEdgeOffsetForCompassDirection(j,sf);
                oa->offsetLeft  = oa->left == src ? offset.first : offset.second;
                oa->offsetRight = oa->left == src ? offset.second : offset.first;
            }
        }
    }
    // Did we find an alignment?
    if (oa) {
        // If so, then complete the OrderedAlignment object.
        int l = oa->left, r = oa->right;
        vpsc::Rectangle *rl = m_rs.at(l), *rr = m_rs.at(r);
        // Determine dimensions.
        vpsc::Dim sepDim   = oa->af == ACAHORIZ ? vpsc::XDIM : vpsc::YDIM;
        vpsc::Dim alignDim = oa->af == ACAHORIZ ? vpsc::YDIM : vpsc::XDIM;
        // Create the separation constraint.
        double sep = oa->af == ACAHORIZ ?
                    (rl->width()+rr->width())/2.0 : (rl->height()+rr->height())/2.0;
        oa->separation = new cola::SeparationConstraint(sepDim,l,r,sep);
        // Create the alignment constraint.
        oa->alignment = new cola::AlignmentConstraint(alignDim);
        oa->alignment->addShape(l,oa->offsetLeft);
        oa->alignment->addShape(r,oa->offsetRight);
    }
    return oa;
}

// Say whether the proposed separation is a bad one.
bool ACALayout::badSeparation(int j, ACASepFlags sf)
{
    // If allowed separations have been set, then check whether sf is allowed.
    if (m_allowedSeps.size()>0) {
        ACASepFlags allowed = m_allowedSeps.at(j);
        if ( (sf&allowed) != sf ) return true; // it is a bad separation
    }
    // If we are /not/ doing aggressive ordering, then
    // first check if sf is ruled out for reversing existing order.
    cola::Edge e = m_es.at(j);
    int src = e.first, tgt = e.second;
    if (!m_aggressiveOrdering) {
        vpsc::Rectangle *rs = m_rs.at(src), *rt = m_rs.at(tgt);
        double dx = rt->getCentreX() - rs->getCentreX();
        double dy = rt->getCentreY() - rs->getCentreY();
        ACASepFlags currPos = vectorToSepFlag(dx,dy);
        bool conflict = propsedSepConflictsWithExistingPosition(sf,currPos);
        if (conflict) return true; // it is a bad separation
    }
    // Finally check if sf conflicts with the separation state table.
    ACASepFlags currConstraint = (ACASepFlags)(*m_separationState)(src,tgt);
    bool conflict = propsedSepConflictsWithExistingSepCo(sf,currConstraint);
    return conflict;
}

bool ACALayout::createsOverlap(int src, int tgt, ACASepFlags sf)
{
    // Determine which shape is low and which is high in the dimension of interest.
    int lowIndex = sf==ACANORTH || sf==ACAWEST ? tgt : src;
    int highIndex = (lowIndex==tgt ? src : tgt);
    ACAFlags af = sepToAlignFlag(sf);
    // Determine the coordinates in the dimension of interest.
    vpsc::Rectangle *rLow = m_rs.at(lowIndex), *rHigh = m_rs.at(highIndex);
    double lowCoord = af==ACAHORIZ ? rLow->getCentreX() : rLow->getCentreY();
    double highCoord = af==ACAHORIZ ? rHigh->getCentreX() : rHigh->getCentreY();
    // Let L and H be the low and high shapes respectively.
    // We consider each node U which is already aligned with either L or H.
    // Any such node must have lower coord than L if it is connected to L, and
    // higher coord than H if it is connected to H. If either of those conditions
    // fails, then we predict overlap.
    bool overlap = false;
    for (int j = 0; j < m_n; j++) {
        if (j==lowIndex || j==highIndex) continue;
        vpsc::Rectangle *r = m_rs.at(j);
        int lj = (*m_alignmentState)(lowIndex, j);
        int hj = (*m_alignmentState)(highIndex, j);
        if (lj&af || hj&af) {
            double z = af==ACAHORIZ ? r->getCentreX() : r->getCentreY();
            // low shape
            if ( lj&ACACONN && lowCoord < z ) {
                overlap = true; break;
            }
            // high shape
            if ( hj&ACACONN && z < highCoord ) {
                overlap = true; break;
            }
        }
    }
    return overlap;
}

/* Compute a score in [0.0, 1.0] measuring how far the "edge" in question E is
 * deflected from horizontal or vertical, depending on the passed alignment flag.
 * The "edge" E is the straight line from the centre of the src rectangle to that
 * of the tgt rectangle, as given by the src and tgt indices, and the vector of
 * Rectangles.
 *
 * If t is the angle that E makes with the positive x-axis, then the score we return
 * is sin^2(t) for horizontal alignments, and cos^2(t) for vertical.
 *
 * So smaller deflection scores mean edges that are closer to axis-aligned.
 */
double ACALayout::deflection(int src, int tgt, ACASepFlags sf)
{
    vpsc::Rectangle *s = m_rs.at(src), *t = m_rs.at(tgt);
    double sx=s->getCentreX(), sy=s->getCentreY(), tx=t->getCentreX(), ty=t->getCentreY();
    double dx = tx-sx, dy = ty-sy;
    double dx2 = dx*dx, dy2 = dy*dy;
    double l = dx2 + dy2;
    double dfl = sf==ACAWEST || sf==ACAEAST ? dy2/l : dx2/l;
    return dfl;
}

double ACALayout::bendPointPenalty(int src, int tgt, ACASepFlags sf)
{
    double penalty = BP_PENALTY;
    ACAFlags af = sepToAlignFlag(sf);
    ACAFlags op = af==ACAHORIZ ? ACAVERT : ACAHORIZ;
    std::set<int> deg2Nodes = m_useNonLeafDegree ? m_nldeg2Nodes : m_deg2Nodes;
    std::multimap<int,int> nbrs = m_useNonLeafDegree ? m_nlnbrs : m_nbrs;
    // First check whether src would be made into a bendpoint.
    if (deg2Nodes.count(src)!=0) {
        // Find neighbour j of src which is different from tgt, by iterating over nbrs of src.
        int j;
        std::pair< std::multimap<int,int>::iterator, std::multimap<int,int>::iterator > range;
        range = nbrs.equal_range(src);
        for (std::multimap<int,int>::iterator it=range.first; it!=range.second; ++it) {
            j = it->second;
            if (j != tgt) break;
        }
        // Now check if they are aligned in the opposite dimension.
        int as = (*m_alignmentState)(src,j);
        if (as & op) return penalty;
    }
    // Now check whether tgt would be made into a bendpoint.
    if (deg2Nodes.count(tgt)!=0) {
        // Find neighbour j of tgt which is different from src, by iterating over nbrs of tgt.
        int j;
        std::pair< std::multimap<int,int>::iterator, std::multimap<int,int>::iterator > range;
        range = nbrs.equal_range(tgt);
        for (std::multimap<int,int>::iterator it=range.first; it!=range.second; ++it) {
            j = it->second;
            if (j != src) break;
        }
        // Now check if they are aligned in the opposite dimension.
        int as = (*m_alignmentState)(tgt,j);
        if (as & op) return penalty;
    }
    return 0;
}

double ACALayout::leafPenalty(int src, int tgt)
{
    double penalty = LEAF_PENALTY;
    return m_leaves.count(src)!=0 || m_leaves.count(tgt)!=0 ? penalty : 0;
}

EdgeOffset ACALayout::getEdgeOffsetForCompassDirection(int j, ACASepFlags sf)
{
    EdgeOffset offset(0,0);
    // If offsets have been specified for this sep flag, then retrieve
    // the offset of the given index.
    std::map<ACASepFlags,EdgeOffsets>::iterator it = m_edgeOffsets.find(sf);
    if (it!=m_edgeOffsets.end()) {
        EdgeOffsets offsets = it->second;
        offset = offsets.at(j);
    }
    return offset;
}

} // namespace cola
