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

#include <algorithm>

#include "libcola/aca.h"

using namespace std;
using namespace vpsc;

namespace cola {

ACASepFlag negateSepFlag(ACASepFlag sf)
{
    unsigned short c = (unsigned short) sf;
    c += 16*c;
    c &= 60; // 00111100
    unsigned short b = c >> 2;
    ACASepFlag nf = (ACASepFlag) b;
    return nf;
}

ACAFlag sepToAlignFlag(ACASepFlag sf)
{
    return sf==ACANORTH || sf==ACASOUTH ? ACAVERT : ACAHORIZ;
}

ACAFlag perpAlignFlag(ACAFlag af)
{
    return af==ACAHORIZ ? ACAVERT : ACAHORIZ;
}

ACASepFlag vectorToSepFlag(double dx, double dy)
{
    int f = 0;
    f |= dx > 0 ? ACAEAST : dx < 0 ? ACAWEST : 0;
    f |= dy > 0 ? ACASOUTH : dy < 0 ? ACANORTH : 0;
    return (ACASepFlag) f;
}

bool propsedSepConflictsWithExistingPosition(ACASepFlag pro, ACASepFlag ex)
{
    // Proposed separation pro conflicts with existing position ex
    // if the union of their bits contains both north and south, or
    // both east and west.
    int u = pro | ex;
    return ( ( (u&5) == 5 ) || ( (u&10) == 10 ) );
}

bool propsedSepConflictsWithExistingSepCo(ACASepFlag pro, ACASepFlag ex)
{
    // Proposed separation pro conflicts with existing separation ex
    // if ex has any of the complementary bits of pro.
    int proComp = 15 - pro;
    return (proComp & ex);
}

bool sortEvents(const AlignedNodes::Event &lhs, const AlignedNodes::Event &rhs)
{
    return lhs.m_pos < rhs.m_pos;
}

AlignedNodes::AlignedNodes(Dim primaryDim, int nodeIndex, Rectangle *nodeRect)  :
    m_primaryDim(primaryDim)
{
    m_secondaryDim = primaryDim==vpsc::XDIM ? vpsc::YDIM : vpsc::XDIM;
    m_nodeIndices.push_back(nodeIndex);
    // Make new copy of rectangle, so existing one will not be modified if and when
    // this AlignedNodes object gets shifted.
    vpsc::Rectangle *S = new vpsc::Rectangle(*nodeRect);
    m_nodeRects.push_back(S);
}

AlignedNodes::~AlignedNodes()
{
    for (unsigned i = 0; i < m_nodeRects.size(); i++) {
        delete m_nodeRects[i];
    }
}

void AlignedNodes::removeNode(int index)
{
    vpsc::Rectangles::iterator rit = m_nodeRects.begin();
    for (std::vector<int>::iterator it=m_nodeIndices.begin(); it!=m_nodeIndices.end(); ++it) {
        if (*it==index) {
            m_nodeIndices.erase(it);
            m_nodeRects.erase(rit);
            break;
        }
        rit++;
    }
}

void AlignedNodes::removeAuxiliaryNodes(int minIndex)
{
    vpsc::Rectangles::iterator rit = m_nodeRects.begin();
    for (std::vector<int>::iterator it=m_nodeIndices.begin(); it!=m_nodeIndices.end(); ++it) {
        if (*it>=minIndex) {
            m_nodeIndices.erase(it);
            m_nodeRects.erase(rit);
        }
        rit++;
    }
}

std::string AlignedNodes::toString(string pad)
{
    std::string s = "";
    char buffer [1000];
    // Nodes
    s += pad+"Nodes:\n";
    for (unsigned i = 0; i < m_nodeIndices.size(); i++) {
        int index = m_nodeIndices.at(i);
        sprintf(buffer,"%d:",index);
        s += pad+std::string(buffer);
        if (m_primaryDim>=0 && m_primaryDim<2 && m_secondaryDim>=0 && m_secondaryDim<2) {
            vpsc::Rectangle *R = m_nodeRects.at(i);
            double u = R->getMinD(m_primaryDim), U = R->getMaxD(m_primaryDim);
            double v = R->getMinD(m_secondaryDim), V = R->getMaxD(m_secondaryDim);
            sprintf(buffer," [%.2f,%.2f] x [%.2f,%.2f]\n",u,U,v,V);
            s += std::string(buffer);
        } else {
            sprintf(buffer," ERROR: bad dimension numbers %d, %d\n",m_primaryDim,m_secondaryDim);
            s += std::string(buffer);
        }
    }
    // Edges
    s += pad+"Edges:\n";
    for (unsigned i = 0; i < m_edgeIndices.size(); i++) {
        int index = m_edgeIndices.at(i);
        double c=m_edgeConstCoords.at(i), l=m_edgeLowerBounds.at(i), u=m_edgeUpperBounds.at(i);
        sprintf(buffer,"%d: at %.2f from %.2f to %.2f\n",index,c,l,u);
        s += pad+std::string(buffer);
    }
    return s;
}

AlignedNodes *AlignedNodes::combineWithEdge(const AlignedNodes &other, int node1, double offset1, int node2, double offset2, int edge)
{
    // Determine the local coordinate of port 1, in this AlignedNodes object.
    unsigned i = 0;
    while (i < m_nodeIndices.size()) {
        if (m_nodeIndices.at(i)==node1) break;
        i++;
    }
    vpsc::Rectangle *R = m_nodeRects.at(i);
    double v1 = R->getCentreD(m_secondaryDim);
    v1 += offset1;
    // Determine the local coordinate of port 2, in the other AlignedNodes object.
    unsigned j = 0;
    while (j < other.m_nodeIndices.size()) {
        if (other.m_nodeIndices.at(j)==node2) break;
        j++;
    }
    vpsc::Rectangle *S = other.m_nodeRects.at(j);
    double v2 = S->getCentreD(m_secondaryDim);
    v2 += offset2;
    // Combine with a shifted copy of other.
    AlignedNodes *a = this->combine(other,v1-v2);
    if (edge >= 0) {
        // Add the new edge.
        double u1 = R->getMinD(m_primaryDim), U1 = R->getMaxD(m_primaryDim);
        double u2 = S->getMinD(m_primaryDim), U2 = S->getMaxD(m_primaryDim);
        double l = U1 < u2 ? U1 : U2;
        double u = U1 < u2 ? u2 : u1;
        a->addEdge(edge,v1,l,u);
    }
    // Return pointer to the new object.
    return a;
}

AlignedNodes *AlignedNodes::combineWithPorts(const AlignedNodes &other, int node1, double offset1, int node2, double offset2)
{
    int edgeIndex = -1;
    return combineWithEdge(other,node1,offset1,node2,offset2,edgeIndex);
}

void AlignedNodes::addEdge(int index, double c, double l, double u)
{
    m_edgeIndices.push_back(index);
    m_edgeConstCoords.push_back(c);
    m_edgeLowerBounds.push_back(l);
    m_edgeUpperBounds.push_back(u);
}

AlignedNodes *AlignedNodes::combine(const AlignedNodes &other, double shift)
{
    AlignedNodes *S = new AlignedNodes(m_primaryDim);
    // node indices
    for (unsigned i = 0; i < m_nodeIndices.size(); i++) {
        S->m_nodeIndices.push_back(m_nodeIndices.at(i));
    }
    for (unsigned i = 0; i < other.m_nodeIndices.size(); i++) {
        S->m_nodeIndices.push_back(other.m_nodeIndices.at(i));
    }
    // node rects
    for (unsigned i = 0; i < m_nodeRects.size(); i++) {
        vpsc::Rectangle *R = new vpsc::Rectangle(*m_nodeRects.at(i));
        S->m_nodeRects.push_back(R);
    }
    for (unsigned i = 0; i < other.m_nodeRects.size(); i++) {
        vpsc::Rectangle *R = new vpsc::Rectangle(*other.m_nodeRects.at(i));
        double z = R->getCentreD(m_secondaryDim);
        R->moveCentreD(m_secondaryDim,z+shift);
        S->m_nodeRects.push_back(R);
    }
    // edge indices
    for (unsigned i = 0; i < m_edgeIndices.size(); i++) {
        S->m_edgeIndices.push_back(m_edgeIndices.at(i));
    }
    for (unsigned i = 0; i < other.m_edgeIndices.size(); i++) {
        S->m_edgeIndices.push_back(other.m_edgeIndices.at(i));
    }
    // edge const coords
    for (unsigned i = 0; i < m_edgeConstCoords.size(); i++) {
        S->m_edgeConstCoords.push_back(m_edgeConstCoords.at(i));
    }
    for (unsigned i = 0; i < other.m_edgeConstCoords.size(); i++) {
        S->m_edgeConstCoords.push_back(other.m_edgeConstCoords.at(i)+shift);
    }
    // edge lower bounds
    for (unsigned i = 0; i < m_edgeLowerBounds.size(); i++) {
        S->m_edgeLowerBounds.push_back(m_edgeLowerBounds.at(i));
    }
    for (unsigned i = 0; i < other.m_edgeLowerBounds.size(); i++) {
        S->m_edgeLowerBounds.push_back(other.m_edgeLowerBounds.at(i));
    }
    // edge upper bounds
    for (unsigned i = 0; i < m_edgeUpperBounds.size(); i++) {
        S->m_edgeUpperBounds.push_back(m_edgeUpperBounds.at(i));
    }
    for (unsigned i = 0; i < other.m_edgeUpperBounds.size(); i++) {
        S->m_edgeUpperBounds.push_back(other.m_edgeUpperBounds.at(i));
    }
    return S;
}

bool AlignedNodes::thereAreOverlaps(void)
{
    // If there are no edges at all, then there are no overlaps.
    if (m_edgeIndices.size() == 0) return false;
    // Will need unique indices for edges, so compute strict upper
    // bound on existing node indices.
    int N = 0;
    // Make events for nodes opening and closing, and for edges.
    std::vector<Event> events;
    for (unsigned i = 0; i < m_nodeIndices.size(); i++) {
        int index = m_nodeIndices.at(i);
        if (index >= N) N = index + 1;
        vpsc::Rectangle *R = m_nodeRects.at(i);
        Event eo(OPEN), ec(CLOSE);
        eo.m_index = index;
        ec.m_index = index;
        double u = R->getMinD(m_primaryDim),   U = R->getMaxD(m_primaryDim);
        double v = R->getMinD(m_secondaryDim), V = R->getMaxD(m_secondaryDim);
        eo.m_pos = v; eo.m_lowerBound = u; eo.m_upperBound = U;
        ec.m_pos = V; ec.m_lowerBound = u; ec.m_upperBound = U;
        events.push_back(eo);
        events.push_back(ec);
    }
    for (unsigned i = 0; i < m_edgeIndices.size(); i++) {
        Event e(LOCAL);
        e.m_index = m_edgeIndices.at(i);
        e.m_pos = m_edgeConstCoords.at(i);
        e.m_lowerBound = m_edgeLowerBounds.at(i);
        e.m_upperBound = m_edgeUpperBounds.at(i);
        events.push_back(e);
    }
    // Sort by m_pos
    std::sort(events.begin(), events.end(), sortEvents);
    // Scan through list looking for overlaps.
    std::vector<Event> openIntervals;
    unsigned i = 0;
    while (i < events.size()) {
        // Get open intervals for next constant coord.
        double v = events.at(i).m_pos;
        std::set<int> closedEvents;
        while (i < events.size() && events.at(i).m_pos == v) {
            Event e = events.at(i);
            if (e.m_type==CLOSE) {
                closedEvents.insert(e.m_index);
            } else {
                openIntervals.push_back(e);
            }
            i++;
        }
        // Now create a vector of the open and close events for the intervals.
        std::vector<Event> intervalEvents;
        for (unsigned j = 0; j < openIntervals.size(); j++) {
            Event e = openIntervals.at(j);
            Event o(OPEN), c(CLOSE);
            int index = e.m_index + (e.m_type==LOCAL ? N : 0);
            o.m_index = c.m_index = index;
            o.m_pos = e.m_lowerBound;
            c.m_pos = e.m_upperBound;
            intervalEvents.push_back(o);
            intervalEvents.push_back(c);
        }
        // Sort
        std::sort(intervalEvents.begin(), intervalEvents.end(), sortEvents);
        // Check for overlapping intervals.
        // No two intervals may be open at the same time, but it is okay if one opens at the same
        // coordinate at which another closes.
        int openIndex = -1;
        for (unsigned k = 0; k < intervalEvents.size(); k++) {
            // If find an overlap, then immediately return true.
            Event e = intervalEvents.at(k);
            if (openIndex < 0) {
                // If there is currently no open interval, then open one.
                openIndex = e.m_index;
            } else {
                // There is currently an open interval.
                // Check the index of the next point.
                int nextIndex = e.m_index;
                if (nextIndex == openIndex) {
                    // If it equals the index of the open interval, then that interval
                    // is now closed.
                    openIndex = -1;
                } else {
                    // Otherwise an interval is open, and we've just encountered the opening end
                    // of a different interval.
                    // This is okay only if the very next event is the closing end for
                    // the first open interval, since at most two intervals can share any
                    // endpoint without there being an overlap.
                    COLA_ASSERT(k+1 < intervalEvents.size());
                    Event f = intervalEvents.at(k+1);
                    if (f.m_index != openIndex) {
                        // There is an overlap.
                        return true;
                    } else {
                        // The first open interval was closed, and the second one is
                        // now to be noted by the openIndex.
                        openIndex = nextIndex;
                        k++;
                    }
                }
            }
        }
        // Remove closed intervals, if any.
        if (closedEvents.size()>0) {
            std::vector<Event> tmp;
            for (unsigned j = 0; j < openIntervals.size(); j++) {
                Event e = openIntervals.at(j);
                if (e.m_type==OPEN && closedEvents.count(e.m_index)!=0) continue;
                tmp.push_back(e);
            }
            openIntervals = tmp;
        }
    }
    return false;
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
      m_aggressiveOrdering(false),
      m_overlapPrevention(ACAOPCENTREALIGN),
      m_alignmentSetsTrackLayout(false),
      m_fdlayout(NULL)
{
    computeDegrees();
    generateVPSCConstraints();
    initStateTables();
    initAlignmentSets(vpsc::XDIM);
    initAlignmentSets(vpsc::YDIM);
}

ACALayout::~ACALayout(void)
{
    delete m_alignmentState;
    delete m_separationState;
    delete m_fdlayout;
    for (unsigned i = 0; i < m_ordAligns.size(); i++) {
        delete m_ordAligns.at(i);
    }
    for (unsigned i = 0; i < m_xvs.size(); i++) {
        delete m_xvs.at(i);
    }
    for (unsigned i = 0; i < m_yvs.size(); i++) {
        delete m_yvs.at(i);
    }
    for (unsigned i = 0; i < m_hSets.size(); i++) {
        delete m_hSets.at(i);
    }
    for (unsigned i = 0; i < m_vSets.size(); i++) {
        delete m_vSets.at(i);
    }
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

void ACALayout::overlapPrevention(ACAOverlapPrevention p)
{
    m_overlapPrevention = p;
}

void ACALayout::setAlignmentOffsetsForCompassDirection(ACASepFlag sf, EdgeOffsets offsets)
{
    COLA_ASSERT(offsets.size()==(size_t)m_m); // There should be one offset for each edge.
    m_edgeOffsets.insert( std::pair<ACASepFlag,EdgeOffsets>(sf,offsets) );
}

void ACALayout::setAllowedSeparations(ACASepFlags seps)
{
    COLA_ASSERT(seps.size()==(size_t)m_m); // There should be one flag for each edge.
    m_allowedSeps = seps;
}

void ACALayout::setAllowedSeparations(ACASepFlagsStruct seps)
{
    setAllowedSeparations(seps.sepFlags);
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

std::string ACALayout::writeAlignmentSets(void)
{
    char buf [100];
    std::string s = "";
    // Horizontal sets
    s += "Horizontal alignment sets:";
    sprintf(buf," (%ld)\n",m_hSets.size());
    s += std::string(buf);
    for (std::map<int,AlignedNodes*>::iterator it=m_hSets.begin(); it!=m_hSets.end(); ++it)
    {
        int index = it->first;
        AlignedNodes *a = it->second;
        sprintf(buf,"  %d:\n",index);
        s += std::string(buf);
        s += a->toString("    ");
    }
    s += "\n";
    // Vertical sets
    s += "Vertical alignment sets";
    sprintf(buf," (%ld)\n",m_vSets.size());
    s += std::string(buf);
    for (std::map<int,AlignedNodes*>::iterator it=m_vSets.begin(); it!=m_vSets.end(); ++it)
    {
        int index = it->first;
        AlignedNodes *a = it->second;
        sprintf(buf,"  %d:\n",index);
        s += std::string(buf);
        s += a->toString("    ");
    }
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

void ACALayout::initAlignmentSets(Dim dim)
{
    // Get the sets and constraints for the named primary dimension.
    std::map<int,AlignedNodes*> &sets = dim==vpsc::XDIM ? m_hSets : m_vSets;
    vpsc::Constraints &cs = dim==vpsc::XDIM ? m_yEqCs : m_xEqCs;
    // Add an alignment for each constraint.
    for (vpsc::Constraints::iterator cit=cs.begin(); cit!=cs.end(); ++cit) {
        vpsc::Constraint c = **cit;
        int l = c.left->id, r = c.right->id;
        double g = c.gap;
        // Get existing set for each variable, or create a new one if it doesn't exist yet.
        AlignedNodes *ls = getAlignmentSet(dim,l);
        AlignedNodes *rs = getAlignmentSet(dim,r);
        // Combine
        AlignedNodes *a = ls->combineWithPorts(*rs,l,0,r,g);
        // Update map.
        for (unsigned i = 0; i < a->m_nodeIndices.size(); i++) {
            int index = a->m_nodeIndices.at(i);
            std::map<int,AlignedNodes*>::iterator it = sets.find(index);
            if (it!=sets.end()) sets.erase(it);
            sets.insert(std::pair<int,AlignedNodes*>(index,a));
        }
        delete ls;
        delete rs;
    }
    // Clean up.
    // Remove auxiliary indices from the map.
    int N = dim==vpsc::XDIM ? m_numExtraYVars : m_numExtraXVars;
    for (int k = 0; k < N; k++) {
        int v = m_n + k;
        sets.erase(sets.find(v));
    }
    // Remove dummy rectangles for auxiliary variables.
    for (std::map<int,AlignedNodes*>::iterator i=sets.begin(); i!=sets.end(); ++i) {
        i->second->removeAuxiliaryNodes(m_n);
    }
}

AlignedNodes *ACALayout::getAlignmentSet(Dim dim, int i)
{
    std::map<int,AlignedNodes*> &sets = dim==vpsc::XDIM ? m_hSets : m_vSets;
    std::map<int,AlignedNodes*>::iterator it = sets.find(i);
    AlignedNodes *a;
    if (it==sets.end()) {
        // create new set
        // If i is the index of an auxiliary variable, then just make a dummy rectangle for it.
        vpsc::Rectangle *R = i < m_n ? m_rs.at(i) : new vpsc::Rectangle(-5,5,-5,5);
        a = new AlignedNodes(dim,i,R);
        // Add it to the map so that we get expected behaviour elsewhere, e.g.
        // when we try to delete existing value from map.
        sets.insert(std::pair<int,AlignedNodes*>(i,a));
        // The AlignedNodes constructor makes a copy of the passed rectangle.
        // Therefore, if we just created the rectangle R here, for an auxiliary
        // variable, then we should delete it now.
        if (i>=m_n) delete R;
    } else {
        // get existing set
        a = it->second;
    }
    return a;
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
            ACASepFlag sf = gap > 0 ? ACAEAST : ACAWEST;
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
            ACASepFlag sf = gap > 0 ? ACASOUTH : ACANORTH;
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

void ACALayout::recordAlignmentWithClosure(int i, int j, ACAFlag af, int numCols)
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
    ACASepFlag sf = af == ACAVERT ? ACAEAST : ACASOUTH;
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
    ACASepFlag nf = negateSepFlag(sf);
    for (std::set<int>::iterator it=L.begin(); it!=L.end(); ++it) {
        for (std::set<int>::iterator jt=U.begin(); jt!=U.end(); ++jt) {
            (*m_separationState)(*it,*jt) |= sf;
            (*m_separationState)(*jt,*it) |= nf;
        }
    }
}

void ACALayout::recordSeparationWithClosure(int i, int j, ACASepFlag sf, int numCols)
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
        ACAFlag af = perpAlignFlag(sepToAlignFlag(sf));
        for (int k = 0; k < numCols; k++) {
            if ( ( (*m_separationState)(k,i) & sf ) || ( (*m_alignmentState)(k,i) & af ) ) L.insert(k);
            if ( ( (*m_separationState)(j,k) & sf ) || ( (*m_alignmentState)(j,k) & af ) ) U.insert(k);
        }
        // Now record that everything in L is west of everything in U.
        // This is the transitive closure of the new separation.
        ACASepFlag nf = negateSepFlag(sf);
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
    m_alignmentSetsTrackLayout = true;
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
    m_alignmentSetsTrackLayout = true;
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
    m_alignmentSetsTrackLayout = false;
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

void ACALayout::updateAlignmentSetRects(vpsc::Dim dim)
{
    std::map<int,AlignedNodes*> &sets = dim==vpsc::XDIM ? m_hSets : m_vSets;
    std::set<int> doneIndices;
    for (std::map<int,AlignedNodes*>::iterator it=sets.begin(); it!=sets.end(); ++it) {
        int index = it->first;
        if (doneIndices.count(index)!=0) continue;
        AlignedNodes *a = it->second;
        for (unsigned i = 0; i < a->m_nodeIndices.size(); i++) {
            int index = a->m_nodeIndices.at(i);
            vpsc::Rectangle *R = m_rs.at(index);
            vpsc::Rectangle *S = a->m_nodeRects.at(i);
            S->moveCentreX(R->getCentreX());
            S->moveCentreY(R->getCentreY());
            doneIndices.insert(index);
        }
    }
}

void ACALayout::updateStateTables(OrderedAlignment *oa)
{
    int l = oa->left, r = oa->right;
    ACASepFlag sf = oa->af==ACAHORIZ ? ACAEAST : ACASOUTH;
    // State tables:
    recordAlignmentWithClosure(l,r,oa->af);
    recordSeparationWithClosure(l,r,sf);
    // Alignment sets:
    if (m_alignmentSetsTrackLayout) {
        updateAlignmentSetRects(vpsc::XDIM);
        updateAlignmentSetRects(vpsc::YDIM);
    }
    vpsc::Dim dim = oa->af==ACAHORIZ ? vpsc::XDIM : vpsc::YDIM;
    AlignedNodes *ls = getAlignmentSet(dim,l);
    AlignedNodes *rs = getAlignmentSet(dim,r);
    AlignedNodes *a = ls->combineWithEdge(*rs,l,oa->offsetLeft,
                                             r,oa->offsetRight,
                                         oa->edgeIndex);
    std::map<int,AlignedNodes*> &sets = dim==vpsc::XDIM ? m_hSets : m_vSets;
    for (unsigned i = 0; i < a->m_nodeIndices.size(); i++) {
        int index = a->m_nodeIndices.at(i);
        std::map<int,AlignedNodes*>::iterator it = sets.find(index);
        if (it!=sets.end()) {
            sets.erase(it);
            sets.insert(std::pair<int,AlignedNodes*>(index,a));
        }
    }
    delete ls;
    delete rs;
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
        ACASepFlag sf;
        int sn = 1;
        while (sn < 16) {
            sf = (ACASepFlag) sn;
            // shift left to next cardinal direction
            // (do it now in case of 'continue's below)
            sn = sn << 1;
            // Check feasibility.
            // First check whether it is a bad separation.
            bool bs = badSeparation(j,sf);
            if (bs) continue;
            // Next check whether the alignment would create
            // edge-edge or edge-node overlaps.
            switch (m_overlapPrevention) {
            case ACAOPCENTREALIGN:
                if (createsOverlap(src,tgt,sf)) continue;
                break;
            case ACAOPWITHOFFSETS: {
                bool co = createsOverlap2(j,sf);
                if (co) continue;
                break;
            }
            default:
                break;
            }
            // Passed feasibility tests.
            // Now compute penalty.
            double p = 0;
            // Deflection:
            p += deflection(src,tgt,sf);
            // Bend points:
            if (m_addBendPointPenalty) p += bendPointPenalty(src,tgt,sf);
            // Leaves:
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
                oa->edgeIndex = j;
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
        // Record oa.
        m_ordAligns.push_back(oa);
    }
    return oa;
}

// Say whether the proposed separation is a bad one.
bool ACALayout::badSeparation(int j, ACASepFlag sf)
{
    // If allowed separations have been set, then check whether sf is allowed.
    if (m_allowedSeps.size()>0) {
        ACASepFlag allowed = m_allowedSeps.at(j);
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
        ACASepFlag currPos = vectorToSepFlag(dx,dy);
        bool conflict = propsedSepConflictsWithExistingPosition(sf,currPos);
        if (conflict) return true; // it is a bad separation
    }
    // Finally check if sf conflicts with the separation state table.
    ACASepFlag currConstraint = (ACASepFlag)(*m_separationState)(src,tgt);
    bool conflict = propsedSepConflictsWithExistingSepCo(sf,currConstraint);
    return conflict;
}

bool ACALayout::createsOverlap(int src, int tgt, ACASepFlag sf)
{
    // Determine which shape is low and which is high in the dimension of interest.
    int lowIndex = sf==ACANORTH || sf==ACAWEST ? tgt : src;
    int highIndex = (lowIndex==tgt ? src : tgt);
    ACAFlag af = sepToAlignFlag(sf);
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

bool ACALayout::createsOverlap2(int j, ACASepFlag sf)
{
    cola::Edge e = m_es.at(j);
    int src = e.first, tgt = e.second;
    EdgeOffset offset = getEdgeOffsetForCompassDirection(j,sf);
    double srcOffset = offset.first, tgtOffset = offset.second;
    vpsc::Dim dim = sepToAlignFlag(sf) == ACAHORIZ ? vpsc::XDIM : vpsc::YDIM;
    AlignedNodes *srcSet = getAlignmentSet(dim,src);
    AlignedNodes *tgtSet = getAlignmentSet(dim,tgt);
    AlignedNodes *a = srcSet->combineWithPorts(*tgtSet,src,srcOffset,tgt,tgtOffset);
    bool ans = a->thereAreOverlaps();
    delete a;
    return ans;
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
double ACALayout::deflection(int src, int tgt, ACASepFlag sf)
{
    vpsc::Rectangle *s = m_rs.at(src), *t = m_rs.at(tgt);
    double sx=s->getCentreX(), sy=s->getCentreY(), tx=t->getCentreX(), ty=t->getCentreY();
    double dx = tx-sx, dy = ty-sy;
    double dx2 = dx*dx, dy2 = dy*dy;
    double l = dx2 + dy2;
    double dfl = sf==ACAWEST || sf==ACAEAST ? dy2/l : dx2/l;
    return dfl;
}

double ACALayout::bendPointPenalty(int src, int tgt, ACASepFlag sf)
{
    double penalty = BP_PENALTY;
    ACAFlag af = sepToAlignFlag(sf);
    ACAFlag op = af==ACAHORIZ ? ACAVERT : ACAHORIZ;
    std::set<int> deg2Nodes = m_useNonLeafDegree ? m_nldeg2Nodes : m_deg2Nodes;
    std::multimap<int,int> nbrs = m_useNonLeafDegree ? m_nlnbrs : m_nbrs;
    // First check whether src would be made into a bendpoint.
    if (deg2Nodes.count(src)!=0) {
        // Find neighbour j of src which is different from tgt, by iterating over nbrs of src.
        int j = 0;
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
        int j = 0;
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

EdgeOffset ACALayout::getEdgeOffsetForCompassDirection(int j, ACASepFlag sf)
{
    EdgeOffset offset(0,0);
    // If offsets have been specified for this sep flag, then retrieve
    // the offset of the given index.
    std::map<ACASepFlag,EdgeOffsets>::iterator it = m_edgeOffsets.find(sf);
    if (it!=m_edgeOffsets.end()) {
        EdgeOffsets offsets = it->second;
        offset = offsets.at(j);
    }
    return offset;
}

} // namespace cola
