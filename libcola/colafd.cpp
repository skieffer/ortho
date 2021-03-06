/*
 * vim: ts=4 sw=4 et tw=0 wm=0
 *
 * libcola - A library providing force-directed network layout using the 
 *           stress-majorization method subject to separation constraints.
 *
 * Copyright (C) 2006-2010  Monash University
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
 *
*/

#include <QtGui>

#include <vector>
#include <cmath>
#include <limits>

#include "libvpsc/solve_VPSC.h"
#include "libvpsc/variable.h"
#include "libvpsc/constraint.h"
#include "libvpsc/rectangle.h"

#include "libcola/commondefs.h"
#include "libcola/cola.h"
#include "libcola/shortest_paths.h"
#include "libcola/straightener.h"
#include "libcola/cc_clustercontainmentconstraints.h"
#include "libcola/cc_nonoverlapconstraints.h"

#ifdef MAKEFEASIBLE_DEBUG
  #include "libcola/output_svg.h"
#endif

// Needs to come last since it will include windows.h on WIN32 and
// may mess up C++ std library include on GCC 4.4
#include "libcola/cola_log.h"

using namespace std;
using vpsc::XDIM;
using vpsc::YDIM;

namespace cola {

template <class T>
void delete_vector(vector<T*> &v) {
    for_each(v.begin(),v.end(),delete_object());
    v.clear();
}
Resizes PreIteration::__resizesNotUsed;
Locks PreIteration::__locksNotUsed;
inline double dotProd(valarray<double> x, valarray<double> y) {
    COLA_ASSERT(x.size()==y.size());
    double dp=0;
    for(unsigned i=0;i<x.size();i++) {
        dp+=x[i]*y[i]; 
    }
    return dp;
}
template <typename T>
void dumpSquareMatrix(unsigned n, T** L) {
    printf("Matrix %dX%d\n{",n,n);
    for(unsigned i=0;i<n;i++) {
        printf("{");
        for(unsigned j=0;j<n;j++) {
            std::cout<<L[i][j];
            char c=j==n-1?'}':',';
            printf("%c",c);
        }
        char c=i==n-1?'}':',';
        printf("%c\n",c);
    }
}

ConstrainedFDLayout::ConstrainedFDLayout(const vpsc::Rectangles& rs,
        const std::vector< Edge >& es, const double idealLength,
        const bool preventOverlaps,
        const bool snapTo, const double snapDistance,
        const double* eLengths,
        TestConvergence *doneTest, PreIteration* preIteration) 
    : n(rs.size()),
      X(valarray<double>(n)),
      Y(valarray<double>(n)),
      done(doneTest),
      using_default_done(false),
      preIteration(preIteration),
      topologyAddon(new TopologyAddonInterface()),
      rungekutta(true),
      desiredPositions(NULL),
      clusterHierarchy(NULL),
      rectClusterBuffer(0),
      m_idealEdgeLength(idealLength),
      m_generateNonOverlapConstraints(preventOverlaps),
      m_solver(NULL),
      m_tentative_constraint_threshold(0.1), // What value should it be?
      m_constraintToReject(NULL),
      m_constraintToRejectForUnsat(NULL),
      m_addSnapStress(snapTo),
      m_addGridSnapStress(false),
      m_snapGridX(100.0),
      m_snapGridY(100.0),
      m_snap_distance(snapDistance),
      m_snap_strength(20.0),
      // snap stress functions:
      //  1: smooth M-stress
      //  2: piecewise linear M-stress
      //  3: dual quadratic
      //  4: inverted quadratic (has long-range forces)
      //  5: quartic
      //  6: smooth V-stress
      //  7: linear V-stress
      //  8: quadraic U-stress
      m_snapStressFunction(8),
      m_debugLineNo(0),
      m_final_stress(0),
      m_addEdgeNodeRepulsion(true)
{
    if (done == NULL)
    {
        done = new TestConvergence();
        using_default_done = true;
    }

    //FILELog::ReportingLevel() = logDEBUG1;
    FILELog::ReportingLevel() = logERROR;
    boundingBoxes = rs;
    done->reset();
    unsigned i=0;
    for(vpsc::Rectangles::const_iterator ri=rs.begin();ri!=rs.end();++ri,++i) {
        X[i]=(*ri)->getCentreX();
        Y[i]=(*ri)->getCentreY();
        FILE_LOG(logDEBUG) << *ri;
    }
    D=new double*[n];
    G=new unsigned short*[n];
    for(unsigned i=0;i<n;i++) {
        D[i]=new double[n];
        G[i]=new unsigned short[n];
    }

    edges = es;

    if(eLengths == NULL) {
        computePathLengths(es,NULL);
    } else {
        valarray<double> eLengthsArray(eLengths,es.size());
        computePathLengths(es,&eLengthsArray);
    }

    switch (m_snapStressFunction) {
    case 1:
        // Smooth M-stress
        m_snapStressAlpha = 1.0;
        m_snapStressBeta = 5.0;
        m_snapStressGamma = 625.0;
        m_snapStressSigma = 25.0;
        m_snapStressRho = 1;
        break;
    case 2:
        // Linear M-stress
        m_snapStressAlpha = 1.0;
        m_snapStressBeta = 5.0;
        m_snapStressGamma = 625.0;
        m_snapStressSigma = 25.0;
        m_snapStressRho = 1;
        break;
    case 3:
        // Dual-quadratic
        m_snapStressAlpha = 0.01;
        m_snapStressBeta = 1.0;
        m_snapStressGamma = 1.0;
        m_snapStressSigma = 50.0;
        m_snapStressRho = 1;
        break;
    case 4:
        // Inverted Quadratic
        m_snapStressAlpha = 1.0;
        m_snapStressBeta = 5.0;
        m_snapStressGamma = 625.0;
        m_snapStressSigma = 25.0;
        m_snapStressRho = 1;
        break;
    case 5:
        // Quartic
        m_snapStressAlpha = 1.0;
        m_snapStressBeta = 5.0;
        m_snapStressGamma = 625.0;
        m_snapStressSigma = 25.0;
        m_snapStressRho = 1;
        break;
    case 6:
        // Smooth V-stress
        m_snapStressAlpha = 1.0;
        m_snapStressBeta = 100.0;
        m_snapStressGamma = 625.0;
        m_snapStressSigma = 50.0;
        m_snapStressRho = 1;
        break;
    case 7:
        // Linear V-stress
        m_snapStressAlpha = 1.0;
        m_snapStressBeta = 100.0;
        m_snapStressGamma = 625.0;
        m_snapStressSigma = 50.0;
        m_snapStressRho = 1;
        break;
    case 8:
        // Quadratic U-stress
        //m_snapStressBeta = 60.0;
        //m_snapStressSigma = 10.0;
        m_snapStressSigma = m_snap_distance;
        //m_snapStressBeta = 2*m_snapStressSigma;
        m_snapStressBeta = m_snap_strength;
    }
}

vpsc::Constraint *ConstrainedFDLayout::getMaxAbsLMConstraint()
{
    return m_solver->max_abs_lm;
}

void dijkstra(const unsigned s, const unsigned n, double* d, 
        const vector<Edge>& es, const double* eLengths)
{
    const valarray<double> eLengthsArray(eLengths, es.size());
    shortest_paths::dijkstra(s,n,d,es,&eLengthsArray);
}

/*
 * Sets up the D and G matrices.  D is the required euclidean distances
 * between pairs of nodes based on the shortest paths between them (using
 * m_idealEdgeLength*eLengths[edge] as the edge length, if eLengths array 
 * is provided otherwise just m_idealEdgeLength).  G is a matrix of 
 * unsigned ints such that G[u][v]=
 *   0 if there are no forces required between u and v 
 *     (for example, if u and v are in unconnected components)
 *   1 if attractive forces are required between u and v
 *     (i.e. if u and v are immediately connected by an edge and there is
 *      no topology route between u and v (for which an attractive force
 *      is computed elsewhere))
 *   2 if no attractive force is required between u and v but there is
 *     a connected path between them.
 */
void ConstrainedFDLayout::computePathLengths(
        const vector<Edge>& es,
        const std::valarray<double>* eLengths) 
{
    /*
    //DEBUG
    QString msg = "============================================\n";
    msg += "ConstrainedFDLayout computing path lengths.\n";
    msg += "Edges:\n";
    for (int i = 0; i < es.size(); i++) {
        Edge e = es.at(i);
        int u = e.first, v = e.second;
        double L = (*eLengths)[i];
        //qDebug() << u << v << L;
        msg += QString("%1").arg(u);
        msg += QString(" %1").arg(v);
        msg += QString(" %1\n").arg(L);
    }
    msg += "============================================\n";
    qDebug() << msg;
    //ENDDEBUG
    */
    shortest_paths::johnsons(n,D,es,eLengths);
    //dumpSquareMatrix<double>(n,D);
    for(unsigned i=0;i<n;i++) {
        for(unsigned j=0;j<n;j++) {
            if(i==j) continue;
            double& d=D[i][j];
            unsigned short& p=G[i][j];
            p=2;
            if(d==DBL_MAX) {
                // i and j are in disconnected subgraphs
                p=0;
            } else if(!eLengths) {
                D[i][j]*=m_idealEdgeLength;
            }
        }
    }
    for(vector<Edge>::const_iterator e=es.begin();e!=es.end();++e) {
        unsigned u=e->first, v=e->second; 
        G[u][v]=G[v][u]=1;
    }
    topologyAddon->computePathLengths(G);
    //dumpSquareMatrix<short>(n,G);
}

typedef valarray<double> Position;
void getPosition(Position& X, Position& Y, Position& pos) {
    unsigned n=X.size();
    COLA_ASSERT(Y.size()==n);
    COLA_ASSERT(pos.size()==2*n);
    for(unsigned i=0;i<n;++i) {
        pos[i]=X[i];
        pos[i+n]=Y[i];
    }
}
/*
 * moves all rectangles to the desired position while respecting
 * constraints.
 * @param pos target positions of both axes
 */
void ConstrainedFDLayout::setPosition(Position& pos) {
    COLA_ASSERT(Y.size()==X.size());
    COLA_ASSERT(pos.size()==2*X.size());
    moveTo(vpsc::HORIZONTAL,pos);
    moveTo(vpsc::VERTICAL,pos);
}
/*
 * Layout is performed by minimizing the P-stress goal function iteratively.
 * At each iteration taking a step in the steepest-descent direction.
 * x0 is the current position, x1 is the x0 - descentvector.
 */
void ConstrainedFDLayout::computeDescentVectorOnBothAxes(
        const bool xAxis, const bool yAxis,
        double stress, Position& x0, Position& x1) {
    setPosition(x0);
    if(xAxis) {
        applyForcesAndConstraints(vpsc::HORIZONTAL,stress);
    }
    if(yAxis) {
        applyForcesAndConstraints(vpsc::VERTICAL,stress);
    }
    getPosition(X,Y,x1);
}

/*
 * run() implements the main layout loop, taking descent steps until
 * stress is no-longer significantly reduced.
 * done is a callback used to check stress but also to report updated 
 * positions.
 */
void ConstrainedFDLayout::run(const bool xAxis, const bool yAxis) 
{
    if (extraConstraints.empty())
    {
        // This generates constraints for non-overlap inside and outside
        // of clusters.  To assign correct variable indexes it requires
        // that vs[] contains elements equal to the number of rectangles.
        vpsc::Variables vs[2];
        vs[0].resize(n);
        vs[1].resize(n);
        generateNonOverlapAndClusterCompoundConstraints(vs);
    }
    FILE_LOG(logDEBUG) << "ConstrainedFDLayout::run...";
    double stress=DBL_MAX;
    do {
        if(preIteration) {
            if(!(*preIteration)()) {
                break;
            }
            //printf("preIteration->changed=%d\n",preIteration->changed);
            if(preIteration->changed) {
                stress=DBL_MAX;
            }
            if(preIteration->resizes.size()>0) {
                FILE_LOG(logDEBUG) << " Resize event!";
                handleResizes(preIteration->resizes);
            }
        }
        unsigned N=2*n;
        Position x0(N),x1(N);
        getPosition(X,Y,x0);
        if(rungekutta) {
            Position a(N),b(N),c(N),d(N),ia(N),ib(N);
            computeDescentVectorOnBothAxes(xAxis,yAxis,stress,x0,a);
            ia=x0+(a-x0)/2.0;
            computeDescentVectorOnBothAxes(xAxis,yAxis,stress,ia,b);
            ib=x0+(b-x0)/2.0;
            computeDescentVectorOnBothAxes(xAxis,yAxis,stress,ib,c);
            computeDescentVectorOnBothAxes(xAxis,yAxis,stress,c,d);
            x1=a+2.0*b+2.0*c+d;
            x1/=6.0;
        } else {
            computeDescentVectorOnBothAxes(xAxis,yAxis,stress,x0,x1);
        }
        setPosition(x1);
        stress=computeStress();
        FILE_LOG(logDEBUG) << "stress="<<stress;
    } while(!(*done)(stress,X,Y));
    for(unsigned i=0;i<n;i++) {
        vpsc::Rectangle *r=boundingBoxes[i];
    FILE_LOG(logDEBUG) << *r;
    }
    FILE_LOG(logDEBUG) << "ConstrainedFDLayout::run done.";
    //m_final_stress = stress;
    most_recent_pure_stress = computePureStress();
}
/*
 * Same as run, but only applies one iteration.  This may be useful
 * where it's too hard to implement a call-back (e.g. in java apps)
 */
void ConstrainedFDLayout::runOnce(const bool xAxis, const bool yAxis) {
    if(n==0) return;
    double stress=DBL_MAX;
    unsigned N=2*n;
    Position x0(N),x1(N);
    getPosition(X,Y,x0);
    if(rungekutta) {
        Position a(N),b(N),c(N),d(N),ia(N),ib(N);
        computeDescentVectorOnBothAxes(xAxis,yAxis,stress,x0,a);
        ia=x0+(a-x0)/2.0;
        computeDescentVectorOnBothAxes(xAxis,yAxis,stress,ia,b);
        ib=x0+(b-x0)/2.0;
        computeDescentVectorOnBothAxes(xAxis,yAxis,stress,ib,c);
        computeDescentVectorOnBothAxes(xAxis,yAxis,stress,c,d);
        x1=a+2.0*b+2.0*c+d;
        x1/=6.0;
    } else {
        computeDescentVectorOnBothAxes(xAxis,yAxis,stress,x0,x1);
    }
}


// Used for sorting the CompoundConstraints from lowest priority to highest. 
static bool cmpCompoundConstraintPriority(const cola::CompoundConstraint *lhs, 
        const cola::CompoundConstraint *rhs)
{
    return lhs->priority() < rhs->priority();
}


void ConstrainedFDLayout::recGenerateClusterVariablesAndConstraints(
        vpsc::Variables (&vars)[2], unsigned int& priority, 
        cola::NonOverlapConstraints *noc, Cluster *cluster, 
        cola::CompoundConstraints& idleConstraints)
{
    for (std::vector<Cluster*>::iterator curr = cluster->clusters.begin();
            curr != cluster->clusters.end(); ++curr)
    {
        // For each of the child clusters, recursively call this function.
        recGenerateClusterVariablesAndConstraints(vars, priority, 
                noc, *curr, idleConstraints);
    }

    if ( (noc == NULL) && (dynamic_cast<RootCluster *> (cluster) == NULL) )
    {
        double freeWeight = 0.00000000001;
        // Then create left and right variables for the boundary of this 
        // cluster.
        vpsc::Variable *variable = NULL;
        cluster->clusterVarId = vars[XDIM].size();
        COLA_ASSERT(vars[XDIM].size() == vars[YDIM].size());
        // Left:
        variable = new vpsc::Variable(vars[XDIM].size(), 
                cluster->bounds.getMinX(), freeWeight);
        vars[XDIM].push_back(variable);
        // Right:
        variable = new vpsc::Variable(vars[XDIM].size(), 
                cluster->bounds.getMaxX(), freeWeight);
        vars[XDIM].push_back(variable);
        // Bottom::
        variable = new vpsc::Variable(vars[YDIM].size(), 
                cluster->bounds.getMinY(), freeWeight);
        vars[YDIM].push_back(variable);
        // Top:
        variable = new vpsc::Variable(vars[YDIM].size(), 
                cluster->bounds.getMaxY(), freeWeight);
        vars[YDIM].push_back(variable);

        RectangularCluster *rc = dynamic_cast<RectangularCluster *> (cluster);
        if (rc)
        {
            rc->generateFixedRectangleConstraints(idleConstraints, 
                    boundingBoxes, vars);
        }

        priority--;
        cola::ClusterContainmentConstraints *ccc = 
                new cola::ClusterContainmentConstraints(cluster, priority,
                        boundingBoxes);
        idleConstraints.push_back(ccc);
    }

    if (noc)
    {
        // Enforce non-overlap between all the shapes and clusters at this 
        // level.
        printf("Cluster #%d non-overlap constraints - nodes %d clusters %d\n",
                (int) cluster->clusterVarId, (int) cluster->nodes.size(), 
                (int) cluster->clusters.size());
        unsigned int group = cluster->clusterVarId;
        for (std::vector<unsigned>::iterator curr = cluster->nodes.begin();
                curr != cluster->nodes.end(); ++curr)
        {
            unsigned id = *curr;
            noc->addShape(id, boundingBoxes[id]->width() / 2,
                    boundingBoxes[id]->height() / 2, group);
        }
        for (std::vector<Cluster*>::iterator curr = cluster->clusters.begin();
                curr != cluster->clusters.end(); ++curr)
        {
            Cluster *cluster = *curr;
            RectangularCluster *rectCluster = 
                    dynamic_cast<RectangularCluster *> (cluster);
            if (rectCluster && rectCluster->clusterIsFromFixedRectangle())
            {
                // Treat it like a shape for non-overlap.
                unsigned id = rectCluster->rectangleIndex();
                noc->addShape(id, boundingBoxes[id]->width() / 2,
                        boundingBoxes[id]->height() / 2, group);
            }
            else
            {
                // Normal cluster.
                noc->addCluster(cluster, group);
            }
        }
    }
}

void ConstrainedFDLayout::generateNonOverlapAndClusterCompoundConstraints(
        vpsc::Variables (&vs)[2])
{
    if (clusterHierarchy && !clusterHierarchy->flat())
    {
        // Add remaining nodes that aren't contained within any clusters
        // as children of the root cluster.
        for (unsigned int i = 0; i < boundingBoxes.size(); ++i)
        {
            if (!clusterHierarchy->containsShape(i))
            {
                clusterHierarchy->nodes.push_back(i);
            }
        }

        // Add non-overlap and containment constraints for all clusters
        // and nodes.
        unsigned int priority = PRIORITY_NONOVERLAP;
        clusterHierarchy->computeBoundingRect(boundingBoxes);
        
        // Generate the containment constraints
        recGenerateClusterVariablesAndConstraints(vs, priority, 
                NULL, clusterHierarchy, extraConstraints);
        
        // Generate non-overlap constraints between all clusters and 
        // all contained nodes.
        if (m_generateNonOverlapConstraints)
        {
            priority--;
            cola::NonOverlapConstraints *noc = 
                    new cola::NonOverlapConstraints(priority);
            if (m_addGridSnapStress) {
                noc->setGridMode(true,m_snapGridX,m_snapGridY);
            }
            recGenerateClusterVariablesAndConstraints(vs, priority, 
                    noc, clusterHierarchy, extraConstraints);
            extraConstraints.push_back(noc);
        }
    }
    else if (m_generateNonOverlapConstraints)
    {
        // Add standard non-overlap constraints between each pair of
        // nodes.
        cola::NonOverlapConstraints *noc = 
                new cola::NonOverlapConstraints();
        if (m_addGridSnapStress) {
            noc->setGridMode(true,m_snapGridX,m_snapGridY);
        }
        for (unsigned int i = 0; i < boundingBoxes.size(); ++i)
        {
            noc->addShape(i, boundingBoxes[i]->width() / 2,
                    boundingBoxes[i]->height() / 2);
        }
        extraConstraints.push_back(noc);
    }
}

void ConstrainedFDLayout::makeFeasible(void)
{
    vpsc::Variables vs[2];
    vpsc::Constraints valid[2];

    vpsc::Rectangle::setXBorder(1);
    vpsc::Rectangle::setYBorder(1);
    
    // Populate all the variables for shapes.
    for (unsigned int dim = 0; dim < 2; ++dim)
    {
        vs[dim] = vpsc::Variables(boundingBoxes.size());
        for (unsigned int i = 0; i < vs[dim].size(); ++i)
        {
            double pos = (dim == 0) ? 
                    boundingBoxes[i]->getCentreX() : 
                    boundingBoxes[i]->getCentreY();
            vs[dim][i] = new vpsc::Variable(i, pos, 1);
        }
    }

    vector<double> priorPos(boundingBoxes.size());

    // Clear extra constraints for cluster containment and non-overlap.
    // We keep a separate list of these since we keep them around for 
    // later solving.
    for_each(extraConstraints.begin(), extraConstraints.end(), delete_object());
    extraConstraints.clear();

    generateNonOverlapAndClusterCompoundConstraints(vs);

    // Make a copy of the compound constraints and sort them by priority.
    cola::CompoundConstraints idleConstraints = ccs;
    // Append extraConstraints to idleConstraints.
    idleConstraints.insert(idleConstraints.end(), 
            extraConstraints.begin(), extraConstraints.end());

    std::sort(idleConstraints.begin(), idleConstraints.end(), 
            cmpCompoundConstraintPriority);

    // Initialise extra variables for compound constraints.
    for (unsigned int dim = 0; dim < 2; ++dim)
    {
        generateVariables(idleConstraints, (vpsc::Dim) dim, vs[dim]);
    }

#ifdef MAKEFEASIBLE_DEBUG
    char filename[200];
    int iteration = 0;
    vector<string> labels(boundingBoxes.size());
    for(unsigned i=0;i<boundingBoxes.size();++i)
    {
        stringstream ss;
        ss << i;
        labels[i]=ss.str();
    }
#endif

    // Main makeFeasible loop.
    while (!idleConstraints.empty())
    {
        // idleConstraints is sorted lowest to highest priority, so the 
        // highest priority constraint will be at the back of the vector.
        cola::CompoundConstraint *cc = idleConstraints.back();
        idleConstraints.pop_back();

#ifdef MAKEFEASIBLE_DEBUG
        // Debugging SVG time slice output.
        std::vector<cola::Edge> es;
        for (unsigned int i = 0; i < boundingBoxes.size(); ++i)
        {
            boundingBoxes[i]->moveCentreX(vs[0][i]->finalPosition);
            boundingBoxes[i]->moveCentreY(vs[1][i]->finalPosition);
        }
        iteration++;
        sprintf(filename, "out/file-%05d.pdf", iteration);

        OutputFile of(boundingBoxes,es,clusterHierarchy,filename,true,false);
        of.setLabels(labels);
        of.generate();
#endif

        cc->markAllSubConstraintsAsInactive();
        bool subConstraintSatisfiable = true;
        
        if (cc->shouldCombineSubConstraints())
        {
            // We are processing a combined set of satisfiable constraints,
            // such as for containment within cluster boundary variables, so
            // we just add all the required constraints and solve in both
            // the X and Y dimension once to set the cluster boundaries to 
            // meaningful values.
            while (cc->subConstraintsRemaining())
            {
                cola::SubConstraintAlternatives alternatives = 
                        cc->getCurrSubConstraintAlternatives(vs);
                // There should be no alternatives, just guaranteed
                // satisfiable constraints.
                COLA_ASSERT(alternatives.size() == 1);
                vpsc::Dim& dim = alternatives.front().dim;
                vpsc::Constraint& constraint = alternatives.front().constraint;
                valid[dim].push_back(new vpsc::Constraint(constraint));
                cc->markCurrSubConstraintAsActive(subConstraintSatisfiable);
            }
            for (size_t dim = 0; dim < 2; ++dim)
            {
                vpsc::IncSolver vpscInstance(vs[dim], valid[dim]);
                vpscInstance.satisfy();
            }
            continue;
        }

        while (cc->subConstraintsRemaining())
        {
            cola::SubConstraintAlternatives alternatives = 
                    cc->getCurrSubConstraintAlternatives(vs);
            alternatives.sort();

            if (alternatives.empty())
            {
                continue;
            }

            while (!alternatives.empty())
            {
                // Reset subConstraintSatisfiable for new solve.
                subConstraintSatisfiable = true;

                vpsc::Dim& dim = alternatives.front().dim;
                vpsc::Constraint& constraint = alternatives.front().constraint;
            
                // Store current values for variables.
                for (unsigned int i = 0; i < priorPos.size(); ++i)
                {
                    priorPos[i] = vs[dim][i]->finalPosition;
                }

                // Some solving...
                try 
                {
                    // Add the constraint from this alternative to the 
                    // valid constraint set.
                    valid[dim].push_back(new vpsc::Constraint(constraint));

                    //fprintf(stderr, ".%d %3d - ", dim, valid[dim].size());
                    // Solve with this constraint set.
                    vpsc::IncSolver vpscInstance(vs[dim], valid[dim]);
                    vpscInstance.satisfy();
                }
                catch (char *str) 
                {
                    subConstraintSatisfiable = false;
                    
                    std::cerr << "++++ IN ERROR BLOCK" << std::endl;
                    std::cerr << str << std::endl;
                    for (vpsc::Rectangles::iterator r = boundingBoxes.begin(); 
                            r != boundingBoxes.end(); ++r) 
                    {
                        std::cerr << **r <<std::endl;
                    }
                }
                for (size_t i = 0; i < valid[dim].size(); ++i) 
                {
                    if (valid[dim][i]->unsatisfiable) 
                    {
                        // It might have made one of the earlier added 
                        // constraints unsatisfiable, so we mark that one 
                        // as okay since we will be reverting the most 
                        // recent one.
                        valid[dim][i]->unsatisfiable = false;
                        
                        subConstraintSatisfiable = false;
                    }
                }

                if (!subConstraintSatisfiable)
                {
                    //fprintf(stderr, "*");

                    // Restore previous values for variables.
                    for (unsigned int i = 0; i < priorPos.size(); ++i)
                    {
                        vs[dim][i]->finalPosition = priorPos[i];
                    }
                    
                    // Delete the newly added (and unsatisfiable) 
                    // constraint from the valid constraint set.
                    delete valid[dim].back();
                    valid[dim].pop_back();
                }
                else
                {
                    break;
                }
                // Move on to the next alternative.
                alternatives.pop_front();
            }
#ifdef MAKEFEASIBLE_DEBUG
            if (true || idleConstraints.size() == 0)
            {
                // Debugging SVG time slice output, but don't show this for
                // constraints that promised satisfiable.
                std::vector<cola::Edge> es;
                for (unsigned int i = 0; i < boundingBoxes.size(); ++i)
                {
                    boundingBoxes[i]->moveCentreX(vs[0][i]->finalPosition);
                    boundingBoxes[i]->moveCentreY(vs[1][i]->finalPosition);
                }
                iteration++;
                sprintf(filename, "out/file-%05d.pdf", iteration);

                OutputFile of(boundingBoxes,es,clusterHierarchy,filename,
                        true,false);
                of.setLabels(labels);
                of.generate();
            }
#endif
            cc->markCurrSubConstraintAsActive(subConstraintSatisfiable);
        }
    }

    // Write positions from solver variables back to Rectangles.
    for (unsigned int i = 0; i < boundingBoxes.size(); ++i)
    {
        boundingBoxes[i]->moveCentreX(vs[0][i]->finalPosition);
        boundingBoxes[i]->moveCentreY(vs[1][i]->finalPosition);
    }

    vpsc::Rectangle::setXBorder(0);
    vpsc::Rectangle::setYBorder(0);

    // Cleanup.
    for (unsigned int dim = 0; dim < 2; ++dim)
    {
        for_each(valid[dim].begin(), valid[dim].end(), delete_object());
        for_each(vs[dim].begin(), vs[dim].end(), delete_object());
    }

    topologyAddon->makeFeasible(m_generateNonOverlapConstraints,
            boundingBoxes, clusterHierarchy);
    
    // Update the X and Y vectors with the new shape positions.
    for (unsigned int i = 0; i < boundingBoxes.size(); ++i)
    {
        X[i] = boundingBoxes[i]->getCentreX();
        Y[i] = boundingBoxes[i]->getCentreY();
    }
}

ConstrainedFDLayout::~ConstrainedFDLayout()
{
    if (using_default_done)
    {
        delete done;
    }

    for (unsigned i = 0; i < n; ++i)
    {
        delete [] G[i];
        delete [] D[i];
    }
    delete [] G;
    delete [] D;
    delete topologyAddon;
}

void ConstrainedFDLayout::freeAssociatedObjects(void)
{
    // Free Rectangles
    for_each(boundingBoxes.begin(), boundingBoxes.end(), delete_object());
    boundingBoxes.clear();
    
    // Free compound constraints
    std::list<CompoundConstraint *> freeList(ccs.begin(), ccs.end());
    freeList.sort();
    freeList.unique();
    if (freeList.size() != ccs.size())
    {
        // The compound constraint list had repeated elements.
        fprintf(stderr, "Warning: CompoundConstraints vector contained %d "
                "duplicates.\n", (int) (ccs.size() - freeList.size()));
    }
    ccs.clear();
    for_each(freeList.begin(), freeList.end(), delete_object());

    if (clusterHierarchy)
    {
        delete clusterHierarchy;
        clusterHierarchy = NULL;
    }
}

void ConstrainedFDLayout::setTopology(TopologyAddonInterface *newTopology)
{
    COLA_ASSERT(topologyAddon);
    delete topologyAddon;
    topologyAddon = newTopology->clone();
}

TopologyAddonInterface *ConstrainedFDLayout::getTopology(void)
{
    return topologyAddon->clone();
}


void setupVarsAndConstraints(unsigned n, const CompoundConstraints& ccs,
        const vpsc::Dim dim, vpsc::Rectangles& boundingBoxes,
        RootCluster *clusterHierarchy,
        vpsc::Variables& vs, vpsc::Constraints& cs, 
        valarray<double> &coords) 
{
    vs.resize(n);
    for (unsigned i = 0; i < n; ++i)
    {
        vs[i] = new vpsc::Variable(i, coords[i]);
    }

    if (clusterHierarchy && !clusterHierarchy->clusters.empty())
    {
        // Create variables for clusters
        clusterHierarchy->computeBoundingRect(boundingBoxes);
        clusterHierarchy->createVars(dim, boundingBoxes, vs);
    }

    for (CompoundConstraints::const_iterator c = ccs.begin();
            c != ccs.end(); ++c) 
    {
        (*c)->generateVariables(dim, vs);
    }
    for (CompoundConstraints::const_iterator c = ccs.begin();
            c != ccs.end(); ++c) 
    {
        (*c)->generateSeparationConstraints(dim, vs, cs, boundingBoxes);
    }
}


static void setupExtraConstraints(const CompoundConstraints& ccs,
        const vpsc::Dim dim, vpsc::Variables& vs, vpsc::Constraints& cs,
        vpsc::Rectangles& boundingBoxes)
{
    for (CompoundConstraints::const_iterator c = ccs.begin();
            c != ccs.end(); ++c) 
    {
        (*c)->generateVariables(dim, vs);
    }
    for (CompoundConstraints::const_iterator c = ccs.begin();
            c != ccs.end(); ++c) 
    {
        (*c)->generateSeparationConstraints(dim, vs, cs, boundingBoxes);
    }
}

void updateCompoundConstraints(const vpsc::Dim dim,
        const CompoundConstraints& ccs) 
{
    for (CompoundConstraints::const_iterator c = ccs.begin();
            c != ccs.end(); ++c) 
    {
        (*c)->updatePosition(dim);
    }
}
void project(vpsc::Variables& vs, vpsc::Constraints& cs, valarray<double>& coords) {
    unsigned n=coords.size();
    vpsc::IncSolver s(vs,cs);
    s.solve();
    for(unsigned i=0;i<n;++i) {
        coords[i]=vs[i]->finalPosition;
    }
}
void setVariableDesiredPositions(vpsc::Variables& vs, vpsc::Constraints& cs,
        const DesiredPositionsInDim& des, valarray<double>& coords)
{
    COLA_UNUSED(cs);

    unsigned n=coords.size();
    COLA_ASSERT(vs.size()>=n);
    for(unsigned i=0;i<n;++i) {
        vpsc::Variable* v=vs[i];
        v->desiredPosition = coords[i];
        v->weight=1;
    }
    for (DesiredPositionsInDim::const_iterator d=des.begin();
            d!=des.end(); ++d) {
        COLA_ASSERT(d->first<vs.size());
        vpsc::Variable* v=vs[d->first];
        v->desiredPosition = d->second;
        v->weight=10000;
    }
}
void checkUnsatisfiable(const vpsc::Constraints& cs, 
        UnsatisfiableConstraintInfos* unsatisfiable) {
    for(vpsc::Constraints::const_iterator c=cs.begin();c!=cs.end();++c) {
        if((*c)->unsatisfiable) {
            UnsatisfiableConstraintInfo* i=new UnsatisfiableConstraintInfo(*c);
            unsatisfiable->push_back(i);
        }
    }
}

void ConstrainedFDLayout::handleResizes(const Resizes& resizeList)
{
    topologyAddon->handleResizes(resizeList, n, X, Y, ccs, boundingBoxes,
            clusterHierarchy);
}
/*
 * move positions of nodes in specified axis while respecting constraints
 * @param dim axis
 * @param target array of desired positions (for both axes)
 */
void ConstrainedFDLayout::moveTo(const vpsc::Dim dim, Position& target) {
    COLA_ASSERT(target.size()==2*n);
    FILE_LOG(logDEBUG) << "ConstrainedFDLayout::moveTo(): dim="<<dim;
    valarray<double> &coords = (dim==vpsc::HORIZONTAL)?X:Y;
    vpsc::Variables vs;
    vpsc::Constraints cs;
    setupVarsAndConstraints(n, ccs, dim, boundingBoxes,
            clusterHierarchy, vs, cs, coords);
    DesiredPositionsInDim des;
    if(preIteration) {
        for(vector<Lock>::iterator l=preIteration->locks.begin();
                l!=preIteration->locks.end();l++) {
            des.push_back(make_pair(l->getID(),l->pos(dim)));
            FILE_LOG(logDEBUG1)<<"desi: v["<<l->getID()<<"]=("<<l->pos(vpsc::HORIZONTAL)
                <<","<<l->pos(vpsc::VERTICAL)<<")";
        }
    }
    for(unsigned i=0, j=(dim==vpsc::HORIZONTAL?0:n);i<n;++i,++j) {
        vpsc::Variable* v=vs[i];
        v->desiredPosition = target[j];
    }
    setVariableDesiredPositions(vs,cs,des,coords);
    if (topologyAddon->useTopologySolver())
    {
        topologyAddon->moveTo(dim, vs, cs, coords, clusterHierarchy);
    } else {
        // Add non-overlap constraints, but not variables again.
        setupExtraConstraints(extraConstraints, dim, vs, cs, boundingBoxes);
        // Projection.
        project(vs,cs,coords);
        moveBoundingBoxes();
    }
    updateCompoundConstraints(dim, ccs);
    for_each(vs.begin(),vs.end(),delete_object());
    for_each(cs.begin(),cs.end(),delete_object());
}

/*
QString ConstrainedFDLayout::writeVADVect(valarray<double> v) const {
    QString s = "";
    for (int i = 0; i < v.size(); i++) {
        s += QString::number(v[i],'g',3)+"  ";
    }
    return s;
}

QString ConstrainedFDLayout::writeSparseMatrix(SparseMatrix M, bool linear) const {
    QString s = "";
    QString sep = linear?" | ":"\n";
    unsigned n = M.rowSize();
    for(unsigned i=0;i<n;i++) {
        for(unsigned j=0;j<n;j++) {
            s += QString::number(M.getIJ(i,j),'g',3)+"  ";
        }
        s += sep;
    }
    return s;
}
*/

/*
 * The following computes an unconstrained solution then uses Projection to
 * make this solution feasible with respect to constraints by moving things as
 * little as possible.  If "meta-constraints" such as avoidOverlaps or edge
 * straightening are required then dummy variables will be generated.
 */
double ConstrainedFDLayout::applyForcesAndConstraints(const vpsc::Dim dim, const double oldStress) {
    FILE_LOG(logDEBUG) << "ConstrainedFDLayout::applyForcesAndConstraints(): dim="<<dim;
    valarray<double> g(n);
    valarray<double> &coords = (dim==vpsc::HORIZONTAL)?X:Y;

    //if (dim==vpsc::HORIZONTAL) qDebug() << "x: " << writeVADVect(coords);

    DesiredPositionsInDim des;
    if(preIteration) {
        for(vector<Lock>::iterator l=preIteration->locks.begin();
                l!=preIteration->locks.end();l++) {
            des.push_back(make_pair(l->getID(),l->pos(dim)));
            FILE_LOG(logDEBUG1)<<"desi: v["<<l->getID()<<"]=("<<l->pos(vpsc::HORIZONTAL)
                <<","<<l->pos(vpsc::VERTICAL)<<")";
        }
    }
    vpsc::Variables vs;
    vpsc::Constraints cs;
    double stress;
    setupVarsAndConstraints(n, ccs, dim, boundingBoxes,
            clusterHierarchy, vs, cs, coords);

    if (topologyAddon->useTopologySolver())
    {
        stress = topologyAddon->applyForcesAndConstraints(this, dim, g, vs, cs, 
                coords, des, oldStress);
    } else {
        // Add non-overlap constraints, but not variables again.
        setupExtraConstraints(extraConstraints, dim, vs, cs, boundingBoxes);
        // Projection.
        SparseMap HMap(n);
        //qDebug() << "Computing Forces ========================================================";

        //qDebug() << m_debugLineNo << " ----------";
        //m_debugLineNo++;

        computeForces(dim,HMap,g);

        // Set dgdv values.
        for (int i = 0; i < n; i++) {
            vs.at(i)->set_dgdv(g[i]);
        }

        //if (dim==vpsc::HORIZONTAL) qDebug() << "g: " << writeVADVect(g);
        SparseMatrix H(HMap);
        //if (dim==vpsc::HORIZONTAL) qDebug() << "H: " << writeSparseMatrix(H,true);
        valarray<double> oldCoords=coords;
        double ss = computeStepSize(H,g,g);
        //if (dim==vpsc::HORIZONTAL) qDebug() << "ss: " << ss;
        applyDescentVector(g,oldCoords,coords,oldStress,ss);
        //if (dim==vpsc::HORIZONTAL) qDebug() << "c1: " << writeVADVect(coords);

        setVariableDesiredPositions(vs,cs,des,coords);
        //if (dim==vpsc::HORIZONTAL) qDebug() << "c2: " << writeVADVect(coords);


        //project(vs,cs,coords);
        // We do the projection here, so we can keep the solver for later use.
        unsigned n=coords.size();
        if (m_solver) {
            delete m_solver;
        }
        m_solver = new vpsc::IncSolver(vs,cs);
        m_solver->solve();
        for(unsigned i=0;i<n;++i) {
            coords[i]=vs[i]->finalPosition;
        }
        // end projection

        // Get the tentative constraint with largest absolute LM (may be null).
        vpsc::Constraint *max_abs_lm = m_solver->max_abs_lm;
        if (max_abs_lm!=NULL) {
            //qDebug() << "HAVE max_abs_lm";
            // If its LM is large enough (in abs. val.), then mark the
            // compound constraint to which it belongs as to-be-rejected.
            //qDebug() << "Max lm in abs val: " << fabs(max_abs_lm->lm);
            //qDebug() << "Threshold: " << m_tentative_constraint_threshold;
            if (fabs(max_abs_lm->lm) >= m_tentative_constraint_threshold) {
                m_constraintToReject = max_abs_lm->compoundOwner;
            }
        } else {
            //qDebug() << "no max_abs_lm";
        }
        // Check if any tentative constraint was marked for rejection.
        vpsc::Constraint *rejected = m_solver->rejected_tentative_constraint;
        if (rejected) {
            m_constraintToRejectForUnsat = rejected->compoundOwner;
            m_alignedShapeVarIndexForUnsat = rejected->alignedShapeVarIndex;
        }

        //if (dim==vpsc::HORIZONTAL) qDebug() << "c3: " << writeVADVect(coords);

        valarray<double> d(n);
        d=oldCoords-coords;
        //if (dim==vpsc::HORIZONTAL) qDebug() << "d: " << writeVADVect(d);
        double stepsize=computeStepSize(H,g,d);
        stepsize=max(0.,min(stepsize,1.));
        //qDebug() << "stepsize: " << stepsize;
        //printf(" dim=%d beta: ",dim);
        stress = applyDescentVector(d,oldCoords,coords,oldStress,stepsize);

        moveBoundingBoxes();
    }
    updateCompoundConstraints(dim, ccs);
    if(unsatisfiable.size()==2) {
        checkUnsatisfiable(cs,unsatisfiable[dim]);
    }
    FILE_LOG(logDEBUG) << "ConstrainedFDLayout::applyForcesAndConstraints... done, stress="<<stress;
    for_each(vs.begin(),vs.end(),delete_object());
    for_each(cs.begin(),cs.end(),delete_object());
    return stress;
}

vpsc::Blocks *ConstrainedFDLayout::getBlocks() {
    if (m_solver==NULL) {
        return NULL;
    } else {
        return m_solver->getBlocks();
    }
}

/*
 * Attempts to set coords=oldCoords-stepsize*d.  If this does not reduce
 * the stress from oldStress then stepsize is halved.  This is repeated
 * until stepsize falls below a threshhold.
 * @param d is a descent vector (a movement vector intended to reduce the
 * stress)
 * @param oldCoords are the previous position vector
 * @param coords will hold the new position after applying d
 * @param stepsize is a scalar multiple of the d to apply
 */
double ConstrainedFDLayout::applyDescentVector(
        valarray<double> const &d,
        valarray<double> const &oldCoords,
        valarray<double> &coords,
        const double oldStress,
        double stepsize
        )
{
    COLA_UNUSED(oldStress);

    COLA_ASSERT(d.size()==oldCoords.size());
    COLA_ASSERT(d.size()==coords.size());
    while(fabs(stepsize)>0.00000000001) {
        coords=oldCoords-stepsize*d;
        double stress=computeStress();
        //printf(" applyDV: oldstress=%f, stress=%f, stepsize=%f\n", oldStress,stress,stepsize);
        //if(oldStress>=stress) {
            return stress;
        //}
        coords=oldCoords;
        stepsize*=0.5;
    }
    return computeStress();
}
        
/*
 * Computes:
 *  - the matrix of second derivatives (the Hessian) H, used in 
 *    calculating stepsize; and
 *  - the vector g, the negative gradient (steepest-descent) direction.
 */
void ConstrainedFDLayout::computeForces(
        const vpsc::Dim dim,
        SparseMap &H,
        valarray<double> &g) {
    if(n==1) return;
    g=0;
#define BASICSTRESS
#ifdef  BASICSTRESS
    // for each node:
    for(unsigned u=0;u<n;u++) {
        // Stress model
        double Huu=0;
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double rx=X[u]-X[v], ry=Y[u]-Y[v];
            double l=sqrt(rx*rx+ry*ry);
            double d=D[u][v];
            if(l>d && p>1) continue; // attractive forces not required
            double d2=d*d;
            /* force apart zero distances */
            if (l < 1e-30) {
                l=0.1;
            }
            double dx=dim==vpsc::HORIZONTAL?rx:ry;
            double dy=dim==vpsc::HORIZONTAL?ry:rx;
            g[u]+=dx*(l-d)/(d2*l);
            Huu-=H(u,v)=(d*dy*dy/(l*l*l)-1)/d2;
        }
        H(u,u)=Huu;
    }
#endif
    if(m_addSnapStress) {
        computeSnapForces(dim,H,g);
    }
    if(m_addGridSnapStress) {
        computeGridSnapForces(dim,H,g);
    }
    if(m_addEdgeNodeRepulsion) {
        computeEdgeNodeRepulsionForces(dim,H,g);
    }
    if(desiredPositions) {
        for(DesiredPositions::const_iterator p=desiredPositions->begin();
            p!=desiredPositions->end();++p) {
            unsigned i = p->id;
            double d=(dim==vpsc::HORIZONTAL)
                ?p->x-X[i]:p->y-Y[i];
            d*=p->weight;
            g[i]-=d;
            H(i,i)+=p->weight;
        }
    }
}

// Force chooser:
void ConstrainedFDLayout::computeSnapForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    switch (m_snapStressFunction) {
    case 1:
        smoothMForces(dim,H,g);
        break;
    case 2:
        linearMForces(dim,H,g);
        break;
    case 3:
        dualQuadraticForces(dim,H,g);
        break;
    case 4:
        invertedQuadraticForces(dim,H,g);
        break;
    case 5:
        quarticForces(dim,H,g);
        break;
    case 6:
        smoothVForces(dim,H,g);
        break;
    case 7:
        linearVForces(dim,H,g);
        break;
    case 8:
#define useTailoredNodeSnapForces
#ifdef  useTailoredNodeSnapForces
        tailoredQuadUForces(dim,H,g);
#else
        quadUForces(dim,H,g);
#endif
        break;
    }
}

// Grid forces; uses quadratic U-stress:
void ConstrainedFDLayout::computeGridSnapForces(const vpsc::Dim dim, SparseMap &H, valarray<double> &g) {
    double b = m_snap_strength;
    double sig = m_snap_distance;
    double k = b/(sig*sig);
    for(unsigned u=0;u<n;u++) {
        double z=dim==vpsc::HORIZONTAL?X[u]:Y[u];
        double D=dim==vpsc::HORIZONTAL?m_snapGridX:m_snapGridY;
        double q,r;
        double d=0,t=0;
        r = fabs(modf(z/D, &q));
// If want to break ties halfway between grid lines so that we always
// favour the side closer to the origin, then define smallerCoordsWin.
// Otherwise zero force will be exerted at the midpoints.
#define smallerCoordsWin
#ifndef smallerCoordsWin
        if (r!=0.5) { // no "tug of war"
            d = r<0.5 ? z-q*D : ( z>0 ? z-(q+1)*D : z-(q-1)*D );
            if (-sig<=d && d<=sig) {
                double dg = k*d;
                double dH = k;
                g[u]+=dg;
                H(u,u)+=dH;
            }
        }
#else
        d = r<=0.5 ? z-q*D : ( z>0 ? z-(q+1)*D : z-(q-1)*D );
        if (-sig<=d && d<=sig) {
            double dg = k*d;
            double dH = k;
            g[u]+=dg;
            H(u,u)+=dH;
        }
#endif
    }
}

void ConstrainedFDLayout::computeEdgeNodeRepulsionForces(const vpsc::Dim dim, SparseMap &H, valarray<double> &g) {
    double b = m_snap_strength;
    double sig = m_snap_distance;
    double k = b/(sig*sig);
    double K = 20*k; // stronger than grid force
    foreach(Edge e, edges) {
        // Is e axis-aligned in the dimension we are handling?
        unsigned u = e.first, v = e.second;
        double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
        double eps = 5;
        if (fabs(d) > eps) continue;
        // If so then proceed.
        double ez, el, eh;
        if (dim==vpsc::HORIZONTAL) {
            ez = (X[u]+X[v])/2.0; el = min(Y[u],Y[v]); eh = max(Y[u],Y[v]);
        } else {
            ez = (Y[u]+Y[v])/2.0; el = min(X[u],X[v]); eh = max(X[u],X[v]);
        }
        // Consider each node.
        for(unsigned w=0;w<n;w++) {
            if (w==u || w==v) continue;
            double z=dim==vpsc::HORIZONTAL?X[w]:Y[w];
            double k=dim==vpsc::HORIZONTAL?Y[w]:X[w];
            if (k<=el || eh<=k) continue;
            double r = z - ez;
            if (fabs(r) > sig) continue;
            r = r>=0?r-sig:r+sig;
            double dg = K*r;
            double dH = K;
            g[u]+=-dg;
            H(u,u)+=dH;
            //fprintf(stderr, "dg=%f, dH=%f\n",dg,dH);
        }
    }
}

// Quadratic U-stress forces tailored to each pair of nodes' dimensions:
void ConstrainedFDLayout::tailoredQuadUForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double b = m_snap_strength;
    foreach(Edge e, edges) {
        unsigned u = e.first, v = e.second;
        double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
        vpsc::Rectangle *ru = boundingBoxes.at(u);;
        vpsc::Rectangle *rv = boundingBoxes.at(v);
        double sig=dim==vpsc::HORIZONTAL?(ru->width()+rv->width())/2:(ru->height()+rv->height())/2;
        double k = b/(sig*sig);
        if (-sig < d && d <= sig) {
            double kd = k*d;
            g[u]+=kd;
            H(u,v)-=k;
            H(u,u)+=k;
            g[v]-=kd;
            H(v,u)-=k;
            H(v,v)+=k;
        }
    }
}

// Quadratic U-stress forces:
void ConstrainedFDLayout::quadUForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double b = m_snapStressBeta;
    double sig = m_snapStressSigma;
    double k = b/(sig*sig);
#define nbrSnap
#ifdef nbrSnap
    foreach(Edge e, edges) {
        unsigned u = e.first, v = e.second;
        double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
        if (-sig < d && d <= sig) {
            double kd = k*d;
            g[u]+=kd;
            H(u,v)-=k;
            H(u,u)+=k;
            g[v]-=kd;
            H(v,u)-=k;
            H(v,v)+=k;
        }
        /*
        unsigned u = e.first, v = e.second;
        double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
        if (-sig <= d && d <= sig) {
            g[u]+=k*d;
            H(u,v)-=k;
            H(u,u)+=k;
        }
        u = e.second, v = e.first;
        d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
        if (-sig <= d && d <= sig) {
            g[u]+=k*d;
            H(u,v)-=k;
            H(u,u)+=k;
        }
        */
    }
#else
    for(unsigned u=0;u<n;u++) {
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            if (-sig <= d && d <= sig) {
                g[u]+=k*d;
                H(u,v)-=k;
                H(u,u)+=k;
            }
            //qDebug() << "H(u,v) -= " << k;
            //qDebug() << "g[u] += " << k*d;
        }
    }
#endif
}

// Smooth V-stress forces:
void ConstrainedFDLayout::smoothVForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double a = m_snapStressAlpha;
    double b = m_snapStressBeta;
    double sig = m_snapStressSigma;
    double l3 = sig*sig*sig;
    for(unsigned u=0;u<n;u++) {
        double Huu=0, gu=0;
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            double r=0, q=1, t=0;
            if (-sig <= d && d < 0) {
                r=-sig; q=0;
                t=-1;
            } else if (0 <= d && d < sig) {
                r=0; q=sig;
                t=1;
            }
            double dmp=d-r, dmq=d-q;
            double c=-3/l3;
            gu+=t*b*c*dmp*dmq;
            //double Huv=-t*b*c*(dmp+dmq);
            double Huv=t*b*c*(dmp+dmq); // !! Try opposite sign, since getting weird results!
            //qDebug() << "Adding " + QString::number(Huv,'g',3) + " to H at " + QString::number(u)+","+QString::number(v);
            H(u,v)+=Huv;
            Huu-=Huv;
        }
        g[u]+=gu;
        H(u,u)+=Huu;
    }
}

// Linear V-stress forces:
void ConstrainedFDLayout::linearVForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double b = m_snapStressBeta;
    double sig = m_snapStressSigma;
    double q = b/sig;
    for(unsigned u=0;u<n;u++) {
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            double r=0;
            if (-sig <= d && d < 0) {
                g[u] -= q;
            } else if (0 <= d && d < sig) {
                g[u] += q;
            }
        }
    }
}

// Smooth M-stress forces:
void ConstrainedFDLayout::smoothMForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double a = m_snapStressAlpha;
    double b = m_snapStressBeta;
    double sig = m_snapStressSigma;
    for(unsigned u=0;u<n;u++) {
        double Huu=0, gu=0;
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            double r=0, q=0, t=1;
            if (-sig-a <= d && d < -sig) {
                r=-sig-a; q=-sig;
                t=1;
            } else if (-sig <= d && d < 0) {
                r=-sig; q=0;
                t=-1;
            } else if (0 <= d && d < sig) {
                r=0; q=sig;
                t=1;
            } else if (sig <= d && d <= sig+a) {
                r=sig; q=sig+a;
                t=-1;
            }
            double dmp=d-r, dmq=d-q, l=q-r;
            double l3=l*l*l;
            double c=-3/l3;
            gu+=t*b*c*dmp*dmq;
            double Huv=-t*b*c*(dmp+dmq);
            H(u,v)+=Huv;
            Huu-=Huv;
        }
        g[u]+=gu;
        H(u,u)+=Huu;
    }
}


// Linear M-stress forces:
void ConstrainedFDLayout::linearMForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double a = m_snapStressAlpha;
    double b = m_snapStressBeta;
    double sig = m_snapStressSigma;
    for(unsigned u=0;u<n;u++) {
        double gu=0;
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            double s = 0;
            if (-sig-a <= d && d < -sig) {
                s = 1.0/a;
            } else if (-sig <= d && d < 0) {
                s = -1.0/sig;
            } else if (0 <= d && d < sig) {
                s = 1.0/sig;
            } else if (sig <= d && d < sig+a) {
                s = -1.0/a;
            }
            gu+=s;
        }
        g[u]+=b*gu;
    }
}

// Summed quadratic bumps snap stress forces
void ConstrainedFDLayout::dualQuadraticForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    //qDebug() << "dual quadratic force";
    double a = m_snapStressAlpha;
    double c = -a*m_snapStressBeta;
    double s = m_snapStressSigma;
    for(unsigned u=0;u<n;u++) {
        double Huu=0, gu=0;
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            double dp=d+s, dn=d-s;
            double dp2=dp*dp, dn2=dn*dn;
            double ep=1+a*dp2, en=1+a*dn2;
            double ep2=ep*ep, en2=en*en;
            double ep3=ep2*ep, en3=en2*en;
            gu+=dp/ep2 + dn/en2;
            double Huv = c*( (3*a*dp2 - 1)/ep3 + (3*a*dn2 - 1)/en3 );
            H(u,v)+=Huv;
            Huu-=Huv;
        }
        g[u]+=c*gu;
        H(u,u)+=Huu;
    }
}


// Snap stress with long-range forces
void ConstrainedFDLayout::invertedQuadraticForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double a = m_snapStressAlpha;
    double k = m_snapStressGamma;
    double c = a*k*m_snapStressBeta;
    for(unsigned u=0;u<n;u++) {
        double Huu=0, gu=0;
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            double d2 = d*d;
            double e = k+a*d2;
            double e2 = e*e;
            double e3 = e2*e;
            gu+=d/e2;
            double Huv = (c/e3)*(3*a*d2 - k);
            H(u,v)+=Huv;
            Huu-=Huv;
        }
        g[u]+=c*gu;
        H(u,u)+=Huu;
    }
}


// Quartic snap stress forces:
void ConstrainedFDLayout::quarticForces(const vpsc::Dim dim, SparseMap &H, std::valarray<double> &g) {
    double a = m_snapStressAlpha;
    double c = 2*a*m_snapStressBeta;
    double sig2 = m_snapStressSigma*m_snapStressSigma;
    for(unsigned u=0;u<n;u++) {
        double Huu=0, gu=0;
        for(unsigned v=0;v<n;v++) {
            if(u==v) continue;
            unsigned short p = G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double d=dim==vpsc::HORIZONTAL?X[u]-X[v]:Y[u]-Y[v];
            double d2 = d*d;
            double f=d*d-sig2;
            double e = 1+a*f*f;
            double e2 = e*e;
            gu-=f*d/e2;
            double Huv = (c/e2)*(2*d2 + f*(1 - 8*a*f*d2/e));
            H(u,v)+=Huv;
            Huu-=Huv;
        }
        g[u]+=c*gu;
        H(u,u)+=Huu;
    }
}


/*
 * Returns the optimal step-size in the direction d, given gradient g and 
 * hessian H.
 */
double ConstrainedFDLayout::computeStepSize(
        SparseMatrix const &H, 
        valarray<double> const &g, 
        valarray<double> const &d) const
{
    //qDebug() << "g: " << writeVADVect(g);
    //qDebug() << "H: " << writeSparseMatrix(H,true);
    //qDebug() << "d: " << writeVADVect(d);
    COLA_ASSERT(g.size()==d.size());
    COLA_ASSERT(g.size()==H.rowSize());
    // stepsize = g'd / (d' H d)
    double numerator = dotProd(g,d);
    valarray<double> Hd(d.size());
    H.rightMultiply(d,Hd);
    double denominator = dotProd(d,Hd);
    //COLA_ASSERT(numerator>=0);
    //COLA_ASSERT(denominator>=0);
    //qDebug() << "Step size numerator: " << numerator;
    //qDebug() << "Step size denominator: " << denominator;
    double epsilon = 0.0001;
    //if(denominator<epsilon) denominator=0;
    //if(denominator!=0) qDebug() << "step size: " << numerator/denominator;
    if(denominator==0) return 0;
    return numerator/denominator;
}
/*
 * Just computes the cost (Stress) at the current X,Y position
 * used to test termination.
 * This method will call preIteration if one is set.
 */
double ConstrainedFDLayout::computeStress() const {
    FILE_LOG(logDEBUG)<<"ConstrainedFDLayout::computeStress()";
    double stress=0;
#ifdef  BASICSTRESS
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double rx=X[u]-X[v], ry=Y[u]-Y[v];
            double l=sqrt(rx*rx+ry*ry);
            double d=D[u][v];
            if(l>d && p>1) continue; // no attractive forces required
            // We add stress for node pair (u,v) if l <= d, i.e.
            // if u and v are "too close", so that such nodes repel;
            // if the nodes are "too distant" then we add stress only
            // if p is 0 or 1, i.e. if there is not a path which will
            // generate an attractive force between these nodes through
            // the topology add-on.
            double d2=d*d;
            double rl=d-l;
            double s=rl*rl/d2;
            stress+=s;
            FILE_LOG(logDEBUG2)<<"s("<<u<<","<<v<<")="<<s;
        }
    }
#endif
    //qDebug() << "stress: " << stress;
    if(m_addSnapStress) {
        double snapstress = computeSnapStress();
        //qDebug() << "ss: " << snapstress;
        stress += snapstress;
    }
    if(m_addGridSnapStress) {
        double gridSnapStress = computeGridSnapStress();
        //qDebug() << "gs: " << gridSnapStress;
        stress += gridSnapStress;
    }
    if(m_addEdgeNodeRepulsion) {
        double edgeNodeRepulsionStress = computeEdgeNodeRepulsionStress();
        stress += edgeNodeRepulsionStress;
    }
    //qDebug() << "plus snap stress: " << stress;
    if(preIteration) {
        if ((*preIteration)()) {
            for(vector<Lock>::iterator l=preIteration->locks.begin();
                    l!=preIteration->locks.end();l++) {
                double dx=l->pos(vpsc::HORIZONTAL)-X[l->getID()], dy=l->pos(vpsc::VERTICAL)-Y[l->getID()];
                double s=10000*(dx*dx+dy*dy);
                stress+=s;
                FILE_LOG(logDEBUG2)<<"d("<<l->getID()<<")="<<s;
            }
        }
    }
    stress += topologyAddon->computeStress();
    if(desiredPositions) {
        for(DesiredPositions::const_iterator p = desiredPositions->begin();
            p!=desiredPositions->end();++p) {
            double dx = X[p->id] - p->x, dy = Y[p->id] - p->y;
            stress+=0.5*p->weight*(dx*dx+dy*dy);
        }
    }
    return stress;
}

double ConstrainedFDLayout::computePureStress() const {
    double stress=0;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double rx=X[u]-X[v], ry=Y[u]-Y[v];
            double l=sqrt(rx*rx+ry*ry);
            double d=D[u][v];
            if(l>d && p>1) continue; // no attractive forces required
            // We add stress for node pair (u,v) if l <= d, i.e.
            // if u and v are "too close", so that such nodes repel;
            // if the nodes are "too distant" then we add stress only
            // if p is 0 or 1, i.e. if there is not a path which will
            // generate an attractive force between these nodes through
            // the topology add-on.
            double d2=d*d;
            double rl=d-l;
            double s=rl*rl/d2;
            stress+=s;
        }
    }
    if(preIteration) {
        if ((*preIteration)()) {
            for(vector<Lock>::iterator l=preIteration->locks.begin();
                    l!=preIteration->locks.end();l++) {
                double dx=l->pos(vpsc::HORIZONTAL)-X[l->getID()], dy=l->pos(vpsc::VERTICAL)-Y[l->getID()];
                double s=10000*(dx*dx+dy*dy);
                stress+=s;
                FILE_LOG(logDEBUG2)<<"d("<<l->getID()<<")="<<s;
            }
        }
    }
    stress += topologyAddon->computeStress();
    if(desiredPositions) {
        for(DesiredPositions::const_iterator p = desiredPositions->begin();
            p!=desiredPositions->end();++p) {
            double dx = X[p->id] - p->x, dy = Y[p->id] - p->y;
            stress+=0.5*p->weight*(dx*dx+dy*dy);
        }
    }
    return stress;
}

// Snap stress chooser:
double ConstrainedFDLayout::computeSnapStress() const {
    switch (m_snapStressFunction) {
    case 1:
        return smoothMStress();
    case 2:
        return linearMStress();
    case 3:
        return dualQuadraticStress();
    case 4:
        return invertedQuadraticStress();
    case 5:
        return quarticStress();
    case 6:
        return smoothVStress();
    case 7:
        return linearVStress();
    case 8:

//#define useTailoredNodeSnapForces
#ifdef  useTailoredNodeSnapForces
        return tailoredQuadUStress();
#else
        return quadUStress();
#endif
    }
}

double ConstrainedFDLayout::computeEdgeNodeRepulsionStress() const {
    double aarStress = 0;
    double sig = m_snap_distance;
    double sig2 = sig*sig;
    foreach(Edge e, edges) {
        // Is e axis-aligned in either dimension?
        unsigned u = e.first, v = e.second;
        double dx=X[u]-X[v], dy=Y[u]-Y[v];
        double eps = 5;
        if (fabs(dx) > eps && fabs(dy) > eps) continue;
        // Let dim be the dimension in which the force acts (so the opposite
        // dimension to that in which e is aligned).
        vpsc::Dim dim = fabs(dx)<=eps?vpsc::HORIZONTAL:vpsc::VERTICAL;
        double ez, el, eh;
        if (dim==vpsc::HORIZONTAL) {
            ez = (X[u]+X[v])/2.0; el = min(Y[u],Y[v]); eh = max(Y[u],Y[v]);
        } else {
            ez = (Y[u]+Y[v])/2.0; el = min(X[u],X[v]); eh = max(X[u],X[v]);
        }
        // Consider each node.
        for(unsigned w=0;w<n;w++) {
            if (w==u || w==v) continue;
            double z=dim==vpsc::HORIZONTAL?X[w]:Y[w];
            double k=dim==vpsc::HORIZONTAL?Y[w]:X[w];
            if (k<=el || eh<=k) continue;
            double r = z - ez;
            if (fabs(r) > sig) continue;
            r = r>=0?r-sig:r+sig;
            aarStress += r*r/sig2;
        }
    }
    aarStress *= 20*m_snap_strength;
    return aarStress;
}

// Grid snap stress; will use Quadratic U-stress:
double ConstrainedFDLayout::computeGridSnapStress() const {
    double stress=0;
    double sig = m_snap_distance;
    double sig2 = sig*sig;
    double W = m_snapGridX, H = m_snapGridY;
    for(unsigned u=0;u<n;u++) {
        double x=X[u], y=Y[u];
        double q,r;
        double dx=0,dy=0;
        double sx=0,sy=0;
#ifndef smallerCoordsWin
        // x-dimension
        r = fabs(modf(x/W, &q));
        if (r!=0.5) { // no "tug of war"
            dx=r<0.5?r*W:(1-r)*W;
            if (dx<=sig) {
                sx = dx*dx/sig2;
            }
        }
        // y-dimension
        r = fabs(modf(y/H, &q));
        if (r!=0.5) {
            dy=r<0.5?r*H:(1-r)*H;
            if (dy<=sig) {
                sy = dy*dy/sig2;
            }
        }
#else
        // x-dimension
        r = fabs(modf(x/W, &q));
        dx=r<=0.5?r*W:(1-r)*W;
        if (dx<=sig) {
            sx = dx*dx/sig2;
        }
        // y-dimension
        r = fabs(modf(y/H, &q));
        dy=r<=0.5?r*H:(1-r)*H;
        if (dy<=sig) {
            sy = dy*dy/sig2;
        }
#endif

        stress+=sx+sy;
    }
    stress = m_snap_strength*stress;
    return stress;
}

// Quadratic U-stress tailored to each pair of nodes' dimensions
double ConstrainedFDLayout::tailoredQuadUStress() const {
    double stress = 0;
    foreach(Edge e, edges) {
        unsigned u = e.first, v = e.second;
        double dx=X[u]-X[v], dy=Y[u]-Y[v];
        vpsc::Rectangle *ru = boundingBoxes.at(u);
        double wu = ru->width(), hu = ru->height();
        vpsc::Rectangle *rv = boundingBoxes.at(v);
        double wv = rv->width(), hv = rv->height();
        double sigX = (wu+wv)/2, sigY = (hu+hv)/2;

        double sx = 0;
        if (-sigX < dx && dx <= sigX) {
            sx = dx*dx/(sigX*sigX);
        }

        double sy = 0;
        if (-sigY < dy && dy <= sigY) {
            sy = dy*dy/(sigY*sigY);
        }

        stress+=sx+sy;
    }
    stress *= m_snap_strength;
    return stress;
}

// Quadratic U-stress:
double ConstrainedFDLayout::quadUStress() const {
    double stress=0;
    double sig = m_snapStressSigma;
    double sig2 = sig*sig;
    double c=0.0; // the constant stress -- experiment with 0 and 1 here!
#ifdef nbrSnap
    foreach(Edge e, edges) {
        unsigned u = e.first, v = e.second;
        double dx=X[u]-X[v], dy=Y[u]-Y[v];

        double sx = 0;
        if (dx <= -sig) {
            sx = c;
        } else if (-sig < dx && dx <= sig) {
            sx = dx*dx/sig2;
        } else { // sig < dx
            sx = c;
        }

        double sy = 0;
        if (dy <= -sig) {
            sy = c;
        } else if (-sig < dy && dy <= sig) {
            sy = dy*dy/sig2;
        } else { // sig < dy
            sy = c;
        }

        stress+=sx+sy;
    }
#else
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double dx=X[u]-X[v], dy=Y[u]-Y[v];

            double sx = 0;
            if (dx < -sig) {
                sx = c;
            } else if (-sig <= dx && dx <= sig) {
                sx = dx*dx/sig2;
            } else { // sig < dx
                sx = c;
            }

            double sy = 0;
            if (dy < -sig) {
                sy = c;
            } else if (-sig <= dy && dy <= sig) {
                sy = dy*dy/sig2;
            } else { // sig < dy
                sy = c;
            }

            stress+=sx+sy;
        }
    }   
#endif
    return m_snapStressBeta*stress;
}

// Smooth V-stress:
// (letter "V" with infinite serif)
// ----    ----
//     \  /
//      \/
double ConstrainedFDLayout::smoothVStress() const {
    // For now we disregard parameter rho.
    double stress=0;
    double a = m_snapStressAlpha;
    double sig = m_snapStressSigma;
    double w=0;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double dx=X[u]-X[v], dy=Y[u]-Y[v];

            double sx = 0;
            if (dx < -sig) {
                sx = 1.0;
            } else if (-sig <= dx && dx < 0) {
                w = (dx+sig)/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sx = 2*w3 - 3*w2 + 1;
            } else if (0 <= dx && dx < sig) {
                w = dx/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sx = -2*w3 + 3*w2;
            } else if (sig <= dx) {
                sx = 1.0;
            }

            double sy = 0;
            if (dy < -sig) {
                sy = 1.0;
            } else if (-sig <= dy && dy < 0) {
                w = (dy+sig)/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sy = 2*w3 - 3*w2 + 1;
            } else if (0 <= dy && dy < sig) {
                w = dy/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sy = -2*w3 + 3*w2;
            } else if (sig <= dy) {
                sy = 1.0;
            }

            stress+=sx+sy;
        }
    }
    return m_snapStressBeta*stress;
}

// Linear V-stress:
// (letter "V" with infinite serif)
// ----    ----
//     \  /
//      \/
double ConstrainedFDLayout::linearVStress() const {
    // For now we disregard parameter rho.
    double stress=0;
    double sig = m_snapStressSigma;
    double w=0;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double dx=X[u]-X[v], dy=Y[u]-Y[v];

            double sx = 0;
            if (dx < -sig) {
                sx = 1.0;
            } else if (-sig <= dx && dx < 0) {
                sx = -dx/sig;
            } else if (0 <= dx && dx < sig) {
                sx = dx/sig;
            } else if (sig <= dx) {
                sx = 1.0;
            }

            double sy = 0;
            if (dy < -sig) {
                sy = 1.0;
            } else if (-sig <= dy && dy < 0) {
                sy = -dy/sig;
            } else if (0 <= dy && dy < sig) {
                sy = dy/sig;
            } else if (sig <= dy) {
                sy = 1.0;
            }

            stress+=sx+sy;
        }
    }
    return m_snapStressBeta*stress;
}

// Smooth M-stress:
double ConstrainedFDLayout::smoothMStress() const {
    // For now we disregard parameter rho.
    double stress=0;
    double a = m_snapStressAlpha;
    double sig = m_snapStressSigma;
    double w=0;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double dx=X[u]-X[v], dy=Y[u]-Y[v];

            double sx = 0;
            if (-sig-a <= dx && dx < -sig) {
                w = (dx - (-sig-a))/a;
                double w2 = w*w;
                double w3 = w2*w;
                sx = -2*w3 + 3*w2;
            } else if (-sig <= dx && dx < 0) {
                w = (dx+sig)/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sx = 2*w3 - 3*w2 + 1;
            } else if (0 <= dx && dx < sig) {
                w = dx/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sx = -2*w3 + 3*w2;
            } else if (sig <= dx && dx <= sig+a) {
                w = (dx-sig)/a;
                double w2 = w*w;
                double w3 = w2*w;
                sx = 2*w3 - 3*w2 + 1;
            }

            double sy = 0;
            if (-sig-a <= dy && dy < -sig) {
                w = (dy - (-sig-a))/a;
                double w2 = w*w;
                double w3 = w2*w;
                sy = -2*w3 + 3*w2;
            } else if (-sig <= dy && dy < 0) {
                w = (dy+sig)/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sy = 2*w3 - 3*w2 + 1;
            } else if (0 <= dy && dy < sig) {
                w = dy/sig;
                double w2 = w*w;
                double w3 = w2*w;
                sy = -2*w3 + 3*w2;
            } else if (sig <= dy && dy <= sig+a) {
                w = (dy-sig)/a;
                double w2 = w*w;
                double w3 = w2*w;
                sy = 2*w3 - 3*w2 + 1;
            }

            stress+=sx+sy;
        }
    }
    return m_snapStressBeta*stress;
}

// Linear M-stress:
double ConstrainedFDLayout::linearMStress() const {
    // For now we disregard parameter rho.
    double stress=0;
    double a = m_snapStressAlpha;
    double sig = m_snapStressSigma;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double dx=X[u]-X[v], dy=Y[u]-Y[v];
            double sx = 0;
            if (-sig-a <= dx && dx < -sig) {
                sx = (dx - (-sig-a))/a;
            } else if (-sig <= dx && dx < 0) {
                sx = -dx/sig;
            } else if (0 <= dx && dx < sig) {
                sx = dx/sig;
            } else if (sig <= dx && dx < sig+a) {
                sx = (sig+a-dx)/a;
            }
            double sy = 0;
            if (-sig-a <= dy && dy < -sig) {
                sy = (dy - (-sig-a))/a;
            } else if (-sig <= dy && dy < 0) {
                sy = -dy/sig;
            } else if (0 <= dy && dy < sig) {
                sy = dy/sig;
            } else if (sig <= dy && dy < sig+a) {
                sy = (sig+a-dy)/a;
            }
            stress+=sx+sy;
        }
    }
    return m_snapStressBeta*stress;
}


// Summed quadratic bumps snap stress:
double ConstrainedFDLayout::dualQuadraticStress() const {
    qDebug() << "dual quadratic stress";
    // For now we disregard parameter rho.
    double stress=0;
    double a = m_snapStressAlpha;
    double sig = m_snapStressSigma;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double dx=X[u]-X[v], dy=Y[u]-Y[v];
            double dxp=dx+sig, dxn=dx-sig;
            double dyp=dy+sig, dyn=dy-sig;
            double exp=1+a*dxp*dxp, exn=1+a*dxn*dxn;
            double eyp=1+a*dyp*dyp, eyn=1+a*dyn*dyn;
            double s = 1/exp + 1/exn + 1/eyp + 1/eyn;
            stress+=s;
        }
    }
    return m_snapStressBeta*stress;
}

// Snap stress with long range force:
double ConstrainedFDLayout::invertedQuadraticStress() const {
    // For now we disregard parameter rho.
    double stress=0;
    double a = m_snapStressAlpha;
    double g = m_snapStressGamma;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double rx=X[u]-X[v], ry=Y[u]-Y[v];
            double rx2=rx*rx, ry2=ry*ry;
            double ex=g+a*rx2, ey=g+a*ry2;
            double s = rx2/ex + ry2/ey;
            stress+=s;
        }
    }
    return a*m_snapStressBeta*stress;
}

// Quartic snap stress:
double ConstrainedFDLayout::quarticStress() const {
    // For now we disregard parameter rho.
    double stress=0;
    double sig2=m_snapStressSigma*m_snapStressSigma;
    for(unsigned u=0;(u + 1)<n;u++) {
        for(unsigned v=u+1;v<n;v++) {
            unsigned short p=G[u][v];
            // no forces between disconnected parts of the graph
            if(p==0) continue;
            double rx=X[u]-X[v], ry=Y[u]-Y[v];
            double rx2=rx*rx, ry2=ry*ry;
            double f=rx2-sig2, g=ry2-sig2;
            double d=1+m_snapStressAlpha*f*f, e=1+m_snapStressAlpha*g*g;
            double s = 1/d + 1/e;
            stress+=s;
        }
    }
    return m_snapStressBeta*stress;
}

void ConstrainedFDLayout::moveBoundingBoxes() {
    for(unsigned i=0;i<n;i++) {
        boundingBoxes[i]->moveCentre(X[i],Y[i]);
    }
}


static const double LIMIT = 100000000;

static void reduceRange(double& val)
{
    val = std::min(val, LIMIT);
    val = std::max(val, -LIMIT);
}

void ConstrainedFDLayout::outputInstanceToSVG(std::string instanceName)
{
    std::string filename;
    if (!instanceName.empty())
    {
        filename = instanceName;
    }
    else
    {
        filename = "libcola-debug";
    }
    filename += ".svg";
    FILE *fp = fopen(filename.c_str(), "w");

    if (fp == NULL)
    {
        return;
    }


    double minX = LIMIT;
    double minY = LIMIT;
    double maxX = -LIMIT;
    double maxY = -LIMIT;

    // Find the bounds of the diagram.
    for (size_t i = 0; i < boundingBoxes.size(); ++i)
    {
        double rMinX = boundingBoxes[i]->getMinX();
        double rMaxX = boundingBoxes[i]->getMaxX();
        double rMinY = boundingBoxes[i]->getMinY();
        double rMaxY = boundingBoxes[i]->getMaxY();
   
        reduceRange(rMinX);
        reduceRange(rMaxX);
        reduceRange(rMinY);
        reduceRange(rMaxY);
        
        if (rMinX > -LIMIT)
        {
            minX = std::min(minX, rMinX);
        }
        if (rMaxX < LIMIT)
        {
            maxX = std::max(maxX,rMaxX);
        }
        if (rMinY > -LIMIT)
        {
            minY = std::min(minY, rMinY);
        }
        if (rMaxY < LIMIT)
        {
            maxY = std::max(maxY, rMaxY);
        }
    }
 
    minX -= 50;
    minY -= 50;
    maxX += 50;
    maxY += 50;

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp, "<svg xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%%\" height=\"100%%\" viewBox=\"%g %g %g %g\">\n", minX, minY, maxX - minX, maxY - minY);

    // Output source code to generate this ConstrainedFDLayout instance.
    fprintf(fp, "<!-- Source code to generate this instance:\n");
    fprintf(fp, "#include <vector>\n");
    fprintf(fp, "#include <utility>\n");
    fprintf(fp, "#include \"libcola/cola.h\"\n");
    fprintf(fp, "using namespace cola;\n");
    fprintf(fp, "int main(void) {\n");
    fprintf(fp, "    CompoundConstraints ccs;\n");
    fprintf(fp, "    std::vector<Edge> es;\n");
    fprintf(fp, "    double defaultEdgeLength=%g;\n", m_idealEdgeLength);
    fprintf(fp, "    std::vector<vpsc::Rectangle*> rs;\n");
    fprintf(fp, "    vpsc::Rectangle *rect = NULL;\n\n");
    for (size_t i = 0; i < boundingBoxes.size(); ++i)
    {
        fprintf(fp, "    rect = new vpsc::Rectangle(%g, %g, %g, %g);\n",
               boundingBoxes[i]->getMinX(), boundingBoxes[i]->getMaxX(),
               boundingBoxes[i]->getMinY(), boundingBoxes[i]->getMaxY());
        fprintf(fp, "    rs.push_back(rect);\n\n");
    }
    
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j =  i + 1; j < n; ++j)
        {
            if (G[i][j] == 1)
            {
                fprintf(fp, "    es.push_back(std::make_pair(%lu, %lu));\n", i, j);
            }
        }
    }
    fprintf(fp, "\n");

    for (cola::CompoundConstraints::iterator c = ccs.begin(); 
            c != ccs.end(); ++c)
    {
        (*c)->printCreationCode(fp);
    }

    fprintf(fp, "    ConstrainedFDLayout alg(rs, es, defaultEdgeLength, %s);\n",
            (m_generateNonOverlapConstraints) ? "true" : "false");
    if (clusterHierarchy)
    {
        clusterHierarchy->printCreationCode(fp);
        fprintf(fp, "    alg.setClusterHierarchy(cluster%llu);\n",
                (unsigned long long) clusterHierarchy);
    }
    fprintf(fp, "    alg.setConstraints(ccs);\n");
    fprintf(fp, "    alg.makeFeasible();\n");
    fprintf(fp, "    alg.run();\n");
    fprintf(fp, "};\n");
    fprintf(fp, "-->\n");

    fprintf(fp, "<g inkscape:groupmode=\"layer\" "
            "inkscape:label=\"Rects\">\n");
    for (size_t i = 0; i < boundingBoxes.size(); ++i)
    {
        double minX = boundingBoxes[i]->getMinX();
        double maxX = boundingBoxes[i]->getMaxX();
        double minY = boundingBoxes[i]->getMinY();
        double maxY = boundingBoxes[i]->getMaxY();
    
        fprintf(fp, "<rect id=\"rect-%u\" x=\"%g\" y=\"%g\" width=\"%g\" "
                "height=\"%g\" style=\"stroke-width: 1px; stroke: black; "
                "fill: blue; fill-opacity: 0.3;\" />\n",
                (unsigned) i, minX, minY, maxX - minX, maxY - minY);
    }
    fprintf(fp, "</g>\n");

    fprintf(fp, "<g inkscape:groupmode=\"layer\" "
            "inkscape:label=\"Edges\">\n");
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j =  i + 1; j < n; ++j)
        {
            if (G[i][j] == 1)
            {
                fprintf(fp, "<path d=\"M %g %g L %g %g\" "
                        "style=\"stroke-width: 1px; stroke: black;\" />\n",
                        boundingBoxes[i]->getCentreX(),
                        boundingBoxes[i]->getCentreY(),
                        boundingBoxes[j]->getCentreX(),
                        boundingBoxes[j]->getCentreY());
            }
        }
    }
    fprintf(fp, "</g>\n");

    fprintf(fp, "</svg>\n");
    fclose(fp);
}


} // namespace cola
