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

#include <QtGui>

#include <assert.h>

#include "dsbnode.h"
#include "dsbclone.h"
#include "dsbspecies.h"
#include "dsbreaction.h"
#include "dsbbranch.h"

namespace dunnart {

bool DSBNode::s_followTransporters = false;

void DSBNode::setBranchHeadNumber(int n)
{
    m_branchHeadNumber = n;
}

void DSBNode::setBranch(DSBBranch *b)
{
    m_branch = b;
}

DSBBranch *DSBNode::getBranch()
{
    return m_branch;
}

void DSBNode::setPathway(DSBPathway *pw)
{
    m_pathway = pw;
}

DSBPathway *DSBNode::getPathway()
{
    return m_pathway;
}

bool DSBNode::isConnectedTo(DSBNode *other)
{
    DSBClone *cl    = dynamic_cast<DSBClone*>(this);
    DSBReaction *re = dynamic_cast<DSBReaction*>(this);
    if (cl)
    {
        re = dynamic_cast<DSBReaction*>(other);
    }
    else if (re)
    {
        cl = dynamic_cast<DSBClone*>(other);
    }
    else
    {
        qDebug() << "ERROR: Found Node that is neither Clone nor Reaction.";
        return false; // this shouldn't happen!
    }
    if (!cl || !re) { return false; } // can only be connected if one clone and one reaction
    return re->hasCloneAsInputOrOutput(cl);
}

#if 0
void DSBNode::addBranch(DSBBranch *branch)
{
    m_branches.append(branch);
}

void DSBNode::addFork(DSBFork *fork)
{
    m_forks.append(fork);
}
#endif

QList<DSBBranch*> DSBNode::findBranches(
        QList<QString> blacklist, bool forward, bool extended)
{
    QList<QString> seen;
    QList<DSBBranch*> branches = findBranchesRec(seen, blacklist, forward);
    if (extended)
    {
        // Throw away branches of length 1.
        QList<DSBBranch*> keep;
        for (int i = 0; i < branches.size(); i++)
        {
            DSBBranch *b = branches.at(i);
            if (b->nodes.size() > 1) { keep.append(b); }
        }
        branches = keep;
    }
    return branches;
}

DSBBranch *DSBNode::findMergeTarget(
        QList<DSBBranch *> branches, QList<QString> blacklist)
{
    // Will find longest linear branch, or, failing that, longest cycle.

    // Separate branches into linear and cycles, and determine the length
    // of the longest linear branch, the the length of the longest cycle.
    QList<DSBBranch*> linearBranches;
    QList<DSBBranch*> cycles;
    int linL = 0, cycL = 0;
    for (int i = 0; i < branches.size(); i++)
    {
        DSBBranch *b = branches.at(i);
        if (b->cycle) {
            cycles.append(b);
            int s = b->nodes.size();
            cycL = (s > cycL ? s : cycL);
        }
        else {
            linearBranches.append(b);
            int s = b->nodes.size();
            linL = (s > linL ? s : linL);
        }
    }
    QList<DSBBranch*> opts = linL > 0 ? linearBranches : cycles;
    int maxL = linL > 0 ? linL : cycL;
    // Now will chose from among all branches of length maxL in list opts.
    // We will prefer a branch that does not begin with a blacklisted species.

    // Grab first branches of maxL length whose lead nodes
    // a) are not on blacklist
    // b) are on blacklist (in case there is none of the former kind)
    // Note that a reaction node always counts as "not on blacklist".
    DSBBranch *white = 0;
    DSBBranch *black = 0;
    for (int i = 0; i < opts.size(); i++)
    {
        DSBBranch *b = opts.at(i);
        if (b->nodes.size() == maxL)
        {
            DSBNode *n = b->nodes.first();
            DSBClone *cl = dynamic_cast<DSBClone*>(n);
            DSBReaction *reac = dynamic_cast<DSBReaction*>(n);
            if (reac)
            {
                white = b;
                break;
            }
            else if (cl)
            {
                QString name = cl->getSpecies()->getName();
                if (blacklist.contains(name) && black == 0)
                {
                    black = b;
                }
                else
                {
                    white = b;
                    break;
                }
            }
            else
            {
                qDebug() << "ERROR: DSBNode in findMergeTarget which is neither Clone nor Reaction.";
            }
        }
    }
    assert( !(white==0 && black==0) ); // Must have been at least one branch.
    DSBBranch *mergeTarget;
    if (white) { mergeTarget = white; }
    else { mergeTarget = black; }
    return mergeTarget;

    /*
    // Find a longest linear branch and longest cycle.
    int lin = 0; DSBBranch *linB = 0;
    int cyc = 0; DSBBranch *cycB = 0;
    for (int i = 0; i < linearBranches.size(); i++)
    {
        DSBBranch *b = linearBranches.at(i);
        if (b->nodes.size() > lin)
        {
            lin = b->nodes.size();
            linB = b;
        }
    }
    for (int i = 0; i < cycles.size(); i++)
    {
        DSBBranch *b = cycles.at(i);
        if (b->nodes.size() > cyc)
        {
            cyc = b->nodes.size();
            cycB = b;
        }
    }
    DSBBranch *mergeTarget;
    if (linB) { mergeTarget = linB; }
    else { mergeTarget = cycB; }
    return mergeTarget;
    */
}

QList<DSBBranch*> DSBNode::mergeSelfWithBranches(
        QList<DSBBranch*> branches, QList<QString> blacklist)
{
    // Were there no branches?
    if (branches.isEmpty())
    {
        // Then return one branch, containing just self.
        DSBBranch *b = new DSBBranch;
        b->nodes.append(this);
        branches.append(b);
    }
    // Or was there exactly one branch?
    else if (branches.size() == 1)
    {
        // Then must merge self with the one branch.
        branches.first()->nodes.prepend(this);
        // And set self as parent.
        branches.first()->parent = this;
    }
    else
    {
        // Otherwise there were two or more branches.
        // Choose merge target.
        DSBBranch *mergeTarget = findMergeTarget(branches, blacklist);
        // Merge.
        mergeTarget->nodes.prepend(this);
        // And set self as parent of that branch.
        mergeTarget->parent = this;
        // And set self as parent of all branches that do not already have one,
        // or do have one but it is the same as their first node.
        for (int i = 0; i < branches.size(); i++)
        {
            DSBBranch *b = branches.at(i);
            if (!b->parent) { b->parent = this; }
            else
            {
                DSBNode *p = b->parent;
                DSBNode *f = b->nodes.first();
                if (p == f) { b->parent = this; }
            }
        }
    }

    return branches;
}

}
















