/*
 * vim: ts=4 sw=4 et tw=0 wm=0
 *
 * libvpsc - A solver for the problem of Variable Placement with 
 *           Separation Constraints.
 *
 * Copyright (C) 2005-2008  Monash University
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
*/


#ifndef VPSC_CONSTRAINT_H
#define VPSC_CONSTRAINT_H

#include <iostream>
#include <vector>

namespace cola {
class CompoundConstraint;
}

namespace vpsc {

class Variable;
typedef std::vector<Variable *> Variables;

//! @brief A constraint determines a minimum or exact spacing required between
//!        two Variable objects.
//!
class Constraint
{
	friend std::ostream& operator <<(std::ostream &os,const Constraint &c);
public:
    //! @brief Constructs a minimum or exact spacing constraint between two
    //!        Variable objects.
    //!
    //! (left + gap < right)  or (left + gap == right)
    //!
    //! @param[in] left      The left Variable.
    //! @param[in] right     The right Variable.
    //! @param[in] gap       The minimum or exact distance to separate the 
    //!                      variables by.
    //! @param[in] equality  Whether the separation is an exact distance or
    //!                      not.  The default is false.
	Constraint(Variable *left, Variable *right, double gap, 
            bool equality = false);
	~Constraint();
	
    //! @brief The left Variable.
    Variable *left;
    //! @brief The right Variable.
	Variable *right;
    //! @brief The minimum or exact distance to separate the variables by.
	double gap;
	double lm;
	double slack() const;
	long timeStamp;
	bool active;
    //! @brief Whether the separation is an exact distance or not.
	const bool equality;
    //! @brief Denote whether this constraint was unsatisifable (once the VPSC 
    //!        instance has been solved or satisfied).
	bool unsatisfiable;
    // For use with tentative constraints:
    bool tentative;
    bool rejected;
    long tentativeTimestamp;
    bool temporarilyUnsatisfiable;
    cola::CompoundConstraint* compoundOwner;
    unsigned alignedShapeVarIndex;
};

class CompareConstraints {
public:
	bool operator() (Constraint *const &l, Constraint *const &r) const;
};

//! @brief A vector of pointers to Constraint objects.
typedef std::vector<Constraint*> Constraints;

/** @brief Given a set of variables and constraints, returns a modified set
 *         of constraints with all redundant equality constraints removed.
 *
 * VPSC doesn't work well with redundant equality constraints, usually showing
 * them as unsatisfiable.  This function looks for cycles of equality 
 * constraints and removes the redundant ones.
 */
extern Constraints constraintsRemovingRedundantEqualities(
        const Variables& vars, const Constraints& constraints);

}

#endif // VPSC_CONSTRAINT_H
