/*
 * vim: ts=4 sw=4 et tw=0 wm=0
 *
 * libcola - A library providing force-directed network layout using the 
 *           stress-majorization method subject to separation constraints.
 *
 * Copyright (C) 2006-2008  Monash University
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
*/

#ifndef CYCLE_DETECTOR_H
#define CYCLE_DETECTOR_H

#include <map>
#include <vector>
#include <stack>

#include "libcola/cola.h"


namespace cycle_detector  {
  typedef std::vector<cola::Edge> Edges;
  typedef unsigned TimeStamp;
  typedef std::vector<bool> CyclicEdges;

  class Node  {
    public:
      enum StatusType { NotVisited, BeingVisited, DoneVisiting };

      unsigned id;
      TimeStamp stamp;
      Node *cyclicAncestor;
      std::vector<unsigned> dests;
      StatusType status;

      Node(unsigned id)  { this->id = id; cyclicAncestor = NULL; status = NotVisited; }
      ~Node()  {}
  };

  class CycleDetector  {
    public:
      CycleDetector(unsigned numVertices, Edges *edges);
      ~CycleDetector();
      std::vector<bool> *detect_cycles();
      void mod_graph(unsigned numVertices, Edges *edges);
      unsigned getV()  { return this->V; }
      Edges *getEdges()  { return this->edges; }
      std::vector<TimeStamp> *getTimeStamps();
      Node *get_highest_ca(Node *n);

    private:
      // attributes
      unsigned V;
      Edges *edges;

      // internally used variables.
      std::vector<Node *> *nodes; // the nodes in the graph
      std::vector<bool> *cyclicEdgesMapping; // the cyclic edges in the graph.
      std::vector<unsigned> traverse; // nodes still left to visit in the graph

      // internally used methods
      void make_matrix();
      void visit(unsigned k);
      bool isSink(Node *node);
      bool find_node(std::vector<Node *> *& list, unsigned k);
      std::pair< bool, std::vector<unsigned>::iterator > find_node(std::vector<unsigned>& list, unsigned k);
  };
}
#endif
