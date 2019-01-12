/****************************************************************************/
/* This file is part of REBATE project.                                     */
/*                                                                          */
/* REBATE is free software: you can redistribute it and/or modify           */
/* it under the terms of the GNU General Public License as published by     */
/* the Free Software Foundation, either version 3 of the License, or        */
/* (at your option) any later version.                                      */
/*                                                                          */
/* REBATE is distributed in the hope that it will be useful,                */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU General Public License for more details.                             */
/*                                                                          */
/* You should have received a copy of the GNU General Public License        */
/* along with REBATE.  If not, see <http://www.gnu.org/licenses/>.          */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*  Author:    Dmitrii Chemodanov, University of Missouri-Columbia          */
/*  Title:     REBATE: a REpulsive-BAsed Traffic Engineering Approach       */
/*                        for Dynamic Scale-Free Networks                   */
/*  Revision:  1.0         1/12/2019                                        */
/****************************************************************************/

package edu.mu.vimanlab.chemodanov.rebate.util;


import edu.mu.vimanlab.chemodanov.rebate.graph.*;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;
import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.Map;

/**
 * Created by dmitriichemodanov on 2/6/17.
 */
public class MinMaxUtil {

    public double solveMinMaxUtil(List<Pair<Node, Node>> srcDstPairs, double[] D, Topology t, boolean allocate) {
        //solve final primal problem to obtain optimal mapping for one-shot VNE
        double minMaxUtil = solvePrimal(srcDstPairs, D, t.getGraph(), allocate);

        return minMaxUtil;
    }

    private double solvePrimal(List<Pair<Node, Node>> srcDstPairs, double[] D, Graph g, boolean allocate) {
        //construct and solve LP-P to get final (near) optimal mapping
        try {
            IloCplex cplex = new IloCplex();

            IloNumVar[][] var = new IloNumVar[3][]; //we have 3 types of vars, i.e., per path and per vnode-pnode mapping
            IloRange[][] rng = new IloRange[2][]; //we have 8 constraints

            constructLP_P(cplex, var, rng, srcDstPairs, D, g); //construct LP-P with int vars and solve for (near) optimal mapping

            cplex.setOut(null); // turn off logging
            if (cplex.solve()) // if cplex found LP-P solution continue, otherwise reject VNR
            {
                double minMaxUtil = cplex.getObjValue();
                if (allocate) {
                    //allocate flow
                    double[] f = cplex.getValues(var[0]); //flows in initial directions
                    double[] fr = cplex.getValues(var[1]); //flows in reverse directions
                    for (int pairI = 0; pairI < srcDstPairs.size(); pairI++)
                        for (Edge pe : g.getEdgesArray()) {
                            float bwToAllocate = (float) (f[pairI * g.getNumEdges() + pe.getID()]
                                    + fr[pairI * g.getNumEdges() + pe.getID()]);
                            if (bwToAllocate > pe.getAvBW())
                                bwToAllocate = pe.getAvBW();
                            pe.allocateResources(bwToAllocate);
                        }
                }
                cplex.end();
                return minMaxUtil;
            } else
                throw new RuntimeException("[MinMaxUtil] failed while trying to solve LP");
        } catch (IloException e) {
            System.err.println("Concert exception caught '" + e + "' caught");
        }
        return 0.0; //never gets here
    }

    private void constructLP_P(IloCplex model,
                               IloNumVar[][] var,
                               IloRange[][] rng, List<Pair<Node, Node>> srcDstPairs, double[] D, Graph g) throws IloException {
        // Vars
        // First define the variables and obj values for each vnode
        int pairEdgeNum = g.getNumEdges() * srcDstPairs.size();
        //first var type is per each path available (known so far)
        double[] f_lb = new double[pairEdgeNum];
        double[] f_ub = new double[pairEdgeNum];
        IloNumVarType[] f_t = new IloNumVarType[pairEdgeNum];
        String[] f_names = new String[pairEdgeNum];
        String[] fr_names = new String[pairEdgeNum];
        for (int i = 0; i < pairEdgeNum; i++) {
            Pair<Node, Node> pair = srcDstPairs.get(i / g.getNumEdges());
            Edge pe = g.getEdgeFromID(i % g.getNumEdges());
            f_lb[i] = 0;
            f_ub[i] = D[i / g.getNumEdges()];
            f_t[i] = IloNumVarType.Float;
        }
        IloNumVar[] f = model.numVarArray(pairEdgeNum, f_lb, f_ub, f_t);
        IloNumVar[] fr = model.numVarArray(pairEdgeNum, f_lb, f_ub, f_t);
        var[0] = f;
        var[1] = fr;

        IloNumVar[] t = new IloNumVar[1];
        t[0] = model.numVar(0, Double.POSITIVE_INFINITY); // var to track maximum link utilization
        var[2] = t;

        // Objective Function:  maximize flow over all src-dst pairs
        model.addMinimize(model.sum(var[2]));

        // Constraints
        // Physical Edge Capacity  Constraints
        rng[0] = new IloRange[g.getNumEdges()];
        int pl = 0;
        for (Edge pedge : g.getEdgesArray()) {
            IloNumExpr[] pairMappings = new IloNumExpr[srcDstPairs.size()];

            for (int pairI = 0; pairI < srcDstPairs.size(); pairI++)
                pairMappings[pairI] = model.sum(model.prod(-1.0 / pedge.getAvBW(), var[0][pairI * g.getNumEdges() + pedge.getID()]),
                        model.prod(-1.0 / pedge.getAvBW(), var[1][pairI * g.getNumEdges() + pedge.getID()]));

            rng[0][pl] = model.addGe(model.sum(model.sum(pairMappings), model.prod(1.0, var[2][0])), 0);
            pl++;
        }
        //Flow Conservation Constraints
        int pairNum = srcDstPairs.size();
        int pairNodeNum = g.getNumNodes() * pairNum;
        rng[1] = new IloRange[pairNodeNum];
        for (int pairI = 0; pairI < srcDstPairs.size(); pairI++) {
            Node pairSrc = srcDstPairs.get(pairI).first();
            Node pairDst = srcDstPairs.get(pairI).second();
            //intermediate node constraints
            for (Node pnode : g.getNodesArray()) {
                //check if pnode is src, dst or none of these
                double flowD = pnode.equals(pairSrc) ? D[pairI] : pnode.equals(pairDst) ? -D[pairI] : 0;
                IloNumExpr[] inFlow = new IloNumExpr[pnode.getNeighbors().size()];
                IloNumExpr[] outFlow = new IloNumExpr[pnode.getNeighbors().size()];
                int nhi = 0;
                for (Map.Entry<Node, Edge> neighbor : pnode.getNeighbors().entrySet()) {
                    Edge pedge = neighbor.getValue();
                    int dir = pedge.getSrc().equals(pnode) ? 0 : 1; // if source than var 0 otherwise var 1
                    inFlow[nhi] = model.prod(-1, var[1 - dir][pairI * g.getNumEdges() + pedge.getID()]);
                    outFlow[nhi] = model.prod(1, var[dir][pairI * g.getNumEdges() + pedge.getID()]);
                    nhi++;
                }
                rng[1][pairI * g.getNumNodes() + pnode.getID()] = model.addEq(model.sum(model.sum(inFlow), model.sum(outFlow)), flowD);
            }
        }
    }
}