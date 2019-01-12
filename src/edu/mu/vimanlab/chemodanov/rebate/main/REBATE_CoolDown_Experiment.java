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

package edu.mu.vimanlab.chemodanov.rebate.main;

import org.apache.commons.math3.distribution.ParetoDistribution;
import edu.mu.vimanlab.chemodanov.rebate.forwarding.*;
import edu.mu.vimanlab.chemodanov.rebate.hyperbolic.GreedyEmbedding;
import edu.mu.vimanlab.chemodanov.rebate.hyperbolic.SpanningTree;
import edu.mu.vimanlab.chemodanov.rebate.graph.Edge;
import edu.mu.vimanlab.chemodanov.rebate.graph.Node;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;
import edu.mu.vimanlab.chemodanov.rebate.util.*;

import java.util.*;

/**
 * Created by dmitriichemodanov on 7/30/17.
 */
public class REBATE_CoolDown_Experiment {

    public static void main(String[] args) {
        String tier1 = "data/tier1_topo.graphml";

        //Experiment settings
        Random r = new Random();
        int N = 286; //num of nodes
        int TTL = 256;//t.getGraph().getNumNodes(); //set packet's ttl policy here
        int[] pairNum = {100}; //
        int trialNum = 100;

        //max/min flow demands
        double minD = 0.01;
        double maxD = 1;

        double[][] RG1_125_HM_util = new double[pairNum.length][trialNum];

        double[][] RL1_050_HM_util = new double[pairNum.length][trialNum];
        double[][] RL1_075_HM_util = new double[pairNum.length][trialNum];
        double[][] RL1_100_HM_util = new double[pairNum.length][trialNum];
        double[][] RL1_125_HM_util = new double[pairNum.length][trialNum];
        double[][] RL1_150_HM_util = new double[pairNum.length][trialNum];

        double[][] RG1_100_GE_util = new double[pairNum.length][trialNum];

        double[][] RL1_050_GE_util = new double[pairNum.length][trialNum];
        double[][] RL1_075_GE_util = new double[pairNum.length][trialNum];
        double[][] RL1_100_GE_util = new double[pairNum.length][trialNum];
        double[][] RL1_125_GE_util = new double[pairNum.length][trialNum];
        double[][] RL1_150_GE_util = new double[pairNum.length][trialNum];

        ///START Experiment
        long startExperiment = System.currentTimeMillis();
        for (int p = 0; p < pairNum.length; p++)
            for (int trial = 0; trial < trialNum; trial++) {
                long startTrial = System.currentTimeMillis();
                System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                System.out.println("+ Experiment with " + pairNum[p] + " src-dst pairs. Trial " + (trial + 1) + " out of " + trialNum);
                System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                System.out.println();

                //recreate rebate.topology and RTs
                Topology t = ImportZooGraph.importZooTopo(tier1);
                System.out.println("+Physical Topology has been generated...");
                t.getGraph().initNeighbors();
                System.out.println("+Neighbors were initialized...");

                //re-assign edge capacities
                //Note comment when bw is stored
                for (Edge e : t.getGraph().getEdgesArray())
                    e.setBW((float) 10);// all physical links have 10 Gbps bandwidth

                //GENERATE PAIRS and FLOW DEMAND FIRST
                List<Pair<Integer, Integer>> pairs = new ArrayList<>();
                List<Pair<Node, Node>> srcDstPairs = new ArrayList<>();
                List<Integer> nodeIds = new ArrayList<>();
                for (Node n : t.getGraph().getNodesArray())
                    if (!n.isFailed() && n.getNeighbors().size() <= 2)
                        nodeIds.add(n.getID());
                int nodesNum = nodeIds.size();
                while (pairs.size() < pairNum[p]) {
                    int src = nodeIds.get(r.nextInt(nodesNum));
                    int dst = nodeIds.get(r.nextInt(nodesNum));
                    if (src != dst) {
                        pairs.add(new Pair<>(src, dst));
                        srcDstPairs.add(new Pair<>(t.getGraph().getNodeFromID(src),
                                t.getGraph().getNodeFromID(dst)));
                    }
                }
                double[] D = new double[srcDstPairs.size()];
                ParetoDistribution pareto = new ParetoDistribution();
                double paretoMax = 0;
                double paretoMin = Double.MAX_VALUE;
                for (int d = 0; d < srcDstPairs.size(); d++) {
                    D[d] = pareto.sample(); //r.nextDouble() * (maxD - minD) + minD;//
                    paretoMax = paretoMax < D[d] ? D[d] : paretoMax;
                    paretoMin = paretoMin > D[d] ? D[d] : paretoMin;
                }
                for (int d = 0; d < srcDstPairs.size(); d++)
                    D[d] = ((D[d] - paretoMin) * (maxD - minD) / (paretoMax - paretoMin) + minD); //(maxD + minD) -
                System.out.println("+" + pairs.size() + " src-dst pairs with pareto demands have been generated...");

                System.out.println("+ Start Weird Traffic Engineering...");
                System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                System.out.println();


                boolean isFlowIncreasing = true;
                double[] currentFlow = new double[D.length];


                /////////////////////////////////////
                ////Hyperbolic Mapping starts here...
                /////////////////////////////////////

                System.out.println("+ Hyperbolic Mapping is done.");
                t.getGraph().initMinMaxHypR();
                System.out.println("+ Max Hyp R=" + t.getGraph().getMaxR() + ", Min Hyp R=" + t.getGraph().getMinR());
                System.out.println();

                isFlowIncreasing = true;
                REBATE_DUAL_Global rebateGlobal = new REBATE_DUAL_Global();
                rebateGlobal.setDeg(1.25);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateGlobal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RG1_125_HM] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RG1_125_HM_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                t.getGraph().releaseAllResources();

                //REBATE Local version 1 goes here...
                isFlowIncreasing = true;
                REBATE_DUAL_Local rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.25);
                rebateLocal.setLambda(0.0);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_050_HM] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_050_HM_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_050_HM_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_050_HM_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.25);
                rebateLocal.setLambda(0.25);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_075_HM] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_075_HM_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_075_HM_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_075_HM_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.25);
                rebateLocal.setLambda(0.5);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_100_HM] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_100_HM_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_100_HM_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_100_HM_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.25);
                rebateLocal.setLambda(0.75);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_125_HM] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_125_HM_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_125_HM_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_125_HM_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.25);
                rebateLocal.setLambda(1.0);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_150_HM] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_150_HM_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_150_HM_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_150_HM_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                System.out.println();

                /////////////////////////
                //Greedy Embedding strategy here!!!
                ////////////////////////
                Node root = t.getGraph().getNodeFromID(32);
                root = SpanningTree.buildSpanningTree(t, root); //MaxUtil
                GreedyEmbedding.greedyEmbedding(t, root);
                System.out.println("+ Greedy Embedding is done.");
                t.getGraph().initMinMaxHypR();
                System.out.println("+Max Hyp R=" + t.getGraph().getMaxR() + ", Min Hyp R=" + t.getGraph().getMinR());
                System.out.println();

                isFlowIncreasing = true;
                rebateGlobal = new REBATE_DUAL_Global();
                rebateGlobal.setDeg(1);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateGlobal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RG1_100_GE] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RG1_100_GE_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RG1_100_GE_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RG1_100_GE_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                //REBATE Local version 1 goes here...
                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.0);
                rebateLocal.setLambda(0.0);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_050_GE] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_050_GE_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_050_GE_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_050_GE_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.0);
                rebateLocal.setLambda(0.25);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_075_GE] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_075_GE_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_075_GE_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_075_GE_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.0);
                rebateLocal.setLambda(0.5);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_100_GE] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_100_GE_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_100_GE_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_100_GE_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.0);
                rebateLocal.setLambda(0.75);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_125_GE] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_125_GE_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_125_GE_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_125_GE_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.0);
                rebateLocal.setLambda(1.0);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                double curFlow = allocateFlowOnPath(rgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[RL1_150_GE] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                RL1_150_GE_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment1/Util/RL1_150_GE_PairNum_" + pairNum[p] + ".txt", t.getGraph());
                Util.exportValue("results/experiment1/Satisfied/RL1_150_GE_PairNum_" + pairNum[p] + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow))/pairNum[p]);
                t.getGraph().releaseAllResources();

                System.out.println();
                System.out.println("+ Total trial=" + trial +" time=" + Double.valueOf(System.currentTimeMillis() - startTrial) / 60000 + " minutes...");
            }
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("+ Total Experiment Time=" + Double.valueOf(System.currentTimeMillis() - startExperiment) / (1000 * 3600) + " hours...");
        System.out.println("+ Experiment Utilization Results...");
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("RG_HM = " + Arrays.deepToString(RG1_125_HM_util) + ";");
        System.out.println();
        System.out.println("RL_000_HM = " + Arrays.deepToString(RL1_050_HM_util) + ";");
        System.out.println("RL_025_HM = " + Arrays.deepToString(RL1_075_HM_util) + ";");
        System.out.println("RL_050_HM = " + Arrays.deepToString(RL1_100_HM_util) + ";");
        System.out.println("RL_075_HM = " + Arrays.deepToString(RL1_125_HM_util) + ";");
        System.out.println("RL_100_HM = " + Arrays.deepToString(RL1_150_HM_util) + ";");
        System.out.println();
        System.out.println("RG_GE = " + Arrays.deepToString(RG1_100_GE_util) + ";");
        System.out.println();
        System.out.println("RL_000_GE = " + Arrays.deepToString(RL1_050_GE_util) + ";");
        System.out.println("RL_025_GE = " + Arrays.deepToString(RL1_075_GE_util) + ";");
        System.out.println("RL_050_GE = " + Arrays.deepToString(RL1_100_GE_util) + ";");
        System.out.println("RL_075_GE = " + Arrays.deepToString(RL1_125_GE_util) + ";");
        System.out.println("RL_100_GE = " + Arrays.deepToString(RL1_150_GE_util) + ";");
    }

    private static int computeNumOfSatisfied(double[] d, double[] cf) {
        int satsfied = 0;
        for (int i = 0; i < d.length; i++)
            satsfied += Math.abs(d[i] - cf[i]) < 0.0000001 ? 1 : 0;
        return satsfied;
    }

    private static double computeMaxUtil(Edge[] edges) {
        double maxUtil = 0;
        for (Edge e : edges) {
            double util = (e.getBW() - e.getAvBW()) / e.getBW();
            if (util > maxUtil)
                maxUtil = util;
        }
        return maxUtil;
    }

    private static int numOfSaturatedEdges(Edge[] edges) {
        int num = 0;
        for (Edge e : edges)
            if ((e.getBW() - e.getAvBW()) / e.getBW() == 1)
                num++;
        return num;
    }

    private static double sum(double[] values) {
        double sum = 0;
        for (double v : values)
            sum += v;
        return sum;
    }

    private static double allocateFlowOnPath(List<Node> path, double d) {
        float minBw = Float.MAX_VALUE;
        Map<Edge, Float> bw = new LinkedHashMap<>();
        Map<Edge, Integer> times = new LinkedHashMap<>();

        for (int i = 0; i < path.size() - 1; i++) {
            Edge e = path.get(i).getNeighbors().get(path.get(i + 1));

            if (bw.containsKey(e)) {
                int k = times.remove(e);
                times.put(e, ++k);
            } else {
                bw.put(e, e.getAvBW());
                times.put(e, 1);
            }
        }

        for (Edge e : bw.keySet()) {
            float avBw = bw.get(e) / times.get(e);
            double bwPacket = 0.01;//Double.valueOf(1500 * 8) / 1000000000; //one packet Bw in Gbps

            if (bwPacket < avBw)
                avBw = (float) bwPacket;

            if (avBw < minBw)
                minBw = avBw;
        }

        //check if a single packet is less than demand;
        minBw = minBw < d ? minBw : (float) d;

        for (int i = 0; i < path.size() - 1; i++) {
            Edge e = path.get(i).getNeighbors().get(path.get(i + 1));
            if (e.getAvBW() - minBw < 0 && Math.abs(e.getAvBW() - minBw) < 0.0000001)
                e.allocateResources(e.getAvBW());
            else if (e.getAvBW() - minBw >= 0)
                e.allocateResources(minBw);
            else
                throw new RuntimeException("Smth wrong with path allocation...");
        }

        return minBw;
    }
}