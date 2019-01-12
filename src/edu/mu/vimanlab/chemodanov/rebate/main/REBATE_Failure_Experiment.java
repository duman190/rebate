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

import edu.mu.vimanlab.chemodanov.rebate.forwarding.*;
import edu.mu.vimanlab.chemodanov.rebate.hyperbolic.GreedyEmbedding;
import edu.mu.vimanlab.chemodanov.rebate.hyperbolic.SpanningTree;
import edu.mu.vimanlab.chemodanov.rebate.graph.Edge;
import edu.mu.vimanlab.chemodanov.rebate.graph.Node;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;
import edu.mu.vimanlab.chemodanov.rebate.util.*;
import org.apache.commons.math3.distribution.ParetoDistribution;

import java.util.*;

/**
 * Created by dmitriichemodanov on 7/30/17.
 */
public class REBATE_Failure_Experiment {

    public static void main(String[] args) {
        String tier1 = "data/tier1_topo.graphml";

        //Experiment settings
        Random r = new Random();
        int N = 286; //num of nodes
        int TTL = 256;//set packet's ttl policy here
        int pairNum = 1000; //1000
        int trialNum = 50; //100
        double[] failure = {0, 0.1, 0.2, 0.3};

        //max/min flow demands
        double minD = 0.01;
        double maxD = 1;

        //store variables
        double[][] GP_HM_util = new double[failure.length][trialNum];
        double[][] GP_GE_util = new double[failure.length][trialNum];
        double[][] RL1_125_HM_util = new double[failure.length][trialNum];
        double[][] RL1_100_GE_util = new double[failure.length][trialNum];

        Set<Integer> failedEdgeIds = new HashSet<>();

        ///START Experiment
        long startExperiment = System.currentTimeMillis();
        for (int trial = 0; trial < trialNum; trial++) {
            //clear previous failures
            failedEdgeIds.clear();

            //GENERATE PAIRS and FLOW DEMAND FIRST
            List<Pair<Integer, Integer>> pairs = new ArrayList<>();
            List<Integer> nodeIds = new ArrayList<>();
            for (int n = 0; n < N; n++)
                nodeIds.add(n);
            int nodesNum = nodeIds.size();
            while (pairs.size() < pairNum) {
                int src = nodeIds.get(r.nextInt(nodesNum));
                int dst = nodeIds.get(r.nextInt(nodesNum));
                if (src != dst) {
                    pairs.add(new Pair<>(src, dst));
                }
            }
            double[] D = new double[pairs.size()];
            ParetoDistribution pareto = new ParetoDistribution();
            double paretoMax = 0;
            double paretoMin = Double.MAX_VALUE;
            for (int d = 0; d < pairs.size(); d++) {
                D[d] = pareto.sample();
                paretoMax = paretoMax < D[d] ? D[d] : paretoMax;
                paretoMin = paretoMin > D[d] ? D[d] : paretoMin;
            }
            for (int d = 0; d < pairs.size(); d++)
                D[d] = ((D[d] - paretoMin) * (maxD - minD) / (paretoMax - paretoMin) + minD); //(maxD + minD) -
            System.out.println("+" + pairs.size() + " src-dst pairs with pareto demands have been generated...");


            for (int p = 0; p < failure.length; p++) {

                long startTrial = System.currentTimeMillis();
                System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                System.out.println("+ Experiment with " + failure[p] * 100 + "% failures. Trial " + (trial + 1) + " out of " + trialNum);
                System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                System.out.println();

                //recreate rebate.topology and RTs
                Topology t = ImportZooGraph.importZooTopo(tier1);
                System.out.println("+Physical Topology has been generated...");
                t.getGraph().initNeighbors();
                System.out.println("+Neighbors were initialized...");

                //generate failures
                while (failedEdgeIds.size() < failure[p] * t.getGraph().getNumEdges())
                    failedEdgeIds.add(r.nextInt(t.getGraph().getNumEdges()));

                //re-assign edge capacities
                //Note comment when bw is stored
                for (Edge e : t.getGraph().getEdgesArray())
                    if (failedEdgeIds.contains(e.getID()))
                        e.setBW(0);
                    else
                        e.setBW((float) 10);// all physical links have 10 Gbps bandwidth

                //reassign src-dst pairs
                List<Pair<Node, Node>> srcDstPairs = new ArrayList<>();
                for (Pair<Integer, Integer> pair : pairs)
                    srcDstPairs.add(new Pair<>(t.getGraph().getNodeFromID(pair.first()),
                            t.getGraph().getNodeFromID(pair.second())));


                System.out.println("+ Start Weird Traffic Engineering...");
                System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                System.out.println();

                /////////////////////////////////////
                ////Hyperbolic Mapping starts here...
                /////////////////////////////////////

                System.out.println("+ Hyperbolic Mapping is done.");
                t.getGraph().initMinMaxHypR();
                System.out.println("+ Max Hyp R=" + t.getGraph().getMaxR() + ", Min Hyp R=" + t.getGraph().getMinR());
                System.out.println();

                //start randomly allocating with gpgf
                List<Integer> pathLengths = new LinkedList<>();
                boolean isFlowIncreasing = true;
                GPGFHypFlow gpgfFlow = new GPGFHypFlow();
                double[] currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> gpgfPath = gpgfFlow.greedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL);
                            if (gpgfPath.get(gpgfPath.size() - 1).equals(srcDst.second())) {
                                pathLengths.add(gpgfPath.size() - 1);
                                double curFlow = allocateFlowOnPath(gpgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[GP_HM] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                GP_HM_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment_failure1/Util/GP_HM_PairNum_" + p + ".txt", t.getGraph());
                Util.exportValue("results/experiment_failure1/Satisfied/GP_HM_PairNum_" + p + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow)) / pairNum);
                Util.exportValues("results/experiment_failure1/Length/GP_HM_PairNum_" + p + ".txt", pathLengths);
                t.getGraph().releaseAllResources();

                pathLengths.clear();
                isFlowIncreasing = true;
                REBATE_DUAL_Local rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1.25);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                pathLengths.add(rgfPath.size() - 1);
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
                Util.exportEdgesUtil("results/experiment_failure1/Util/RL1_125_HM_PairNum_" + p + ".txt", t.getGraph());
                Util.exportValue("results/experiment_failure1/Satisfied/RL1_125_HM_PairNum_" + p + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow)) / pairNum);
                Util.exportValues("results/experiment_failure1/Length/RL1_125_HM_PairNum_" + p + ".txt", pathLengths);
                t.getGraph().releaseAllResources();

                /////////////////////////
                //Greedy Embedding strategy here!!!
                ////////////////////////
                Node root = t.getGraph().getNodeFromID(32);
                root = SpanningTree.buildSpanningTree(t, root); //MaxUtil
                GreedyEmbedding.greedyEmbedding(t, root);
                System.out.println();
                System.out.println("+ Greedy Embedding is done.");
                t.getGraph().initMinMaxHypR();
                System.out.println("+Max Hyp R=" + t.getGraph().getMaxR() + ", Min Hyp R=" + t.getGraph().getMinR());
                System.out.println();

                //start randomly allocating with gpgf
                pathLengths.clear();
                isFlowIncreasing = true;
                gpgfFlow = new GPGFHypFlow();
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> gpgfPath = gpgfFlow.greedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL);
                            if (gpgfPath.get(gpgfPath.size() - 1).equals(srcDst.second())) {
                                pathLengths.add(gpgfPath.size() - 1);
                                double curFlow = allocateFlowOnPath(gpgfPath, D[i] - currentFlow[i]);
                                if (curFlow > 0)
                                    isFlowIncreasing = true;
                                currentFlow[i] += curFlow;
                            }
                        }
                }
                System.out.println("+[GP_GE] min of the max edge utilization of " + computeMaxUtil(t.getGraph().getEdgesArray())
                        + " has been achieved; demands of " + computeNumOfSatisfied(D, currentFlow)
                        + " pairs have been satisfied in total. Num of fully utilized edges=" + numOfSaturatedEdges(t.getGraph().getEdgesArray())
                        + ". Total throughput=" + sum(currentFlow));
                GP_GE_util[p][trial] = computeMaxUtil(t.getGraph().getEdgesArray());
                Util.exportEdgesUtil("results/experiment_failure1/Util/GP_GE_PairNum_" + p + ".txt", t.getGraph());
                Util.exportValue("results/experiment_failure1/Satisfied/GP_GE_PairNum_" + p + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow)) / pairNum);
                Util.exportValues("results/experiment_failure1/Length/GP_GE_PairNum_" + p + ".txt", pathLengths);
                t.getGraph().releaseAllResources();

                pathLengths.clear();
                isFlowIncreasing = true;
                rebateLocal = new REBATE_DUAL_Local(N);
                rebateLocal.setDeg(1);
                currentFlow = new double[D.length];
                while (isFlowIncreasing) {
                    isFlowIncreasing = false;
                    for (int i = 0; i < D.length; i++)
                        if (currentFlow[i] < D[i]) {
                            Pair<Node, Node> srcDst = srcDstPairs.get(i);
                            List<Node> rgfPath = rebateLocal.potentialGreedyForwarding(srcDst.first().getID(), srcDst.second().getID(), t, TTL, D[i]);
                            if (rgfPath.get(rgfPath.size() - 1).equals(srcDst.second())) {
                                pathLengths.add(rgfPath.size() - 1);
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
                Util.exportEdgesUtil("results/experiment_failure1/Util/RL1_100_GE_PairNum_" + p + ".txt", t.getGraph());
                Util.exportValue("results/experiment_failure1/Satisfied/RL1_100_GE_PairNum_" + p + ".txt", Double.valueOf(computeNumOfSatisfied(D, currentFlow)) / pairNum);
                Util.exportValues("results/experiment_failure1/Length/RL1_100_GE_PairNum_" + p + ".txt", pathLengths);
                t.getGraph().releaseAllResources();

                System.out.println();
                System.out.println("+ Total trial=" + trial + " time=" + Double.valueOf(System.currentTimeMillis() - startTrial) / 60000 + " minutes...");
            }
        }
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("+ Total Experiment Time=" + Double.valueOf(System.currentTimeMillis() - startExperiment) / (1000 * 3600) + " hours...");
        System.out.println("+ Experiment Utilization Results...");
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("GP_HM=" + matlabDeepToString(GP_HM_util));
        System.out.println();
        System.out.println("RL1_125_HM=" + matlabDeepToString(RL1_125_HM_util));
        System.out.println();
        System.out.println("GP_GE=" + matlabDeepToString(GP_GE_util));
        System.out.println();
        System.out.println("RL1_100_GE=" + matlabDeepToString(RL1_100_GE_util));
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

    public static String matlabDeepToString(double[][] a) {
        StringBuilder str = new StringBuilder("[");
        for (int i = 0; i < a.length; i++) {
            str.append("[");
            for (int j = 0; j < a[i].length - 1; j++)
                str.append(a[i][j]).append(" ");
            str.append(a[i][a[i].length - 1]).append("];");
        }
        str.append("];");

        return str.toString();
    }
}