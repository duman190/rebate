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
import edu.mu.vimanlab.chemodanov.rebate.graph.Node;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;
import edu.mu.vimanlab.chemodanov.rebate.util.*;

import java.util.*;

/**
 * Created by dmitriichemodanov on 7/30/17.
 */
public class REBATE_TTL_Experiment {
    public static void main(String[] args) {
        String tier1 = "data/tier1_topo.graphml";

        //Experiment settings
        Random r = new Random();
        int N = 286; //num of nodes
        int[] TTL = {16, 32, 64, 128, 256};

        //Generate pairs
        List<Pair<Integer, Integer>> pairs = new ArrayList<>();
        for (int i = 0; i < N; i++)
            for (int k = i + 1; k < N; k++)
                pairs.add(new Pair<>(i, k));

        //store variables
        double[] RL1_125_HM_delivery = new double[TTL.length];
        int[][] RL1_125_HM_header = new int[TTL.length][pairs.size()];
        double[][] RL1_125_HM_row = new double[TTL.length][pairs.size()];
        double[][] RL1_125_HM_row2 = new double[TTL.length][pairs.size()];
        double[][] RL1_125_HM_row_col = new double[TTL.length][pairs.size()];

        //recreate rebate.topology and RTs
        Topology t = ImportZooGraph.importZooTopo(tier1);
        System.out.println("+Physical Topology has been generated...");
        t.getGraph().initNeighbors();
        System.out.println("+Neighbors were initialized...");
        System.out.println();

        /////////////////////////////////////
        ////Hyperbolic Mapping starts here...
        /////////////////////////////////////
        t.getGraph().initMinMaxHypR();

        ///START Experiment
        long startExperiment = System.currentTimeMillis();
        for (int p = 0; p < TTL.length; p++) {
            long startTrial = System.currentTimeMillis();
            System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            System.out.println("+ Experiment with TTL=" + TTL[p] + " policy for " + pairs.size() + " src-dst pairs has been started.");
            System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

            int deliveryCnt = 0;
            List<Integer> header = new ArrayList<>();

            //REBATE goes here..
            REBATE_DUAL_Local rebateLocal = new REBATE_DUAL_Local(N);
            rebateLocal.setDeg(1.25);
            deliveryCnt = 0;
            for (int i = 0; i < pairs.size(); i++) {
                header = new ArrayList<>();
                Pair<Integer, Integer> srcDst = pairs.get(i);
                List<Node> path = rebateLocal.potentialGreedyForwarding(srcDst.first(), srcDst.second(), t, TTL[p], 0, header);
                if (path.get(path.size() - 1).getID() == srcDst.second()) {
                    deliveryCnt++;
                    RL1_125_HM_header[p][i] = header.get(0);
                    if (RL1_125_HM_header[p][i] > 0) {
                        RL1_125_HM_row[p][i] = Double.valueOf(header.get(1));
                        RL1_125_HM_row2[p][i] = Double.valueOf(header.get(2));
                        RL1_125_HM_row_col[p][i] = Double.valueOf(header.get(3));
                    }
                }
            }
            RL1_125_HM_delivery[p] = Double.valueOf(deliveryCnt) / pairs.size();
            System.out.println("+[RL_HM] TTL =" + TTL[p] + "; SR=" + RL1_125_HM_delivery[p]
                    + "; Avg. HS=" + avg(RL1_125_HM_header[p]) + "; Max. HS=" + max(RL1_125_HM_header[p])
                    + "; Avg. RC=" + avg(RL1_125_HM_row[p]) + "; Max. RC=" + max(RL1_125_HM_row[p])
                    + "; Avg. RC2=" + avg(RL1_125_HM_row2[p]) + "; Max. RC2=" + max(RL1_125_HM_row2[p])
                    + "; Avg. RCC=" + avg(RL1_125_HM_row_col[p]) + "; Max. RCC=" + max(RL1_125_HM_row_col[p]));

            System.out.println();
            System.out.println("+ Total scenario time=" + Double.valueOf(System.currentTimeMillis() - startTrial) / 60000 + " minutes...");

        }
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("+ Total Experiment Time=" + Double.valueOf(System.currentTimeMillis() - startExperiment) / (1000 * 3600) + " hours...");
        System.out.println("+ Experiment TTL/Header Results...");
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        Util.exportMatrix("results/experiment1/RL_HM_header.txt", RL1_125_HM_header);
        Util.exportValues("results/experiment1/RL_HM_8bit.txt", RL1_125_HM_row[TTL.length-1]);
        Util.exportValues("results/experiment1/RL_HM_12bit.txt", RL1_125_HM_row2[TTL.length-1]);
        Util.exportValues("results/experiment1/RL_HM_16bit.txt", RL1_125_HM_row_col[TTL.length-1]);
        System.out.println("RL_HM_delivery=" + Arrays.toString(RL1_125_HM_delivery));
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    }

    private static double avg(int[] n) {
        int sum = 0;
        for (int i = 0; i < n.length; i++)
            sum += n[i];
        return Double.valueOf(sum) / n.length;
    }

    private static double max(int[] n) {
        double max = 0;
        for (int i = 0; i < n.length; i++)
            max = max > n[i] ? max : n[i];
        return max;
    }

    private static double avg(double[] n) {
        int sum = 0;
        for (int i = 0; i < n.length; i++)
            sum += n[i];
        return Double.valueOf(sum) / n.length;
    }

    private static double max(double[] n) {
        double max = 0;
        for (int i = 0; i < n.length; i++)
            max = max > n[i] ? max : n[i];
        return max;
    }
}