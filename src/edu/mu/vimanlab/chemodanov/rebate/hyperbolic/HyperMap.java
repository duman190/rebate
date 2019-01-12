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

package edu.mu.vimanlab.chemodanov.rebate.hyperbolic;


import edu.mu.vimanlab.chemodanov.rebate.graph.Node;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;
import org.apache.commons.math3.analysis.function.Acosh;
import org.apache.commons.math3.complex.Complex;

import java.util.*;

/**
 * Created by dmitriichemodanov on 1/12/17.
 */
public class HyperMap {
    public static void hyperbolicMap(Topology t, double m, boolean correctionStep) {
        //sort nodes by their node degree
        Node[] nodes = t.getGraph().getNodesArray();
        Arrays.sort(nodes, Node.DegreeComparator);

        Set<Integer> correctionSteps = new LinkedHashSet<>();
        correctionSteps.add((int) Math.round(0.02 * nodes.length));
        correctionSteps.add((int) Math.round(0.05 * nodes.length));
        correctionSteps.add((int) Math.round(0.2 * nodes.length));
        correctionSteps.add((int) Math.round(0.5 * nodes.length));
        if (correctionStep)
            System.out.println("[HyperMap]: Correction steps are enabled and are:" + correctionSteps);

        if (nodes.length > 0) {
            //compute avd node degree
            int total = 0;
            for (Node n : nodes)
                total += n.getNeighbors().size();

            //other constant
            double avgDeg = Double.valueOf(total) / Double.valueOf(nodes.length);
            double LL = (avgDeg - 2 * m) < 0 ? 0 : (avgDeg - 2 * m)/ 2;
            double yu = 2.03;//.15;//.15;//2.15;// //degree distribution
            double betta = 1.0 / (yu - 1);
            double T = 0.65;//0.6;

            System.out.println("[HyperMap]: Constants are avgDeg=" + avgDeg + " m=" + m
                    + " L=" + LL + " yu=" + yu + " betta=" + betta + " T=" + T);

            //assign (0,0) polar coordinates to root
            nodes[0].ro = 0;
            nodes[0].phi = 0;

            //rebate.main algorithm
            for (int i = 1; i < nodes.length; i++) {
                //step 5 - assign radial coordinate
                nodes[i].ro = 2 * Math.log(i + 1);

                //step 6 - increase radial coordinate of every existing node j<i
                for (int j = 0; j < i; j++)
                    nodes[j].ro = betta * nodes[j].ro + (1 - betta) * nodes[i].ro;

                //step 7 - assign angular coordinate by maximizing L2i
                double angularStep = 1.0 / (i + 1);
                double bestAngle = 0;
                double maxL = -1;
                for (double angle = 0.0; angle <= Math.PI * 2; angle = angle + angularStep) {
                    double L = 1;
                    for (int j = 0; j < i; j++)
                        L *= computeL(nodes[j], nodes[i], i + 1, angle, T, betta, m, LL);
                    if (L > maxL) {
                        maxL = L;
                        bestAngle = angle;
                    }
                }
                nodes[i].phi = bestAngle;

                //step 8 - correction step
                if (correctionStep && correctionSteps.contains(i + 1)) {
                    for (int j = 0; j <= i; j++) {
                        double angularStep2 = 1.0 / (i + 1);
                        double bestAngle2 = 0;
                        double maxL2 = -1;
                        for (double angle = 0.0; angle <= Math.PI * 2; angle = angle + angularStep2) {
                            double L = 1;
                            for (int l = 0; l <= i; l++)
                                if (l != j)
                                    L *= computeL(nodes[l], nodes[j], j + 1, angle, T, betta, m, LL);
                            if (L > maxL2) {
                                maxL2 = L;
                                bestAngle2 = angle;
                            }
                        }
                        nodes[j].phi = bestAngle2;
                    }

                System.out.println("[HyperMap]: Finished with node[" + (i + 1) + "] id=" + nodes[i].getID() + " deg=" + nodes[i].getNeighbors().size());
                }
            }


            for (int i = 0; i < nodes.length; i++) {
                //double norm = maxR + 1;
                nodes[i].z2 = new Complex(nodes[i].ro * Math.cos(nodes[i].phi), nodes[i].ro * Math.sin(nodes[i].phi));//before normalization
                nodes[i].ro /= 2 * Math.log(nodes.length);// + 1;
                nodes[i].z1 = new Complex(nodes[i].ro * Math.cos(nodes[i].phi), nodes[i].ro * Math.sin(nodes[i].phi));
            }
        }
    }

    private static double computeL(Node j, Node i, int ii, double angle, double T, double betta, double m, double LL) {
        double p = computeP(j, i, ii, angle, T, betta, m, LL);
        if (j.getNeighbors().containsKey(i)) {
            return p;
        } else {
            return 1.0 - p;
        }
    }

    private static double computeP(Node j, Node i, int ii, double angle, double T, double betta, double m, double LL) {
        double Qij = Math.PI - Math.abs(Math.PI - Math.abs(j.phi - angle));//angle
        //double xij = i.ro + j.ro + 2 * Math.log(Qij / 2);//
        double xij_new = new Acosh().value(Math.cosh(i.ro) * Math.cosh(j.ro)
                - Math.sinh(i.ro) * Math.sinh(j.ro) * Math.cos(Qij));
        double Ii = (1 - Math.pow(ii, betta - 1)) / (1 - betta);
        double Li = (2 * LL * (1 - betta) * (Math.pow(Double.valueOf(ii + 1)
                / Double.valueOf(ii), 2 * betta - 1) - 1) * (1 - Math.pow(ii, betta - 1)))
                / (Math.pow((1 - Math.pow(ii + 1, betta - 1)), 2) * (2 * betta - 1));
        double mi = m + Li;
        double Ri = i.ro - 2 * Math.log((2 * T * Ii) / (mi * Math.sin(T) * Math.PI));
        double p = 1.0 / (1 + Math.exp((xij_new - Ri) / (2 * T)));
        return p;
    }
}
