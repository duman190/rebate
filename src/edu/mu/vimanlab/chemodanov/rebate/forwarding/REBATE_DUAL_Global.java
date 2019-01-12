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

package edu.mu.vimanlab.chemodanov.rebate.forwarding;


import edu.mu.vimanlab.chemodanov.rebate.graph.*;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;
import edu.mu.vimanlab.chemodanov.rebate.util.Dist;
import org.apache.commons.math3.complex.Complex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by chemo_000 on 12/28/2015.
 */
public class REBATE_DUAL_Global {
    private double deg = 1.0;

    public void setDeg (double deg)
    {
        this.deg = deg;
    }

    public List<Node> potentialGreedyForwarding(int src, int dst, Topology t, int maxLength, double demand) {
        return potentialGreedyForwarding(src, dst, t, maxLength, demand, new ArrayList<Integer>(1));
    }

    public List<Node> potentialGreedyForwarding(int src, int dst, Topology t, int maxLength, double demand, List<Integer> headerSize) {
        Node[] nodes = t.getGraph().getNodesArray();
        List<Node> path = new ArrayList<>();
        path.add(nodes[src]);
        Node dstN = nodes[dst];
        //information stored on the packet
        double lastPRepulsion = Double.MAX_VALUE;
        double lastPAttraction = Double.MAX_VALUE;
        Map<Node, Integer> visits = new HashMap<>(); //keep previously found min hops

        double maxUtil = computeMaxUtil(t);

        //start Forwarding
        while (path.get(path.size() - 1).getID() != dst && path.size() < maxLength) {
            Node n = path.get(path.size() - 1);
            Node next = null;

            ////compute n potential
            double oPotential = computeNodePotential(n, dstN, t.getGraph().getMinR(), demand, maxUtil, t);
            //if last Repulsion potential is greater than for current hop and in Repulsion zone, proceed safely with Repulsion mode
            if (oPotential < lastPRepulsion) {
                //learn potential of current Repulsion mode
                lastPRepulsion = oPotential;
                //compute neighbors potential in Repulsion mode
                double minPotential = oPotential;
                for (Node nh : n.getNeighbors().keySet())
                    if (n.getNeighbors().get(nh).getAvBW() > 0) {
                        double potential = computeNodePotential(nh, dstN, t.getGraph().getMinR(), demand, maxUtil, t);
                        if (potential < minPotential) {
                            next = nh;
                            minPotential = potential;
                        }
                    }
            }
            //compute current dist
            double nDist = Dist.poincare(n, dstN);
            if (next == null && (-1 / nDist) < lastPAttraction) { //Repulsion local minimum, proceed in Attraction mode (i.e., GF mode)
                //learn potential of current Attraction mode
                lastPAttraction = -1 / nDist;
                //compute neighbors potential in Attraction mode
                double minPotential = -1 / nDist;
                for (Node nh : n.getNeighbors().keySet())
                    if (n.getNeighbors().get(nh).getAvBW() > 0) {
                        double potential = -1 / Dist.poincare(nh, dstN);
                        if (potential < minPotential) {
                            next = nh;
                            minPotential = potential;
                        }
                    }
            }
            if (next == null) { //RGF in both Repulsion and Attraction modes is unavailable to route packet proceed in Pressure mode
                int minVisits = Integer.MAX_VALUE;
                //find min visited next hop candidates
                for (Node neighbor : n.getNeighbors().keySet())
                    if (n.getNeighbors().get(neighbor).getAvBW() > 0)
                        if (!visits.containsKey(neighbor)) {
                            minVisits = 0;
                            break;
                        } else if (visits.get(neighbor) < minVisits)
                            minVisits = visits.get(neighbor);
                List<Node> candidates = new ArrayList<>();
                for (Node neighbor : n.getNeighbors().keySet())
                    if (n.getNeighbors().get(neighbor).getAvBW() > 0)
                        if (!visits.containsKey(neighbor) || visits.get(neighbor) == minVisits)
                            candidates.add(neighbor);
                //compute candidates potential in Pressure mode
                double min = Double.MAX_VALUE;
                for (Node c : candidates) {
                    double potential = computeNodePotential(c, dstN, t.getGraph().getMinR(), demand, maxUtil, t);;//-1 / Dist.poincare(c, dstN);
                    if (potential < min) {
                        next = c;
                        min = potential;
                    }
                }
                //increment visit of next node
                if (next != null && !visits.containsKey(next))
                    visits.put(next, 1);
                else if (next != null)
                    visits.put(next, visits.remove(next) + 1);
                else {
//                    System.out.println("EPGF-LOCAL-HYP DID NOT FOUND NEXT HOP!!!:" + path.size() + " Path " + src + "->" + dst +
//                            " Last Repulsion Proximity=" + lastPRepulsion + "Last Attraction Proximity=" + lastPAttraction);
                    headerSize.add(visits.size());
                    return path;
                }
            }

            path.add(next);
        }

//        if (path.size() >= maxLength)
//            System.out.println("EPGF-LOCAL-HYP DID NOT FOUND DST!!!");//:" + path.size() + " Path " + src + "->" + dst +
//                    " Last Repulsion Proximity=" + lastPRepulsion + "Last Attraction Proximity=" + lastPAttraction);

        headerSize.add(visits.size());

        return path;
    }

    private double computeMaxUtil(Topology t) {
        double maxUtil = 0;

        for (Edge e : t.getGraph().getEdgesArray())
            if (maxUtil < (e.getBW() - e.getAvBW()) / e.getBW())
                maxUtil = (e.getBW() - e.getAvBW()) / e.getBW();

        return maxUtil;
    }

/*
    private double computeNodePotential(Node n, Node dst, double R_min, double c_min, double maxUtil) {
        Complex o = new Complex(0, 0); //origin of the poincare model;

        //Compute equipotential circumference radius first
        double d_ot = Dist.poincare(o, dst.z1);
        double y_max = 1 - maxUtil; //1 -
        double B = c_min * (d_ot - R_min) / (c_min + y_max * (d_ot - R_min));
        double R = 0.5 * (B + Math.sqrt(B * B + 4 * d_ot * R_min));


        //Compute XYZ potential - TODO: think on proper name for this part of potential
        double d_io = Dist.poincare(n.z1, o);
        double d_it = Dist.poincare(n, dst);
        double gCos = (d_io * d_io + d_ot * d_ot - d_it - d_it) / (2 * d_io * d_ot);
        double xi = R / Math.sqrt(Math.pow(d_io, 2) * Math.pow(d_ot, 2) + Math.pow(R, 4)
                - 2 * Math.pow(R, 2) * d_io * d_ot * gCos);

        //return node potential
        return -1/d_it + xi;
    }
    */
private double computeNodePotential(Node n, Node dst, double D, double R_min, double maxUtil, Topology t) {
    Complex o = new Complex(0, 0); //origin of the poincare model;

    //Compute equipotential circumference radius first
    double R = computeObstacleR2(t);
    double d_ot = Dist.poincare(o, dst.z1);
    //double R_max = d_ot; //t.getGraph().getMaxR(); //



    double gamma = maxUtil; //1 -
    //double B = c_min * (d_ot - R_min) / (c_min + gamma * (d_ot - R_min));
    //double R = 0.5 * (B + Math.sqrt(B * B + 4 * d_ot * R_min));
    //double a = d_ot + Math.sqrt(Math.abs(d_ot * d_ot - R * R));
    //Compute XYZ potential - TODO: think on proper name for this part of potential
    double d_io = Dist.poincare(n.z1, o);
    //double TwoCos = (d_io * d_io + d_ot * d_ot - d_it * d_it) / (d_io * d_ot);
    //double xi = (2 * a) / (Math.pow(d_io, 2) + Math.pow(a, 2) - d_io * a * TwoCos);
    double Q = R == R_min ? 0 :
            (Math.pow(R_min, deg) * Math.pow(R, deg) * (gamma * (R_min + d_ot) * (R + d_ot) + (R - R_min)))
                    / ((Math.pow(R, deg) - Math.pow(R_min, deg)) * (R_min + d_ot) * (R + d_ot));
    double xi = Q / Math.pow(d_io, deg);

    //return node potential
    double d_it = Dist.poincare(n, dst);
    return -1 / d_it + xi;
}

    private double computeObstacleR2(Topology t) {
        double maxUtil = 0;

        for (Edge e : t.getGraph().getEdgesArray())
            if (maxUtil < (e.getBW() - e.getAvBW()) / e.getBW())
                maxUtil = (e.getBW() - e.getAvBW()) / e.getBW();

        return (t.getGraph().getMaxR() - t.getGraph().getMinR()) * maxUtil + t.getGraph().getMinR();
    }
}