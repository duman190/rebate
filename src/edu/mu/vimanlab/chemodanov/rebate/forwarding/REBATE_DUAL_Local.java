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
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by chemo_000 on 12/28/2015.
 */
public class REBATE_DUAL_Local {
    private double deg = 1;
    private double lambda = 0;
    //private double sigma = 2.16;

    //information stored on the src-dst flow packets
    private double[][] knownUtil;
    private double[][] bidUtil;

    public REBATE_DUAL_Local(int N) {
        knownUtil = new double[N][N];
        bidUtil = new double[N][N];
    }

    public void setDeg(double deg) {
        this.deg = deg;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
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
        //decide what max utilization use to send this packet
        double maxUtil = bidUtil[src][dst] >= knownUtil[src][dst] ? bidUtil[src][dst] : lambda * bidUtil[src][dst] + (1 - lambda) * knownUtil[src][dst];
        knownUtil[src][dst] = maxUtil;
        bidUtil[src][dst] = 0;

        //start Forwarding
        while (path.get(path.size() - 1).getID() != dst && path.size() < maxLength) {
            Node n = path.get(path.size() - 1);
            Node next = null;

            //bid node n max adjacent rebate.util
            double nodeUtilBid = computeAdjacentMaxUtil(n);
            bidUtil[src][dst] = nodeUtilBid > bidUtil[src][dst] ? nodeUtilBid : bidUtil[src][dst];

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
                    double potential = computeNodePotential(c, dstN, t.getGraph().getMinR(), demand, maxUtil, t);
                    ;//-1 / Dist.poincare(c, dstN);
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
                    //headerSize.add(row_collision(visits, t.getGraph().getMinR(), t.getGraph().getMaxR()));
                    headerSize.add(col_row_collision(visits, t.getGraph().getMinR(), t.getGraph().getMaxR(), 16));
                    headerSize.add(col_row_collision(visits, t.getGraph().getMinR(), t.getGraph().getMaxR(), 256));
                    return path;
                }
            }


            path.add(next);
        }

//        if (path.size() >= maxLength)
//            System.out.println("EPGF-LOCAL-HYP DID NOT FOUND DST!!!");//:" + path.size() + " Path " + src + "->" + dst +
//                    " Last Repulsion Proximity=" + lastPRepulsion + "Last Attraction Proximity=" + lastPAttraction);

        headerSize.add(visits.size());
        //headerSize.add(row_collision(visits, t.getGraph().getMinR(), t.getGraph().getMaxR()));
        headerSize.add(col_row_collision(visits, t.getGraph().getMinR(), t.getGraph().getMaxR(), 16));
        headerSize.add(col_row_collision(visits, t.getGraph().getMinR(), t.getGraph().getMaxR(), 64));
        headerSize.add(col_row_collision(visits, t.getGraph().getMinR(), t.getGraph().getMaxR(), 256));
        return path;
    }

    private int col_row_collision(Map<Node, Integer> visits, double minR, double maxR, int bucketNum) {
        int[][] hash_table = new int[bucketNum][bucketNum];
        double r_dif = maxR - minR;
        for (Node n : visits.keySet()) {
            Double r = Dist.poincare(new Complex(0, 0), n.z1) - minR;
            int row = (int) Math.round(r * (bucketNum - 1) / r_dif);
            Double phi = Math.atan2(n.z1.getImaginary(), n.z1.getReal());
            if (phi < 0)
                phi += Math.PI * 2;
            int col = (int) Math.round(phi * (bucketNum - 1) / (2 * Math.PI));
            if (row < 0 || row >= bucketNum || col < 0 || col >= bucketNum )
                System.out.print("Exception");
            hash_table[row][col]++;
        }
        //compute collisions
        int sum = 0;
        for (int i = 0; i < bucketNum; i++)
            for (int j = 0; i < bucketNum; i++)
                sum += hash_table[i][j] > 1 ? (hash_table[i][j] - 1) : 0;
        return sum;
    }

    private int row_collision(Map<Node, Integer> visits, double minR, double maxR) {
        int[] hash_table = new int[256];
        double r_dif = maxR - minR;
        for (Node n : visits.keySet()) {
            Double r = n.z1.abs();//;Dist.poincare(new Complex(0, 0), n.z1) - minR;
            int b = (int) Math.round(r * 7) + 1; //find number of buckets /// r_dif


            Double phi = Math.atan2(n.z1.getImaginary(), n.z1.getReal());
            if (phi < 0)
                phi += Math.PI * 2;

            int key = 0;
            if(b > 1) {
                for (int bi = 1; bi < b; bi++)
                    key += Math.pow(2, bi - 1);
                key += (int) Math.round(phi * Math.pow(2, b - 1) / (2 * Math.PI));
            }
            if (key < 0 || key >= 256 )
                System.out.print("Exception");
            hash_table[key]++;
        }
        //compute collisions
        int sum = 0;
        for (int i = 0; i < 256; i++)
                sum += hash_table[i] > 1 ? (hash_table[i] - 1) : 0;
        return sum;
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
        double sigma = 0.3 * (t.getGraph().getMaxR() - R_min);
        NormalDistribution norm = new NormalDistribution(R_min, sigma);

        double R = norm.inverseCumulativeProbability((1.0+maxUtil)/2);//(t.getGraph().getMaxR() - t.getGraph().getMinR()) * maxUtil + t.getGraph().getMinR();//computeObstacleR2(t);
        if (R > t.getGraph().getMaxR())
            R = t.getGraph().getMaxR();
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

    private double computeAdjacentMaxUtil(Node n) {
        double maxUtil = 0;

        for (Edge e : n.getNeighbors().values()) {
            double util = (e.getBW() - e.getAvBW()) / e.getBW();
            if (util > maxUtil)
                maxUtil = util;
        }

        return maxUtil;
    }
}