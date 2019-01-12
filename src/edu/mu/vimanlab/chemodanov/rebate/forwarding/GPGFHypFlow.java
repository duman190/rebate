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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by chemo_000 on 1/9/2016.
 */
public class GPGFHypFlow {
    public List<Node> greedyForwarding(int src, int dst, Topology t, int maxLength) {
        return greedyForwarding(src, dst, t, maxLength, new ArrayList<Integer>(1));
    }

    public List<Node> greedyForwarding(int src, int dst, Topology t, int maxLength, List<Integer> headerSize) {
        Node[] nodes = t.getGraph().getNodesArray();
        List<Node> path = new ArrayList<>();
        path.add(nodes[src]);
        Node dstN = nodes[dst];

        //information stored on packet
        Map<Node, Integer> visits = new HashMap<>();
        double lastProximity = Double.MAX_VALUE;

        while (path.get(path.size() - 1).getID() != dst && path.size() < maxLength) {
            Node n = path.get(path.size() - 1);
            Node next = null;
            //compute current dist
            double nDist = Dist.poincare(n.z1, dstN.z1);
            if (nDist < lastProximity) { //Proceed in Greedy Forwarding mode
                //learn potential of current Attraction mode
                lastProximity = nDist;
                //compute neighbors potential in Attraction mode
                double minDist = nDist;
                for (Node neighbor : n.getNeighbors().keySet())
                    if (n.getNeighbors().get(neighbor).getAvBW() > 0)
                    {
                        double dist = Dist.poincare(neighbor.z1, dstN.z1);
                        if (dist < minDist) {
                            next = neighbor;
                            minDist = dist;
                        }
                    }
            }
            if (next == null) { //GPGF in Greedy mode is unavailable to route packet proceed in Pressure mode
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
                //compute candidates distance in Pressure mode
                double min = Double.MAX_VALUE;
                for (Node candidate : candidates) {
                    double dist = Dist.poincare(candidate.z1, dstN.z1);
                    if (dist < min) {
                        next = candidate;
                        min = dist;
                    }
                }
//                if (path.size()>64)
//                    System.out.print("Do smth");
                //increment visit of next node
                if (next != null && !visits.containsKey(next))
                    visits.put(next, 1);
                else if (next != null)
                    visits.put(next, visits.remove(next) + 1);
                else {
//                    System.out.println("GPGF-Hyp DID NOT FOUND NEXT HOP!!!:" + path.size() + " Path " + src + "->" + dst +
//                            " Last Geo Proximity=" + lastProximity);
                    headerSize.add(visits.size());
                    return path;
                }
            }

//            if(next == null)
//                return path;

            path.add(next);
        }

//        if (path.size() >= maxLength)
//            System.out.println("GPGF-Hyp DID NOT FOUND DST!!!");
        headerSize.add(visits.size());

        return path;

    }
}
