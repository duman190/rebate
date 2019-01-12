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


import edu.mu.vimanlab.chemodanov.rebate.graph.*;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;

import java.util.*;

/**
 * Created by dmitriichemodanov on 1/9/17.
 */
public class SpanningTree {
    /**
     * Spanning tree based on topological order -> minimizing number of hops
     * @param t
     */
    public static Node buildSpanningTree(Topology t, Node root)
    {
        //clear tree info
        for (Node n : t.getGraph().getNodesArray())
        {
            n.depth = 0;
            n.parent = null;
            n.inTree = false;
        }

        root.inTree = true;
        Queue<Node> q = new LinkedList<>();
        q.add(root);
        int maxRootDepth = 0;
        while (!q.isEmpty())
        {
            Node u = q.poll();
            if(u.depth > maxRootDepth)
                maxRootDepth = u.depth;
            for (Map.Entry<Node, Edge> e : u.getNeighbors().entrySet())
            {
                Node v = e.getKey();
                if (!v.inTree)
                {
                    v.parent=u;
                    v.depth = u.depth + 1;
                    v.inTree=true;
                    q.add(v);
                }
            }
        }
        return root;
    }
}
