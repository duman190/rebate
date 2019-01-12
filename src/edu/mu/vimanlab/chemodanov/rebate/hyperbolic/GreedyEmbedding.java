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
import org.apache.commons.math3.complex.Complex;
import org.apfloat.Apcomplex;
import org.apfloat.ApcomplexMath;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;

import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;

/**
 * Created by dmitriichemodanov on 1/9/17.
 */
public class GreedyEmbedding {

    public static void greedyEmbedding(Topology t, Node root) {
        int accuracy = 24;
        Apfloat two = new Apfloat(2);
        Node[] nodes = t.getGraph().getNodesArray();
        Queue<Node> q = new LinkedList<>();
        q.add(root);
        root.z = new Apcomplex(new Apfloat("-0.5", accuracy), new Apfloat("-0.1", accuracy)); //root coordinate in Poincare Disk Model
        root.alpha = ApfloatMath.pi(accuracy);
        root.betta = ApfloatMath.pi(accuracy).multiply(two);

        while (!q.isEmpty()) {
            Node parent = q.poll();
//            System.out.println(String.format("Node id=" + parent.getID() + " p="
//                    + (parent.parent == null ? "null" : parent.parent.getID()) + " z=" + parent.z));
            for (Map.Entry<Node, Edge> e : parent.getNeighbors().entrySet()) {
                Node child = e.getKey();
                if (child.parent == parent) {
                    //decide on angles
                    child.alpha = parent.alpha;
                    child.betta = parent.alpha.add(parent.betta);
                    child.betta = child.betta.divide(two);
                    parent.alpha = child.betta;
                    //compute coordinate and update angle
                    child.z = computeComplexCoordinate(child, parent, two);
                    //convert to polar coordinates
                    child.ro0 = ApcomplexMath.abs(child.z);
                    child.phi0 = ApfloatMath.atan2(child.z.imag(), child.z.real());

                    child.alpha = child.alpha.add(child.betta);
                    child.alpha = child.alpha.divide(two);
                    q.add(child);
                }
            }
        }

        //correction step

        Apfloat maxR = Apfloat.ZERO;
        int wrongCnt = 0;
        for (Node n : nodes) {
            if (ApcomplexMath.abs(n.z).doubleValue() > maxR.doubleValue())
                maxR = ApcomplexMath.abs(n.z);
            if (ApcomplexMath.abs(n.z).doubleValue() > 1)
                wrongCnt++;
        }
//        System.out.println("+ total of " + wrongCnt +" nodes have wrong precise hyper coordinates! MaxR=" + maxR);

        //convert to apache complex
        for(Node n : t.getGraph().getNodesArray())
            n.z1 = new Complex(n.z.real().doubleValue(), n.z.imag().doubleValue());

        double maxR1 = 0;
        double wrongCnt1 = 0;
        for (Node n : t.getGraph().getNodesArray()){
            if (n.z1.abs() > maxR1)
                maxR1 = n.z1.abs();
            if (n.z1.abs() >= 1)
                wrongCnt1++;
        }
//        System.out.println("+ total of " + wrongCnt1 +" nodes have wrong double hyper coordinates! MaxR=" + maxR1);
        if (maxR1 >= 1) {
            double norm = maxR1 + 0.0000000000000001;//Double.MIN_NORMAL; //0.000000001;//
            for (Node n : t.getGraph().getNodesArray()) {
                double newR = ApcomplexMath.abs(n.z).divide(new Apfloat(String.valueOf(norm), 100)).doubleValue();
                Apfloat newPhi = ApfloatMath.atan2(n.z.imag(), n.z.real());
                n.z1 = new Complex(newR * ApfloatMath.cos(newPhi).doubleValue(), newR * ApfloatMath.sin(newPhi).doubleValue());
            }
//            System.out.println("+Wrong double coordinates are fixed...");
        }
    }

    private static Apcomplex computeComplexCoordinate(Node child, Node parent, Apfloat two) {
        Apcomplex a = new Apcomplex(ApfloatMath.cos(child.alpha), ApfloatMath.sin(child.alpha));
        Apcomplex b = new Apcomplex(ApfloatMath.cos(child.betta), ApfloatMath.sin(child.betta));
        Apcomplex m = a.add(b);
        m = m.divide(two);
        Apcomplex c = new Apcomplex(Apfloat.ONE, Apfloat.ZERO).divide(m.conj());
        Apfloat r = ApfloatMath.pow(ApfloatMath.pow(ApcomplexMath.abs(m), two), new Apfloat(-1));
        r = r.subtract(Apfloat.ONE);

        Apcomplex denominator = parent.z.conj();
        denominator = denominator.subtract(c.conj());
        Apcomplex z = new Apcomplex(r, Apfloat.ZERO).divide(denominator);
        z = z.add(c);
        return z;
    }
}
