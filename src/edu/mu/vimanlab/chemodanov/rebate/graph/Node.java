/****************************************************************************/
/*                  Copyright 2001, Trustees of Boston University.          */
/*                               All Rights Reserved.                       */
/*                                                                          */
/* Permission to use, copy, or modify this software and its documentation   */
/* for educational and research purposes only and without fee is hereby     */
/* granted, provided that this copyright notice appear on all copies and    */
/* supporting documentation.  For any other uses of this software, in       */
/* original or modified form, including but not limited to distribution in  */
/* whole or in part, specific prior permission must be obtained from Boston */
/* University.  These programs shall not be used, rewritten, or adapted as  */
/* the basis of a commercial software or hardware product without first     */
/* obtaining appropriate licenses from Boston University.  Boston University*/
/* and the author(s) make no representations about the suitability of this  */
/* software for any purpose.  It is provided "as is" without express or     */
/* implied warranty.                                                        */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*  Author:     Alberto Medina                                              */
/*              Anukool Lakhina                                             */
/*  Title:     BRITE: Boston university Representative Topology gEnerator   */
/*  Revision:  2.0         4/02/2001                                        */
/****************************************************************************/

package edu.mu.vimanlab.chemodanov.rebate.graph;

import org.apache.commons.math3.complex.Complex;
import org.apfloat.Apcomplex;
import org.apfloat.Apfloat;
import org.apfloat.ApfloatMath;

import java.util.*;


/**
 * The Node class for our Graph contains only the base
 * minimum member variables.  Environment specific semantics (such as Router
 * or AS etc) that are often attached to a node are divorced from the
 * Node class and stored in a decorator class, NodeConf (short for
 * NodeConfiguration).  As such you can add/remove attributes to the Node at
 * run-time without having to change the Node representation.
 * <p/>
 * We use the decorator pattern as the relationship between a Node and
 * its configuration, NodeConf.  This pattern is described in Design Patterns:
 * Elements of Reusable Object-Oriented Software by Gamma et al.
 * ISBN#: 0-201-63361-2.
 * <p/>
 * Unique NodeIDs (used in Graph) are computed by maintaining a static
 * variable that is incremented each time a new nodes is created.
 */
public final class Node implements Comparable<Node> {

    int id;
    int indegree;
    int outdegree;
    int color;
    private Map<Node, Edge> neighbors;
    public Map<Node, List<List<Node>>> rt;
    public double minDistance = Double.POSITIVE_INFINITY;
    public Map<Integer, Node> nm_predecessor = null;
    Map<List<Node>, List<Double>> dominantPaths = null;
    private boolean virtual = false;
    private boolean isFailed = false;
    private boolean failureMode = false;
    private Map<Node, Edge> neighborsInFailureMode;

    //GEmbedding settings
    public Node parent = null;
    public boolean inTree = false;
    public int depth = 0;
    public Apcomplex z;
    public Apfloat ro0; //is used for normalization purpose
    public Apfloat phi0;
    public Apfloat alpha = ApfloatMath.pi(24);
    public Apfloat betta = ApfloatMath.pi(24).multiply(new Apfloat(2));

    //rebate.hyperbolic mapping
    public double ro;
    public double phi;
    public Complex z1; //coordinate in Poincare Disk Model;
    public Complex z2;

    //geo-data + outage risk
    public String label = "";
    public Double lat = null;
    public Double lon = null;
    private double risk = 0.0;

    public static int nodeCount = -1;

    /**
     * provides a comparator to sort nodes by their IDs
     */
    public static NodeIDComparator IDcomparator = new NodeIDComparator();

    /**
     * provides a comparator to sort nodes by their degree ()
     */
    public static NodeDegreeComparator DegreeComparator = new NodeDegreeComparator();


    /**
     * Class Constructor1.  Assigns the node a unique ID (maintained
     * by a static var) and does other initialization.
     */
    public Node() {
        this.id = ++nodeCount;
        this.color = GraphConstants.WHITE;
        indegree = 0;
        outdegree = 0;
        neighbors = new LinkedHashMap<>();
        neighborsInFailureMode = new LinkedHashMap<>();
        rt = new LinkedHashMap<>();
        nm_predecessor = new LinkedHashMap<>();
        dominantPaths = new LinkedHashMap<>();
    }


    /**
     * Constructor2: same as Constructor 1 but allows convenience of
     * specifying indegree and outdegree of the nodes.  Useful when
     * importing other topologies.
     */
    public Node(int inDeg, int outDeg) {
        this.id = ++nodeCount;
        this.color = GraphConstants.WHITE;
        this.indegree = inDeg;
        this.outdegree = outDeg;

    }

    /**
     * Constructor3: used for VNode creation
     *
     * @param id
     */
    public Node(int id, double target) {
        this.id = id;
        this.color = GraphConstants.WHITE;
        indegree = 0;
        outdegree = 0;
        neighbors = new LinkedHashMap<>();
        neighborsInFailureMode = new LinkedHashMap<>();
        rt = new LinkedHashMap<>();
        nm_predecessor = new LinkedHashMap<>();
        dominantPaths = new LinkedHashMap<>();
        this.virtual = true;
    }

    /**
     * Class Constructor4.  Assigns the node a unique ID (maintained
     * by a static var) and does other initialization.
     */
    public Node(int id) {
        this.id = id;
        this.color = GraphConstants.WHITE;
        indegree = 0;
        outdegree = 0;
        neighbors = new LinkedHashMap<>();
        neighborsInFailureMode = new LinkedHashMap<>();
        rt = new LinkedHashMap<>();
        nm_predecessor = new LinkedHashMap<>();
        //nm_minDist = new LinkedHashMap<>();
        dominantPaths = new LinkedHashMap<>();
    }

    public int compareTo(Node other) {
        return Double.compare(minDistance, other.minDistance);
    }

    /*get methods*/
    public int getID() {
        return id;
    }

    public int getOutDegree() {
        return outdegree;
    }

    public int getColor() {
        return color;
    }

    /*set methods*/
    public void setColor(int c) {
        this.color = c;
    }

    public void incrementInDegree() {
        ++this.indegree;
    }

    public void incrementOutDegree() {
        ++this.outdegree;
    }

    public void addNeighbor(Node n, Edge e) {
        neighbors.put(n, e);
    }

    public void clearNeighbors() {
        neighbors.clear();
    }

    public void setNMPredecessor(int nh, Node n) {
        if (nm_predecessor != null)
            this.nm_predecessor.put(nh, n);
        else
            throw new RuntimeException("Cannot set predecessor for " + nh + " neighborhood!");
    }

    public Map<Node, Edge> getNeighbors() {
        if (failureMode)
            return this.neighborsInFailureMode;
        else
            return this.neighbors;
    }

    public List<List<Node>> getRTPaths(Node dst) {
        if (rt.containsKey(dst))
            return rt.get(dst);
        else
            return null;
    }

    /**
     * @return maximum number of path hops stored in RT. Used to estimate net diameter.
     */
    public int maxRTPathLength() {
        int maxHops = 0;

        for (Map.Entry<Node, List<List<Node>>> entry : rt.entrySet()) {
            int minDstHops = Integer.MAX_VALUE; //find min hop path among k-shortest for a particular dst
            for (List<Node> pathToDst : entry.getValue())
                if (pathToDst.size() < minDstHops)
                    minDstHops = pathToDst.size();
            if (minDstHops > maxHops) //find max length path among min hop paths to dsts
                maxHops = minDstHops;
        }

        return maxHops;
    }

    public boolean isFailed() {
        return this.isFailed;
    }

    public void setFailed(boolean isFailed) {
        this.isFailed = isFailed;
    }

    public void setFailureMode(boolean failureMode) {
        this.failureMode = failureMode;
        if (failureMode)
            initNeighborsInFailureMode();
        else
            neighborsInFailureMode.clear();
    }

    private void initNeighborsInFailureMode() {
        if (!this.isFailed())
            for (Map.Entry<Node, Edge> nhEntry : this.neighbors.entrySet())
                if (!nhEntry.getKey().isFailed() && !nhEntry.getValue().isFailed())
                    neighborsInFailureMode.put(nhEntry.getKey(), nhEntry.getValue());
    }

    public String toString() {
        StringBuilder str = new StringBuilder();
        if (this.virtual)
            str.append("vnode[");
        else
            str.append("pnode[");

        str.append(this.getID()).append("]");

        return str.toString();
    }

    public int hashCode() {
        return this.toString().hashCode();
    }
}


/**
 * NodeID comparator provides a comparator to compare Node IDs.  You
 * can follow this template and trivially write your own comparator if
 * you need for instance, to sort the nodes in another fashion, eg
 * indegrees.  We use this comparator to sort nodes when printing them
 * to a file.
 */
class NodeIDComparator implements Comparator {

    public int compare(Object n1, Object n2) {
        int n1id = ((Node) n1).getID();
        int n2id = ((Node) n2).getID();

        if (n1id < n2id) return -1;
        if (n1id == n2id) return 0;
        if (n1id > n2id) return 1;

        return 1;  //never gets here but javac complains
    }
}

/**
 * NodeDegree comparator provides a comparator to compare Node degrees.
 */
class NodeDegreeComparator implements Comparator {

    public int compare(Object n1, Object n2) {
        double n1Cost = ((Node) n1).getNeighbors().size();
        double n2Cost = ((Node) n2).getNeighbors().size();

        if (n1Cost > n2Cost) return -1;
        if (n1Cost == n2Cost) return 0;
        if (n1Cost < n2Cost) return 1;

        return 1;  //never gets here but javac complains
    }
}

