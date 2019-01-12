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

import java.lang.*;
import java.util.Comparator;
import java.util.Random;


/**
 * The Edge class for our Graph contains only the base minimum member
 * variables.  Environment specific semantics (such as RouterEdge or
 * bandwidth etc) that are often attached to an edge are divorced from
 * the this class and stored in a decorator class, EdgeConf (short for
 * EdgeConfiguration).  As such you can add/remove attributes to the Edge at
 * run-time without having to change the Edge representation.  <p> We
 * use the decorator pattern as the relationship between a Edge and
 * its configuration, EdgeConf.  This pattern is described in Design Patterns:
 * Elements of Reusable Object-Oriented Software by Gamma et al.
 * ISBN#: 0-201-63361-2.  <p>
 * <p/>
 * Like NodeIDs, unique edge ids are determined by using a static int
 * and incrementing it each time a new edge is created.  This is the
 * default method.  However, for improved lookup performance, we also
 * provide alternative methods of computing EdgeIDs which embed the
 * ids of the source and destination nodes of this edge in the EdgeID.
 * As such given the source and destination node, one can compute the
 * EdgeID in constant time. See Edge.computeID(..) and
 * Edge.computeDirectedID(..) methods for more on how this is done.  <p>
 * <p/>
 * NOTE: The direction of the edge can be either DIRECTED or UNDIRECTED.
 * This allows for graphs that contain a hybrid of directed and
 * undirected edges.
 */

public final class Edge {
    Node src;
    int direction;
    Node dst;
    int id;
    int color;
    float BW;
    float avBW;
    double delay = 0;
    double cost = 0;
    public double delayDijkstra = 0;
    public double costDijkstra = 0;
    public double hopDijkstra = 1;
    public int popularity = 0; //edge appearance in shortest paths between all nodes in the network

    //K-Shortest Path Settings
    public boolean isDeleted = false;

    float dist = -1;
    final int C = 1;

    public static int edgeCount = -1;
    public static EdgeIDComparator IDcomparator = new EdgeIDComparator();
    public static EdgePopularityComparator PopularityComparator = new EdgePopularityComparator();
    public static EdgeCostComparator CostComparator = new EdgeCostComparator();
    private boolean isFailed = false;

    /**
     * Class Constructor.
     *
     * @param src The source node of this edge
     * @param dst The destination node of this edge
     */
    public Edge(Node src, Node dst) {
        this.src = src;
        this.dst = dst;
        this.delay = 0;
        this.delayDijkstra = this.delay;
        this.cost = new Random().nextDouble() * 8 + 1;
        this.costDijkstra = this.cost;
        this.id = ++edgeCount;
        this.direction = GraphConstants.UNDIRECTED;  //by default we select this
    }

    /**
     * Constructor2: use for VLink creation.
     *
     * @param src The source node of this edge
     * @param dst The destination node of this edge
     */
    public Edge(int id, Node src, Node dst, double bw, double delayC) {
        this.src = src;
        this.dst = dst;
        this.BW = (float) bw;
        this.avBW = this.BW;
        this.delay = delayC;
        this.delayDijkstra = this.delay;
        this.cost = 0;
        this.costDijkstra = this.cost;
        this.id = id;
        this.direction = GraphConstants.UNDIRECTED;  //by default we select this
    }

    /**
     * Constructor3: use for custom rebate.graph.
     *
     * @param src The source node of this edge
     * @param dst The destination node of this edge
     */
    public Edge(int id, Node src, Node dst, double delayC) {
        this.src = src;
        this.dst = dst;
        this.BW = 0; //randomly generate bw between 10 and 0 Gbps
        this.avBW = this.BW;
        this.delay = delayC;
        this.delayDijkstra = this.delay;
        this.cost = 0;
        this.costDijkstra = this.cost;
        this.id = id;
        this.direction = GraphConstants.UNDIRECTED;  //by default we select this
    }

    /**
     * Set the direction of this edge to either GraphConstants.DIRECTED or GraphConstants.UNDIRECTED
     *
     * @param d the direction of the rebate.graph, one of the possible values specified in class GraphConstants
     */
    public void setDirection(int d) {
        direction = d;
    }

    /**
     * returns direction of this edge. (either GraphConstants.DIRECTED or GraphConstants.UNDIRECTED)
     *
     * @return int The Direction of the edge.  Either GraphConstants.DIRECTED or GraphConstants.UNDIRECTED
     */
    public int getDirection() {
        return direction;
    }

    public double getDelay() {
        return delay;
    }

    public double getCost() {
        return cost;
    }

    public void setCost(double cost) {
        this.cost = cost;
    }


    /**
     * Computes a unique EdgeID These IDs have the property that
     * id(src,dest) = id(dest,src) and so should be used with undirected
     * graphs only.  See the computeDirectedID method for computing
     * edgeIDs for directed graphs.
     * <p/>
     * An EdgeID embeds the IDs of the srouce and destination nodes
     * of this edge in it.  This is done by simply concatenating the
     * bit represenation of the source and destination.  If the
     * concatenated represenation is larger than an int, -1 is
     * returned.  The caller of the function should check for this
     * condition and if a -1 is returned, computeLongID should be
     * called instead.
     *
     * @param srcID generally, the id of the source-node
     * @param dstID generally, the id of the dest-node
     * @return int returns an int which is srcID concattenated with
     * dstID. -1 if concattenated result overflows int.
     */
    public static int computeID(int srcID, int dstID) {
        //WARNING: only works for undirected  (a,b) has same id as (b,a)

        int d = dstID >> 16;
        int s = srcID >> 16;
        if (d == 0 && s == 0) {
            if (srcID < dstID)
                return ((srcID << 16) | dstID);   //this gurantees (s,d) to have sameid as (d,s)
            else return ((dstID << 16) | srcID);
        }
        //System.out.println("DEBUG:  need long edgeid for "+srcID + " and " + dstID);
        return -1;

    }

    /**
     * Similar to computeID(src,dest) except returns a long id.  this
     * should be used if the srcID and destID are too large to yield
     * an EdgeID which can fit in an int.  Gurantees that id(src,dst)
     * == id(dst,src)
     *
     * @param srcID the node-id of the srouce node
     * @param dstID the node-id of the dest node
     * @return long  this edgeID is a long repr. of srcID concattenated with dstID
     */
    public static long computeLongID(int srcID, int dstID) {
        long lo;
        if (srcID < dstID) {
            lo = ((long) srcID << 32 | dstID);
        } else {
            lo = ((long) dstID << 32 | srcID);
        }
        //System.out.println("** lo="+lo);
        return lo;

    }


    /*Analagous to computeID() above but this computes IDs for directed rebate.graph.
      That is, it  computeDirectedID(a,b) != computeDirectedID(b,a)

      @param srcID  the source node id
      @param dstID  the desitnation node id
      @return int  the result of concatenating srcID with dstID, or -1 if the result overflows an int.

    */
    public static int computeDirectedID(int srcID, int dstID) {
        int d = dstID >> 16;
        int s = srcID >> 16;

        if (d == 0 && s == 0)
            return ((srcID << 16) | dstID);
        return -1;

    }

    /**
     * Analagous to computeLongID above but computes ID for directed rebate.graph.  That is, ids returned by
     * this method gurantee that id(a,b) ! = id(b,a).
     *
     * @param srcID the source node id
     * @param dstID the dest node id
     * @return long  the result of concatenating srcID with dstID
     */
    public static long computeDirectedLongID(int srcID, int dstID) {
        long lo = ((long) srcID << 32 | dstID);
        return lo;
    }



    /*get methods*/
    public Node getSrc() {
        return src;
    }

    public Node getDst() {
        return dst;
    }

    public int getID() {
        return this.id;
    }

    public static int getEdgeCount() {
        return edgeCount;
    }

    public int getColor() {
        return color;
    }

    public float getBW() {
        return this.BW;
    }

    ;

    public float getAvBW() {
        return this.avBW;
    }

    ;

    /*set methods*/
    public void setSrc(Node src) {
        this.src = src;
    }

    public void setDst(Node dst) {
        this.dst = dst;
    }

    public void setColor(int c) {
        color = c;
    }

    public void setBW(float bw) {
        this.BW = bw;
        this.avBW = bw;
//        this.delay = C*bw/calculateDist(src, dst);
    }

    public void releaseResources() {
        this.avBW = this.BW;
    }

    public void releaseSpecifiedResources(float bw) {
        this.avBW += bw;
    }

    public void allocateResources(float bw) {
        if (bw <= avBW)
            this.avBW = this.avBW - bw;
        else
            throw new RuntimeException("Resources of " + bw + " for link:" + src.getID() + "->" + dst.getID()
                    + " with avBw=" + avBW + " can't be allocated!");
    }

    public void setEuclideanDist(float d) {
        this.dist = d;
    }

    public boolean isFailed() {
        return this.isFailed;
    }

    public void setFailed(boolean isFailed) {
        this.isFailed = isFailed;
    }
}

/**
 * EdgeID comparator provides a comparator to compare Edge IDs.  You
 * can follow this template and trivially write your own comparator if
 * you need for instance, to sort the edges in another fashion, eg
 * source node-ids.  We use this comparator to sort edges when printing them
 * to a file.
 */
class EdgeIDComparator implements Comparator {

    public int compare(Object e1, Object e2) {
        int e1id = ((Edge) e1).getID();
        int e2id = ((Edge) e2).getID();
    
    /*if e1 < e2 then return -1*/
        if (e1id < e2id) return -1;
    /*if e1==e2, then return 0*/
        if (e1id == e2id) return 0;
    /*if e1> e2 then return 1*/
        if (e1id > e2id) return 1;

    /*should never get here*/
        return 1;

    }
}

/**
 * EdgePopularity comparator provides a comparator to compare Edge popularities.
 */
class EdgePopularityComparator implements Comparator {
    public int compare(Object e1, Object e2) {
        double popularity1 = ((Edge) e1).popularity;
        double popularity2 = ((Edge) e2).popularity;

        if (popularity1 > popularity2) return -1;
        if (popularity1 == popularity2) return 0;
        if (popularity1 < popularity2) return 1;

        return 1;  //never gets here but javac complains
    }
}


/**
 * EdgePopularity comparator provides a comparator to compare Edge popularities.
 */
class EdgeCostComparator implements Comparator {
    public int compare(Object e1, Object e2) {
        double cost1 = ((Edge) e1).costDijkstra;
        double cost2 = ((Edge) e2).costDijkstra;

        if (cost1 < cost2) return -1;
        if (cost1 == cost2) return 0;
        if (cost1 > cost2) return 1;

        return 0;
    }
}





