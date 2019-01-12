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

import edu.mu.vimanlab.chemodanov.rebate.forwarding.model.Obstacle;
import edu.mu.vimanlab.chemodanov.rebate.util.*;
import org.apache.commons.math3.complex.Complex;

import java.util.*;


/**
 * An Internet rebate.topology is represented as rebate.graph with the nodes
 * representing routers (or ASs, depending on the rebate.topology level) and
 * the edges representing links between them.   <p> We use an adjacency list
 * representation as our Graph implementation. Nodes are represented
 * by their own class, Node and similary, edges by a class, Edge.  The
 * Graph implementation is independent from any semantics attached to
 * a Graph, its nodes, and its edges.  As such, if you decide to use a
 * different representation of a rebate.graph, you can do so by simply
 * replacing this representation with your own and not affecting
 * BRITE's generation (as long as you expose a similar interface).
 * The idea is to allow for the ability to "plug and play" different
 * rebate.graph representations depending on your need.  (Some rebate.graph
 * representations might be faster than others for specific analysis
 * tasks etc) <p>
 * <p/>
 * We now provide implementation details of this Graph representation:
 * <br> This Graph representation has three HashMaps:
 * <code>Nodes</code>, <code>Edges</code> and
 * <code>adjList</code>. The <code>Nodes</code> HashMap contains as
 * keys the node-ids and as values, the <code>Node</code> objects
 * themselves. Similarly, the <code>Edges</code> HashMap contains
 * Edge-IDs as keys and <code>Edge</code> objects as values.  The
 * <code>adjList</code> HashMap contains as keys node-ids of source
 * nodes and destination node-ids as values.  (How NodeIDs and EdgeIDs
 * are computed can be found in the documentation for the Node and
 * Edge class respectively.)
 */
public class Graph implements GraphConstants {
    protected int numNodes;
    protected int numEdges;

    protected HashMap Nodes;  /*this is our repository of nodes & edges*/
    protected HashMap Edges;
    protected Map<Integer, Edge> IdToEdge;

    protected HashMap adjList;   /*this is adjaceny list representation of the
                         rebate.graph with nodeIDs */

    private int netDiam = 0;
    private double maxR = 0;
    private double minR = 0;

    private List<Obstacle> obstacles = new ArrayList<>();
    private List<Obstacle> obstaclesHyp = new ArrayList<>();

    /**
     * Create a rebate.graph with default initial capacity. The rebate.graph can grow beyond
     * the initial capacity as you add more nodes.
     */
    public Graph() {
        Nodes = new HashMap();
        Edges = new HashMap();
        adjList = new HashMap();
        IdToEdge = new HashMap<>();
    }

    /**
     * Create a rebate.graph with specified intial capacity.  The rebate.graph can
     * growth beyond this but its initial size is set.  If you know the
     * approximate size of the rebate.graph, use this constructor as it helps
     * performance.
     *
     * @param numNodes The initial number of nodes in the rebate.graph
     */
    public Graph(int numNodes) {
        Nodes = new HashMap(numNodes);
        Edges = new HashMap(2 * numNodes);
        IdToEdge = new HashMap<>(2 * numNodes);
        adjList = new HashMap(numNodes);
    }

    /**
     * Given a vector of graphs, produce a single rebate.graph with
     * disconnected components.  Useful for flattening/combining
     * different graphs.  This may be a memory-expensive operation since
     * you are creating a copy of all the nodes and edges of all the
     * graphs.
     */
    public Graph(ArrayList graphs) {
        //assumption: all graphs in vector are same size, but if not,
        //thats ok because datastructure automatically  increases
        int N = graphs.size() * (((Graph) graphs.get(0)).getNumNodes());
        ;
        Nodes = new HashMap(N);
        Edges = new HashMap(N);
        adjList = new HashMap(N);
        int size = graphs.size();
        for (int i = 0; i < size; ++i) {
            Graph g = (Graph) graphs.get(i);
            Nodes.putAll(g.getNodes());
            Edges.putAll(g.getEdges());
            //   System.out.println("Dumping Edges of rebate.graph: " + i+ "\n===========\n");
            //System.out.println(Edges.toString());
            //System.out.println("\n\n");
            adjList.putAll(g.getAdjList());
        }
        //System.out.println(Edges);
        numNodes = Nodes.size();
        numEdges = Edges.size();

        for (Object edge : Edges.values()) {
            Edge e = (Edge) edge;
            IdToEdge.put(e.getID(), e);
        }
    }

    /**
     * The addEdge method checks if the Edge is directed or
     * undirected.  If it is, it calls addDirectedEdge.  Otherwise, it
     * calls addUndirectedEdge.
     *
     * @param e The edge to be added
     */
    public void addEdge(Edge e) {
        if (e.getDirection() == GraphConstants.DIRECTED)
            addDirectedEdge(e);
        else addUndirectedEdge(e);
        IdToEdge.put(e.getID(), e);
    }


    /**
     * The addUndirectedEdge method adds the given edge's source node-id and
     * destination node-id to the adjList twice: once with the source
     * node-id as the key and next with the desitnation node-id as the
     * key.  We only increment the numEdges count once however.  Both
     * the indegree and outdegree of both nodes is incremented.
     * Finally, the edge itself is added to the Edges HashMap.
     *
     * @param e The edge to be added to our rebate.graph
     */
    public void addUndirectedEdge(Edge e) {
        ++numEdges;

	/*since this is an undirected rebate.graph, add edge from src to dest
      and another from dest to src*/
        Node src = e.getSrc();
        Node dst = e.getDst();

        if (src == null) {
            dumpToOutput();
            Util.ERR("src is null! in addEdge() of UndirectedGraph ");

        }
        if (dst == null) {

            dumpToOutput();
            Util.ERR("dst is null! in addEdge() of UndirectedGraph");
        }

        int intID = Edge.computeID(src.getID(), dst.getID());
        if (intID == -1) {
            long longID = Edge.computeLongID(src.getID(), dst.getID());
            Edges.put(new Long(longID), e);
        } else Edges.put(new Integer(intID), e);

	/*increment in/out degree of both nodes*/
        src.incrementInDegree();
        src.incrementOutDegree();
        dst.incrementInDegree();
        dst.incrementOutDegree();

	/*add to adj list*/
        Integer srcID = new Integer(src.getID());
        Integer dstID = new Integer(dst.getID());

        if (!adjList.containsKey(srcID)) {
            Nodes.put(srcID, src); //new node
            ++numNodes;

            ArrayList vect = new ArrayList();
            vect.add(dstID);
            adjList.put(srcID, vect);
        } else {
            ArrayList vect = (ArrayList) adjList.get(srcID);
            vect.add(dstID);
            adjList.put(srcID, vect);
        }

	/*now add reverse direction edge*/
        if (!adjList.containsKey(dstID)) {
            Nodes.put(dstID, dst); //new node
            ++numNodes;

            ArrayList vect = new ArrayList();
            vect.add(srcID);
            adjList.put(dstID, vect);

        } else {
            ArrayList vect = (ArrayList) adjList.get(dstID);
            vect.add(srcID);
            adjList.put(dstID, vect);
        }

    }


    /**
     * The addDirectedEdge method of a Graph takes the source node of
     * the edge to be added and adds it as a key in the adjList.  The
     * destination node-id is added as one of the values of this
     * source node-id.  The indegree of the destination node and the
     * outdegree of the source node is incremented here.  Next, the
     * edge itself is added to the Edges hashmap. Finally, the number
     * of edges counter is incremented.
     *
     * @param e the Edge to be added to our rebate.graph
     */
    public void addDirectedEdge(Edge e) {
        ++numEdges;

        Node src = e.getSrc();
        Node dst = e.getDst();

        if (src == null)
            Util.ERR("src is null! in addEdge() DirectedGraph");
        if (dst == null)
            Util.ERR("dst is null! in addEdge() DirectedGraph");
        Edge a;
        int intID = Edge.computeDirectedID(src.getID(), dst.getID());
        if (intID == -1) {
            long longID = Edge.computeDirectedLongID(src.getID(), dst.getID());
            a = (Edge) Edges.put(new Long(longID), e);
        } else a = (Edge) Edges.put(new Integer(intID), e);

        if (a != null)
            --numEdges;  //this was a repeat edge and therefore is thrown out.

	/*increment inDegree of dst and outdegree of src*/
        dst.incrementInDegree();
        src.incrementOutDegree();

	/*add to adj list*/
        Integer srcID = new Integer(src.getID());
        Integer dstID = new Integer(dst.getID());
        if (!adjList.containsKey(srcID)) {
            Nodes.put(srcID, src); //new node
            ++numNodes;
            ArrayList vect = new ArrayList();
            vect.add(dstID);
            adjList.put(srcID, vect);
        } else {
            ArrayList vect = (ArrayList) adjList.get(srcID);
            vect.add(dstID);
            adjList.put(srcID, vect);
        }

    }

    /**
     * given IDs of source and destination nodes, returns true if
     * an edge between the corresponding nodes exists.
     *
     * @param srcID the node id of the source node
     * @param dstID the node if of the destination node
     */
    public boolean hasEdge(int srcID, int dstID) {

        int intID = Edge.computeID(srcID, dstID);
        if (intID == -1) {
            long longID = Edge.computeLongID(srcID, dstID);
            if (Edges.containsKey(new Long(longID)))
                return true;
        } else {
            if (Edges.containsKey(new Integer(intID)))
                return true;
        }
        return false;
    }

    /**
     * add a node to our adjacency list.
     *
     * @param n the node to be added
     */
    public void addNode(Node n) {
        Integer nID = new Integer(n.getID());
        Nodes.put(nID, n);
        if (!adjList.containsKey(nID)) {
            adjList.put(nID, new ArrayList());
            ++numNodes;
        }
    }

    /**
     * a useful debug routine.  dumps the rebate.graph (nodes and edges) in NLANR ASConnlist-list format.
     * to the standard output stream (usually the console)
     * The output format looks like:
     * <pre>
     * from ->  to1, to2, to3
     * from2 -> to7, to8, to9
     * ...
     * NumEdges = <numEdges>
     * NumNodes = <numNodes>
     * </pre>
     */
    public void dumpToOutput() {
        Iterator kI = adjList.keySet().iterator();

        while (kI.hasNext()) {
            Integer n = (Integer) kI.next();
            ArrayList v = (ArrayList) adjList.get(n);
            if (v.size() == 0)
                continue;
            System.out.print(n + " -> ");
            int size = v.size();
            for (int i = 0; i < size - 1; ++i) {
                Integer ni = (Integer) v.get(i);
                System.out.print(ni + ", ");
            }
            Integer lastNode = (Integer) v.get(v.size() - 1);
            System.out.println(lastNode + ".");
        }
        System.out.println("NumEdges = " + numEdges);
        System.out.println("NumNodes = " + numNodes);
    }

    /**
     * Returns an Array representation (copy) of the ndoes in this rebate.graph
     */
    public Node[] getNodesArray() {
        return (Node[]) Nodes.values().toArray(new Node[numNodes]);
    }

    /**
     * Returns an ArrayList representation (copy) of all the edges in this rebate.graph
     */
    public ArrayList getEdgesVector() {
        return new ArrayList(Edges.values());
    }

    /**
     * returns an Array representation (copy) of all the edges in this rebate.graph
     */
    public Edge[] getEdgesArray() {
        return (Edge[]) Edges.values().toArray(new Edge[Edges.size()]);
    }

    /**
     * get number of nodes in this rebate.graph
     */
    public int getNumNodes() {
        return numNodes;
    }

    /**
     * get number of edges in this rebate.graph
     */
    public int getNumEdges() {
        return numEdges;
    }

    /**
     * returns a HashMap representation of the Nodes in the rebate.graph. THe
     * keys are the Node-IDs and the values are the Node object
     * references themselves.  WARNING: This is not a copy of the
     * rebate.graph, so changes here will affect the rebate.graph structure!
     */
    public HashMap getNodes() {
        return Nodes;
    }

    /**
     * returns a HashMap representation of the Edges in the rebate.graph. THe
     * keys are the Edge-IDs (as computed by the Edge.ComputeID() method)
     * and the values are the Node object references themselves.
     * WARNING: This is not a copy of the rebate.graph, so changes here will
     * affect the rebate.graph structure!
     */
    public HashMap getEdges() {
        return Edges;
    }

    /**
     * returns a HashMap representation of the adjacency list. The
     * keys are node-ids of the source nodes and the values are an
     * Arraylist of node-ids of destination nodes.
     * WARNING: This is not a copy of the rebate.graph, so changes here will
     * affect the rebate.graph structure!
     */
    public HashMap getAdjList() {
        return adjList;
    }

    /**
     * Given node n, returns an array of Node objects that are n's neighbors.
     *
     * @return Node[] Array of nodes
     */
    public Node[] getNeighborsOf(Node n) {
        Integer nID = new Integer(n.getID());
        ArrayList neighborID = (ArrayList) adjList.get(nID);
        int size = neighborID.size();

        Node[] neighbors = new Node[size];
        for (int i = 0; i < size; ++i) {
            int id = ((Integer) neighborID.get(i)).intValue();
            neighbors[i] = (Node) getNodeFromID(id);

        }
        return neighbors;
    }

    /**
     * a lookup function to get a node from its id.  assumes id is valid.  returns null if invalid id
     */
    public Node getNodeFromID(int id) {
        return (Node) Nodes.get(new Integer(id));
    }

    /**
     * a lookup function to get an edge from its id.  assumes id is valid.  returns null if invalid id
     */
    public Edge getEdgeFromID(int id) {
        return IdToEdge.get(id);
    }

    public void initNeighbors() {
        cleanNeighbors();

        for (Object o : getEdgesVector()) {
            Edge e = (Edge) o;
            Node src = e.getSrc();
            Node dst = e.getDst();
            src.addNeighbor(dst, e);
            dst.addNeighbor(src, e);
        }
    }

    public int getNetDiameter() {
        return this.netDiam;
    }

    public void initMinMaxHypR() {
        this.maxR = 0;
        this.minR = Double.MAX_VALUE;
        Complex origin = new Complex(0, 0);
        for (Node n : getNodesArray()) {
            double dist = Dist.poincare(origin, n.z1);
            if(dist > this.maxR)
                this.maxR = dist;
            if(dist < this.minR)
                this.minR = dist;
//            System.out.println("[Topo] RT and diameter were calculated for node id=" + n.getID() +
//                    " with estimated so far diam=" + netDiam + ". [DEBUG] current max hop path in RT: " + currentDiam);
        }
    }

    public double getMaxR()
    {
        return this.maxR;
    }

    public double getMinR()
    {
        return this.minR;
    }

    public void cleanNeighbors() {
        for (Node n : getNodesArray())
            n.clearNeighbors();
    }

    public void initNeighborsWithFailures() {
        for (Node n : getNodesArray())
            n.setFailureMode(true);
    }

    public void releaseAllResources() {
        for (Object o : getEdgesVector()) {
            Edge e = (Edge) o;
            e.releaseResources();
        }
    }
}

















