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

package edu.mu.vimanlab.chemodanov.rebate.util;

import edu.mu.vimanlab.chemodanov.rebate.graph.*;
import edu.mu.vimanlab.chemodanov.rebate.topology.Topology;
import org.apache.commons.math3.complex.Complex;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;
import java.util.*;

/**
 * Created by dmitriichemodanov on 6/1/17.
 */
public class ImportZooGraph implements InternetZooConstants {
    public static Topology importZooTopo(String input) {
        Map<Pair<String, String>, Integer> gpsToNodeId = new LinkedHashMap<>(); //latitude first, longitude second
        Map<Integer, Complex> nodeIdToVCoordinate = new HashMap<>();
        Map<Integer, String> nodeIdToLabel = new LinkedHashMap<>(); //latitude first, longitude second
        Map<Pair<Integer, Integer>, Integer> srcDstToEdgeId = new LinkedHashMap<>();
        Map<Pair<Integer, Integer>, Float> srcDstToEdgeBw = new LinkedHashMap<>();

        System.out.print("Import " + input + "... ");
        try {
            File inputFile = new File(input); //

            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(inputFile);

            Map<Integer, Integer> localToGlobalNodeId = new HashMap<>();
            Map<Integer, String> localNodeIdToLabel = new HashMap<>();

            NodeList nList = doc.getElementsByTagName("node");
            for (int temp = 0; temp < nList.getLength(); temp++) {
                org.w3c.dom.Node nNode = nList.item(temp);
                if (nNode != null && nNode.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE) {
                    Element eElement = (Element) nNode;
                    int localNodeId = Integer.valueOf(eElement.getAttribute("id"));
                    NodeList nodeData = eElement.getElementsByTagName("data");
                    String lat = null;
                    String lon = null;
                    String label = null;
                    Double real = 0.0;
                    Double img = 0.0;
                    for (int ni = 0; ni < nodeData.getLength(); ni++) {
                        org.w3c.dom.Node dataI = nodeData.item(ni);
                        Element elementI = (Element) dataI;
                        if (elementI.getAttribute("key").equals(D29))
                            lat = dataI.getTextContent();
                        else if (elementI.getAttribute("key").equals(D32))
                            lon = dataI.getTextContent();
                        else if (elementI.getAttribute("key").equals(D33))
                            label = dataI.getTextContent();
                        else if (elementI.getAttribute("key").equals(D98))
                            real = Double.valueOf(dataI.getTextContent());
                        else if (elementI.getAttribute("key").equals(D99))
                            img = Double.valueOf(dataI.getTextContent());
                    }
                    if (label != null)
                        localNodeIdToLabel.put(localNodeId, label);
                    if (lon != null && lat != null) {
                        Pair<String, String> nodeGeoKey = new Pair<>(lat, lon);
                        if (gpsToNodeId.containsKey(nodeGeoKey))
                            localToGlobalNodeId.put(localNodeId, gpsToNodeId.get(nodeGeoKey)); // we need this set to merge edges
                        else {
                            localToGlobalNodeId.put(localNodeId, gpsToNodeId.size()); // we need this set to merge edges
                            if (label != null)
                                nodeIdToLabel.put(gpsToNodeId.size(), label);
                            gpsToNodeId.put(nodeGeoKey, gpsToNodeId.size());
                        }
                        nodeIdToVCoordinate.put(gpsToNodeId.get(nodeGeoKey), new Complex(real, img));
                    }
                }
            }
            NodeList eList = doc.getElementsByTagName("edge");
            for (int temp = 0; temp < eList.getLength(); temp++) {
                org.w3c.dom.Node eNode = eList.item(temp);
                if (eNode != null && eNode.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE) {
                    Element eElement = (Element) eNode;
                    int localSrcId = Integer.valueOf(eElement.getAttribute("source"));
                    int localDstId = Integer.valueOf(eElement.getAttribute("target"));
                    NodeList edgeData = eElement.getElementsByTagName("data");
                    String edgeId = null;
                    float edgeBw = 0;
                    for (int ni = 0; ni < edgeData.getLength(); ni++) {
                        org.w3c.dom.Node dataI = edgeData.item(ni);
                        Element elementI = (Element) dataI;
                        if (elementI.getAttribute("key").equals(D37))
                            edgeId = dataI.getTextContent();
                        if (elementI.getAttribute("key").equals(D36))
                            edgeBw = Float.valueOf(dataI.getTextContent());
                    }
                    //System.out.println("Edge old id=" + edgeId + " new id=" + srcDstToEdgeId.size());

                    Integer globalSrcId = localToGlobalNodeId.get(localSrcId);
                    Integer globalDstId = localToGlobalNodeId.get(localDstId);
                    if (!localToGlobalNodeId.containsKey(localSrcId) && localNodeIdToLabel.containsKey(localSrcId)) {
                        String srcLabel = localNodeIdToLabel.get(localSrcId);
                        for (Map.Entry<Integer, String> entry : nodeIdToLabel.entrySet())
                            if (entry.getValue().equals(srcLabel))
                                globalSrcId = entry.getKey();
                    }
                    if (!localToGlobalNodeId.containsKey(localDstId) && localNodeIdToLabel.containsKey(localDstId)) {
                        String srcLabel = localNodeIdToLabel.get(localDstId);
                        for (Map.Entry<Integer, String> entry : localNodeIdToLabel.entrySet())
                            if (entry.getValue().equals(srcLabel) && localToGlobalNodeId.containsKey(entry.getKey()))
                                globalDstId = entry.getKey();
                    }
                    if (globalSrcId != null && globalDstId != null) {
                        Pair<Integer, Integer> srcDstEdgeKey = new Pair<>(globalSrcId, globalDstId);
                        Pair<Integer, Integer> srcDstEdgeKeyRev = new Pair<>(globalDstId, globalSrcId);
                        if (!srcDstToEdgeId.containsKey(srcDstEdgeKey) && !srcDstToEdgeId.containsKey(srcDstEdgeKeyRev)) {
                            srcDstToEdgeId.put(srcDstEdgeKey, srcDstToEdgeId.size()); //we omit storing of the reverse key by assuming undirected rebate.graph
                            srcDstToEdgeBw.put(srcDstEdgeKey, edgeBw); //we omit storing of the reverse key by assuming undirected rebate.graph
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        cleanData(nodeIdToVCoordinate, gpsToNodeId, srcDstToEdgeId);
        System.out.println("Graph with total " + gpsToNodeId.size() + " nodes and " + srcDstToEdgeId.size() + " edges is parsed!");

        System.out.print("Normalizing... ");
        normilizeNodeIds(nodeIdToVCoordinate, gpsToNodeId, srcDstToEdgeId, srcDstToEdgeBw);
        System.out.println("Graph with total " + gpsToNodeId.size() + " nodes and " + srcDstToEdgeId.size() + " edges is normilized!");

        System.out.print("Convert to topo... ");
        Topology t = convertToTopo(nodeIdToVCoordinate, gpsToNodeId, srcDstToEdgeId, srcDstToEdgeBw, nodeIdToLabel);
        System.out.println("Done: rebate.graph with " + t.getGraph().getNumNodes() + " nodes and " + t.getGraph().getNumEdges() + " edges.");

        return t;
    }

    private static void cleanData(Map<Integer, Complex> nodeIdToVCoordinate,
                                  Map<Pair<String, String>, Integer> gpsToNodeId,
                                  Map<Pair<Integer, Integer>, Integer> srcDstToEdgeId) {
        Set<Integer> connectedNodes = new HashSet<>();
        for (Map.Entry<Pair<Integer, Integer>, Integer> entry : srcDstToEdgeId.entrySet()) {
            connectedNodes.add(entry.getKey().first());
            connectedNodes.add(entry.getKey().second());
        }
        Set<Pair<String, String>> keysToDelete = new HashSet<>();
        for (Map.Entry<Pair<String, String>, Integer> entry : gpsToNodeId.entrySet())
            if (!connectedNodes.contains(entry.getValue()))
                keysToDelete.add(entry.getKey());
        for (Pair<String, String> key : keysToDelete) {
            nodeIdToVCoordinate.remove(gpsToNodeId.remove(key));
        }
    }

    private static void normilizeNodeIds(Map<Integer, Complex> nodeIdToVCoordinate,
                                         Map<Pair<String, String>, Integer> gpsToNodeId,
                                         Map<Pair<Integer, Integer>, Integer> srcDstToEdgeId,
                                         Map<Pair<Integer, Integer>, Float> srcDstToEdgeBw) {
        Map<Integer, Complex> nodeIdToVCoordinateNew = new LinkedHashMap<>();
        Map<Pair<String, String>, Integer> gpsToNodeIdNew = new LinkedHashMap<>(gpsToNodeId.size());
        Map<Pair<Integer, Integer>, Integer> srcDstToEdgeIdNew = new LinkedHashMap<>(srcDstToEdgeId.size());
        Map<Pair<Integer, Integer>, Float> srcDstToEdgeBwNew = new LinkedHashMap<>(srcDstToEdgeBw.size());
        Map<Integer, Integer> oldToNewId = new HashMap<>(gpsToNodeId.size());

        for (Map.Entry<Pair<String, String>, Integer> entry : gpsToNodeId.entrySet()) {
            oldToNewId.put(entry.getValue(), gpsToNodeIdNew.size());
            gpsToNodeIdNew.put(entry.getKey(), gpsToNodeIdNew.size());
        }

        for (Map.Entry<Integer, Integer> idEntry : oldToNewId.entrySet())
            nodeIdToVCoordinateNew.put(idEntry.getValue(), nodeIdToVCoordinate.get(idEntry.getKey()));

        for (Map.Entry<Pair<Integer, Integer>, Integer> entry : srcDstToEdgeId.entrySet())
        {
            Integer newSrc = oldToNewId.get(entry.getKey().first());
            Integer newDst = oldToNewId.get(entry.getKey().second());
            if(newSrc != null && newDst != null) {
                Pair<Integer, Integer> edgeKey = new Pair<>(newSrc, newDst);
                srcDstToEdgeIdNew.put(edgeKey, srcDstToEdgeIdNew.size());
                srcDstToEdgeBwNew.put(edgeKey, srcDstToEdgeBw.get(entry.getKey()));
            }
        }

        //update data sets
        gpsToNodeId.clear();
        gpsToNodeId.putAll(gpsToNodeIdNew);
        nodeIdToVCoordinate.clear();
        nodeIdToVCoordinate.putAll(nodeIdToVCoordinateNew);
        srcDstToEdgeId.clear();
        srcDstToEdgeId.putAll(srcDstToEdgeIdNew);
        srcDstToEdgeBw.clear();
        srcDstToEdgeBw.putAll(srcDstToEdgeBwNew);
    }

    private static Topology convertToTopo(Map<Integer, Complex> nodeIdToVCoordinate,
                                          Map<Pair<String, String>, Integer> gpsToNodeId,
                                          Map<Pair<Integer, Integer>, Integer> srcDstToEdgeId,
                                          Map<Pair<Integer, Integer>, Float> srcDstToEdgeBw,
                                          Map<Integer, String> nodeIdToLabel) {
        Graph g = new Graph(gpsToNodeId.size());

        for (Map.Entry<Pair<String, String>, Integer> entry : gpsToNodeId.entrySet()) {
            Node n = new Node(entry.getValue());
            n.label = nodeIdToLabel.containsKey(entry.getValue()) ? nodeIdToLabel.get(entry.getValue()) : "";
            n.lat = Double.valueOf(entry.getKey().first());
            n.lon = Double.valueOf(entry.getKey().second());
            n.z1 = nodeIdToVCoordinate.get(entry.getValue());
            g.addNode(n);
        }

        for (Map.Entry<Pair<Integer, Integer>, Integer> entry : srcDstToEdgeId.entrySet())
        {
            Node src = g.getNodeFromID(entry.getKey().first());
            Node dst = g.getNodeFromID(entry.getKey().second());
            Edge e = new Edge(entry.getValue(), src, dst, Util.distKmToDelayMs(Util.gpsDistKm(src.lat, src.lon, dst.lat, dst.lon)));
            e.setCost(1);
            e.setBW(srcDstToEdgeBw.get(entry.getKey()));
            g.addEdge(e);
        }

        Topology t = new Topology(g);
        return t;
    }
}
