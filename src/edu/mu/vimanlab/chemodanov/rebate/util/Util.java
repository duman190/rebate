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
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

/**
 * A repository of miscellaneous utility functions.
 */
public final class Util {

    /**
     * Print an error formatted string and exit app
     */
    public static void ERR(String err) {
        System.out.println("[ERROR]  : " + err);
        System.exit(0);
    }

    public static double gpsDistKm(double lat1, double lon1, double lat2, double lon2) {
        int earthRadiusKm = 6371;
        double dLat = (lat2 - lat1) * Math.PI / 180;
        double dLon = (lon2 - lon1) * Math.PI / 180;
        lat1 = lat1 * Math.PI / 180;
        lat2 = lat2 * Math.PI / 180;
        double a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
                Math.sin(dLon / 2) * Math.sin(dLon / 2) * Math.cos(lat1) * Math.cos(lat2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return earthRadiusKm * c;
    }

    public static double distKmToDelayMs(double dist) {
        double lightSpeedFiber = 200; // km per ms
        return dist / lightSpeedFiber;
    }

    public static void exportEdgesUtil(String output, Graph g) {
        try {
            FileWriter fw = new FileWriter(output, true);
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter writer = new PrintWriter(bw);
            for (Edge e : g.getEdgesArray()) {
                double allocatedBw = e.getBW() - e.getAvBW();
                double utilization = allocatedBw / Double.valueOf(e.getBW());
                writer.println(utilization);
            }
            writer.close();
        } catch (IOException e) {
        }
    }

    public static void exportValue(String output, double d) {
        try {
            FileWriter fw = new FileWriter(output, true);
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter writer = new PrintWriter(bw);
            writer.println(d);
            writer.close();
//            System.out.print("File saved! ");
        } catch (IOException e) {
            // do something
        }
    }

    public static void exportValues(String output, double[] d) {
        try {
            FileWriter fw = new FileWriter(output, true);
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter writer = new PrintWriter(bw);
            for (double i : d)
                writer.println(i);
            writer.close();
        } catch (IOException e) {
        }
    }

    public static void exportValues(String output, List d) {
        try {
            FileWriter fw = new FileWriter(output, true);
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter writer = new PrintWriter(bw);
            for (Object i : d)
                writer.println(i);
            writer.close();
        } catch (IOException e) {
        }
    }

    public static void exportMatrix(String output, int[][] m) {
        try {
            FileWriter fw = new FileWriter(output, true);
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter writer = new PrintWriter(bw);
            for (int[] i : m) {
                for (int j : i)
                    writer.print(j + " ");
                writer.println();
            }
            writer.close();
        } catch (IOException e) {
        }
    }
}








