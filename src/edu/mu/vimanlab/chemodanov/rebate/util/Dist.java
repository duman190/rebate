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

import edu.mu.vimanlab.chemodanov.rebate.graph.Node;
import org.apache.commons.math3.analysis.function.Acosh;
import org.apache.commons.math3.complex.Complex;

/**
 * Created by chemo_000 on 12/28/2015.
 */
public class Dist {
    public static double poincare(Complex z1, Complex z2) {
        double d = Math.pow(z1.subtract(z2).abs(), 2) /
                ((1 - Math.pow(z1.abs(), 2)) * (1 - Math.pow(z2.abs(), 2)));
        return new Acosh().value(1 + 2 * d);
    }

    public static double poincare(Node src, Node dst) {
        return poincare(src.z1, dst.z1);
    }
}
