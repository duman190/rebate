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

package edu.mu.vimanlab.chemodanov.rebate.main;

/**
 * Created by dmitriichemodanov on 1/12/19.
 */
public class Main {
    public static void main (String[] args)
    {
        if (args.length >= 1 && Integer.valueOf(args[0]) >= 0 && Integer.valueOf(args[0]) <=3)
            switch(Integer.valueOf(args[0])){
                case 0: REBATE_Demand_Experiment.main(null);
                    break;
                case 1: REBATE_CoolDown_Experiment.main(null);
                    break;
                case 2: REBATE_TTL_Experiment.main(null);
                    break;
                case 3: REBATE_Failure_Experiment.main(null);
                    break;
            }
        else
            System.out.println("Not enough or wrong input arguments. Please specify at a valid scenario number from 0 to 3!");
    }
}
