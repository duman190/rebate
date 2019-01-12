REBATE: a REpulsive-BAsed Traffic Engineering
Approach for Dynamic Scale-Free Networks  
=================
In this project, we develop the first to our knowledge REpulsive-BAsed Traffic Engineering (REBATE) approach for dynamic scale-free networks. REBATE is built upon dual principles of the demand-aware TE and fundamentals properties of hyperbolic spaces. Using trace-driven numerical simulations from this repository, we show how REBATE can reduce the maximum link utilization up to 25% when compared to a common geometric routing-based traffic steering. However, REBATE can be still subject to practically large optimality gaps compared to both demands-aware and oblivious TE. Thus, our work should pave the way for more efficient TE in the nextgeneration
dynamic scale-free networks.

Authors
=================
2019 VIMAN laboratory, Computer Science Department, University of Missouri-Columbia.

```
Updated January 12, 2019 by Dmitrii Chemodanov
```

All feedback appreciated to dycbt4@mail.missouri.edu 

License
=================
This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details


What is inside?
================
The source code of the REBATE numerical simulations. Java simulation is used to evaluate our novel REBATE algorithm with alternative geometric hyperbolic routing solutions in dynamic scale-free networks using Tier-1 core network topology.

Distribution
================
The distribution tree contains: 

* README

	- this file
    
* LICENSE

	- GNU license file

* build.xml (build file for ant)    
    
* lib/      (java library dependencies)
    
* results/   (folder to store simulation results of experiments in a form of txt files)
    
* src/      (java source files)
        
        ```
        edu/mu/vimanlab/chemodanov/rebate/main/Main (java sim main file)
        ```
        

Compilation and run
============
Compiling and run of this software requires Ant, Java 1.8 and IBM ILOG CPLEX Optimization Studio (v.12.7 or higher) installed. These can be downloaded respectively from:  
http://jakarta.apache.org/ant/index.html 
http://java.sun.com/j2se/
http://www.ibm.com/products/ilog-cplex-optimization-studio

## Run
* Specify path to cplex in build.xml file

    ```
    <property name="cplex_path" location="/Users/user1/Applications/IBM/ILOG/CPLEX_Studio127/cplex/bin/x86-64_osx"/> 
    ```

* clean

    ```
    ant clean 
    ```
    
* compile

    ```
    ant
    ```

* run any of 4 experiments (0: the general experiment comparing optimal, geometric hyperbolic routing and rebate in terms of their optimization; experiment 1: a cold-down propery evaluation experiment; 2: a TTL experiment to access the rebate dynamic packet header size; 3: an experiment with network dynamics/failures)
    
    ```
    ant -Dexperiment=(from 0 to 3) run 
    
    ```