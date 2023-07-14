# Permutation-Based Approximate Multiplier with High Accuracy

This repository contains:

- software: the MATLAB code of the optimization method.
- multipliers: Verilog models of reproduced multipliers and generated multipliers.
- scripts: the scripts to synthesize multipliers with [Arizona State Predictive PDK (ASAP) 7nm process library](https://github.com/The-OpenROAD-Project/asap7) in Synopsys Design Compiler (DC).

## software

The compressed partial products of all possible optimizations are first generated through MATLAB, followed by the GA to give the estimated parameters of the hardware cost and the results of the final optimization. We feed the obtained Verilog code into Design Compiler (DC) for Static Timing Analysis (STA) to get the PDA results; the obtained C file is co-compiled with the test code to get the accuracy results.

The 'src' folder contains the MATLAB code of the method. Please follow the steps to generate multipliers:

- step1:  select the unsigned multiplier or the signed multiplier: sign = 0 or 1 in 'Exchange.m'.
- step2: decide the number of rows of the partial products $l$ to be compressed.
- step3: run 'Exchange.m' to generate '.mat'.
- step4:  find a a control parameter $\lambda$ for a given desired percent reduction of area $R$ by 'findLambEX.m'.
- step5: run 'GAEX.m' to solve the optimization objective and **directly** generate Verilog and C models of multipliers.

## multipliers

We provide unsigned8b multipliers of different l values and different $\lambda$ values in folder unsigned8b multipliers.

## results

We present two comparison charts that illustrate the relationship between hardware overhead and accuracy of our approximate multiplier compared to other existing approximate multipliers.

![MSE](.\result_pic\MSE.png)

![MAE](.\result_pic\MAE.png)

## scripts

The files 'constraints_comb.tcl' and 'constraints_seq.tcl' work for synthesis of multipliers and accelerators respectively.