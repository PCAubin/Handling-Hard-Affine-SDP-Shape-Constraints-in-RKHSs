# Handling-Hard-Affine-SDP-Shape-Constraints-in-RKHSs

Matlab code to reproduce the numerical results "Handling Hard Affine SDP Shape Constraints in RKHSs", Pierre-Cyril Aubin-Frankowski and Zoltan Szabo, JMLR 2022, https://arxiv.org/abs/2101.01519. One can find on https://pcaubin.github.io/ several presentations related to the topic (or freely contact the author).

Catenary corresponds to Experiment 6.1. Control corresponds to Experiment 6.2. Robotic arm corresponds to Experiment 6.3. Economy corresponds to Experiment 6.4. 

The key message is that by simply adding an extra second-order cone term in finitely-many affine constraints (with derivatives), one can guarantee that the affine constraints will hold for ALL the points in a compact set. This provides a solvable finite-dimensional problem even though the original problem was over infinite dimensions with infinitely-many constraints.

The source files in this repository can of course be used for implementing problems. They could otherwise be used as a fruitful inspiration since several files need to be manually adapted to fit the user's specific problem.

#### Abstract
Shape constraints, such as non-negativity, monotonicity, convexity or supermodularity, play a key role in various applications of machine learning and statistics. However, incorporating this side information into predictive models in a hard way (for example at all points of an interval) for rich function classes is a notoriously challenging problem. We propose a unified and modular convex optimization framework, relying on second-order cone (SOC) tightening, to encode hard affine SDP constraints on function derivatives, for models belonging to vector-valued reproducing kernel Hilbert spaces (vRKHSs). The modular nature of the proposed approach allows to simultaneously handle multiple shape constraints, and to tighten an infinite number of constraints into finitely many. We prove the convergence of the proposed scheme and that of its adaptive variant, leveraging geometric properties of vRKHSs. Due to the covering-based construction of the tightening, the method is particularly well-suited to tasks with small to moderate input dimensions. The efficiency of the approach is illustrated in the context of shape optimization, safety-critical control, robotics and econometrics. 

## Reference

Arxiv link: https://arxiv.org/abs/2101.01519

JMLR: https://jmlr.org/papers/v23/21-0007.html

## Requirements

The code is written in Matlab and uses CVX http://cvxr.com/cvx/ to solve the QPQC and SOCP problems through MOSEK (https://www.mosek.com/downloads/). We also use YALMIP (https://yalmip.github.io/) for the robotic arm.
Type the following code with directories replaced by your own to have CVX+MOSEK work
% addpath C:\Users\pierr\Mosek\9.2\toolbox\R2015aom
% cd C:\Users\pierr\Desktop\cvx
% cvx_setup

## Download

```
git clone https://github.com/PCAubin/Handling-Hard-Affine-SDP-Shape-Constraints-in-RKHSs
cd Handling-Hard-Affine-SDP-Shape-Constraints-in-RKHSs
```
Or download from the "Download ZIP" button and unzip it.

Pierre-Cyril Aubin-Frankowski, 01/01/2023
