# Arbitrage with Power Factor Correction using Energy Storage

Co-optimization formulation with energy arbitrage and power factor correction is discussed in our paper:
Cite: Hashmi, Md Umar, Deepjyoti Deka, Ana Busic, Lucas Pereira, and Scott Backhaus. "Arbitrage with power factor correction using energy storage." IEEE Transactions on Power Systems (2020).

INRIA, DI ENS, Ecole Normale Sup\'erieure, CNRS, PSL Research University, Paris, France

Contact: umar.hashmi123@gmail.com

This paper proposes co-optimization formulations and also compares them with benchmarks. Five models are provided:

## Model 1: Only arbitrage
The storage is considered to be only performing arbitrage under NEM 1.0 (Net-energy metering with equal buying and selling price of electricity)

## Model 2: Arbitrage with hard power factor constraint using McCormick Relaxations
This model performs arbitrage and PFC. However, provides an infesible output if any time instant PF cannot be corrected due to storage constraint or instantaneous load values.

## Model 3: Receding horizon PFC with Arbitrage
In order to demonstrate PFC being myopic we use a model which myopically performs PFC and lookahead is used for arbitrage. It is observed that the results converge with Model 4. (please refer to the paper for details)

## Model 4: Arbitrage with linear penalty on PF violation


## Model 5: Arbitrage with linear penalty on PF violation & cost imposed on converter usage



![alt text](https://github.com/umar-hashmi/powerfactorcorrection/blob/master/modelsComp.png)

