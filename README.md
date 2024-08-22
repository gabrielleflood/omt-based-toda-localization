# Optimal Mass Transport Based TDOA Localization
Repository for Matlab code for the paper Multi-Source Localization and Data Association for Time-Difference of Arrival Measurements[^1]. The problem considered is that of localizing multiple signal sources based on time-difference of arrival (TDOA) measurements. The source positions and signals are unknown, while the receiver positions are assumed to be known as well as the TDOA measurements (which e.g. can be achieved by correlating received signals).

This method performs joint localization and data association by means of an optimal mass transport formulation. The method can be divided into two steps. First, a candidate source position set is constructed using a multilateration solver[^2] based on minimal sets of receiver pairs. Then, the best candidate positions are found together with the association between these positions and the measured TDOA values. Finally, as the association is found, local optimization can be performed, minimizing the error between the measured TDOAs and the relative distances between the found source positions and the receivers.

## Run the code
To run a simple experiment, where TDOAs are measured from a number of *R* simulated receiver positions and *S* simulated source positions, run the script 
```
simple_experiment.m
```
The TDOA measurements *tdoas_measured* and the receiver positions *r* can also be provided if a true experiment has been run. The method then computes and outputs the source positions.

If the true source positions *s* and the true TDOA values are known, the parameter *exists_gt* can be set to 1, and then the Euclidean error for the result as well as the CRLB value for the setup are also provided. This is default for the simulated setting.

The 

### Experiment from the paper
The script
```
exp_err_wrt_tdoapert_2024.m
```
reproduces the experiment in the paper where the error of the proposed method is evaluated for different levels of noise in the TDOA measurements. Here, the results for the method SMTL[^3] is also included, for comparison. Note that this is our implementation of that method. This does require the Matlab package CVX[^4]. For information on how to install this, please visit their [homepage](https://cvxr.com/cvx). Un disable the SMTL computations, comment out line 174
```
SMTL_mean_distance = get_SMTL_mean_distance(ggrid,r,s,S,tdoas_measured, SMTL_lambda);
```

## Citation 
If you use this in your research, please cite
```
@inproceedings{flood2024multi,
  title={Multi-Source Localization and Data Association for Time-Difference of Arrival Measurements},
  author={Flood, Gabrielle and Elvander, Filip},
  booktitle={2024 32th European Signal Processing Conference (EUSIPCO)},
  year={2024},
  organization={IEEE}
}
```

[^1]: Flood, G. and Elvander, F., 2024. Multi-Source Localization and Data Association for Time-Difference of Arrival Measurements. In *Proceedings of European Signal Processing Conference (EUSIPCO)*.
[^2]: Åström, K., Larsson, M., Flood, G. and Oskarsson, M., 2021.Extension of time-difference-of-arrival self calibration solutions using robust multilateration. In *Proceedings of European Signal Processing Conference (EUSIPCO)*, pp. 870–874.
[^3]: Jamali-Rad, H. and Leus, G., 2013. Sparsity-aware multi-source TDOA localization. *IEEE Transactions on Signal Processing*, vol. 61, no. 19, pp. 4874–4887.
[^4]: Grant, M. and Boyd, S. 2014. CVX: Matlab Software for Disciplined Convex Programming, version 2.1. https://cvxr.com/cvx
