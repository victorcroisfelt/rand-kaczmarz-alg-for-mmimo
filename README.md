# Randomized Kaczmarz Algorithm for Massive MIMO Systems with Channel Estimation and Spatial Correlation

This is a research-oriented code package that is primarily intended to allow readers to replicate the results of the articles mentioned below and also encourage and accelerate further research on this topic:

Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão, "Kaczmarz Precoding and Detection for Massive MIMO Systems,"
IEEE Wireless Communications and Networking Conference,, pp. 1–6, Marrakech, Morroco, 2019. Publication is expected soon.

Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão, "[Randomized Kaczmarz Algorithm for Massive MIMO Systems with Channel Estimation and Spatial Correlation](https://doi.org/10.1002/dac.4158)", Wiley International Journal of Communication Systems, 2019;e4158. Also available on arXiv: https://arxiv.org/abs/1904.04376.

The package is based on the Matlab language and can, in fact, reproduce all the numerical results and figures discussed in the article. To contextualize, in the sequel, we present the abstract of the article and other important information.

I hope this content helps in your research and contributes to building the precepts behind open science. Remarkably, in order to boost the idea of open science and further drive the evolution of science, we also motivate you to share your published results to the public.

If you have any questions and if you have encountered any inconsistency, please do not hesitate to contact me via victorcroisfelt@gmail.com.

## Abstract
To exploit the benefits of massive multiple-input multiple-output (M-MIMO) technology in scenarios where base stations (BSs) need to be cheap and equipped with simple hardware, the computational complexity of classical signal processing schemes for spatial multiplexing of users shall be reduced. This calls for suboptimal designs that perform well the combining/precoding steps and simultaneously achieve low computational complexities. An approach based on the iterative Kaczmarz algorithm (KA) has been recently investigated, assuring well execution without the knowledge of second order moments of the wireless channels in the BS, and with easiness, since no tuning parameters, besides the number of iterations, are required. In fact, the randomized version of KA (rKA) has been used in this context due to global convergence properties. Herein, modifications are proposed on this first rKA-based attempt, aiming to improve its performance–complexity trade-off solution for M-MIMO systems. We observe that long-term channel effects degrade the rate of convergence of the rKA-based schemes. This issue is then tackled herein by means of a hybrid rKA initialization proposal that lands within the region of convexity of the algorithm and assures fairness to the communication system. The effectiveness of our proposal is illustrated through numerical results which bring more realistic system conditions in terms of channel estimation and spatial correlation than those used so far. We also characterize the computational complexity of the proposed rKA scheme, deriving upper bounds for the number of iterations. A case study focused on a dense urban application scenario is used to gather new insights on the feasibility of the proposed scheme to cope with the inserted BS constraints.

## Content
The codes provided herein can be used to simulate the Figs. 1 to 5 contained in the full article (Figs 1 and 2 of the conference article). This is done by running the scripts that have "simulation" in their names, while those with "function" in their names are called by the main scripts. Further details about each file can be found inside them.

## Acknowledgments
This research was supported in part by the Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) of Brazil, under grant 371462/2019-3.

## Citing this Repository and License
This code is subject to the MIT license. If you use any part of this repository for research, please consider to cite our aforementioned work.

```bibtex
@article{https://doi.org/10.1002/dac.4158,
author = {Rodrigues, Victor Croisfelt and Marinello Filho, José Carlos and Abrão, Taufik},
title = {Randomized Kaczmarz algorithm for massive MIMO systems with channel estimation and spatial correlation},
journal = {International Journal of Communication Systems},
volume = {32},
number = {18},
pages = {e4158},
keywords = {combining, computational complexity, Kaczmarz algorithm, massive MIMO, precoding},
doi = {https://doi.org/10.1002/dac.4158},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/dac.4158},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/dac.4158},
note = {e4158 dac.4158},
year = {2019}
}
