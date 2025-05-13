# SPP-SBL: Space-Power Prior Sparse Bayesian Learning for Block Sparse Recovery
This set of Matlab (Vision R2021b) functions contain the core code to reproduce the results of the above paper (https://arxiv.org/).



> ### `SPP_SBL.m`  
>*The main function for the Space-Power Prior Sparse Bayesian Learning algorithm.*
>#### Input:
>- `y`: measurements  
>- `A`: sensing matrix  
>- `sigma`: initial value of the noise standard deviation  
>- `beta`: initial value of the coupling parameter
>#### Output:
>- `x_new`: estimated sparse signal  
>- `eta`: estimated coupling parameters in each iteration  

---

## Demos

A set of demonstration scripts for evaluating the performance of the SPP-SBL algorithm on different types of data.

- `demo_synthetic_heteroscedastic.m`: Demonstrates SPP-SBL on heteroscedastic data from [1].
- `demo_synthetic_chain.m`: Tests the algorithm on synthetic data with chain-type signal [2].
- `demo_synthetic_multipattern.m`: Evaluates SPP-SBL on synthetic multi-pattern signals [3].
- `demo_audioset.m`: Applies the algorithm to audio signal recovery using a subset of real-world audio data [4].
- `demo_image_cameraman.m`: Applies SPP-SBL to image reconstruction using the standard *Cameraman* image.





##

**For bug reports, please contact me at email: yanhaozhang@buaa.edu.cn.**


Authors: Yanhao Zhang, Zhihan Zhu, Yong Xia.

Beihang University,  May, 12, 2025.

##

### References
[1] Y. Zhang, Z. Zhu, and Y. Xia, “Block sparse Bayesian learning: A diversified scheme,” in Advances in Neural Information Processing Systems, vol. 37, pp.129 988–130 017, 2024.

[2] M. Korki, J. Zhang, C. Zhang, and H. Zayyani, “Iterative Bayesian reconstruction of non-iid block-sparse signals,” IEEE Transactions on Signal Processing, vol. 64, no. 13, pp. 3297–3307, 2016.

[3] A. Sant, M. Leinonen, and B. D. Rao, “Block-sparse signal recovery via general total variation regularized sparse Bayesian learning,” IEEE Transactions on Signal Processing, vol. 70, pp. 1056–1071, 2022.

[4] Gemmeke, Jort F., et al. "Audio set: An ontology and human-labeled dataset for audio events." 2017 IEEE international conference on acoustics, speech and signal processing (ICASSP). IEEE, 2017.


