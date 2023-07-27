## Abstract:
Predicting the eventual volume of tumor cells, that might proliferate from a given tumor, can help in cancer early detection and medical procedures planning to prevent their migration to other organs. In this work, a new statistical framework is proposed using Bayesian techniques for detecting the eventual volume of cells expected to proliferate from a Glioblastoma (GBM) tumor. Specifically, the tumor region is first extracted using a parallel image segmentation algorithm. Once the tumor region is determined, we are interested in the number of cells that can proliferate from this tumor until its survival time. For this, we construct the posterior distribution of the tumor cell numbers based on proposed likelihood function and a certain prior. We also determine the corresponding probability that no tumor cell goes undetected when we find the ultimate eventual volume. The model so developed gives excellent results on our dataset. It is expected that the model will work on any dataset where the changes are not measured with regards time. This is particularly crucial for predicting tumor volumes for a deadly tumor like glioblastoma. Furthermore, we extend the detection model and conduct a Bayesian regression analysis by incorporating radiomic features to discover %that their inclusion enhances the chances of those non-tumor cells remaining undetected. The main focus of the study was to develop a time-independent prediction model that can reliably predict the ultimate volume a malignant tumor of the fourth-grade severity can attain, and can also determine if incorporation of radiomic properties of the tumor enhances the chances of no malignant cells remaining undetected.

## Code Usage:
### (a) Running the proposed method:
**Simulation Run:**
> Total number of tests done = 10
> 
> Frequency of Malignant cells detected at each test:
>
> ```ruby
> aj=c(61,9,34,83,14,21,5,11,8,5,4) ## malignant cells detected at each test
> ```
>
> 
> Frequency of Overall tumor cells detected at each test:
>
> ```ruby
> nj=c(226,107,119,307,74,21,33,225,268,125) ## overall tumor cells detected at each test
> ```
>
> Detecting the 

