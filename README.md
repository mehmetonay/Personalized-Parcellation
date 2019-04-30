# fMRI Guided Personalization of Cortical Parcellation Maps

M S Onay, BSC; U Kucukaslan, BSC; C Ulasoglu Yildiz, PhD; B Acar, PhD

## Abstract

### PURPOSE
To personalize the 3D cortical parcellation of T1w MRI volumes to increase intra-parcel functional homogeneity by means of co-registered fMRI guidance. 
 
### METHODS AND MATERIALS
Inspired by the Demons Algorithm (Med. Image Anal., 2(3), pp 243-260,1998), we propose an iterative 3D deformation of the cortical parcellation map of a given subject, constrained to the cortex, to maximize the pairwise intra-parcel correlation coefficients (CCs) of fMRI BOLD signals (B). The co-registered T1w MRI, fMRI and (148 parcel, P) Destrieux cortical parcellation map are the inputs. The iterations are driven by a vector-field, F(rij), normal to boundary (rij) of Pi and Pj, and zero elsewhere. F is computed according to a local cost that takes into account the change in intra-parcel functional homogeneity and parcel volume change (ΔV). The former uses the CCs between the local signal, B(rij), and the parcel representative signals (Bi, Bj) defined as the in-parcel B that is closest to the left principal singular vector of the matrix composed of parcel's B's. The latter is added with a multiplicative factor λ. We used FreeSurfer and FSL for all co-registration and parcellation steps. We tested the algorithm on 5 normal subjects for λ=0 and 0.25, by measuring the mean of pairwise CCs between in-parcel B's and the mean B per parcel, before and after the algorithm, for each subject. Statistical significance is assessed with paired t-test of CCs per parcel per subject and mean values over subjects is reported. The mean of subject-wise minimum and maximum parcel volumes (Vmin, Vmax) are reported together with the percentage of parcels with ΔV>20%.
 
### RESULTS
The initial mean intra-parcel CC of 0.43±0.04 was increased to 0.52±0.04 (0.51±0.04) with mean p-value 1.21e-5 (4.50e-5) for λ=0 (0.25). The mean [Vmin/Vmax] (mm3) was changed from [229.5/15373.2] to [482.7/10762.9] ([286.9/13797]) with 63% (10%) of parcels having ΔV>20% for λ=0 (0.25).
 
### CONCLUSION
The proposed functionally homogeneizing personalization of cortical parcellations showed a statistically significant improvement in intra-parcel correlations accompanied by a decrease in parcel volume variation. Volume constraint with λ=0.25 does not affect the final intra-parcel CCs while limiting the volume change.
 
### CLINICAL RELEVANCE/APPLICATION
fMRI guided personalization of 3D cortical parcellation maps significantly improves intra-parcel functional homogeneity, and can potentially improve brain network models.

## Installation

In order to use the code, you need the following libraries installed on your computer
- ITK
- Eigen

```bash
pip install numpy
```

Python code can be shown like
```python
import os

print("hey man")
```

## Usage
In order to use our code, 

C++ code can be shown as 
```c++
include <iostream>
main(){
  print("hey man");
}
```

## Header 2
Something about header 2

## License
[MIT](https://choosealicense.com/licenses/mit/)
