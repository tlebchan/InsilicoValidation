# Mutation calling optimisation in bulk WES using in silico approach
## Aim
Explore the use of in silico methods to fine-tune and validate a pipeline for identifying of somatic mutations in Whole Exome Sequencing samples
## Objectives
* Develop pipelines for sample dilution and mutation introduction
* Compare reproducibility of mutation filtration approaches on artificially diluted samples
* Introduce mutations in real samples and compare performance of different mutation approaches
## Methods
In order to evaluate BostonGene mutation calling pipeline performance two artificial approaches were utilised:

1. **Artificial dilution**

![](images/artificial_dilution.png)

2. **Mutation introduction**

![](images/mutation_introduction.png)
## Results

1. Somatic mutation VAF for patients with downsampled purity

![](images/VAF.png)

2. Reproducibility against purity for several mutation filter approaches

![](images/reproducability.png)

3. Recall and precision on in silico generated data

![](images/recall.png)

![](images/precision_corrected.png)
## Conclusions
* BostonGene internal filters perform much better across the range of purity available, compared to basic filtrations
* The current filter works much better at high purity, while the new filter copes well with low purity samples
* While both filters do well in calling precision, the new filter has a much better performance in terms of recall especially in the low purity range
## References
