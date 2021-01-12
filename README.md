# Structural templates for imaging EEG cortical sources in infants

This repository contains the code used for the paper entitled "Structural templates for imaging EEG cortical sources in infants" from NeuroImage ([open access](https://www.sciencedirect.com/science/article/pii/S1053811920311678)). If you use these templates, please cite:

1. C. O’Reilly, E. Larson, J. E. Richards, and M. Elsabbagh, “Structural templates for imaging EEG cortical sources in infants,” NeuroImage, vol. 227, p. 117682, Feb. 2021, doi: [10.1016/j.neuroimage.2020.117682](http://dx.doi.org/10.1016/j.neuroimage.2020.117682).
2. J. E. Richards, C. Sanchez, M. Phillips-Meek, and W. Xie, “A database of age-appropriate average MRI templates,” NeuroImage, vol. 124, pp. 1254–1259, Jan. 2016, doi: [10.1016/j.neuroimage.2015.04.055](http://dx.doi.org/10.1016/j.neuroimage.2015.04.055).


<details>
<summary>BibTex entries</summary>

```Bibtex
@article{oreilly_structural_2021,
	title = {Structural templates for imaging {EEG} cortical sources in infants},
	volume = {227},
	issn = {1053-8119},
	url = {http://www.sciencedirect.com/science/article/pii/S1053811920311678},
	doi = {10.1016/j.neuroimage.2020.117682},
	language = {en},
	urldate = {2021-01-12},
	journal = {NeuroImage},
	author = {O'Reilly, Christian and Larson, Eric and Richards, John E. and Elsabbagh, Mayada},
	month = feb,
	year = {2021},
	keywords = {Electroencephalography, Forward model, Infant, Neurodevelopment, Population template, Source reconstruction},
	pages = {117682}
}
@article{richards_database_2016,
	series = {Sharing the wealth: {Brain} {Imaging} {Repositories} in 2015},
	title = {A database of age-appropriate average {MRI} templates},
	volume = {124},
	issn = {1053-8119},
	url = {http://www.sciencedirect.com/science/article/pii/S1053811915003559},
	doi = {10.1016/j.neuroimage.2015.04.055},
	language = {en},
	journal = {NeuroImage},
	author = {Richards, John E. and Sanchez, Carmen and Phillips-Meek, Michelle and Xie, Wanze},
	month = jan,
	year = {2016},
	keywords = {Average MRI templates, Brain development, Lifespan MRI, Neurodevelopmental MRI Database},
	pages = {1254--1259}
}
```

</details>

## Abstract

Electroencephalographic (EEG) source reconstruction is a powerful approach that allows anatomical localization of electrophysiological brain activity. Algorithms used to estimate cortical sources require an anatomical model of the head and the brain, generally reconstructed using magnetic resonance imaging (MRI). When such scans are unavailable, a population average can be used for adults, but no average surface template is available for cortical source imaging in infants. To address this issue, we introduce a new series of 13 anatomical models for subjects between zero and 24 months of age. These templates are built from MRI averages and boundary element method (BEM) segmentation of head tissues available as part of the Neurodevelopmental MRI Database. Surfaces separating the pia mater, the gray matter, and the white matter were estimated using the Infant FreeSurfer pipeline. The surface of the skin as well as the outer and inner skull surfaces were extracted using a cube marching algorithm followed by Laplacian smoothing and mesh decimation. We post-processed these meshes to correct topological errors and ensure watertight meshes. Source reconstruction with these templates is demonstrated and validated using 100 high-density EEG recordings from 7-month-old infants. Hopefully, these templates will support future studies on EEG-based neuroimaging and functional connectivity in healthy infants as well as in clinical pediatric populations.
