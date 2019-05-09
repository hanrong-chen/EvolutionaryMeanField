EvolutionaryMeanField

This repository contains all codes, written in Julia, that produced the data and figures in the manuscript
```
[Mean-field computational approach to HIV dynamics on a fitness landscape](https://www.biorxiv.org/content/10.1101/518704v2)
```
authored by Hanrong Chen and Mehran Kardar.

# Computational pipeline

1. LANL HIV sequence database
We downloaded HIV-1 group M subtype B gag protein sequences from the [LANL HIV sequence database](https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html).
The sequences are provided in the ```LANL data``` folder.

2. Data pre-processing
In ```process_p24_sequences.ipynb```, we processed the sequences to obtain the one and two-point mutational correlations of p24 gag, which are provided in the ```p24 landscape``` folder.

3. Inference of p24 prevalence landscape
We used [available software](https://bartonlab.ucr.edu/projects/ACE/) implementing the [adaptive cluster expansion algorithm](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.106.090601) to infer the prevalence landscape of p24.
The results, for two different regularization strengths, are provided in the ```p24 landscape``` folder. See the end of ```process_p24_sequences.ipynb``` for a comparison between the obtained fields and couplings. We used ```gh-0.001``` in our subsequent analyses.

4. Patient data
Data for the example patient used in the text (taken from [1]) are provided in the ```patient data``` folder.

5. Run EMF method
The evolutionary mean-field (EMF) method is a computationally tractable, high-recombination-rate model of HIV dynamics on a fitness landscape. In ```EMF.ipynb```, we applied this method to the dynamics of p24 in the example patient described in step 4, whose infecting strain and T cell responses against the virus are known. All simulations and figures in the main text were produced by this code.

6. Comparison with Wright–Fisher method
In ```Wright-Fisher.ipynb```, we performed fixed-population-size stochastic simulations (similar to [2]) as a comparison with the results obtained in step 5.

References:
[1] Liu, M. K. P., Hawkins, N., Ritchie, A. J., Ganusov, V. V., Whale, V., Brackenridge, S., ... & Goonetilleke, N. (2013). Vertical T cell immunodominance and epitope entropy determine HIV-1 escape. Journal of Clinical Investigation, 123(1), 380–393.
[2] Barton, J. P., Goonetilleke, N., Butler, T. C., Walker, B. D., McMichael, A. J., & Chakraborty, A. K. (2016). Relative rate and location of intra-host HIV evolution to evade cellular immunity are predictable. Nature Communications, 7, 11660.