Todo: upload large files:

results/raw_data/20200904_IMP_diff_models_all_eval.csv ~38 MB

data/edgestats/tw-inststatsedge_rt_by_ot.csv ~50MB


[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# The impact of passive social media viewers in influence maximization

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

This repository contains codes (written in Julia) and results that were used in the research 
reported on in the paper: Kahr Michael, Leitner Markus, Ljubic Ivana. The impact of passive social media viewers in influence maximization, 2024,
(https://doi.org/10.1287/ijoc.2023.0047). 



**Important: This repository will not receive any updates. Moreover, the code provided here is unlikely to work on newer machines because it uses outdated Julia package and CPLEX versions (see subsection Requirements below). An updated version (and potential future updates) can be found here:
https://github.com/m-kahr/mnlcimp**.

## Cite

To cite the contents of this repository, please cite both the paper and this repository, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0047 

https://doi.org/10.1287/ijoc.2023.0047.cd 

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{kahr2024imp,
  author =        {Kahr Michael, Leitner Markus, Ljubic Ivana},
  publisher =     {INFORMS Journal on Computing},
  title =         {{The impact of passive social media viewers in influence maximization}},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0047.cd},
  url =           {https://github.com/INFORMSJoC/2023.0047},
}  
```

## Sources Codes

All codes are written in Julia 1.4, and tested on Ubuntu 18.04 (and Julia 1.1 and 1.4). The used packages and dependencies are provided in the file `Project.toml` in the [src](src/) folder (where the main code is located). The code for data evaluation is located in the [src/evaluation](src/evaluation/) folder (where all scripts for evaluation are located).

### Requirements (main code)
- A UNIX-based machine, because the code accesses the file `/proc/self/stat` to keep track of the used memory during the optimization procedure. 
- IBM CPLEX 12.9 or 12.10. Newer versions will not work because they are not supported by `CPLEX.jl v0.6.6` which is used in the code. Moreover, legacy callbacks are implmented (not the generic ones).

### Requirements (evaluation)
- a set up `Python` environment for use with `PyCall` and `PyPlot`
- Python packages: `matplotlib 3.3.4`, `seaborn 0.11.1`, `numpy 1.20.2`, and `tikzplotlib 0.9.8`. Newer versions of the latter packages might work, but are not tested.

### Usage (main code)
Set the path to your CPLEX installation directory at the beginning of file `mainLogit.jl` in the [src](src/) folder. Excecute `mainLogit.jl` in your favorite editor (e.g., VScode). You can set the instance parameters in `params.jl`, and instance graphs in `misc.jl` in function `setInstancePathAndName!( params::Dict )`. Alternatively, you can start the code from the shell with the `--console false` option, e.g.,
```
/path/to/julia-1.4.0/bin/julia /path/to/mainLogit.jl --ipath /path/to/instances/ --ifile TW-orms-20200807-Y2020-anonymized --precompile true --writeoutput true --nscenarios 100 --Nscenarios 1000 --kL 0 --kF 5 --nSSAItr 10 --console false 
```
See `params.jl` for a description of the parameters and other available parameters (not listed here). Note some of the latter are unused (and marked as unused in that file).

The results are written in csv format in the [results](folder/) folder. Each run two files are written. The first one contains the headlines of for the data in the second one. 

### Usage (evaluation)
The number of experiments we conducted evolved over time. Thus there exist several different raw data result files in [results/rawdata](results/raw_data/). To excecute a specific analyses choose one of the (commented) types in file `analyze_results.jl`. For instance, `rtype = "ANAL_IMP"` analyizes IMP variants (see the presice description next to the available types). 

## Data
The [data](data/) folder contains three subfolders:

### edgestats
This folder contains the file `tw-inststatsedge_rt_by_ot.csv` from which forwarding probabilities (extracted from twitter data) are sampled when instance graphs from the literature are used, namely `msg-college.im`, `msg-email-eu.im`,`soc-advogato.im`, and `soc-anybeat.im`.

### instances
This folder contains three types of instances:
- `D-n25-k4-b0.3-beta2.0-18.0-i1` an artifical instance only used for precompiling functions in a dummy run to avoid measuring compilation times.
- `msg-*`, `soc-*` are instance graphs from the literature. These graphs only contain edge lists, with the number of nodes and arcs reporterd on top.
- `TW-*` are instances we extracted from X (former Twitter). The user IDs and the user names are anonymized. The instances are organized in user profile data and user relations 
   
   - **user profile data:**
   - `ind`: node index
   - `tid`: Twitter ID (anonymized)
   - `tname`: user profile name (anonymized)
   - `friends`: # of friends
   - `followers`: # of followers
   - `otweets`: # of original tweets
   - `retweets/quotes`: # of retweets and quotes user profile does
   - `replies`: # of replies user profile does
   - `likes`: # of likes user profile receives
   - **user relations**
   - `ind`: user relation (or arc) index
   - `influencing_user`: node index of influencing user
   - `influenced_user`: node index of influenced user
   - `retweets|quotes`: # of influencing retweets and quotes (e.g., influenced user retweets content from influencing user)
   - `mentions`: # of mentions (e.g., influenced user mentions influencing user)
   - `answers`: # of answers (e.g., influenced user answers to influencing user)

### leader_seed_sets
This folder contains leader seed sets computed with all studied different methods F, O, T, R25, R50, R75 and $|\Omega'| = 100$, using the (generalized) Benders decomposition approach. The first number corresponds to the seed set cardinality, and the remaining numbers to the node ID's in the seed set.

## Results

The [results](results/) contains the raw results from our computational experiments, and the plots we derived from the latter data and which are used in the paper. The raw results are used in the [src/evaluation](src/evaluation/) scripts.


## Scripts

Folder [scripts](scripts/) contains bash scripts we used to start our computational experiments on a high-performance cluster via `qsub`.



## Contact & Support
See file `AUTHORS`.
