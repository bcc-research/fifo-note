# Code for A Note on the Welfare Gap in Fair Ordering
This repository contains the code used to generate the figures in the paper
[A Note on the Welfare Gap in Fair Ordering](https://angeris.github.io/papers/note-on-fifo.pdf).

## Running the script
Clone the repository, navigate to the folder in your terminal, and (assuming Julia is installed), simply run

```bash
julia sims.jl
```

The script will activate the environment specified by `Project.toml`, 
install all required packages, and then run the program.
Figures are output to the `figs/` folder.

## Changing parameters
The parameters are set to those used in the paper. These can be changed in the `sims.jl` file.