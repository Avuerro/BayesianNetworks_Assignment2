# Bayesian Networks Assignment 2

This repo contains the code developed for the second assignment of the Bayesian Networks & Causal Inference course at the Radboud University (NWI-IMC012).

## Brief Introduction

In recent years, forest fires have been occurring more frequently and more intensely.
Not only does this result in increased environmental damage, but these forest fires also become harder to fight.
Therefore, understanding the factors that influence such forest fires is important, 
and these insights might also be used to predict the severity of forest fires.
In the first assignment we set out to create a model based on prior knowledge that best captures the information in the data.
We used the model to predict the size of the burned area.
In the current assignment the goal was to learn the model by applying two different algorithms, the SI-HITON-PC algorithm and the Tabu Search Algorithm.

## Structure of this Repo

This repo contains the following folders:

- `imgs` contains image of the Bayesian networks developed, and an overview of the Fire Weather Index (FWI).
- `plts` contains some examples of plots that were generated in this project.

This repo also contains the following important files:

 - `assignment2.R` contains the fundamental code regarding the project: Building the Bayesian network, fitting it, performing predictions and computing path coefficients.
 - `data_exploration.R` contains all the code that we used for exploring and preprocessing the dataset. 
 - `forestfires.csv` is the original dataset, as it is publicly available on the [UCI Machine Learning Repository](http://archive.ics.uci.edu/ml/datasets/Forest+Fires). This data is used in `data_exploration.R`.
 - `explored_forestfires.csv` contains the preprocessed data, after executing `data_exploration.R`. This file is used in `assignment1.R`.

## Collaboration

This project has been done in a team of three people:
 - [Avuerro](https://github.com/Avuerro)
 - [ElFilosofoX](https://github.com/ElFilosofoX)
 - [JordAI](https://github.com/jordai)

