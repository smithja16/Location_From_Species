# Location_From_Species

This code accompanies the article: "Smith J.A. (2026) Methods for estimating location from spatial patterns in species composition: a fishing location case study. Ecological Modelling."
https://www.sciencedirect.com/science/article/pii/S0304380026000669. 
Pre-print here: https://www.biorxiv.org/content/10.64898/2025.11.30.691468v1.abstract.

The article and this code evaluate a few methods for predicting survey location from patterns in species composition.

It simulates count data with specified properties (such as the strength of the spatial signal in those counts)
with both accurate and inaccurate location data. Three approaches are used to refine those locations given the
observed species compositions for each observation.

This is a problem for some fisheries catch data, so that is the context of the code.

The top script that requires running is 'Estimating Location from Species - TOP SCRIPT.R'.

![Picture1](https://github.com/user-attachments/assets/fc73d30b-958c-4a30-a5aa-262dd56d9de3)
