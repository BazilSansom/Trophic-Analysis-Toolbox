# Python module for trophic-analysis and network visualisation

*trophic_tools* is a Python module that provides tools for tophic analysis as introduced by [1] 
and also provides some tools developed for visualising networks using these methods.

Functions included are:

- *trophic_levels*       : returns trophic levels as per [1]

- *trophic_incoherence*  : returns trophic incoherence as per [1] (where trophic coherence is 1-incoherence)

- *trophic_layout*       : returns a layout where y-possitions are given by trophic levels, and x-possitions based on a modified force-dircted graph drawing algorithm to spread nodes out on the x-axis

- *trophic_plot*         : plots network according to *trophic_layout*

All functions take networkx graph object as input. A full list of dependencies is provided in *requirements*.

A small demo of the use of these functions is provided in *demo_trophic_analysis*.

Refs:

- [1] MacKay, Johnson & Sansom (2020), "How directed is a directed network", *Royal Society Open Science*, vol. 7, issue 9. doi: https://doi.org/10.1098/rsos.201138

