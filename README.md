# KEGGconnect
### Connects interactions between genes across KEGG Pathways
This tool is still WIP, but can currently produce a Cytoscape sif file of relationships between genes across KEGG pathways

### Table of Contents
* [Installation](#in)
* [Running](#run)
  * [Arguments](#arg)
  * [Example](#eg)
* [To-do](#todo)

<br><br>

<a name="in"/>

## Installation
Install the latest version of [R](https://cran.rstudio.com).

This script also requires the [KEGGREST](https://bioconductor.org/packages/release/bioc/html/KEGGREST.html) package and [KEGGgraph](https://www.bioconductor.org/packages/release/bioc/html/KEGGgraph.html) package.

Then clone the repository and run the script using RScript (see [Running](#run)).

<a name="run">

## Running
Use Rscript kconnect.r to run, or pass command-line arguments as outlined in [Arguments](#arg).

<a name="arg"/>

### Arguments
The script is run according to the format [option] [organism] [filename] [outputdir]

<a name="eg"/>

### Example
The tool can be tested with the example file provided. Simply run the command below in the main directory to create a Cytoscape sif file in the directory "test".
Rscript kconnect.r cg hsa ./exampledata/genelist.txt test
         

<a name="todo"/>

## To-do
- Improve runtime
- Collapse edges and simply note the corresponding pathways in an output file (currently gives an edge for every pathway between pairs of genes)
- Special formatting of nodes and edges
