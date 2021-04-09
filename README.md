
## Watershed PeakCaller

The Watershed Algorithm is adapted as a peak caller to visualise periodicity of cell cycle-regulated genes of the yeast Saccharomyces cerevisiae. 

Gene_expression data: https://pubmed.ncbi.nlm.nih.gov/9843569/

  Run the program with:

```python3 Watershed_01.py yeast_cell_cycle.txt 3 YER111C YMR043W YLR131C```



 1. Adjacency Threshold - specifies the neighborhood size of the current datapoint.
    Preferably, value between 1 to 5.
    
 2. Find systematic name of the gene/s of interest from yeastgenome.org 
    *Eg.  SWI4 - YER111C* 
