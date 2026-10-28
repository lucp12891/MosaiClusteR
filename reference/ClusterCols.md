# Helper that colours dendrogram leaves by cluster membership

Internal helper applied via `dendrapply` that colours the leaves and
edges of a dendrogram according to the cluster a leaf belongs to, or to
a specific selection of objects. Used by `ClusterPlot`.

## Usage

``` r
ClusterCols(x, Data, nrclusters = NULL, cols = NULL, colorComps = NULL)
```

## Arguments

- x:

  A node of a dendrogram (passed by `dendrapply`).

- Data:

  The clustering result (an `agnes`/`hclust` object) used to determine
  cluster membership.

- nrclusters:

  Optional. The number of clusters to cut the dendrogram in. Default is
  NULL.

- cols:

  The colours for the clusters if nrclusters is specified. Default is
  NULL.

- colorComps:

  Optional. A character vector of objects to highlight in red. Default
  is NULL.

## Value

The dendrogram node `x` with coloured `nodePar` and `edgePar`
attributes.
