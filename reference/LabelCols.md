# Helper that colours specific dendrogram leaves

Internal helper applied via `dendrapply` that colours the leaves and
edges of a dendrogram for one or two selections of objects. Used by
`LabelPlot`.

## Usage

``` r
LabelCols(x, Sel1, Sel2 = NULL, col1 = NULL, col2 = NULL)
```

## Arguments

- x:

  A node of a dendrogram (passed by `dendrapply`).

- Sel1:

  The selection of objects to be colored.

- Sel2:

  An optional second selection to be colored. Default is NULL.

- col1:

  The color for the first selection. Default is NULL.

- col2:

  The color for the optional second selection. Default is NULL.

## Value

The dendrogram node `x` with coloured `nodePar` and `edgePar`
attributes.
