# Create a color palette to be used in the plots

In order to facilitate the visualization of the influence of the
different methods on the clustering of the objects, colours can be used.
The function `ColorPalette` is able to pick out as many colours as there
are clusters. This is done with the help of the `ColorRampPalette`
function of the grDevices package

## Usage

``` r
ColorPalette(colors = c("red", "green"), ncols = 5)
```

## Arguments

- colors:

  A vector containing the colors of choice

- ncols:

  The number of colors to be specified. If higher than the number of
  colors, it specifies colors in the region between the given colors.

## Value

A vector containing the hex codes of the chosen colors.

## Examples

``` r

Colors1<-ColorPalette(c("cadetblue2","chocolate","firebrick2",
"darkgoldenrod2", "darkgreen","blue2","darkorchid3","deeppink2"), ncols=8)
```
