# Find an element in a data structure

The function `FindElement` is used internally in the `PreparePathway`
function but might come in handy for other uses as well. Given the name
of an object, the function searches for that object in the data
structure and extracts it. When multiple objects have the same name, all
are extracted.

## Usage

``` r
FindElement(what = NULL, object = NULL, element = list())
```

## Arguments

- what:

  A character string indicating which object to look for. Default is
  NULL.

- object:

  The data structure to look into. Only the classes data frame and list
  are supported. Default is NULL.

- element:

  Not to be specified by the user.

## Value

The returned value is a list with an element for each object found. The
element contains everything the object contained in the original data
structure.

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- list(a = 1, b = list(TopDE = 42, c = list(TopDE = 7)))
FindElement(what = "TopDE", object = obj)
} # }
```
