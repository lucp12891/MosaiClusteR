test_that("pathway/gene-set functions are exported functions", {
	expect_true(is.function(DiffGenes))
	expect_true(is.function(DiffGenesSelection))
	expect_true(is.function(FindGenes))
	expect_true(is.function(Geneset.intersect))
	expect_true(is.function(Geneset.intersectSelection))
	expect_true(is.function(PathwayAnalysis))
	expect_true(is.function(Pathways))
	expect_true(is.function(PathwaysIter))
	expect_true(is.function(PathwaysSelection))
	expect_true(is.function(PreparePathway))
	expect_true(is.function(SharedGenesPathsFeat))
	expect_true(is.function(SharedSelection))
	expect_true(is.function(SharedSelectionLimma))
	expect_true(is.function(SharedSelectionMLP))
})

test_that("DiffGenes errors on missing limma", {
	skip_if(requireNamespace("limma", quietly = TRUE))
	expect_error(DiffGenes(list(), geneExpr = NULL), "limma")
})

test_that("DiffGenesSelection errors on missing limma", {
	skip_if(requireNamespace("limma", quietly = TRUE))
	expect_error(DiffGenesSelection(list(), Selection = "a"), "limma")
})

test_that("Geneset.intersect errors on missing plyr", {
	skip_if(requireNamespace("plyr", quietly = TRUE))
	expect_error(Geneset.intersect(list()), "plyr")
})

test_that("Geneset.intersectSelection errors on missing plyr", {
	skip_if(requireNamespace("plyr", quietly = TRUE))
	expect_error(Geneset.intersectSelection(list()), "plyr")
})

test_that("PathwayAnalysis errors on missing MLP", {
	skip_if(requireNamespace("MLP", quietly = TRUE))
	expect_error(PathwayAnalysis(list()), "MLP")
})

test_that("Pathways errors on missing MLP", {
	skip_if(requireNamespace("MLP", quietly = TRUE))
	expect_error(Pathways(list()), "MLP")
})

test_that("PathwaysIter errors on missing MLP", {
	skip_if(requireNamespace("MLP", quietly = TRUE))
	expect_error(PathwaysIter(list()), "MLP")
})

test_that("PathwaysSelection errors on missing MLP", {
	skip_if(requireNamespace("MLP", quietly = TRUE))
	expect_error(PathwaysSelection(Selection = "a"), "MLP")
})
