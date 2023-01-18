test_that("forming pseudobulk works", {
  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(100))
  colnames(Y) <- paste0("cell_", seq_len(50))
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = Y, logcounts = log(Y + 1)),
                                                    colData  = data)
  expect_error(pseudobulk(sce))
  expect_error(pseudobulk(sce, NULL))
  expect_error(pseudobulk(sce, vars(1:2)))
  expect_error(pseudobulk(sce, group_by = vars(city), age = diff(age)))

  pseudobulk(sce, group_by = vars(city), age = mean(age), head(fav_food, n = 1))

  tmp <- pseudobulk(sce, group_by = vars(city, fav_food))
  cd <- SummarizedExperiment::colData(tmp)
  expect_equal(rownames(cd), as.character(paste0(cd$city, ".", cd$fav_food)))

  tmp2 <- pseudobulk(sce[,1], group_by = vars(city), age = mean(age), fav_food = head(fav_food, n = 1))
  colnames(tmp2) <- "cell_1"
  SummarizedExperiment::colData(tmp2)$city <- as.character(SummarizedExperiment::colData(tmp2)$city)
  expect_equal(SummarizedExperiment::colData(sce[,1])[,c("city", "age", "fav_food")],
               SummarizedExperiment::colData(tmp2))

  SummarizedExperiment::colData(sce)$fact <- factor(sample(letters[1:3], 50, replace = TRUE),
                                                    levels = letters[1:4])

  tmp3 <- pseudobulk(sce, group_by = vars(fact), age = mean(age))
  expect_equal(levels(SummarizedExperiment::colData(tmp3)$fact), letters[1:4])
  expect_equal(tmp3$age[1], mean(sce$age[sce$fact == tmp3$fact[1]]))

  tmp4 <- pseudobulk(sce, group_by = vars(fact), aggregation_functions = list(counts = matrixStats::rowMins))
  expect_equal(unname(SummarizedExperiment::assay(tmp4, "counts")[,"a"]),
               matrixStats::rowMins(SummarizedExperiment::assay(sce[,sce$fact == "a"], "counts")))
  expect_equal(unname(SummarizedExperiment::assay(tmp4, "logcounts")[,"b"]),
               matrixStats::rowMeans2(SummarizedExperiment::assay(sce[,sce$fact == "b"], "logcounts")))

  pca <- stats::prcomp(t(SummarizedExperiment::assay(sce,"logcounts")), rank. = 2)
  SingleCellExperiment::reducedDim(sce, "PCA") <- pca$x
  SingleCellExperiment::reducedDim(sce, "PCA2") <- SingleCellExperiment::LinearEmbeddingMatrix(pca$x, pca$rotation)
  tmp5 <- pseudobulk(sce, group_by = vars(fav_food))
  expect_equal(dim(SingleCellExperiment::reducedDim(tmp5, "PCA")), c(3, 2))
  expect_equal(dim(SingleCellExperiment::reducedDim(tmp5, "PCA2")), c(3, 2))

  # Try advanced metaprogramming features
  fav_food <- "test"
  pseudobulk(sce, group_by = vars(city), age = mean(.data$age), .env$fav_food)

  f <- function(arg, arg2){
    pseudobulk(sce, group_by = vars({{arg}}), mean({{arg2}}))
  }
  f(city, arg2 = age)

  pseudobulk(sce, group_by = vars(city), aggregation_functions = list(.default = Matrix::rowSums))
  # Try context specific functions
  pseudobulk(sce, group_by = vars(city), n = n())
})


test_that("function labelling works", {
  res1 <- get_aggregation_function("test", aggregation_functions = list(.default = "rowSums2", test = "rowMeans2"))
  expect_equal(res1$fnc, MatrixGenerics::rowMeans2)
  expect_equal(res1$label, "rowMeans2")

  res1 <- get_aggregation_function("test", aggregation_functions = list(.default = "rowSums2", test = rowMeans))
  expect_equal(res1$fnc, rowMeans)
  # I don't know if there is anyway to get anything more helpful
  expect_equal(res1$label, "custom function")

})

