test_that("forming pseudobulk works", {
  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(100))
  colnames(Y) <- paste0("cell_", seq_len(50))
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = Y, logcounts = log(Y + 1)),
                                                    colData  = data)
  expect_error(pseudobulk_sce(sce))
  expect_error(pseudobulk_sce(sce, NULL))
  expect_error(pseudobulk_sce(sce, 1:2))
  expect_error(pseudobulk_sce(sce, pseudobulk_by = city, age = diff(age)))

  pseudobulk_sce(sce, pseudobulk_by = city, age = mean(age), head(fav_food, n = 1))

  tmp <- pseudobulk_sce(sce, pseudobulk_by = paste0(city, "-", fav_food))
  cd <- SummarizedExperiment::colData(tmp)
  expect_equal(rownames(cd), as.character(cd$`paste0(city, "-", fav_food)`))

  tmp2 <- pseudobulk_sce(sce[,1], pseudobulk_by = city, age = mean(age), fav_food = head(fav_food, n = 1))
  colnames(tmp2) <- "cell_1"
  SummarizedExperiment::colData(tmp2)$city <- as.character(SummarizedExperiment::colData(tmp2)$city)
  expect_equal(SummarizedExperiment::colData(sce[,1])[,c("city", "age", "fav_food")],
               SummarizedExperiment::colData(tmp2))

  SummarizedExperiment::colData(sce)$fact <- factor(sample(letters[1:3], 50, replace = TRUE),
                                                    levels = letters[1:4])

  tmp3 <- pseudobulk_sce(sce, pseudobulk_by = fact, age = mean(age))
  expect_equal(levels(SummarizedExperiment::colData(tmp3)$fact), letters[1:4])

  tmp4 <- pseudobulk_sce(sce, pseudobulk_by = fact, aggregation_functions = list(counts = matrixStats::rowMins))
  expect_equal(unname(SummarizedExperiment::assay(tmp4, "counts")[,"a"]),
               matrixStats::rowMins(SummarizedExperiment::assay(sce[,sce$fact == "a"], "counts")))
  expect_equal(unname(SummarizedExperiment::assay(tmp4, "logcounts")[,"b"]),
               matrixStats::rowMeans2(SummarizedExperiment::assay(sce[,sce$fact == "b"], "logcounts")))
})



