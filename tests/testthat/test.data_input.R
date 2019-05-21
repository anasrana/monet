context("Reading csv data files")


gene_path <- "gene_exp.csv"
gene_path2 <- "gene_exp_v2.csv"
atac_path <- "atac_seq_tf.csv"
atac_wrong <- "atac_seq.csv"

gene_dt <- monet:::prepGeneExp(gene_path, "gene")
gene_dt2 <- monet:::prepGeneExp(gene_path2, "V1")
atac_dt <- monet:::prepGeneExp(atac_path, "Gene")
atac_dt_w <- monet:::prepGeneExp(atac_wrong, "Gene")

test_that("reading gene expression", {

    expect_that(data.table::key(gene_dt), equals("gene"))
    expect_that(nrow(gene_dt), equals(40))
    expect_that(ncol(gene_dt), equals(6))

    expect_that(data.table::key(gene_dt2), equals("gene"))
    expect_that(nrow(gene_dt2), equals(20))
    expect_that(ncol(gene_dt2), equals(6))

    expect_true(all.equal(gene_dt[gene_dt2$gene, ], gene_dt2))
})

test_that("reading ATACseq", {

    expect_that(data.table::key(atac_dt), equals("gene"))
    expect_that(nrow(atac_dt), equals(50))
    expect_that(ncol(atac_dt), equals(5))

})

test_that("Fails with wrong gene col", {
    expect_error(monet:::prepGeneExp(gene_path, "Gene"),
                 regexp = "No column 'Gene'")

    expect_error(
        monetData_init(gene_path, atac_path, "gende", "asinh", "Gene"),
        regexp = "No column 'gende' found in input")

    expect_error(monet:::prepGeneExp(gene_dt, "Gene"),
                 regexp = "No column 'Gene' found in input data provided.")
})

context("monet data prep")

test_that("testing with path", {
    monet_data <- monetData_init(gene_path, atac_path, "gene", "asinh", "Gene")

    expect_is(monet_data, "monetDataRaw")
    expect_equal(monet_data@no_tpts, 4)
    expect_equal(monet_data@no_genes, 40)
    expect_equal(monet_data@data_test, "Different number of rows")
    expect_equal(monet_data@gene_exp_path, "gene_exp.csv")
})

context("Failing with bad data")

test_that("mismatch in time-points file",
    expect_error(
         monetData_init(gene_path, atac_wrong, "gene", "asinh", "Gene"))
          )

test_that("mismatch in time-points data.table",
  expect_error(
           monetData_init(gene_dt, atac_dt_w, "gene", "none"))
          )

