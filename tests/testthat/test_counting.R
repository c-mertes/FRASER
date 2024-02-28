context("Test counting")

test_that("Count junctions", {
    out <- capture.output(attach(test_generate_count_example()))

    expect_is(test_fdsSample3, "FraserDataSet")

    # test how many ranges we found
    expect_equal(length(test_rangeFDS), 3)
    expect_equal(length(nonSplicedReads(test_rangeFDS)), 5)

    # test the manually counted positions
    expect_equal(as.vector(counts(test_rangeFDS, type="j")), test_rawCountsJ)
    # expect_equal(as.vector(counts(test_rangeFDS, type="ss")), test_rawCountsSS)
})

test_that("Strand spcific counting", {
    suppressMessages(attach(test_generate_strand_specific_count_example()))
    
    expect_is(test_fdsSample3_stranded, "FraserDataSet")
    
    # test how many ranges we found
    expect_equal(length(test_rangeFDS_stranded), 5)
    expect_equal(length(nonSplicedReads(test_rangeFDS_stranded)), 9)
    
    # test the manually counted positions
    expect_equal(as.vector(counts(test_rangeFDS_stranded, type="j")), 
                    test_rawCountsJ_stranded)
    expect_equal(as.vector(counts(test_rangeFDS_stranded, type="ss")),
                    test_rawCountsSS_stranded)
    
    # check for non empty chromosome but no split reads present
    fds <- createTestFraserSettings()
    strandSpecific(fds) <- TRUE
    ans <- countSplitReadsPerChromosome("chrUn_gl000218", bamFile(fds)[1], 
            pairedEnd=TRUE, strandMode=strandSpecific(fds)[1], genome=NULL,
            scanBamParam=scanBamParam(fds))
    expect_equivalent(ans, GRanges())
    
})

test_that("test minAnchor", {
    fds <- createTestFraserSettings()
    features <- makeGRangesFromDataFrame(data.table(GeneID=1:3, Chr="chr19",
            Start=c(7592515, 7594599, 7594599),
            End=c(7592749, 7595171, 7595320), Strand="*"))
    
    out <- capture.output({ctnNS5 <- as.matrix(countNonSplicedReads(
            "sample3", features, fds, minAnchor=5, recount=TRUE)) })
    out <- capture.output({ctnNS25 <- as.matrix(countNonSplicedReads(
        "sample3", features, fds, minAnchor=25, recount=TRUE)) })
    
    expect_equivalent(c(7, 8, 0, 0, 7), ctnNS5[,1])
    expect_equivalent(c(5, 8, 0, 0, 7), ctnNS25[,1])
})

test_that("Test psi values", {
    attach(test_generate_count_example())

    expect_equal(as.vector(counts(test_rangeFDS, type="psi3")), test_rawCountsJ)
    expect_equal(test_p3rawOCounts,
        as.vector(counts(test_rangeFDS, type="psi3", side="other"))
    )

    expect_equal(as.vector(counts(test_rangeFDS, type="psi5")), test_rawCountsJ)
    expect_equal(test_p5rawOCounts,
        as.vector(counts(test_rangeFDS, type="psi5", side="other"))
    )

    #expect_equal(as.vector(counts(test_rangeFDS, type="theta")), test_rawCountsSS)
    expect_equal(test_pSrawOCounts,
        as.vector(counts(test_rangeFDS, type="theta", side="other"))
    )

    expect_equal(as.vector(assays(test_rangeFDS)[["psi3"]]),
        test_rawCountsJ / (test_rawCountsJ + test_p3rawOCounts)
    )

    expect_equal(as.vector(assays(test_rangeFDS)[["psi5"]]),
        test_rawCountsJ / (test_rawCountsJ + test_p5rawOCounts)
    )

    #expect_equal(as.vector(assays(test_rangeFDS)[["theta"]]),
    #    test_rawCountsSS / (test_rawCountsSS + test_pSrawOCounts)
    #)
})
