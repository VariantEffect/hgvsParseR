library(hgvsParseR)

context("Building HGVS strings")



test_that("genomic strings are built correctly", {
	builder <- new.hgvs.builder.g()

	expected <- "g.123A>C"
	result <- builder$substitution(123,"A","C")
	expect_equal(expected, result)
	print(result)

	expected <- "g.123_125del"
	result <- builder$deletion(123,125)
	expect_equal(expected, result)
	print(result)

	expected <- "g.123_125inv"
	result <- builder$inversion(123,125)
	expect_equal(expected, result)
	print(result)
	
	expected <- "g.123_125dup"
	result <- builder$duplication(123,125)
	expect_equal(expected, result)
	print(result)
	
	expected <- "g.123_124insACTTG"
	result <- builder$insertion(123,"ACTTG")
	expect_equal(expected, result)
	print(result)

	expected <- "g.123_125delinsACTTG"
	result <- builder$delins(123,125,"ACTTG")
	expect_equal(expected, result)
	print(result)
	
})


test_that("coding strings are built correctly", {
	builder <- new.hgvs.builder.c()

	expected <- "c.123A>C"
	result <- builder$substitution(123,"A","C")
	expect_equal(expected, result)
	print(result)

	expected <- "c.123+2A>C"
	result <- builder$substitution(123,"A","C",2)
	expect_equal(expected, result)
	print(result)

	expected <- "c.1-2A>C"
	result <- builder$substitution(1,"A","C",-2)
	expect_equal(expected, result)
	print(result)

	expected <- "c.123_125del"
	result <- builder$deletion(123,125)
	expect_equal(expected, result)
	print(result)

	expected <- "c.123_125inv"
	result <- builder$inversion(123,125)
	expect_equal(expected, result)
	print(result)
	
	expected <- "c.123_125dup"
	result <- builder$duplication(123,125)
	expect_equal(expected, result)
	print(result)
	
	expected <- "c.123_124insACTTG"
	result <- builder$insertion(123,"ACTTG")
	expect_equal(expected, result)
	print(result)

	expected <- "c.123+2_123+3insACTTG"
	result <- builder$insertion(123,"ACTTG",2)
	expect_equal(expected, result)
	print(result)

	expected <- "c.123_125delinsACTTG"
	result <- builder$delins(123,125,"ACTTG")
	expect_equal(expected, result)
	print(result)
	
})



test_that("protein strings are built correctly", {
	builder <- new.hgvs.builder.p(aacode=1)

	expected <- "p.R123K"
	result <- builder$substitution(123,"Arg","Lys")
	expect_equal(expected, result)
	print(result)
	
	expected <- "p.R123_L152del"
	result <- builder$deletion(123,"Arg",152,"Leu")
	expect_equal(expected, result)
	print(result)

	expected <- "p.R123_L152dup"
	result <- builder$duplication(123,"Arg",152,"Leu")
	expect_equal(expected, result)
	print(result)

	expected <- "p.R123_M124insLMI"
	result <- builder$insertion(123,"Arg","Met",c("Leu","Met","Ile"))
	expect_equal(expected, result)
	print(result)

	expected <- "p.R123_L152delinsKWS"
	result <- builder$delins(123,"Arg",152,"Leu",c("Lys","Trp","Ser"))
	expect_equal(expected, result)
	print(result)

	expected <- "p.R123fs"
	result <- builder$frameshift(123,"Arg")
	expect_equal(expected, result)
	print(result)

	expected <- "p.R123Kfs*201"
	result <- builder$frameshift(123,"Arg",variantAA="Lys",newStop=201)
	expect_equal(expected, result)
	print(result)
})
