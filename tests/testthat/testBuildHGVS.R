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

	# expected <- "c.123_125del"
	# result <- builder$deletion(123,125)
	# expect_equal(expected, result)
	# print(result)

	# expected <- "c.123_125inv"
	# result <- builder$inversion(123,125)
	# expect_equal(expected, result)
	# print(result)
	
	# expected <- "c.123_125dup"
	# result <- builder$duplication(123,125)
	# expect_equal(expected, result)
	# print(result)
	
	# expected <- "c.123_124insACTTG"
	# result <- builder$insertion(123,"ACTTG")
	# expect_equal(expected, result)
	# print(result)

	# expected <- "c.123_125delinsACTTG"
	# result <- builder$delins(123,125,"ACTTG")
	# expect_equal(expected, result)
	# print(result)
	
})