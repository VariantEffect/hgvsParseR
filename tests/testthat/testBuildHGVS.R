# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of hgvsParseR.
#
# hgvsParseR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hgvsParseR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hgvsParseR.  If not, see <https://www.gnu.org/licenses/>.

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

	expected <- "g.[123A>C;231G>A]"
	result <- with(builder,cis(substitution(123,"A","C"),substitution(231,"G","A")))
	expect_equal(expected,result)
	print(result)
	
	expected <- "g.[123A>C];[231G>A];[311_320del]"
	result <- with(builder,trans(substitution(123,"A","C"),substitution(231,"G","A"),deletion(311,320)))
	expect_equal(expected,result)
	print(result)

	expected <- "g.[123A>C(;)231G>A]"
	result <- with(builder,nophase(substitution(123,"A","C"),substitution(231,"G","A")))
	expect_equal(expected,result)
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
	
	expected <- "c.[123A>C;231G>A]"
	result <- with(builder,cis(substitution(123,"A","C"),substitution(231,"G","A")))
	expect_equal(expected,result)
	print(result)
	
	expected <- "c.[123A>C];[231G>A];[311_320del]"
	result <- with(builder,trans(substitution(123,"A","C"),substitution(231,"G","A"),deletion(311,320)))
	expect_equal(expected,result)
	print(result)

	expected <- "c.[123A>C(;)231G>A]"
	result <- with(builder,nophase(substitution(123,"A","C"),substitution(231,"G","A")))
	expect_equal(expected,result)
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

	expected <- "p.[R123K;S125_L152del]"
	result <- with(builder,cis(substitution(123,"R","K"),deletion(125,"S",152,"L")))
	expect_equal(expected,result)
	print(result)

	expected <- "p.R123="
	result <- builder$synonymous(123,"Arg")
	expect_equal(expected,result)
	print(result)
	
	expected <- "p.="
	result <- builder$synonymous()
	expect_equal(expected,result)
	print(result)

})
