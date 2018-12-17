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

context("Parsing HGVS strings")


test_that("common HGVS strings are parsed correctly", {

	common.hgvs.strings <- c(
		"g.1318G>T","g.3661_3706del","g.495_499inv","g.3661_3706dup",
		"g.7339_7340insTAGG","g.112_117delinsTG",
		"c.1318G>T","c.3661_3706del","c.495_499inv","c.3661_3706dup",
		"c.7339_7340insTAGG","c.112_117delinsTG",
		"p.Arg490Ser","p.R490S","p.Trp87Ter","p.Trp78*","p.W87*","p.Asp388del",
		"p.Asp388_Gln393del","p.Asp388_Gln393dup","p.Ala228_Val229insTrpPro",
		"p.Ala228_Val229insLys*","p.L7_H8delinsWQQFRTG","p.Arg98fs","p.Arg98Profs*23",
		"p.=","p.Arg12=","p.M15="
	)

	result <- parseHGVS(common.hgvs.strings)

	cat("\n")
	print(result)

	# expect_equal(str_length("a"), 1)
	# expect_equal(str_length("ab"), 2)
	# expect_equal(str_length("abc"), 3)
})

test_that("multi-mutant HGVS strings are parsed correctly", {

	multi.hgvs.strings <- c(
		"c.[76C>T;283G>C]","c.[76C>T];[283G>C]","c.[76C>T(;)283G>C]",
		"c.[76C>T;133_134insTAGG]"
	)

	result <- parseHGVS(multi.hgvs.strings)

	cat("\n")
	print(result)

})

test_that("multi-mutants are distinguished correctly from single mutants", {

	mixed.hgvs.strings <- c(
		"c.[76C>T;283G>C]",
		"c.495_499inv",#mix in a single mutant to test
		"c.[76C>T;133_134insTAGG]"
	)

	result <- parseHGVS(mixed.hgvs.strings)

	cat("\n")
	print(result)

})



test_that("conversion between AA codes works", {

	hgvs.strings <- c(
		"p.Arg490Ser","p.R490S","p.Trp87Ter","p.Trp78*","p.W87*"
	)

	result <- parseHGVS(hgvs.strings,aacode=1)
	cat("\n")
	print(result)

	result <- parseHGVS(hgvs.strings,aacode=3)
	cat("\n")
	print(result)

})


test_that("invalid strings get marked correctly", {

	hgvs.strings <- c(
		"p.Arg490Ser","x.R490S","pTrp87Ter","p.Trp7&8*","r.1318g>u"
	)

	result <- parseHGVS(hgvs.strings,aacode=1)
	cat("\n")
	print(result)

	result <- parseHGVS(hgvs.strings,aacode=3)
	cat("\n")
	print(result)

})


# test_that("exotic HGVS strings are parsed correctly", {

# 	exotic.hgvs.strings <- c(
# 		"g.1318G>T","g.3661_3706del","g.495_499inv","g.3661_3706dup",
# 		"g.7339_7340insTAGG","g.333_590con1844_2101","g.112_117delinsTG",
# 		"c.1318G>T","c.3661_3706del","c.495_499inv","c.3661_3706dup",
# 		"c.7339_7340insTAGG","c.333_590con1844_2101","c.112_117delinsTG",
# 		"c.7309+1160T[22]","c.31+59093TAA[9]","c.31+59093_31+59095[9]",
# 		"c.[76C>T;283G>C]","c.[76C>T];[283G>C]","c.[76C>T(;)283G>C]",
# 		"r.1318g>u","r.3661_3706del","r.3661_3706dup","r.7339_7340instagg",
# 		"r.456_457ins456+87_456+121","r.3334_3350inv","r.112_117delinsug",
# 		"r.[76c>u;283g>c]","r.[76c>u];[283g>c]","r.[76a>c,70_77del]",
# 		"r.spl","r.?",
# 		"p.(Arg490Ser)","p.Arg490Ser","p.R490S","p.Trp87Ter","p.Trp78*","p.W87*",
# 		"p.Asp388_Gln393del","p.Asp388_Gln393dup","p.Ala228_Val229insTrpPro",
# 		"p.Ala228_Val229insLys*","p.L7_H8delinsWQQFRTG","p.(Arg98fs)","p.Arg98Profs*23",
# 		"p.Met1ValextMet-12","p.Ter110GlnextTer17","p.Gln34[22]","p.Ser7_Ala9[6]",
# 		"p.[Trp78*;Arg490Ser]","p.[Trp78*];[Arg490Ser]","p.[Trp78*(;)Arg490Ser]","p.?",
# 		"m.1318G>T","m.3661_3706del","n.1318G>T","n.3661_3706del"
# 	)

# 	result <- parseHGVS(exotic.hgvs.strings)

# 	cat("\n")
# 	print(result)

# 	# expect_equal(str_length("a"), 1)
# 	# expect_equal(str_length("ab"), 2)
# 	# expect_equal(str_length("abc"), 3)
# })

