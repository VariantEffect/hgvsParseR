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

options(stringsAsFactors=FALSE)
library(hgvsParseR)

context("Translating HGVS strings")

test_that("HGVS strings are translated correctly", {

	cds <- "ATGTTGTCACCACCCTGA"

	cbuilder <- new.hgvs.builder.c()
	hgvs <- cbuilder$substitution(4,"T","A")
	result <- translateHGVS(hgvs,cds)
	expect_equal("p.Leu2Met",result[["hgvsp"]])

})
