# File for running bfast0n on a datacube provided by rpy2
# Dependencies:
# devtools::install_github("bfast2/strucchange@devel")
# devtools::install_github("bfast2/bfast@dev")

library(pbapply)

suppressWarnings(source("bfast0n.r"))
suppressWarnings(source("bfastmon.r"))

source("fast-bfast.r")
EnableFastBfast()

# rpy2 will put an object called "DataCube" into the environment

#ExampleTS = c(897, NA, 1052, 785, 749, 712, 828, 755, 731, 623, 578, 459, 594, 533, 485, 411, 425, 413, 369, 353, 280, 286, 319, 262, 290, 267, 251, 268, 254, 266, 276, 275, 278, 301, 288, 349, 406, 414, 452, 477, 506, 596, 588, 540, 541, 579, NA, 520, 708, 564, 564, 455, 558, 504, 437, 431, 383, 386, 487, 463, 593, 458, 500, 464, 483, 424, 397, 303, 302, 277, 266, 256, 242, 271, 247, 248, 261, 275, 276, 284, 298, 299, 315, 300, 312, 329, 387, NA, 416, 368, 571, 598, NA, NA, 591, 543, 482, 601, 493, 477, 471, 476, 607, 501, 423, 467, 404, 420, 377, 329, 309, 298, 281, 267, 269, 252, 241, 250, 250, 241, 264, 243, 240, 229, 241, 241, 235, 240, 251, 269, 298, 321, 324, 399, 405, 552, 597, 554, 554, 779, 531, 409, 413, NA, 453, 466, 369, 367, 396, 388, 400, 377, 326, 324, 304, 297, 272, 290, 271, 272, 269, 257, 251, 247, 244, 249, 239, 251, 245, 231, 238, 273, 244, 257, 257, 304, 389, 460, 374, 535, 586, 588, 557, 754, 762, 587, 751, 547, 629, 690, 572, 501, 466, 432, 422, 395, 345, 341, 316, 350, 325, 315, 323, 301, 298, 286, 282, 274, 252, 264, 259, 254, 247, 254, 247, 248, 243, 251, 252, 283, 283, 300, 283, 308, 254, 293, 275, 328, 388, 459, 509, 555, 567, 549, NA, 557, 481, 461, NA, 603, 790, 613, 792, 594, 550, 536, 535, 516, 489, 432, 408, 371, 368, 363, 323, 334, 310, 311, 302, 326, 297, 319, 292, 312, 295, 315, 290, 304, NA, 308, 279, 330, 331, 381, 336, 402, 402, 337, 422, 348, 373, 329, 328, 316, 330, 305, 320, 302, 333, 297, 340, 303, 300, 284, 272, 265, 263, 262, 267, 241, 244, 239, 240, 226, 235, 234, 241, 262, 247, 262, 254, 245, 239, 234, 259, 235, 244, 297, 344, 402, 369, 336, 317, 371, 421, 455, 417, 393, 372, 499, 448, 646, 523, 403, 350, 341, 268, 328, 303, 321, 287, 274, 261, 248, 259, 252, 251, 247, 250, 242, 249, 237, 180, 239, 231, 231, 248, 223, 227, 256, 283, 333, 346, 425, 421, NA, 419, NA, NA, 807, 806, 963, 840, 722, NA, 765, 749, 786, 687, 653, 737, 484, 448, 436, 416, 386, 375, 368, 351, 340, 340, 305, 278, 276, 287, 270, 271, 288, 276, 260, 275, 299, 264, 250, 279, 274, 264, 258, 243, 246, 318, 388, NA, 328, 328, 305, 305, 364, 374, 304, NA, 577, 427, 480, 527, 562, 536, 608, 520, 498, 472, 441, 450, 399, 373, 362, 342, 388, 278, 289, 265, 268, 259, 262, 247, 263, 245, 266, 246, 255, 255, 254, 270, 254, 248, 265, 244, 260, 334, NA)
#DataCube = array(ExampleTS, c(1, 1, length(ExampleTS)))
#DataCube = abind::abind(DataCube, DataCube, DataCube, DataCube, along=1)
#DataCube = abind::abind(DataCube, DataCube, DataCube, along=2)

# Function to apply the BFAST0NBreaks function onto a datacube
# The result doesn't make sense for R (run aperm() on it to make it have sense)
# but it does make sense in Python.
# Using pbapply, so you get a progress bar by default. Disable it with pboptions(type="none")
runBFAST = function(DataCube, BFFunction=BFAST0NBreaks, ..., MissingValue=NA)
{
    # Replace all the missing values with NA
    if (!is.na(MissingValue))
        DataCube[DataCube == MissingValue] = NA
    
    BFAST_result = pbapply(DataCube, c(1,2), BFFunction, ...)
    
    # Replace all the NAs with the missing value
    if (!is.na(MissingValue))
        BFAST_result[is.na(BFAST_result)] = MissingValue
    
    return(BFAST_result)
}

# For max sensitivity BFAST0N, use runBFAST(DataCube, BFFunction=BFAST0NBreaks, breaks="BIC", scsig=NA)

# To run BFAST Monitor with trend (conservative) model:
# runBFAST(DataCube, BFFunction=BFMBreaks)
# With the seasonal (ultraconservative) model, devtools::install_github("bfast2/bfast@dev") and then:
# runBFAST(DataCube, BFFunction=BFMBreaks, formula=response~trend+season, sbins=4)
