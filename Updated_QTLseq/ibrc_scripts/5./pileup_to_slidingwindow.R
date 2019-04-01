##################################################
## データ読み込み（tab区切り、headerあり）
## 引数：
##      infilename      ファイル名
## 戻り値：
##      data.frame型
##################################################
load_data_with_header <- function(filename){
        ##  Read in tab-delimited file as a table
        ##  Use the first line as column headers
        data <- read.table(filename, head=TRUE, as.is=TRUE, quote="", comment.char="", sep="\t")
        return (data)
}

##################################################
## データ読み込み（tab区切り、headerなし）
## 引数：
##      infilename      ファイル名
## 戻り値：
##      data.frame型
##################################################
load_data_without_header <- function(filename){
        ##  Read in tab-delimited file as a table
        ##  Use the first line as column headers
        data <- read.table(filename, head=FALSE, as.is=TRUE, quote="", comment.char="", sep="\t")
        return (data)
}

##################################################
## 
##################################################
load_pileup <- function(filename){
	data <- load_data_without_header(filename)
	return (data)
}

##################################################
## 
##################################################
get_chromosome_name_size <- function(filename){
	data <- load_data_without_header(filename)
	ChrNum <<- nrow(data)
	for (ii in 1:ChrNum){
		ChrName <<- data[, 1]
		ChrSize <<- data[, 2]
	}
}


#	##################################################
#	## 出力ファイル名を決める
#	## コマンドラインの第3引数args[3]があるときは、それを採用し、
#	## ないときは、第2引数(=入力ファイル名)の後ろにdescript（".png"）を付けたものにする
#	## 引数：
#	##      args    コマンドラインオプション
#	## 戻り値：
#	##      出力用ファイル名
#	##################################################
#	set_outputname <- function(args, descript=".png"){
#			##  Apply user-defined output filename, if it was specified
#			if (length(args) > 2){
#					outfile <- args[3]
#			} else {
#					outfile <- paste(args[2], descript, sep="")
#			}
#			return(outfile)
#	}

##################################################
## 主処理
##################################################
make_data <- function(data, chrnum, chrname, chrsize, winwidth, winshift, mincount, outpath, outfile){

	outfile2 <- paste("mask", mincount, "_", outfile, sep="")

	outfile  <- paste(outpath, "/", outfile, sep="")
	outfile2 <- paste(outpath, "/", outfile2, sep="")

	if(file.exists(outfile)){
		file.remove(outfile)
	}
	if(file.exists(outfile2)){
		file.remove(outfile2)
	}

	myncol <- ncol(data)
	whichcol <- NULL
	if(AutoCalcNumericClmn){
		for (cc in 1:myncol){
			whichcol <- c(whichcol, is.numeric(data[1, cc]))
		}
	}
	else{
		whichcol <- rep(FALSE, myncol)
		for (cc in WhichClmnYouCalc){
			whichcol[cc] <- TRUE	
		}
	}
	print(whichcol)
	print(which(whichcol==TRUE))
	kk <- length(which(whichcol==TRUE))
	print(kk)

	for (nn in 1:chrnum){

		chrdata <- data[data[ClnmChrName] == chrname[nn],]
			#- 2014.06.05 kikuchi
		if(nrow(chrdata) < 1) next
			#- 
		steps <- ceiling(chrsize[nn] / winshift)
		for (ss in 1:steps){
			pos0 <- winshift * (ss - 1)
			pos1 <- pos0 - winwidth / 2
			pos2 <- pos0 + winwidth / 2
			if(pos1 < 0)		pos1 <- 0
			if(pos2 > chrsize[nn])	pos2 <- chrsize[nn]
			actwidth <- pos2 - pos1 + 1
			chrposdata <- chrdata[(chrdata[ClnmChrPos]>=pos1 & chrdata[ClnmChrPos]<=pos2), ]
			
			# # chrposdata_mean <- apply(chrposdata[5:9], 2, mean)
			# # chrposdata_max  <- apply(chrposdata[5:9], 2, max)
			# # chrposdata_min  <- apply(chrposdata[5:9], 2, min)
			
			# chrposdata_mean <- NULL
			# for (cc in 1:myncol){
			# 	if(whichcol[cc]){
			# 		chrposdata_mean <- c(chrposdata_mean, apply(chrposdata[cc], 2, mean))
			# 	}
			# 	else{
			# 		chrposdata_mean <- c(chrposdata_mean,"e")
			# 	}
			# }
			
			chrposdata_mean <- apply(chrposdata[whichcol], 2, mean)
			chrposdata_cnt  <- nrow(chrposdata)
			
			newdata <- c(chrname[nn], pos0, actwidth, chrposdata_mean, chrposdata_cnt)
			write(newdata, file = outfile, ncolumns = myncol + 10, append = TRUE, sep = OutputSep)
		
			if(mincount > 0){	
				if(chrposdata_cnt < mincount){
					chrposdata_null <- rep(NaN, kk)
					newdata <- c(chrname[nn], pos0, actwidth, chrposdata_null, chrposdata_cnt)
				}
				write(newdata, file = outfile2, ncolumns = myncol + 10, append = TRUE, sep = OutputSep)
			}
		}
 	}
	warnings()
}



##################################################
## グローバル変数に格納
## 引数:
##      なし
## 戻り値：
##      なし
##################################################
set_gloval_vars <- function(){
	# グローバル変数PercentProbs
	OutputSep <<- "\t"
	PercentProbs <<- c(90.0, 95.0, 99.0)
		
	# 永続付値
	ChrNum  <<- 0
	ChrName <<- NULL
	ChrSize <<- NULL

	ClnmChrName <<- 1
	ClnmChrPos  <<- 2
	
	# AutoCalcNumericClmn <<- TRUE
	AutoCalcNumericClmn <<- FALSE
	WhichClmnYouCalc <<- c(8:9, 19:20, 23:24, 30:39)
		
	# 1		chromosome
	# 2 	pos
	# 3 	reference base
	# 4		consensus base
	# 5		q1
	# 6 	q2
	# 7 	q3
	# 8 	depth
	# 9 	SNP-index
	# 10	bases
	# 11	qualities
	# 12	chromosome
	# 13 	pos
	# 14 	reference base
	# 15	consensus base
	# 16	q1
	# 17 	q2
	# 18 	q3
	# 19 	depth
	# 20 	SNP-index
	# 21	bases
	# 22	qualities
	# 23	min_depth
	# 24	delta(SNP-index)
	# 25	H90%
	# 26	L90%
	# 27	delta(SNP-index) < L90%
	# 28	L90% <= delta(SNP-index) <= H90%
	# 29	H90% < delta(SNP-index)
	# 30	H95%
	# 31	L95%
	# 32	delta(SNP-index) < L95%
	# 33	L95% <= delta(SNP-index) <= H95%
	# 34	H95% < delta(SNP-index)
	# 35	H99%
	# 36	L99%
	# 37	delta(SNP-index) < L99%
	# 38	L99% <= delta(SNP-index) <= H99%
	# 39	H99% < delta(SNP-index)
}



##################################################
## Main routine
##################################################

# ==================================================
# set default for global variants
# ==================================================
set_gloval_vars()

# ==================================================
# Capture command line arguments
# ==================================================
args <- (commandArgs(TRUE))

# ==================================================
# read data
# ==================================================
data <- load_pileup(args[1])

# ==================================================
# read fai -> chromosome size
# ==================================================
get_chromosome_name_size(args[2])
print(ChrNum)
print(ChrName)
print(ChrSize)

# ==================================================
# parameters for sliding window
# ==================================================
winwidth <- as.numeric(args[3])
winshift <- as.numeric(args[4])
mincount <- as.numeric(args[5])	# masked if the counts in a window is less than this value

# ==================================================
# decide output file name
# ==================================================
# outfile <- set_outputname(args, ".png")
outpath <- args[6]
outfile <- args[7]

# ==================================================
# run process
# ==================================================

make_data(data, ChrNum, ChrName, ChrSize, winwidth, winshift, mincount, outpath, outfile)

