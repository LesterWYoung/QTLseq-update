##################################################
## データ読み込み（tab区切り、headerあり）； 戻り値：data.frame型
##################################################
load_data_with_header <- function(filename){
	##  Read in tab-delimited file as a table
	##  Use the first line as column headers
 	data <- read.table(filename, head=TRUE, as.is=TRUE, quote="", comment.char="", sep="\t")
	return (data)
}

##################################################
## データ読み込み（tab区切り、headerなし）；　戻り値：data.frame型
##################################################
load_data_without_header <- function(filename){
	##  Read in tab-delimited file as a table
	##  Not use the first line as column headers
	data <- read.table(filename, head=FALSE, as.is=TRUE, quote="", comment.char="", sep="\t")
	return (data)
}

##################################################
# Prototype of how to read INI files to process olfactometer data
# efg, 13 June 2007
# Thanks to Gabor Grothendieck for helpful suggestions in the R-Help
# mailing list on how to parse the INI file.
# 
# (https://stat.ethz.ch/pipermail/r-help/2007-June/134115.html)
##################################################
Parse.INI <- function(INI.filename){
	connection <- file(INI.filename)
	Lines  <- readLines(connection)
	close(connection)
	
	Lines <- chartr("[]", "==", Lines)  # change section headers

	connection <- textConnection(Lines)
	d <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE)
	close(connection)

	L <- d$V1 == ""                    # location of section breaks
	d <- subset(transform(d, V3 = V2[which(L)[cumsum(L)]])[1:3], V1 != "")

	ToParse  <- paste("INI.list$", d$V3, "$",  d$V1, " <- '", d$V2, "'", sep="")

	INI.list <- list()
	eval(parse(text=ToParse))

	return(INI.list)
}

##################################################
# configファイルの型変換
##################################################
set_config <- function(myconfig){

	cfg <- myconfig
	
	cfg$Common$chrcolumn_file1 <- as.numeric(cfg$Common$chrcolumn_file1)
	cfg$Common$chrcolumn_file2 <- as.numeric(cfg$Common$chrcolumn_file2)
	
	cfg$Common$mfcol_nr <- as.numeric(cfg$Common$mfcol_nr)

	cfg$Common$xmaxauto <- as.logical(cfg$Common$xmaxauto)
	cfg$Common$xminauto <- as.logical(cfg$Common$xminauto)
	cfg$Common$ymaxauto <- as.logical(cfg$Common$ymaxauto)
	cfg$Common$yminauto <- as.logical(cfg$Common$yminauto)

	cfg$Common$xminfixed <- as.numeric(cfg$Common$xminfixed)
	cfg$Common$xmaxfixed <- as.numeric(cfg$Common$xmaxfixed)
	cfg$Common$yminfixed <- as.numeric(cfg$Common$yminfixed)
	cfg$Common$ymaxfixed <- as.numeric(cfg$Common$ymaxfixed)
	
	cfg$Common$showablineh1 <- as.logical(cfg$Common$showablineh1)
	cfg$Common$ablineh1     <- as.numeric(cfg$Common$ablineh1)
	cfg$Common$showablineh2 <- as.logical(cfg$Common$showablineh2)
	cfg$Common$ablineh2     <- as.numeric(cfg$Common$ablineh2)
	cfg$Common$showablinev1 <- as.logical(cfg$Common$showablinev1)
	cfg$Common$ablinev1     <- as.numeric(cfg$Common$ablinev1)
	cfg$Common$showablinev2 <- as.logical(cfg$Common$showablinev2)
	cfg$Common$ablinev2     <- as.numeric(cfg$Common$ablinev2)
	
	for(id in 1:8){
		if(id == 1)     {mycfg <- cfg$Data1}
		else if(id == 2){mycfg <- cfg$Data2}
		else if(id == 3){mycfg <- cfg$Data3}
		else if(id == 4){mycfg <- cfg$Data4}
		else if(id == 5){mycfg <- cfg$Data5}
		else if(id == 6){mycfg <- cfg$Data6}
		else if(id == 7){mycfg <- cfg$Data7}
		else if(id == 8){mycfg <- cfg$Data8}
    
		mycfg$enabled <- as.logical  (mycfg$enabled)
		mycfg$file    <- as.numeric  (mycfg$file   )
		mycfg$xcolumn <- as.numeric  (mycfg$xcolumn)
		mycfg$ycolumn <- as.numeric  (mycfg$ycolumn)
		mycfg$type    <- as.character(mycfg$type   )
		mycfg$color   <- as.character(mycfg$color  )
		mycfg$pch     <- as.numeric  (mycfg$pch)
		mycfg$ps      <- as.numeric  (mycfg$ps)
		mycfg$lty     <- as.numeric  (mycfg$lty)
		mycfg$lwd     <- as.numeric  (mycfg$lwd)
		
		if(id == 1)     {cfg$Data1 <- mycfg}
		else if(id == 2){cfg$Data2 <- mycfg}
		else if(id == 3){cfg$Data3 <- mycfg}
		else if(id == 4){cfg$Data4 <- mycfg}
		else if(id == 5){cfg$Data5 <- mycfg}
		else if(id == 6){cfg$Data6 <- mycfg}
		else if(id == 7){cfg$Data7 <- mycfg}
		else if(id == 8){cfg$Data8 <- mycfg}
	}
	
	return(cfg)
	
}



##################################################
## ここでは使用いない
##################################################
## ##################################################
## ## 1列目を行名にして読みだす
## ##################################################
## load_configs <- function(filename){
##  	data <- read.table(filename, head=TRUE, as.is=TRUE, quote="", comment.char="", sep="\t", row.names=1)
## 	return (data)
## }
## 
## ##################################################
## ## 出力ファイル名を決める
## ## コマンドラインの第3引数args[3]があるときは、それを採用し、
## ## ないときは、第2引数(=入力ファイル名)の後ろにdescript（".png"）を付けたものにする
## ## 引数：
## ##	args	コマンドラインオプション
## ## 戻り値：
## ## 	出力用ファイル名
## ##################################################
## set_outputname <- function(args, descript=".png"){
## 	##  Apply user-defined output filename, if it was specified
## 	if (length(args) > 2){
## 		outfile <- args[3]
## 	} else {
## 		outfile <- paste(args[2], descript, sep="")
## 	}
## 	return(outfile)
## }


##################################################
## abline用データ組を抽出
##################################################
set_abline <- function(h_or_v, myconfig){
	
	if(h_or_v == "h"){
		showabline1 <- as.logical(myconfig$Common$showablineh1)
		abline1     <- as.numeric(myconfig$Common$ablineh1)
		showabline2 <- as.logical(myconfig$Common$showablineh2)
		abline2     <- as.numeric(myconfig$Common$ablineh2)
	}
	else{
		showabline1 <- as.logical(myconfig$Common$showablinev1)
		abline1     <- as.numeric(myconfig$Common$ablinev1)
		showabline2 <- as.logical(myconfig$Common$showablinev2)
		abline2     <- as.numeric(myconfig$Common$ablinev2)		
	}
	
	return (c(showabline1, abline1, showabline2, abline2))
}

##################################################
## show abline
##################################################
show_ablines <- function(myablineh, myablinev){
	if(myablineh[1]){
		abline(h=myablineh[2], lty=3, lwd=3)
	}
	if(myablineh[3]){
		abline(h=myablineh[4], lty=3, lwd=3)
	}

	if(myablinev[1]){
		myablinevMb <- myablinev[2] / OneMega
		abline(v=myablinevMb, lty=3, lwd=3)
	}
	if(myablinev[3]){
		myablinevMb <- myablinev[4] / OneMega
		abline(v=myablinevMb, lty=3, lwd=3)
	}
}

##################################################
## ラベル用データを抽出
##################################################
set_label <- function(x_or_y, myconfig){
	if(x_or_y == "x")	value <- myconfig$Common$xlab	
	else				value <- myconfig$Common$ylab		
	return (value)
}
	

##################################################
## Data1~8に切り分ける
##################################################
select_1_8 <- function(myconfig, id){
	if(id == 1){mycfg <- myconfig$Data1}
	else if(id == 2){mycfg <- myconfig$Data2}
	else if(id == 3){mycfg <- myconfig$Data3}
	else if(id == 4){mycfg <- myconfig$Data4}
	else if(id == 5){mycfg <- myconfig$Data5}
	else if(id == 6){mycfg <- myconfig$Data6}
	else if(id == 7){mycfg <- myconfig$Data7}
	else if(id == 8){mycfg <- myconfig$Data8}
	
	return (mycfg)
}

##################################################
## 軸の範囲を決めるため、最大最小値を探す
##################################################	
search_maxminauto <- function(d1, d2, myconfig, x_or_y){
	mycfg <- NULL
	foundmax <- -Inf
	foundmin <- Inf
	for(id in 1:8){
		mycfg <- select_1_8(myconfig, id)
		if(x_or_y == "x"){ xycol <- mycfg$xcol }
		if(x_or_y == "y"){ xycol <- mycfg$ycol }
		if(mycfg$enabled){
			if(mycfg$file == 1){
				dmax <- max(d1[xycol], na.rm = TRUE)
				dmin <- min(d1[xycol], na.rm = TRUE)
			}
			else{
				dmax <- max(d2[xycol], na.rm = TRUE)
				dmin <- min(d2[xycol], na.rm = TRUE)
			}
			if(foundmax < dmax){ foundmax <- dmax }
			if(foundmin > dmin){ foundmin <- dmin }
		}
	}
	if(foundmax > 0){ foundmax <- (foundmax * 1.1) }
	else            { foundmax <- (foundmax * 0.9) }

	if(foundmin > 0){ foundmin <- (foundmin * 0.9) }
	else            { foundmin <- (foundmin * 1.1) }

	foundmaxmin <- list(foundmax, foundmin)	
	return(foundmaxmin)
}

##################################################
## 軸の範囲を決める 
##################################################
set_axis_lim <- function(d1, d2, myconfig, x_or_y){
	if(x_or_y == "x"){
		minauto  <- myconfig$Common$xminauto
		maxauto  <- myconfig$Common$xmaxauto
		minfixed <- myconfig$Common$xminfixed
		maxfixed <- myconfig$Common$xmaxfixed
	}
	if(x_or_y == "y"){
		minauto  <- myconfig$Common$yminauto
		maxauto  <- myconfig$Common$ymaxauto
		minfixed <- myconfig$Common$yminfixed
		maxfixed <- myconfig$Common$ymaxfixed
	}
	if(maxauto || minauto){
		foundmaxmin <- search_maxminauto(d1, d2, myconfig, x_or_y)
		foundmax <- as.numeric(foundmaxmin[1])
		foundmin <- as.numeric(foundmaxmin[2])
		print(foundmax)
		print(foundmin)
	}
	if(! minauto){
		foundmin <- minfixed
	}
	if(! maxauto){
		foundmax <- maxfixed
	}
	mylim <- c(foundmin, foundmax)
	
	return (mylim)
}

##################################################
## mfcol_ncを決める 
##################################################
set_mfcol_nc <- function(chr_length, mfcol_nr){

	chr_quotient <- chr_length %/% mfcol_nr
	chr_reminder <- chr_length %% mfcol_nr
	if(chr_reminder == 0){
		mfcol_nc <- chr_quotient
	}else{
		mfcol_nc <- chr_quotient + 1
	}
	return (mfcol_nc)
}

##################################################
## 主処理
##################################################
make_graph <- function(d1, d2, myconfig, outfile){
	
	# --------------------
	# 染色体名の列番号
	# --------------------
	ccol1 <- myconfig$Common$chrcolumn_file1 
	ccol2 <- myconfig$Common$chrcolumn_file2 

	# --------------------
	# xy軸の範囲
	# --------------------
	myxlim <- set_axis_lim(d1, d2, myconfig, "x")
	myxlim <- myxlim / OneMega
	myylim <- set_axis_lim(d1, d2, myconfig, "y")

	# --------------------
	# 軸のラベル
	# --------------------
	myxlab <- set_label("x", myconfig)
	myylab <- set_label("y", myconfig)

	# --------------------
	# abline(h/v)
	# --------------------
	myablinehs <- set_abline("h", myconfig)
	myablinevs <- set_abline("v", myconfig)

	# --------------------
	# get unique chromosome name
	# --------------------
	# chr_names <- unique(d2[, ccol2])
	chr_names  <- unique(c(d1[, ccol1], d2[, ccol2]))
	chr_length <- length(chr_names)
	
	# --------------------
	# set mfcol_nr, mfcol_nc
	# --------------------
	mymfcol_nr <- myconfig$Common$mfcol_nr
	mymfcol_nc <- set_mfcol_nc(chr_length, mymfcol_nr)

	png(outfile, width=800*mymfcol_nc, height=400*mymfcol_nr)
	par(mfcol=c(mymfcol_nr, mymfcol_nc))
	par(ps=30, lwd=4, bty="l", tck=-0.02)
	par(xaxs="i", yaxs="i")
	par(mar=c(6, 7.0, 4, 8))
	
	for(mychr in chr_names){
		
		search <- mychr
		title  <- mychr
		d1_chr <- d1[d1[[ccol1]]==search,]
		d2_chr <- d2[d2[[ccol2]]==search,]
		
		plot1st <- FALSE
		for(id in 1:8){
			mycfg <- select_1_8(myconfig, id)
			if(mycfg$enabled){
				
				if(mycfg$file == 1){ dd_chr <- d1_chr }
				else			   { dd_chr <- d2_chr }
				xcol    <- mycfg$xcolumn
				ycol    <- mycfg$ycolumn
				mycolor <- mycfg$color
				mytype  <- mycfg$type
				mypch   <- mycfg$pch
				myps    <- mycfg$ps
				mylty   <- mycfg$lty
				mylwd   <- mycfg$lwd
				if(plot1st == FALSE){
				 	plot(dd_chr[[xcol]]/OneMega, dd_chr[[ycol]], col=mycolor, type=mytype, xlim=myxlim, ylim=myylim, xlab=myxlab, ylab=myylab, main=title, cex.axis=1, cex.lab=1, mgp=c(4,1.5,0), pch=mypch, ps=myps, lty=mylty, lwd=mylwd)
					show_ablines(myablinehs, myablinevs)
					plot1st <- TRUE
				}
				else{
					par(new=T)
					plot(dd_chr[[xcol]]/OneMega, dd_chr[[ycol]], col=mycolor, type=mytype, xlim=myxlim, ylim=myylim, xlab="", ylab="", xaxt="n", yaxt="n", pch=mypch, ps=myps, lty=mylty, lwd=mylwd)
				}
			}
		}
	}
	dev.off()
	warnings()
}


##################################################
## グローバル変数に格納
##################################################
set_gloval_vars <- function(){
	# グローバル変数
	OutputSep <<- "\t"
	OneMega <<- 1000000
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
# load configulation file
# ==================================================
myconfig <- Parse.INI(args[1])
myconfig <- set_config(myconfig)
summary(myconfig)

# ==================================================
# read data
# ==================================================
data <- load_data_without_header(args[2])

# ==================================================
# read database
# ==================================================
datasw <- load_data_without_header(args[3])


# ==================================================
# decide output file name
# ==================================================
# outfile <- set_outputname(args, ".png")
outfile <- args[4]

# ==================================================
# run process
# ==================================================

make_graph(data, datasw, myconfig, outfile)

