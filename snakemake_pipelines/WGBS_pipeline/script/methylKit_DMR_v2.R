## R --slave --args 100 treat1.txt,treat2.txt ctrl1.txt,ctrl2.txt rice CpG H19-04_vs_WT /home/yuming/data6/Methydata/project/rice/80-196880609/DMR/output < methylKit_DMR_v1.R &

library(methylKit)
Args <- commandArgs(TRUE)
bin_size <- as.integer(Args[1])  ### bin size 
min_cov <- as.integer(Args[2])  ### 4，depth for each cytosine
treat_files <- as.character(Args[3])   ### treat1.txt,treat2.txt
ctrl_files <- as.character(Args[4])    ### ctrl1.txt,ctrl2.txt
assembly <- as.character(Args[5])   ## rice, human, tair
context <- as.character(Args[6])   ###context :'CpG', 'CHG' or 'CHH'
diff <- as.integer(Args[7])      ### 70 50 10
output_obj <- as.character(Args[8])  ### blank_vs_H01-03-1-2.cg.chr1.4DMR.obj.RData
output_fil_obj <- as.character(Args[9])  ### blank_vs_H01-03-1-2.cg.chr1.4DMR.fil.obj.RData
output_Diff <- as.character(Args[10])  ### blank_vs_H01-03-1-2.cg.chr1.4DMR.200.Diff70p.RData

qval <- 0.01
if(context %in% "cg"){
	context <- "CpG"
}else if(context %in% "chg"){
	context <- "CHG"
}else if(context %in% "chh"){
	context <- "CHH"
}

treat_files <- as.list(strsplit(treat_files,"\\,")[[1]])
ctrl_files <- as.list(strsplit(ctrl_files,"\\,")[[1]])

file.list <- c(treat_files,ctrl_files)

myobj <- methRead(file.list,sample.id = as.list(c(paste0("T",c(1:length(treat_files))), paste0("C",c(1:length(ctrl_files))))),
assembly = assembly,treatment = c(rep(1, length(treat_files)), rep(0, length(ctrl_files))),context = context,resolution = "base",
header = T,mincov = 0)

filtered.myobj <- filterByCoverage(myobj,lo.count = min_cov,   ### filter < 4 (filter coverage <=3), not including 4, minimum cov is 4
lo.perc = NULL,hi.count = NULL,hi.perc = NULL)

result1 <- tryCatch({
###等待排错的语句
region <- tileMethylCounts(filtered.myobj, win.size = bin_size, step.size = bin_size)

meth <- unite(region, destrand = FALSE)

myDiff <- calculateDiffMeth(meth, adjust = "fdr")

# get hyper methylated bases
#myDiff25p.hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")
# get hypo methylated bases
#myDiff25p.hypo <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hypo")
#
#
# get all differentially methylated bases
getmyDiff <- getMethylDiff(myDiff, difference = diff, qvalue = qval)
#class(getmyDiff)=="methylDiff"
#[1] TRUE
### meth.diff = methyratio of treament 0 - methyratio of treament 1

#myDiff25p.hyper <- data.frame(myDiff25p.hyper, stringsAsFactors=F)
#myDiff25p.hypo <- data.frame(myDiff25p.hypo, stringsAsFactors=F)
#myDiff25p <- data.frame(myDiff25p, stringsAsFactors=F)
}, warning = function(w) {
###捕获警告（警告仅仅善意提醒，不会导致程序中断，
#属于非致命异常，通常以warning开头）
getmyDiff <<- NULL
#print('warning')
}, error   = function(e) { 
###捕获错误（错误是比较严重的故障，倘若不捕获并处理，
###则会通过编辑器抛出错误信息并中断程序运行，
#因而属于致命异常，是我们重点处理对象）
#if(class(filtered.myobj)=="methylRawList"){
getmyDiff <<- NULL
#print('error')
#getmyDiff <- filtered.myobj
#class(filtered.myobj)=="methylRawList"
#[1] TRUE
#}else{
#getmyDiff <- NULL
#class(a)=="NULL"
#[1] TRUE
#}

}, finally = {
###finally属于无论错误与否都会执行的必须语句，
#这一点与Python中的try/expect中的finally语句用法相同
#print('finally')
save(myobj, file = output_obj)
save(filtered.myobj, file = output_fil_obj)
#print(getmyDiff)
save(getmyDiff, file = output_Diff)

})
