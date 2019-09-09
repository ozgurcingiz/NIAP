#'@export
literatureOverlapAnalysis <- function(ppi,anno,mypath="") # evaluateNet <- function(net,ppi,anno)
{
  #mypath=paste(getwd(),"/",mypath,sep="");
  mypath=choose.dir(default = "", caption = "Select folder")
  #print(mypath)
  
  #library("GAnet");
  
  
  #   file_list <- list.files(path =mypath );
  
  file_list=list.files(path = mypath, pattern = NULL, all.files = TRUE,full.names = TRUE, recursive = TRUE,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE);
  fileCount=length(file_list);
  #print(fileCount)
  
  wholeInfo=matrix(data=NA,nrow=fileCount,ncol=8);
  
  counter=1;
  
  for (file in file_list)    {
    
    
    
    #temp_dataset <-read.table(file, header=TRUE, quote="\""); # temp_dataset <-read.table(file, header=TRUE, sep="\t")
    
    ################# for rda files###############
    
    #file_name=basename(file);
    
    
    #full_name=paste(mypath,"/",file_name,sep="")
    
    #jName=substr(file_name,1,nchar(file_name)-4)
    
    #load(file=full_name) 
    
    #temp_dataset=get(jName)
    
    ####################### for text files ##########
    
    temp_dataset <-read.table(file, header=TRUE, quote="\""); # temp_dataset <-read.table(file, header=TRUE, sep="\t")
    
    
    if(nrow(temp_dataset)<2) {
      next ;
    }
    
    
    net=as.matrix(temp_dataset[1:nrow(temp_dataset),1:2]); # hata verirse sadece bu folder'da binary iliskili txt dosyalari olmali
    net=ganet.UniqNetSimp(net);
    
    
    
    
    temp=as.vector((anno[,2])); # rna'da 1 olacak, gexp'de 2 olacak, deneylerde de 1
    temp=unique(temp);
    
    tempL=length(temp);
    nUniverse <- tempL*(tempL- 1)/2;
    nUniverse <- round(nUniverse/2);
    
    
    
    Validated <-ganet.ComLinks(netlist=as.matrix(net), netdata=as.matrix(ppi));
    
    nFocusedSet <- length( intersect(which(ppi[,1] %in% temp),which(ppi[,2]%in% temp)) );
    
    nOverlap=nrow(Validated);nPredicted=nrow(net); nFocusedSet=nFocusedSet;nUniverse=nUniverse;
    
    
    
    fe1 <- nrow(Validated) 
    tp=fe1;# TP
    fe2 <- nPredicted - nOverlap # FP
    fe3 <- nFocusedSet - nOverlap # FN
    fe4 <- nUniverse - nPredicted - nFocusedSet + nOverlap 
    
    recall= fe1/(fe1+fe3);
    r=as.numeric(recall);
    prec= fe1/(fe1+fe2);
    p=as.numeric(prec);
    
    fmeasure= (2*r*p)/ (r+p);
    
   # print(nrow(Validated))
    #print(nrow(net))
    #print(nFocusedSet)
    #print(nUniverse)
    
    
    
     fin <- ganet.FEtest(nOverlap=nrow(Validated),nPredicted=nrow(net), nFocusedSet=nFocusedSet,nUniverse=nUniverse)
    
    
    
    res <- new.env();
    
    assign("precision", fe1/(fe1+fe2), envir = res); #assign("precision", fin$precision, envir = res);
    assign("tp", fe1, envir = res);
    assign("fp", fe2, envir = res);
    assign("fn", fe3, envir = res);
    #assign("pval", fin$stats$p.value, envir = res);
    assign("recall", r, envir = res);
    assign("fmeasure", fmeasure, envir = res);
    
    wholeInfo[counter,1]=basename(file);
    wholeInfo[counter,2]=fe1/(fe1+fe2);  #wholeInfo[counter,2]=fin$precision;
    wholeInfo[counter,3]=fe1;
    wholeInfo[counter,4]=fe2;
    wholeInfo[counter,5]=fe3;
    wholeInfo[counter,6]=fin$stats$p.value;
    wholeInfo[counter,7]=r;
    wholeInfo[counter,8]=fmeasure;
    
    counter=counter +1;
  }
  
  colnames(wholeInfo)=c("name","precision","tp","fp","fn","pval","recall","fmeasure");
  
  #rm(list = setdiff(ls(), lsf.str()))
  
  wholeInfo
  
  
  
  
  # res
}