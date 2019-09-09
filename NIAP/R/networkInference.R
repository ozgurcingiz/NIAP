#'@export
networkInference <- function(dataset,anno,gnimethod)
{
  #library("GAnet"); 
  #library("minet");
  #library("c3net");
 
  if(gnimethod=="c3net")
  {
    net <- c3net(dataset, network=FALSE);
  }
  if(gnimethod=="aracne")
  {
    mim <- build.mim(t(dataset), estimator="spearman");
    #mim  <- knnmi.all(dataset);
    net=aracne(mim);
  }
  if(gnimethod=="clr")
  {
    #mim <- build.mim(t(dataset), estimator="spearman");
    mim  <- knnmi.all(dataset);
    net=clr(mim);
  }
  if(gnimethod=="mrnet")
  {
    #mim <- build.mim(t(dataset), estimator="spearman");
    mim  <- knnmi.all(dataset);
    net=mrnet(mim);
  }
  if(gnimethod=="wgcna")
  {
    #thrs=pickSoftThreshold(t(dataset), dataIsExpr = TRUE, RsquaredCut = 0.9,  powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),  removeFirst = FALSE, nBreaks = 10, blockSize = NULL,  corFnc = cor, corOptions = list(use = 'p'),  networkType = "unsigned", moreNetworkConcepts = FALSE, verbose = 0, indent = 0);
    
    net=adjacency(datExpr=t(dataset),type = "unsigned",power=7) ; # gexp icin power 7, r.sq degeri 0.9'dan buyuk ilk deger
    net=adjacency(datExpr=t(dataset),type = "unsigned",power=8) ; # rna icin power 8, r.sq degeri 0.9'dan buyuk ilk deger
  }
  
  indices=which(net>0,arr.ind=TRUE); #indices=which(net>0,arr.ind=TRUE); indices=which(mim>0,arr.ind=TRUE)
  
  
  indices=ganet.UniqNetSimp(indices);
  class(indices)="numeric";
  #rownamesAll=as.matrix(rownames(pn_priors));
  #firstColumn=as.matrix(rownamesAll[indices[,1],1]);
  firstColumn=as.matrix(anno[indices[,1],2]); # rna'da 1,gene exp'de 2 olacak rna-seq 1 olacak
  #secondColumn=as.matrix(rownamesAll[indices[,2],1]);
  secondColumn=as.matrix(anno[indices[,2],2]);  #rna'da 1,  gene exp'de 2 olacak rna-seq 1 olacak
  interactions=cbind(firstColumn,secondColumn);
  interactions[,1]=toupper(interactions[,1]);
  interactions[,2]=toupper(interactions[,2]);
  ind=which(as.matrix(interactions[,1])!=as.matrix(interactions[,2]));
  interactions=interactions[ind,];
  interactions=ganet.UniqNetSimp(interactions);
  dosya=paste(gnimethod,".txt");
  #write.table(interactions,dosya); 
  
  # rm(list = setdiff(ls(), lsf.str()))
  
  interactions
  
  
  
}
