#'@export
topologicalBiologicalEvaluation <- function (mypath="",omim)
{
  #library(ProNet);
  #library(igraph);
  #library(gProfileR);
  
  #mypath=paste(getwd(),"/",mypath,sep="");
  
  mypath=choose.dir(default = "", caption = "Select folder")
  
  res <- new.env();
  file_list=list.files(path = mypath, pattern = NULL, all.files = TRUE,full.names = TRUE, recursive = TRUE,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE);
  fileCount=length(file_list);
  fault=FALSE;
  ################################## BIOMART OMIM GENE LIST #####################3
  omim=as.vector(as.matrix(omim));
  
  
  # ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org");
  ###############################
  #   ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl");
  #   gtruthOmim=getBM(attributes=c("goslim_goa_accession","goslim_goa_description"),filters="hgnc_symbol",values=omim, mart=ensembl);
  #   tempIndice=which(as.matrix(gtruthOmim[,1])!="" & as.matrix(gtruthOmim[,2])!=""); # 86 TANE GO SLIM CIKARIYOR
  #   gtruthOmim=gtruthOmim[tempIndice,];
  #   allGoTermsNo=nrow(gtruthOmim);
  ####################################################### gProfileR based GSEA ###########################
  gtruthOmim=gprofiler(omim, organism = "hsapiens",correction_method="bonferroni");
  #gtruthOmim=gprofiler(omim, organism = "hsapiens",correction_method="bonferroni");
  gtruthOmim= as.matrix(gtruthOmim[,9]);
  allGoTermsNo=nrow(gtruthOmim);
  
  ######################################################################################
  #wholeInfo=matrix(data=NA,nrow=50000,ncol=11);
  wholeInfo=matrix(data=NA,nrow=50000,ncol=9);
  topologicInfo=matrix(data=NA,nrow=fileCount,ncol=7);
  counter=0;
  subnetCount=1;
  
  for (file in file_list)    
  {
    
    counter=counter+1;
    
    #temp_dataset <-read.table(file, header=TRUE, quote="\""); # temp_dataset <-read.table(file, header=TRUE, sep="\t");
    
    ############################ RDA FILES READING ###############
    #file_name=basename(file);
    
    #full_name=paste(mypath,"/",file_name,sep="")
    #jName=substr(file_name,1,nchar(file_name)-4)
    
    #load(file=full_name) 
    
    
    #temp_dataset=get(jName)
    
    ############################ text FILES READING ###############
    
    temp_dataset <-read.table(file, header=TRUE, quote="\""); # temp_dataset <-read.table(file, header=TRUE, sep="\t");
    
    #############
    
    if(nrow(temp_dataset)==0) next;
    file_name=basename(file);
    mypath=dirname(file);
   
    
    #file_data=as.vector(as.matrix(temp_dataset));
    G <- graph.data.frame(temp_dataset,directed=FALSE);
    # A<-as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE,attr=NULL); # ONEMLIIIII
    
    ##########################################TOPOLOGICAL FEATURE EXTRACTION #####################
    
    topVal=topology_simple(G);
    
    numNode=topVal[1][[1]];
    numEdge=topVal[2][[1]];
    avgNumNeigh=topVal[6][[1]];
    diamet= topVal[8][[1]];
    density= topVal[9][[1]];
    cluscoef=topVal[10][[1]];
    
    topologicInfo[counter,1]= file_name;
    topologicInfo[counter,2]= numNode;
    topologicInfo[counter,3]= numEdge;
    topologicInfo[counter,4]= avgNumNeigh;
    topologicInfo[counter,5]= diamet;
    topologicInfo[counter,6]= density;
    topologicInfo[counter,7]= cluscoef;
    
    colnames(topologicInfo)=c("file_name","numNode","numEdge","avgNumNeigh","diamet","density","cluscoef");
    
    
    
    # myDegree<-topology_degree(G,power.law=TRUE);
    
    
    
    ##############################################################################################
    
    # myFlag=TRUE;
    
    result<- tryCatch(mcode(G,vwp=0.005,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE),finally=print(''), error= function(err) NULL)
    
    #        result <- mcode(G,vwp=0.1,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE);
    #        print("okk11111");
    #     }, error= function(err)
    #     {
    #       myFlag=FALSE;
    #       #next;
    #  });
    if(is.null(result))
    {
      next;
    }
    
    #print("ooooookkk");
    
    #result <- mcode(G,vwp=0.1,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE);
    temp=result$COMPLEX;
    number_of_subnet=length(temp);
    
    
    tempCounter=0;
    
    
    
    for(i in 1:number_of_subnet)
    {
      
      
      #tempCounter=tempCounter+1;
      cluster1<-induced.subgraph(G,result$COMPLEX[[i]]); 
      mygenes=V(cluster1)$name;
      file_data=as.vector(mygenes);
      file_data <- file_data[grep("^[A-Za-z0-9]+$", file_data, perl=TRUE)];
      
     # print(length(file_data));
      if(length(file_data)<3 | length(file_data)>500 ) next;
      
      
      
      tempCounter=tempCounter+1;
      
      tryCatch(  #fault = tryCatch(
{  
  ############# biomart ############
  #       gtruthFile=getBM(attributes=c("goslim_goa_accession","goslim_goa_description"),filters="hgnc_symbol",values=file_data, mart=ensembl);
  #       tempIndice=which(as.matrix(gtruthFile[,1])!="" & as.matrix(gtruthFile[,2])!="");
  #       gtruthFile=gtruthFile[tempIndice,]; # inferred go terms of file
  #       allEstimations=nrow(gtruthFile);
  
  ##############gProfileR GSEA #############
  gtruthFile=gprofiler(file_data, organism = "hsapiens",correction_method="bonferroni");
  #gtruthFile=gprofiler(file_data, organism = "hsapiens",correction_method="bonferroni");
  gtruthFile= as.matrix(gtruthFile[,9]);
  allEstimations=nrow(gtruthFile);
  overlapTerms=which(as.matrix(gtruthFile[,1]) %in% as.matrix(gtruthOmim[,1]));
  
  
  
  tp=length(overlapTerms); # number of TP
  fp=allEstimations-tp;
  fn= allGoTermsNo- tp;
  #print(tp)
  #print("okkkkkkkkkkkkkkkkk")
  
  prec= tp / allEstimations;
  prec=as.numeric(prec);
  rec= tp / (tp+fn);
  rec=as.numeric(rec);
  fmeasure= (2*rec*prec)/ (rec+prec);
  
  wholeInfo[subnetCount,1]=file_name;
  wholeInfo[subnetCount,2]=mypath;
  # wholeInfo[counter,2]=fin$precision;
  #wholeInfo[subnetCount,3]=length(file_data); # shows how many genes in which group contains
  #wholeInfo[subnetCount,4]= counter;
  #wholeInfo[subnetCount,5]= tempCounter;
  wholeInfo[subnetCount,3]=tp;
  wholeInfo[subnetCount,4]=fp;
  wholeInfo[subnetCount,5]=fn;
  wholeInfo[subnetCount,6]=prec;
  wholeInfo[subnetCount,7]=rec;
  wholeInfo[subnetCount,8]=fmeasure; # 11_subat_coomentli_ilk_hali
  wholeInfo[subnetCount,9]= length(file_data);
  
  subnetCount=subnetCount+1;
  # fault=FALSE;
  # return (fault);
  
}, error= function(err)
{
  #return (res);
  assign("biolEval", wholeInfo, envir = res);
  assign("topological", topologicInfo, envir = res);
  return (res);
  #fault=TRUE;
  
}
      )

#print(fmeasure);

#       if(fault==TRUE) 
#       {
#         assign("biolEval", wholeInfo, envir = res);
#         assign("topological", topologicInfo, envir = res);
#         return (res);
#       }
#       

    }    

#print("#####");
#print(counter);

  }
colnames(wholeInfo)=c("file_name","path","tp","fp","fn","prec","rec","f_measure","genesInGroup");
ind=which(is.na(wholeInfo[,1])==FALSE);
wholeInfo=wholeInfo[ind,];
#res <- new.env();
#   if(fault==FALSE)
#   {
assign("biolEval", wholeInfo, envir = res);
assign("topological", topologicInfo, envir = res);

#rm(list = setdiff(ls(), lsf.str()))

return (res);  
#}



}