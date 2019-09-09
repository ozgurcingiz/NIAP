#'@export
calculateTopologicalFitness <- function (mypath="")
{
  #mypath=paste(getwd(),"/",mypath,sep="");
  #library(igraph);
  
  mypath=choose.dir(default = "", caption = "Select folder")
  
  file_list=list.files(path = mypath, pattern = NULL, all.files = TRUE,full.names = TRUE, recursive = TRUE,ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE);
  fileCount=length(file_list);
  
  
  #wholeInfo=matrix(data=NA,nrow=50000,ncol=10);
  topologicInfo=matrix(data=NA,nrow=fileCount,ncol=4);
  counter=0;
  
  
  for (file in file_list)    
  {
    
    counter=counter+1;
    
    #temp_dataset <-read.table(file, header=TRUE, quote="\""); # temp_dataset <-read.table(file, header=TRUE, sep="\t");
    
    ########################
    #file_name=basename(file);
  
    #full_name=paste(mypath,"/",file_name,sep="")
    #jName=substr(file_name,1,nchar(file_name)-4)
    
    #load(file=full_name) 
    
    #temp_dataset=get(jName)
    
    ######################################
    
    
    temp_dataset <-read.table(file, header=TRUE, quote="\""); # temp_dataset <-read.table(file, header=TRUE, sep="\t");
    
    
    if(nrow(temp_dataset)==0) next;
    file_name=basename(file);
    #print(file_name);
    #file_data=as.vector(as.matrix(temp_dataset));
    G <- graph.data.frame(temp_dataset,directed=FALSE);
    # A<-as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE,attr=NULL); # ONEMLIIIII
    
    ##########################################TOPOLOGICAL FEATURE EXTRACTION #####################
    
    #     topVal=topology_simple(G);
    #     
    #     numNode=topVal[1][[1]];
    #     numEdge=topVal[2][[1]];
    #     avgNumNeigh=topVal[6][[1]];
    #     diamet= topVal[8][[1]];
    #     density= topVal[9][[1]];
    #     cluscoef=topVal[10][[1]];
    
    #mydegree<-degree(G);
    d<-topology_degree(G); #d<-topology_degree(G,power.law=TRUE);
    mypower<-power.law.fit(d[,2]);
    
    rsquare=calculateRSquare (G);
    
    #mypower<-power.law.fit(mydegree);
    
    
    topologicInfo[counter,1]= file_name;
    topologicInfo[counter,2]= mypower$KS.stat;
    topologicInfo[counter,3]= mypower$KS.p;
    topologicInfo[counter,4]= rsquare;
    #     topologicInfo[counter,5]= diamet;
    #     topologicInfo[counter,6]= density;
    #     topologicInfo[counter,7]= cluscoef;
    colnames(topologicInfo)=c("file_name","KS.stat","KS.p","rsquare");
    
    
  }
  
  return (topologicInfo);  
  
  
  
  
}