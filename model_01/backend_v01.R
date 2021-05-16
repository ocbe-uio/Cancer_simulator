
####---------

# load libs
require(truncnorm)
require(colorspace)
require(scales)

# --

#' Initialize Blood vessels
#'
#' Generate a random map of blood vessels to be used in the modeling of cell movement with cellmigRation
#'
#' @param mar integer, Minimum distance between the blood vessels and the sides of the tumor files
#' @param BloVesN integer, Number of blood vessels
#' @param Rrows integer, A number between 1 and 5 determining number of red boxes in the red zone. Specifically:
#'   \item{Rrows= 1 that means that there are 9 boxes in total (8 red + .lumen) }
#'   \item{Rrows= 2 that means that there are 25 boxes in total (24 red + .lumen) }
#'   \item{Rrows= 3 that means that there are 49 boxes in total (48 red + .lumen) }
#'   \item{Rrows= 4 that means that there are 81 boxes in total (80 red + .lumen) }
#'   \item{Rrows= 5 that means that there are 121 boxes in total (120 red + .lumen) }
#' @param minBloVesSep integer, Minimum distance between blood vessels'
#'
#' @param VesselsXY2CSVFile filename used for saving the VesselsXY results (CSV file)
#' @param LatDF2CSVFile filename used for saving the LatDF results (CSV file)
#'
#' @param plot logical, shall a plot be generated to visualize the distribution of vessels on the map
#' @param verbose logical, shall messages be printed to inform user about progress
#' @param seed integer, sets the seed for the random operation
#'
#' @return list including three elements:
#'
#'   \item{VesselsXY, data.frame including all vessel positions on the map }
#'   \item{LatDF, data.frame with 2 columns, i.e. lattice and LatRed }
#'   \item{vesselParams, list including the values of all parameters used for building the vessels map }
#'
#' @export
initializeBloodVessels <- function(mar= 100, BloVesN= 36,  Rrows= 4,
                                   minBloVesSep= 400,
                                   VesselsXY2CSVFile = NULL,
                                   LatDF2CSVFile = NULL,
                                   plot = TRUE, verbose = TRUE,
                                   seed = 12345){

  if (verbose)
    message("Generating vessels")
  # inner function
  inside.circle = function(center, radius, p, plot = FALSE){
    if ((length(center)) !=2){stop ("First argument must be a vector of length two.")}
    if ( radius<=0 ) {stop ("Second argument is not a positive number. Radius must be positive.")}
    if ( (length(radius)) !=1 ) {stop ("Second argument must contain only 1 element.") }

    p.1 = length(p[,1])
    x.coord = numeric(p.1)
    y.coord = numeric(p.1)
    for (i in 1:p.1){
      if ((p[i,1] - center[1])^2 + (p[i,2] - center[2])^2 <= radius^2){
        x.coord[i] = p[i,1]
        y.coord[i] = p[i,2]
      }
    }
    x.coord = x.coord[x.coord != 0]
    y.coord = y.coord[y.coord != 0]
    coordinates = cbind(x.coord,y.coord)
    theta = seq(0,2*pi,length = 2000)
    a = center[1]
    b = center[2]
    x = a + (radius*cos(theta))
    y = b + (radius*sin(theta))

    #Plot cirlce and add points
    if (plot) {
      plot(x,y, type = 'l', xlim = c(0,4000),
           ylim = c(0,4000), ylab = "Y", xlab = "X",col = "red")
      points(p[,1], p[,2], col = "blue", cex=1, pch=19)
      #The points inside the cirlce are red.
      points(x.coord,y.coord, col = "black",cex=1, pch=19)
    }
    return(coordinates)
  }


  # setSeed
  if (is.null(seed)) {
    seed <- sample(x = 1:999999, size = 1)
  } else if (is.na(seed)) {
    seed <- sample(x = 1:999999, size = 1)
  }
  set.seed(seed = seed)

  #initialize DF
  df<-data.frame()
  df[1,1]= runif(1,mar,4000-mar)
  df[1,2]= runif(1,mar,4000-mar)
  for (i in 2:BloVesN){
    if (verbose)
      message(".", appendLF = FALSE)
    df[i,1]= runif(1,mar,4000-mar)
    df[i,2]= runif(1,mar,4000-mar)
    inV<-c()
    for(j in 1:(i-1)){
      p=matrix(c(df[i,1],df[i,2]),1,2)
      cen=c(df[j,1],df[j,2])
      r= minBloVesSep + (Rrows * 20 )
      pp<-inside.circle(cen,r,p)
      inV[j]=pp[1]
    }
    while(length(inV[!is.na(inV)])>0){
      df[i,1]= runif(1,mar,4000-mar)
      df[i,2]= runif(1,mar,4000-mar)
      inV<-c()
      for(j in 1:(i-1)){
        p=matrix(c(df[i,1],df[i,2]),1,2)
        cen=c(df[j,1],df[j,2])
        r= minBloVesSep + (Rrows * 20 )
        pp<-inside.circle(cen,r,p)
        inV[j]=pp[1]
      }
    }
  }

  # convert the X and  Y into boxes
  for (i in 1:BloVesN){
    df[i,3]=(floor(df[i,2]/20)*200)+ceiling(df[i,1]/20)
  }

  colnames(df)<-c("x","y", ".lumen")
  #write.csv(df,file="VesselsXY.csv")
  #save(df,file="VesselsXY.Rdata")



  # Drawing squares
  xleft<-df[,1]- (10 + (Rrows * 20))
  ybottom<-df[,2]-(10 + (Rrows * 20))
  xright<-df[,1]+ (10 + (Rrows * 20))
  ytop<-df[,2]+(10 + (Rrows * 20))

  xleftb<-df[,1]-10
  ybottomb<-df[,2]-10
  xrightb<-df[,1]+10
  ytopb<-df[,2]+10

  # plot
  if (plot) {
    plot(0,0,type = "n", xlim = c(0,4000), ylim = c(0,4000),xlab="X",ylab="Y")
    for (i in 1:BloVesN){
      rect(xleft, ybottom, xright, ytop,col="red",border="red")
    }
    for (i in 1:BloVesN){
      rect(xleftb, ybottomb, xrightb, ytopb,col="black",border="black")
    }
  }


  # making a lattice
  VESdf<-df
  lattice<-c(rep(0,40000))
  length(lattice)
  LatVes<-c()
  for (i in 1:BloVesN){
    LatVes[i]= (floor(VESdf[i,2]/20)*200)+ceiling(VESdf[i,1]/20)
  }
  # print (LatVes)

  for (i in LatVes){    # Now the lattice has the .lumen of the vessels
    lattice[i]=1
  }

  # making a lattice for the red areas
  Red=c()
  TotalRow=c(1:Rrows)
  if (length (TotalRow)<1 | length (TotalRow)> 5)  {
    stop( "Rrows has to be a number between 1 and 5")
  }


  floors1<-c(-200,0,200,NA,NA,NA,NA,NA,NA,NA,NA)
  floors2<-c(-400,-200,0,200,400,NA,NA,NA,NA,NA,NA)
  floors3<-c(-600,-400,-200,0,200,400,600,NA,NA,NA,NA)
  floors4<-c(-800,-600,-400,-200,0,200,400,600,800,NA,NA)
  floors5<-c(-1000,-800,-600,-400,-200,0,200,400,600,800,1000)
  floors<- data.frame(floors1,floors2,floors3,floors4,floors5)

  cDF=floors[,Rrows]
  Floors=cDF[complete.cases(cDF)]

  for (i in LatVes){
    for (n in Floors){
      if (Rrows==1) {
        columns<- c(i+n,i+n+1,i+n-1)
      } else if (Rrows==2) {
        columns<- c(i+n,i+n+1,i+n+2,i+n-1,i+n-2)
      } else if (Rrows==3) {
        columns<- c(i+n,i+n+1,i+n+2,i+n+3,i+n-1,i+n-2,i+n-3)
      } else if (Rrows==4) {
        columns<- c(i+n,i+n+1,i+n+2,i+n+3,i+n+4,i+n-1,i+n-2,i+n-3,i+n-4)
      } else if (Rrows==5) {
        columns<- c(i+n,i+n+1,i+n+2,i+n+3,i+n+4,i+n+5,i+n-1,i+n-2,i+n-3,i+n-4,i+n-5)
      }

      Red=c(Red,columns)
    }
  }

  # removing the .lumena from the red zones
  del<-c()
  for (i in LatVes){
    del<-c(del,which(Red==i))
  }
  RED<-Red[-del]

  #length(Red)
  #length(RED)

  # labeling the boxes of the red zones in the LatRed with "R"
  # generating a lattice with red boxes
  LatRed<-c(rep(0,40000))
  for (i in RED){
    LatRed[i]="R"
  }

  # Making a data frame that has lattice and LatRed
  LatDF<-data.frame(lattice = lattice,LatRed = factor(LatRed))

  if (!is.null(VesselsXY2CSVFile)) {
    write.csv(df,file=VesselsXY2CSVFile, row.names = FALSE)
  }
  if (!is.null(LatDF2CSVFile)) {
    write.csv(LatDF,file=LatDF2CSVFile, row.names = FALSE)
  }


  if (verbose)
    message("Done!")

  return(list(VesselsXY = df, LatDF = LatDF,
              vesselParams = list(mar = mar,
                                  BloVesN = BloVesN,
                                  Rrows = Rrows,
                                  minBloVesSep = minBloVesSep,
                                  seed = seed)))

}


#' Load Blood vessels
#'
#' Load a map of blood vessels to be used in the modeling of cell movement with cellmigRation
#'
#' @param VesselsXY_CSVFile filename including the VesselsXY results (CSV file)
#' @param LatDF_CSVFile filename including the LatDF results (CSV file)
#' @param plot logical, shall a plot be generated to visualize the distribution of vessels on the map
#'
#' @return list including three elements:
#'
#'   \item{VesselsXY, data.frame including all vessel positions on the map }
#'   \item{LatDF, data.frame with 2 columns, i.e. lattice and LatRed }
#'   \item{vesselParams, list including the values of all parameters used for building the vessels map }
#'
#' @export
loadBloodVessels <- function(VesselsXY_CSVFile = NULL,
                             LatDF_CSVFile = NULL,
                             plot = TRUE) {

  # Check
  stopifnot(
    file.exists(VesselsXY_CSVFile),
    file.exists(LatDF_CSVFile))

  # Load files
  Ves <- read.csv(VesselsXY_CSVFile, as.is = TRUE, header = TRUE)
  Lat <- read.csv(LatDF_CSVFile, as.is = TRUE, header = TRUE)

  # COmpute ori params
  VNum <- as.numeric(table(Lat$LatRed))[2]
  NCel <- as.numeric(table(Lat$lattice))[2]
  MCoef <- as.numeric(VNum / NCel)

  # pre-def
  tmp <- data.frame(R = c(1,  2,  3,  4,   5),
                    V = c(8, 24, 48, 80, 120))
  idx <- which.min(abs(tmp$V - MCoef))
  rrow <- tmp[idx, 'R']


  # Plot
  # Drawing squares
  df <- Ves
  Rrows <- rrow
  xleft<-df[,1]- (10 + (Rrows * 20))
  ybottom<-df[,2]-(10 + (Rrows * 20))
  xright<-df[,1]+ (10 + (Rrows * 20))
  ytop<-df[,2]+(10 + (Rrows * 20))

  xleftb<-df[,1]-10
  ybottomb<-df[,2]-10
  xrightb<-df[,1]+10
  ytopb<-df[,2]+10

  # plot
  if (plot) {
    plot(0,0,type = "n", xlim = c(0,4000), ylim = c(0,4000),xlab="X",ylab="Y")
    for (i in 1:VNum){
      rect(xleft, ybottom, xright, ytop,col="red",border="red")
    }
    for (i in 1:VNum){
      rect(xleftb, ybottomb, xrightb, ytopb,col="black",border="black")
    }
  }

  # return
  return(list(VesselsXY = Ves, LatDF = Lat,
              vesselParams = list(mar = NA,
                                  BloVesN = NCel,
                                  Rrows = rrow,
                                  minBloVesSep = NA,
                                  seed = NA)))
}



# --
# Modeling Functions


eval_function <- function(mean_sd){
  mean <- mean_sd[1]
  sd <- mean_sd[2]
  sample <- rtruncnorm(n = DT, a = minSpe, b = maxSpe, mean = mean, sd = sd)
  mean_diff <-abs((desired_meanSpe - mean(sample))/desired_meanSpe)
  sd_diff <- abs((desired_sdSpe - sd(sample))/desired_sdSpe)
  mean_diff + sd_diff
}

Density<-function(B,Lat){
  Den<-c()
  bot<-c(2:199)
  top<-c(39802:39999)
  left<-c()
  right<-c()
  for (r in 1:198){ left<-c(left,1+(r*200))}
  for (r in 1:198){ right<-c(right,200+r*200)}
  loc<-c()

  if (B %in% bot){loc<-c(B+199,B+200,B+201,B-1,B+1,B+39799,B+39800,B+39801)}
  else if  (B %in% top){loc<-c(B-1,B+1,B-201,B-200,B-199,B-39799,B-39800,B-39801)}
  else if  (B %in% left){loc<-c(B+200,B+201,B+1,B-200,B-199,B+199,B+399,B-1)}
  else if  (B %in% right){loc<-c(B+199,B+200,B-1,B-201,B-200,B+1,B-199,B-399)}

  else if  (B ==1){loc<-c(201,202,2,39801,39802,200,400,40000)}
  else if  (B ==200){loc<-c(400,399,199,201,1,40000,39999,39801)}
  else if  (B ==39801){loc<-c(39802,1,2,200,39601,39602,40000,39800)}
  else if  (B ==40000){loc<-c(39999,39799,39800,39801,39601,200,199,1)
  }else{loc=c(B+200,B+201,B+199,B-200,B-201,B-199,B+1,B-1)}
  #for (p in loc){PP<-c(PP,Lat[p])}
  Dens=sum(Lat[loc])
  if (Dens<=2){ Den=2}
  if (Dens<=4 & Dens>2){ Den=1}
  if (Dens<=6 & Dens>4){ Den=0}
  if (Dens==7){ Den=-1}
  if (Dens==8){ Den=-2}

  DenPos<-Lat[loc]     # to know which boxes are empty
  return (list(Den,loc[DenPos==0],DenPos,loc))
}

ForMovement<-function(x1,y1,x2,y2,speed1,speed2,B,Lat){
  alp2=pi/2
  alp1= acos(speed1/sqrt(speed1^2 +  speed2^2))          # initial data
  u=x2-x1
  v=y2-y1
  if (u==0 && v==0){
    u=1
    v=1
  }
  a3=sqrt(u^2+v^2)
  alp3=pi-(alp1+alp2)
  a2=a3*sin(alp2)/sin(alp3)
  RHS1=x1*u+y1*v+a2*a3*cos(alp1)
  RHS2=y2*u-x2*v+a2*a3*sin(alp1)
  x3=(1/a3^2)*(u*RHS1-v*RHS2)
  y3=(1/a3^2)*(v*RHS1+u*RHS2)
  Xr= x2
  Yr= y2
  Xs= x3
  Ys= y3
  An=seq(184,356, by=4)
  AngVx<-c()
  AngVy<-c()
  AvX<-c()
  AvY<-c()
  EmpB<-c()
  AvB<-c()

  for(i in An){
    Nx=Xr + cos(i*pi/180)*(Xs - Xr) - sin(i*pi/180)*(Ys-Yr)
    AngVx[match(i,An)]=Nx
    Ny=Yr + cos(i*pi/180)*(Ys - Yr) + sin(i*pi/180)*(Xs-Xr)
    AngVy[match(i,An)]=Ny
  }
  AngVx<- na.omit(AngVx)
  for (i in 1:length(AngVx)){AngVx[i]<-margin(AngVx[i])}
  AngVy<- na.omit(AngVy)
  for (i in 1:length(AngVy)){AngVy[i]<-margin(AngVy[i])}

  for (i in 1:length(An)){
    EmpB[i]<-(floor(AngVy[i]/20)*200)+ceiling (AngVx[i]/20)
  }
  for (z in 1:length(An)){
    if (EmpB[z]== B){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])

    }
    if (EmpB[z]!= B & Lat[EmpB[z]]== 0){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])
    }
  }
  return(list(v1=AvB,v2=AvX,v3=AvY))
}

BackMovement<-function(x1,y1,x2,y2,speed1,speed2,B,Lat){
  alp2=pi/2
  alp1= acos(speed1/sqrt(speed1^2 +  speed2^2))          # initial data
  u=x2-x1
  v=y2-y1
  if (u==0 && v==0){
    u=1
    v=1
  }
  a3=sqrt(u^2+v^2)
  alp3=pi-(alp1+alp2)
  a2=a3*sin(alp2)/sin(alp3)
  RHS1=x1*u+y1*v+a2*a3*cos(alp1)
  RHS2=y2*u-x2*v+a2*a3*sin(alp1)
  x3=(1/a3^2)*(u*RHS1-v*RHS2)
  y3=(1/a3^2)*(v*RHS1+u*RHS2)
  Xr= x2
  Yr= y2
  Xs= x3
  Ys= y3
  An=-1*seq(180,360, by=4)
  AngVx<-c()
  AngVy<-c()
  AvX<-c()
  AvY<-c()
  EmpB<-c()
  AvB<-c()

  for(i in An){
    Nx=Xr + cos(i*pi/180)*(Xs - Xr) - sin(i*pi/180)*(Ys-Yr)
    AngVx[match(i,An)]=Nx
    Ny=Yr + cos(i*pi/180)*(Ys - Yr) + sin(i*pi/180)*(Xs-Xr)
    AngVy[match(i,An)]=Ny
  }
  AngVx<- na.omit(AngVx)
  for (i in 1:length(AngVx)){AngVx[i]<-margin(AngVx[i])}
  AngVy<- na.omit(AngVy)
  for (i in 1:length(AngVy)){AngVy[i]<-margin(AngVy[i])}

  for (i in 1:length(An)){
    EmpB[i]<-(floor(AngVy[i]/20)*200)+ceiling(AngVx[i]/20)
  }
  for (z in 1:length(An)){
    if (EmpB[z]== B){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])

    }
    if (EmpB[z]!= B & Lat[EmpB[z]]== 0){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])
    }
  }
  return(list(v1=AvB,v2=AvX,v3=AvY))
}

marginXY<-function(NX,NY){
  if(NX>4000){
    NNX=NX-4000
  }else if(NX<0){
    NNX=4000-abs(NX)
  }else{
    NNX=NX
  }

  if(NY>4000){
    NNY=NY-4000
  }else if(NY<0){
    NNY=4000-abs(NY)
  }else{NNY=NY
  }
  return(c(NNX,NNY))
}

margin<-function(object){
  if(object>4000){
    object=object-4000
  }else if(object<0){
    object=4000-abs(object)
  }else{
    object=object
  }
  return(object)
}


NewDauLoc<-function(object){
  xy<-c()
  y=((ceiling(object/200))*20)-10
  x=((abs(object - (floor(object/200))*200))* 20 )-10
  if (x<=0){x=3990}
  xy<-c(x,y)
  return(xy)
}


PhylogenticTree<- function (phylogenicDF,TargetCell){
  steps<-c(1:TargetCell)
  RevSteps<-rev(steps)
  PHYLOGENIC<-c()
  for (i in RevSteps){
    if (TargetCell %in% phylogenicDF[i,]){
      PHYLOGENIC=phylogenicDF[i,]
      PHYLOGENIC=PHYLOGENIC[!is.na(PHYLOGENIC)]
      break
    }
  }
  return(PHYLOGENIC)
}

SimSum<-function(CellTables,cdDF,mainDF,SimTime,AA,CellDeath,RCellNum,
                 LatticeTable,CellSenescence,FinalAveSpe,MCells,V,PER,
                 G1Division,DividingCells,NotMotileCells,RanNum){
  simsumDF<-data.frame()
  Days<- seq(1,SimTime,by=24)
  totalSteps<-c(Days,SimTime)
  RCellNumPerStep<-c()
  CellNumPerStep<-c()
  totalnumCells<-c()
  CellDeathDay<-c()
  CellSenescenceDay<-c()
  Metastasis<-c()

  TTR=cumsum(RCellNum)
  TTIRZ<-TTR[totalSteps]
  for (i in totalSteps){              ### computing the number of cells in red per step
    main<-mainDF[i,]
    m1200<-main[!is.na(main)]
    NRm1200<-c()
    for (n in m1200){
      if (AA[n]==0){
        NRm1200<-c(NRm1200,n)
      }
    }
    RCellNumPerStep<- c(RCellNumPerStep,(length(m1200)- length(NRm1200)))
    CellNumPerStep<- c(CellNumPerStep,length(m1200))
    CellDeathDay<-c(CellDeathDay,CellDeath[i])
    CellSenescenceDay<-c(CellSenescenceDay,CellSenescence[i])
    Metastasis<-c(Metastasis,MCells[i])

    len1200<-cdDF[i,]         ########## computing the ID of the last cell
    FinalnumCell<-length(len1200[!is.na(len1200)])
    totalnumCell<-len1200[!is.na(len1200)]
    nameDF<-rbind(len1200,names(len1200))
    rownames(nameDF)<-NULL
    colnames(nameDF)<-nameDF[2,]
    w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
    totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
    totalnumCells<-c(totalnumCells,totalnumCell)
  }
  visits<-c()
  MaxCell<-c(1:totalnumCell)    ##### to check every cell
  for(d in MaxCell){
    targetCell<-CellTables[[d]][,7]
    if (length(targetCell[!is.na(targetCell)])>0 && "1" %in% targetCell[!is.na(targetCell)]){ visits[d]=d}else{visits[d]=NA}
  }
  visits=visits[!is.na(visits)]

  Day<-totalSteps        ############## computing uniformity index
  UIvalues<-c()
  vDays<- Day[-c(1:12)]                  # Not looking at the UI for the first 12 days since only very few cels are there
  for (i in vDays){
    TS=i
    UI<-UniformityIndex(TS,LatticeTable)
    UIvalues<-c(UIvalues,UI)
  }
  UIVALUES<-c(rep(0, 12),UIvalues)

  simsumDF[1:length(totalSteps),1]<-totalSteps
  simsumDF[1:length(totalSteps),2]<-CellNumPerStep
  simsumDF[1:length(totalSteps),3]<-RCellNumPerStep
  simsumDF[1:length(totalSteps),4]<-round(RCellNumPerStep/CellNumPerStep,digits=2)
  simsumDF[1:length(totalSteps),5]<- totalnumCells
  simsumDF[1:length(totalSteps),6]<- CellDeathDay
  simsumDF[1:length(totalSteps),7]<- CellSenescenceDay

  simsumDF[1:length(totalSteps),8]<-PER
  simsumDF[1:length(totalSteps),9]<-V
  simsumDF[1:length(totalSteps),10]<-TTIRZ
  simsumDF[1:length(totalSteps),11]<- UIVALUES
  simsumDF[1:length(totalSteps),12]<- c(rep(NA,(length(totalSteps)-1)),length(visits))
  simsumDF[1:length(totalSteps),13]<- Metastasis
  simsumDF[1:length(totalSteps),14]<- c(rep(NA,(length(totalSteps)-1)),FinalAveSpe)
  simsumDF[1:length(totalSteps),15]<- c(rep(NA,(length(totalSteps)-1)),length(G1Division))
  simsumDF[1:length(totalSteps),16]<- c(rep(NA,(length(totalSteps)-1)),NotMotileCells)
  simsumDF[1:length(totalSteps),17]<- c(rep(NA,(length(totalSteps)-1)),DividingCells)

  colnames(simsumDF)<-c("Steps","Number of cells","Number of cells in RedZone","RedZone cells ratio","MaxCellID",
                        "Cumulative Dead Cells","Cumulative Senescent cells","Persistence","Speed(um/h)","Cumulative time in RedZone",
                        "Uniformity Index","Total Number of RedZone Visitors","Cumulative Metastatic cells ","Average speed of all cells","Cumulative G1B cells",
                        "L5S Immotile Cells","L5S Dividing Cells")
  write.csv(simsumDF,file=paste0("simsumDF","R",RanNum,".csv"))
  return(simsumDF)
}


UniformityIndex<-function(TS,LatticeTable){
  lat1<-LatticeTable[TS,]
  smallBoxes<- c(0:49)

  oneOne<-seq(1,9801,by=200)
  oneOneBoxes<-c()
  for (i in oneOne){
    for (n in smallBoxes){
      oneOneBoxes<-c(oneOneBoxes,oneOne[match(i,oneOne)] + n)
    }}
  twoOneBoxes<- oneOneBoxes + 50
  threeOneBoxes<- oneOneBoxes + 100
  fourOneBoxes<- oneOneBoxes + 150
  oneTwo<-seq(10001,19801,by=200)
  oneTwoBoxes<-c()
  for (i in oneTwo){
    for (n in smallBoxes){
      oneTwoBoxes<-c(oneTwoBoxes,oneTwo[match(i,oneTwo)] + n)
    }}
  twoTwoBoxes<- oneTwoBoxes + 50
  threeTwoBoxes<- oneTwoBoxes + 100
  fourTwoBoxes<- oneTwoBoxes+ 150

  oneThree<-seq(20001,29801,by=200)
  oneThreeBoxes<-c()
  for (i in oneThree){
    for (n in smallBoxes){
      oneThreeBoxes<-c(oneThreeBoxes,oneThree[match(i,oneThree)] + n)
    }}

  twoThreeBoxes<- oneThreeBoxes + 50
  threeThreeBoxes<- oneThreeBoxes + 100
  fourThreeBoxes<- oneThreeBoxes + 150

  oneFour<-seq(30001,39801,by=200)
  oneFourBoxes<-c()
  for (i in oneFour){
    for (n in smallBoxes){
      oneFourBoxes<-c(oneFourBoxes,oneFour[match(i,oneFour)] + n)
    }}
  twoFourBoxes<- oneFourBoxes + 50
  threeFourBoxes<- oneFourBoxes + 100
  fourFourBoxes<- oneFourBoxes + 150

  LAT11<-lat1[oneOneBoxes]
  LAT11<-LAT11[LAT11==1]
  LAT11<-length(LAT11)

  LAT21<-lat1[twoOneBoxes]
  LAT21<-LAT21[LAT21==1]
  LAT21<-length(LAT21)

  LAT31<-lat1[threeOneBoxes]
  LAT31<-LAT31[LAT31==1]
  LAT31<-length(LAT31)

  LAT41<-lat1[fourOneBoxes]
  LAT41<-LAT41[LAT41==1]
  LAT41<-length(LAT41)



  LAT12<-lat1[oneTwoBoxes]
  LAT12<-LAT12[LAT12==1]
  LAT12<-length(LAT12)

  LAT22<-lat1[twoTwoBoxes]
  LAT22<-LAT22[LAT22==1]
  LAT22<-length(LAT22)

  LAT32<-lat1[threeTwoBoxes]
  LAT32<-LAT32[LAT32==1]
  LAT32<-length(LAT32)

  LAT42<-lat1[fourTwoBoxes]
  LAT42<-LAT42[LAT42==1]
  LAT42<-length(LAT42)

  LAT13<-lat1[oneThreeBoxes]
  LAT13<-LAT13[LAT13==1]
  LAT13<-length(LAT13)

  LAT23<-lat1[twoThreeBoxes]
  LAT23<-LAT23[LAT23==1]
  LAT23<-length(LAT23)

  LAT33<-lat1[threeThreeBoxes]
  LAT33<-LAT33[LAT33==1]
  LAT33<-length(LAT33)

  LAT43<-lat1[fourThreeBoxes]
  LAT43<-LAT43[LAT43==1]
  LAT43<-length(LAT43)

  LAT14<-lat1[oneFourBoxes]
  LAT14<-LAT14[LAT14==1]
  LAT14<-length(LAT14)

  LAT24<-lat1[twoFourBoxes]
  LAT24<-LAT24[LAT24==1]
  LAT24<-length(LAT24)

  LAT34<-lat1[threeFourBoxes]
  LAT34<-LAT34[LAT34==1]
  LAT34<-length(LAT34)

  LAT44<-lat1[fourFourBoxes]
  LAT44<-LAT44[LAT44==1]
  LAT44<-length(LAT44)
  allLATs<-c(LAT11,LAT21,LAT31,LAT41,
             LAT12,LAT22,LAT32,LAT42,
             LAT13,LAT23,LAT33,LAT43,
             LAT14,LAT24,LAT34,LAT44)

  meanallLATs<-mean(allLATs)
  sdallLATs<-sd(allLATs)
  CV=c()
  CV=round(sdallLATs/meanallLATs,digits=2)
  if (CV<=1){
    UI=(1 -CV )*100
  }else{
    UI=0
  }
  return(UI)
}


UIplot<-function(TS,LatticeTable,RanNum,V,PER){
  lat1<-LatticeTable[TS,]
  smallBoxes<- c(0:49)

  oneOne<-seq(1,9801,by=200)
  oneOneBoxes<-c()
  for (i in oneOne){
    for (n in smallBoxes){
      oneOneBoxes<-c(oneOneBoxes,oneOne[match(i,oneOne)] + n)
    }}
  twoOneBoxes<- oneOneBoxes + 50
  threeOneBoxes<- oneOneBoxes + 100
  fourOneBoxes<- oneOneBoxes + 150
  oneTwo<-seq(10001,19801,by=200)
  oneTwoBoxes<-c()
  for (i in oneTwo){
    for (n in smallBoxes){
      oneTwoBoxes<-c(oneTwoBoxes,oneTwo[match(i,oneTwo)] + n)
    }}
  twoTwoBoxes<- oneTwoBoxes + 50
  threeTwoBoxes<- oneTwoBoxes + 100
  fourTwoBoxes<- oneTwoBoxes+ 150

  oneThree<-seq(20001,29801,by=200)
  oneThreeBoxes<-c()
  for (i in oneThree){
    for (n in smallBoxes){
      oneThreeBoxes<-c(oneThreeBoxes,oneThree[match(i,oneThree)] + n)
    }}

  twoThreeBoxes<- oneThreeBoxes + 50
  threeThreeBoxes<- oneThreeBoxes + 100
  fourThreeBoxes<- oneThreeBoxes + 150

  oneFour<-seq(30001,39801,by=200)
  oneFourBoxes<-c()
  for (i in oneFour){
    for (n in smallBoxes){
      oneFourBoxes<-c(oneFourBoxes,oneFour[match(i,oneFour)] + n)
    }}
  twoFourBoxes<- oneFourBoxes + 50
  threeFourBoxes<- oneFourBoxes + 100
  fourFourBoxes<- oneFourBoxes + 150

  LAT11<-lat1[oneOneBoxes]
  LAT11<-LAT11[LAT11==1]
  LAT11<-length(LAT11)

  LAT21<-lat1[twoOneBoxes]
  LAT21<-LAT21[LAT21==1]
  LAT21<-length(LAT21)

  LAT31<-lat1[threeOneBoxes]
  LAT31<-LAT31[LAT31==1]
  LAT31<-length(LAT31)

  LAT41<-lat1[fourOneBoxes]
  LAT41<-LAT41[LAT41==1]
  LAT41<-length(LAT41)



  LAT12<-lat1[oneTwoBoxes]
  LAT12<-LAT12[LAT12==1]
  LAT12<-length(LAT12)

  LAT22<-lat1[twoTwoBoxes]
  LAT22<-LAT22[LAT22==1]
  LAT22<-length(LAT22)

  LAT32<-lat1[threeTwoBoxes]
  LAT32<-LAT32[LAT32==1]
  LAT32<-length(LAT32)

  LAT42<-lat1[fourTwoBoxes]
  LAT42<-LAT42[LAT42==1]
  LAT42<-length(LAT42)

  LAT13<-lat1[oneThreeBoxes]
  LAT13<-LAT13[LAT13==1]
  LAT13<-length(LAT13)

  LAT23<-lat1[twoThreeBoxes]
  LAT23<-LAT23[LAT23==1]
  LAT23<-length(LAT23)

  LAT33<-lat1[threeThreeBoxes]
  LAT33<-LAT33[LAT33==1]
  LAT33<-length(LAT33)

  LAT43<-lat1[fourThreeBoxes]
  LAT43<-LAT43[LAT43==1]
  LAT43<-length(LAT43)

  LAT14<-lat1[oneFourBoxes]
  LAT14<-LAT14[LAT14==1]
  LAT14<-length(LAT14)

  LAT24<-lat1[twoFourBoxes]
  LAT24<-LAT24[LAT24==1]
  LAT24<-length(LAT24)

  LAT34<-lat1[threeFourBoxes]
  LAT34<-LAT34[LAT34==1]
  LAT34<-length(LAT34)

  LAT44<-lat1[fourFourBoxes]
  LAT44<-LAT44[LAT44==1]
  LAT44<-length(LAT44)
  allLATs<-c(LAT11,LAT21,LAT31,LAT41,
             LAT12,LAT22,LAT32,LAT42,
             LAT13,LAT23,LAT33,LAT43,
             LAT14,LAT24,LAT34,LAT44)

  meanallLATs<-mean(allLATs)
  sdallLATs<-sd(allLATs)
  CV<- round(sdallLATs/meanallLATs,digits=2)
  if (CV<=1){ UI=(1 - CV)*100
  }else{
    UI=0
  }

  jpeg(paste0("UniformityPlot-","PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),
       width = 6, height = 6, units = 'in', res = 500)
  barplot(allLATs, main=paste0("UniformityPlot-","PR=",PER,"velecity=",V,"R",RanNum),
          xlab="Square Quadrat",ylab="Number of cells per Square Quadrat",
          sub=(paste0("Uniformity Index =",UI,"%"))  )
  dev.off()

}


PlottingCellsADLastStep<- function(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows){

  BloVesN <- nrow(VesselsXY)

  len1200<-cdDF[SimTime,]         ########## computing the ID of the last cell
  FinalnumCell<-length(len1200[!is.na(len1200)])
  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
  totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
  color<-c()
  nickColors <- function(n, h = c(120,480), l = c(.60,.70), s = c(.8,1), alpha = 1){
    return (alpha(hex(HLS(seq(h[1],h[2],length.out = n), seq(l[1],l[2],length.out = n), seq(s[1],s[2],length.out=n))), alpha))
  }
  n=totalnumCell
  color = nickColors(n)

  s=SimTime
  len1200<-mainDF[s,]         ########## computing the ID of the cell ID in each step
  FinalnumCell<-length(len1200[!is.na(len1200)])
  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  CellID<-as.numeric(nameDF[2,][!is.na(nameDF[1,])])

  lenID<-mainDF[s,]
  x=c(rep(NA,length(lenID[!is.na(lenID)])))
  y=c(rep(NA,length(lenID[!is.na(lenID)])))

  X=c()
  Y=c()
  boxXY<-data.frame(CellID,lenID[!is.na(lenID)],x,y)
  for (t in 1:length(CellID)){
    Nloc<- NewDauLoc(boxXY[t,2])              # the x and y of each cell
    boxXY[t,3]<-Nloc[1]
    boxXY[t,4]<-Nloc[2]
  }

  X=boxXY[,3]
  Y=boxXY[,4]
  black<-VesselsXY[,3]

  blackTable<-data.frame()
  blackTableX<-c()
  blackTableY<-c()
  for (i in black){
    df<- NewDauLoc(i)              # the x and y of each lumin (black boxes)
    blackTableX<-c(blackTableX,df[1])
    blackTableY<-c(blackTableY,df[2])
  }
  blackTable[1:BloVesN,1]<-blackTableX
  blackTable[1:BloVesN,2]<-blackTableY

  xleft<-blackTableX-(10 + (Rrows * 20))
  ybottom<-blackTableY-(10 + (Rrows * 20))
  xright<-blackTableX+(10 + (Rrows * 20))
  ytop<-blackTableY+(10 + (Rrows * 20))

  xleftb<-blackTableX- 10
  ybottomb<-blackTableY-10
  xrightb<-blackTableX+10
  ytopb<-blackTableY+10


  jpeg(paste0("Track plot-",s,"PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
  plot (0,0, type = "n",xlim=c(0,4000),ylim=c(0,4000),xlab="X",ylab="Y",main=paste0(s," Steps   ",length(X), "cells","  Velecity=",V,"   PR=",PER))

  for (i in 1:BloVesN){
    rect(xleft, ybottom, xright, ytop,col=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3),border=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3))
  }
  for (i in 1:BloVesN){
    rect(xleftb, ybottomb, xrightb, ytopb,col="black",border=NULL)
  }
  points(X,Y, col=color[boxXY[,1]],pch=20,cex=0.3)

  dev.off()

}

PlottingCellsAD<- function(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows){

  BloVesN <- nrow(VesselsXY)

  len1200<-mainDF[SimTime,]         ########## computing the ID of the last cell
  FinalnumCell<-length(len1200[!is.na(len1200)])
  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
  totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
  color<-c()
  nickColors <- function(n, h = c(120,480), l = c(.60,.70), s = c(.8,1), alpha = 1){
    return (alpha(hex(HLS(seq(h[1],h[2],length.out = n), seq(l[1],l[2],length.out = n), seq(s[1],s[2],length.out=n))), alpha))
  }
  n=totalnumCell
  color = nickColors(n)

  for (s in 1:SimTime){
    len1200<-mainDF[s,]         ########## computing the ID of the cell ID in each step
    FinalnumCell<-length(len1200[!is.na(len1200)])
    totalnumCell<-len1200[!is.na(len1200)]
    nameDF<-rbind(len1200,names(len1200))
    rownames(nameDF)<-NULL
    colnames(nameDF)<-nameDF[2,]
    CellID<-as.numeric(nameDF[2,][!is.na(nameDF[1,])])

    lenID<-mainDF[s,]
    x=c(rep(NA,length(lenID[!is.na(lenID)])))
    y=c(rep(NA,length(lenID[!is.na(lenID)])))

    X=c()
    Y=c()
    boxXY<-data.frame(CellID,lenID[!is.na(lenID)],x,y)
    for (t in 1:length(CellID)){
      Nloc<- NewDauLoc(boxXY[t,2])              # the x and y of each cell
      boxXY[t,3]<-Nloc[1]
      boxXY[t,4]<-Nloc[2]
    }

    X=boxXY[,3]
    Y=boxXY[,4]
    black<-VesselsXY[,3]


    blackTable<-data.frame()
    blackTableX<-c()
    blackTableY<-c()
    for (i in black){
      df<- NewDauLoc(i)              # the x and y of each lumin (black boxes)
      blackTableX<-c(blackTableX,df[1])
      blackTableY<-c(blackTableY,df[2])
    }
    blackTable[1:BloVesN,1]<-blackTableX
    blackTable[1:BloVesN,2]<-blackTableY

    xleft<-blackTableX-(10 + (Rrows * 20))
    ybottom<-blackTableY-(10 + (Rrows * 20))
    xright<-blackTableX+(10 + (Rrows * 20))
    ytop<-blackTableY+(10 + (Rrows * 20))


    xleftb<-blackTableX- 10
    ybottomb<-blackTableY-10
    xrightb<-blackTableX+10
    ytopb<-blackTableY+10

    jpeg(paste0("Track plot-",s,"PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
    plot (0,0, type = "n",xlim=c(0,4000),ylim=c(0,4000),xlab="X",ylab="Y",main=paste0(s," Steps   ",length(X), "cells","  Velecity=",V,"   PR=",PER))

    for (i in 1:BloVesN){
      rect(xleft, ybottom, xright, ytop,col=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3),border=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3))
    }
    for (i in 1:BloVesN){
      rect(xleftb, ybottomb, xrightb, ytopb,col="black",border=NULL)
    }
    points(X,Y, col=color[boxXY[,1]],pch=20,cex=0.3)

    dev.off()
  }
}

SimulationScript<- function(TimeInterval,DT,DTinR,VesselsXY,PER,LatDF,
                            PeriodForSenescence,PeriodForDeath,LeavingRzone,IntravasationProb,
                            mainDF,cdDF,phylogenicDF,CellTables,speedTable,MetastatisTable,
                            V,SimTime,exportSP=FALSE,exportRdata=FALSE,Rrows){
  AA<-LatDF[,2]
  LLC1<-as.numeric(LatDF[,1])         #for the first cell
  LatC1=LLC1
  LatticeTableC1<-c()
  for (i in 1:DT){
    LatticeTableC1<-rbind(LatticeTableC1,LatC1)
  }
  dim(LatticeTableC1)

  ####################################
  ###########  First cell  ###########
  ####################################

  #¤¤¤¤¤¤¤¤¤¤¤¤¤ First step
  CellTables[[1]][1,1]=1
  CellTables[[1]][1,2]=1
  Rx= runif(1,100,3900)
  Ry= runif(1,100,3900)
  B=(floor(Ry/20)*200)+ceiling(Rx/20)
  #############a while loop to avoid having the first cell starting from a lumin
  if (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
    while (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
      Rx= floor(runif(1,100,3900))
      Ry= floor(runif(1,100,3900))
      B=(floor(Ry/20)*200)+ceiling(Rx/20)
    }
  }

  CellTables[[1]][1,3]=Rx
  CellTables[[1]][1,4]=Ry
  CellTables[[1]][1,5]= (floor(Ry/20)*200)+ceiling(Rx/20)
  CellTables[[1]][1,6]=2
  if (AA[CellTables[[1]][1,5]]==0){
    CellTables[[1]][1,7]=0
  }else{
    CellTables[[1]][1,7]=1
  }
  CellTables[[1]][1,8]=1
  if (CellTables[[1]][1,7]==0){
    CellTables[[1]][1,9]=0
  }else{
    CellTables[[1]][1,9]=DT-DTinR -1
  }
  CellTables[[1]][1,10]=sum(CellTables[[1]][1,6],CellTables[[1]][1,7],CellTables[[1]][1,8],CellTables[[1]][1,9])
  mainDF[1,1]=CellTables[[1]][1,5]

  mainDF[1:5,1:5]
  CellTables[[1]][1:5,]
  cdDF[1,1]=1



  #¤¤¤¤¤¤¤¤¤¤¤¤¤ Second step ############################
  CellTables[[1]][2,1]=1
  CellTables[[1]][2,2]=2
  x=CellTables[[1]][1,3]
  y=CellTables[[1]][1,4]
  ran=round(runif(1,0,1),digits=2)
  dist=CellTables[[1]][2,13] * TimeInterval
  NX<-x + (cos(ran * 2* pi) *dist)
  NY<-y + (sin(ran * 2* pi) *dist)

  mar=marginXY(NX,NY)
  NX=mar[1]
  NY=mar[2]

  if (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
    while (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
      x=CellTables[[1]][1,3]
      y=CellTables[[1]][1,4]
      ran=round(runif(1,0,1),digits=2)
      dist=CellTables[[1]][2,13] * TimeInterval
      NX<-x + (cos(ran * 2* pi) *dist)
      NY<-y + (sin(ran * 2* pi) *dist)
      mar=marginXY(NX,NY)
      NX=mar[1]
      NY=mar[2]
    }
  }

  CellTables[[1]][2,3]= NX
  CellTables[[1]][2,4]= NY
  CellTables[[1]][2,5]= (floor(Ry/20)*200)+ceiling(Rx/20)
  CellTables[[1]][2,6]=2
  if (AA[CellTables[[1]][2,5]]==0){
    CellTables[[1]][2,7]=0
  }else{
    CellTables[[1]][2,7]=1
  }
  CellTables[[1]][2,8]=2
  if (CellTables[[1]][2,7]==0){
    CellTables[[1]][2,9]=0
  }else{
    CellTables[[1]][2,9]=DT-DTinR -1
  }
  CellTables[[1]][2,10]=sum(CellTables[[1]][2,6],CellTables[[1]][2,7],CellTables[[1]][2,8],CellTables[[1]][2,9])
  CellTables[[1]][2,14]= round(sqrt((CellTables[[1]][2,3]- CellTables[[1]][1,3])^2 + (CellTables[[1]][2,4]- CellTables[[1]][1,4])^2),digits=2)

  mainDF[2,1]=CellTables[[1]][2,5]
  cdDF[2,1]=1
  G1Division<-c()


  #¤¤¤¤¤¤¤¤¤¤¤¤¤remaining steps
  for (s in 3:DT){
    Lat=LatticeTableC1[s,]
    CellTables[[1]][s,1]=1
    CellTables[[1]][s,2]=s

    x1=CellTables[[1]][s-2,3]
    y1=CellTables[[1]][s-2,4]

    x2=CellTables[[1]][s-1,3]
    y2=CellTables[[1]][s-1,4]

    speed1=CellTables[[1]][s-1,14]
    if (speed1==0){speed1=1.4}
    speed2=CellTables[[1]][s,13] * TimeInterval

    B=CellTables[[1]][s-1,5]
    Lat[B]=1
    ran=runif(1,0,1)

    if (ran<=(1-PER)){
      BM<-BackMovement(x1,y1,x2,y2,speed1,speed2,B,Lat)
      if (length(BM[[1]])>0){
        BMran<-ceiling(runif(1,0,length(BM[[1]])))
        BB<-BM[[1]][BMran]
        NX<-BM[[2]][BMran]
        NY<-BM[[3]][BMran]
      }else{
        BB<-B
        NX<-CellTables[[1]][s-1,3]
        NY<-CellTables[[1]][s-1,4]
      }
    }

    if (ran>(1-PER)){
      FM<-ForMovement(x1,y1,x2,y2,speed1,speed2,B,Lat)
      if (length(FM[[1]])>0){
        FMran<-ceiling(runif(1,0,length(FM[[1]])))
        BB<-FM[[1]][FMran]
        NX<-FM[[2]][FMran]
        NY<-FM[[3]][FMran]
      }else{
        BB<-B
        NX<-CellTables[[1]][s-1,3]
        NY<-CellTables[[1]][s-1,4]
      }
    }
    CellTables[[1]][s,3]= NX
    CellTables[[1]][s,4]= NY
    CellTables[[1]][s,5]= BB
    CellTables[[1]][s,6]=2
    if (AA[CellTables[[1]][s,5]]==0){
      CellTables[[1]][s,7]=0
    }else{
      CellTables[[1]][s,7]=1
    }

    ranRed<-runif(1,0,1)                                                     ################# Reducing the chance that a cell will leave the red zone
    if (CellTables[[1]][s,7]==0 && CellTables[[1]][s-1,7]==1 && ranRed<=LeavingRzone){
      CellTables[[1]][s,3]= CellTables[[1]][s-1,3]
      CellTables[[1]][s,4]= CellTables[[1]][s-1,4]
      CellTables[[1]][s,5]= CellTables[[1]][s-1,5]
      CellTables[[1]][s,6]= CellTables[[1]][s-1,6]
      CellTables[[1]][s,7]= CellTables[[1]][s-1,7]
    }

    CellTables[[1]][s,8]=s
    if (CellTables[[1]][s,7]==0){
      CellTables[[1]][s,9]=0
    }else{
      CellTables[[1]][s,9]= DT-DTinR - 1
    }

    if (sum(CellTables[[1]][2:s,7])>=G1 && CellTables[[1]][s,7]==0  &&  CellTables[[1]][s,8]>=DTinR  &&  s <DT){
      WWW=which(CellTables[[1]][,7]==1)
      stepID<-CellTables[[1]][,8]
      stepID=stepID[!is.na(stepID)]
      if (length(WWW)>0){
        ones<-stepID[WWW]
        for (i in ones){
          if(i+((G1)-1)>s){
            break
          }else{
            if (sum(CellTables[[1]][i:(i+((G1)-1)),7])==G1){
              CellTables[[1]][s,9]=DT-DTinR
              G1Division<-c(G1Division,1)
              break
            }
          }
        }
      }
    }

    CellTables[[1]][s,10]=sum(CellTables[[1]][s,6],CellTables[[1]][s,7],CellTables[[1]][s,8],CellTables[[1]][s,9])
    CellTables[[1]][s,11]=round(ran,digits=2)
    CellTables[[1]][s,12]=round(ranRed,digits=2)
    CellTables[[1]][s,14]= round(sqrt((CellTables[[1]][s,3]- CellTables[[1]][s-1,3])^2 + (CellTables[[1]][s,4]- CellTables[[1]][s-1,4])^2),digits=2)

    if (ran<=(1-PER)){
      CellTables[[1]][s,15]=length(BM[[1]])
    }else{
      CellTables[[1]][s,15]=length(FM[[1]])
    }

    PRAN<-runif(1,0,1)
    if (CellTables[[1]][s,10]>= DT & PRAN<=0.999){
      CellTables[[1]][s,16]=2
    }
    else if (CellTables[[1]][s,10]< DT && CellTables[[1]][s,10]> G1 && PRAN<=0.001){
      CellTables[[1]][s,16]=2
    }else{
      CellTables[[1]][s,16]=1
    }
    CellTables[[1]][s,17]=round(PRAN,digits=4)
    cdDF[s,1]=CellTables[[1]][s,16]
    mainDF[s,1]=CellTables[[1]][s,5]
    if (CellTables[[1]][s,16]==2){break}
  }
  phylogenicDF[1,1]=1


  ##################################################################
  ####################################### All other cells ##########
  ##################################################################

  target<-CellTables[[1]][,5]
  LenTar<-length(target[!is.na(target)])    # step of the first cell division
  RCellNumPerStep<-c()
  RCellNumPerStep=CellTables[[1]][1:LenTar,7]
  RCellNumPerStep<-c(RCellNumPerStep,rep(NA,(SimTime-LenTar)))

  CellDeath<-c()
  CellDeath<-c(rep(0,LenTar))
  CellDeath<-c(CellDeath,rep(NA,(SimTime-LenTar)))

  CellSenescence<-c()
  CellSenescence<-c(rep(0,LenTar))

  MetastaticCells<-c()
  cat("The first cell devided after ",LenTar, " steps","\n")

  ###########################################  Generating the LatticeTable ##################
  num=round(SimTime/5)
  black<-VesselsXY[,3]
  LL<-as.numeric(LatDF[,1])
  LL[black]=0
  Lat=LL
  LatticeTable<-c()
  for (i in 1:num){
    LatticeTable<-rbind(LatticeTable,Lat)
  }
  LatticeTable2<-LatticeTable
  for (i in 1:4){
    LatticeTable2<-rbind(LatticeTable2,LatticeTable)
  }
  LatticeTable<-rbind(LatticeTableC1[1:(LenTar-1),],LatticeTable2)

  for (step in (LenTar+1):SimTime){
    print(step)
    Lat=LatticeTable[step-1,]
    Lat1<-Lat
    TargetCDrow<-cdDF[step-1,]
    NetTargetCDrow<-which(!is.na(TargetCDrow))
    CellsPerStep<-as.numeric(names(TargetCDrow[NetTargetCDrow]))                        # to know the ID of the cell
    CellState<-cdDF[step-1,]
    CellNum<-length(CellsPerStep)
    HowManyCells<-ceiling(sum(CellState,na.rm = T))      ################## to convert the 0.9999 of a metastatic cell into 1
    MAX<-CellsPerStep[which.max(CellsPerStep)]
    MIN<-CellsPerStep[which.min(CellsPerStep)]

    NewCells<-MAX + (HowManyCells-CellNum)*2                ##### to know how many new cells are there
    NewCellsID<-c(MAX :NewCells)
    if (length(NewCellsID) >1){
      NewCellsID<-NewCellsID[-1]
      cellIDs<-c(CellState[which(CellState!=2)],NewCellsID)
      gone<-as.numeric(names(CellState[which(CellState==2)]))
      go<-c()
      ifold<-c()
      NewOld<-c()                                      # to save the cells that will divide but they are not able anymore to divide
      for(g in gone){go=c(go,rep(g,2))}

      for(g in gone){
        doneSteps<-c()
        doneSteps<-CellTables[[g]][,2]
        doneSteps<-doneSteps[!is.na(doneSteps)]
        doneSteps<-c(doneSteps,step)
        Ma<-match(step,doneSteps)                    # to match the simulation steps with the cell steps

        Dens<-Density(CellTables[[g]][Ma-1,5],Lat1)        # checking the density
        if (Dens[[1]]==-2  || (Dens[[1]]==-1 && Dens[[2]]  %in% black) ){
          NewOld<-c(NewOld,g)
          go<- go[!go %in% g]
          next
        }else{
          options<-which(go==g)      # for every g there are 2 options    options[[1]] and options[[2]]

          ######################## first daughter cell
          mainDF[step,MAX+ options[1]]= CellTables[[g]][Ma-1,5]
          CellTables[[MAX+ options[1]]][1,1]=MAX+ options[1]
          CellTables[[MAX+ options[1]]][1,2]=step-1
          CellTables[[MAX+ options[1]]][1,3]= CellTables[[g]][Ma-1,3]
          CellTables[[MAX+ options[1]]][1,4]= CellTables[[g]][Ma-1,4]

          CellTables[[MAX+ options[1]]][2,1]=MAX+ options[1]
          CellTables[[MAX+ options[1]]][2,2]=step
          Denran<-ceiling(runif(1,0,length(Dens[[2]])))
          Nloc1<-NewDauLoc(CellTables[[g]][Ma-1,5])                       # the x and y of the first daughter cell is
          CellTables[[MAX+ options[1]]][2,3]= Nloc1[1]
          CellTables[[MAX+ options[1]]][2,4]= Nloc1[2]
          CellTables[[MAX+ options[1]]][2,5]= CellTables[[g]][Ma-1,5]
          CellTables[[MAX+ options[1]]][2,6]= Dens[[1]]
          CellTables[[MAX+ options[1]]][2,7]= CellTables[[g]][Ma-1,7]
          CellTables[[MAX+ options[1]]][2,9]= CellTables[[g]][Ma-1,9]
          CellTables[[MAX+ options[1]]][2,10]= sum(Dens[[1]],CellTables[[g]][Ma-1,7],CellTables[[MAX+ options[1]]][2,8],CellTables[[g]][Ma-1,9])
          CellTables[[MAX+ options[1]]][2,11]= NA
          CellTables[[MAX+ options[1]]][2,12]= NA
          CellTables[[MAX+ options[1]]][2,13]= NA
          CellTables[[MAX+ options[1]]][2,14]= round(sqrt((CellTables[[MAX+ options[1]]][2,3]- CellTables[[g]][Ma-1,3])^2 + (CellTables[[MAX+ options[1]]][2,4]- CellTables[[g]][Ma-1,4])^2),digits=2)
          CellTables[[MAX+ options[1]]][2,15]= NA

          CellTables[[MAX+ options[1]]][2,16]=1
          cdDF[step,MAX+ options[1]]=CellTables[[MAX+ options[1]]][2,16]
          phylog<-phylogenicDF[g,]
          len.phylogenic<- length(phylog[!is.na(phylog)])
          num.parents<- len.phylogenic +1
          phylogenicDF[MAX+ options[1],1:num.parents]<-c(phylog[!is.na(phylog)],MAX+ options[1])

          ####################### second daughter cell
          mainDF[step,MAX+ options[2]]= Dens[[2]][Denran]
          CellTables[[MAX+ options[2]]][1,1]=MAX+ options[2]
          CellTables[[MAX+ options[2]]][1,2]=step-1
          CellTables[[MAX+ options[2]]][1,3]= CellTables[[g]][Ma-1,3]
          CellTables[[MAX+ options[2]]][1,4]= CellTables[[g]][Ma-1,4]

          CellTables[[MAX+ options[2]]][2,1]=MAX+ options[2]
          CellTables[[MAX+ options[2]]][2,2]=step
          Nloc2<- NewDauLoc(Dens[[2]][Denran])              # the x and y of the second daughter cell is
          CellTables[[MAX+ options[2]]][2,3]= Nloc2[1]
          CellTables[[MAX+ options[2]]][2,4]= Nloc2[2]

          BSecondDC=Dens[[2]][Denran]  ############# to get the box of the second daughter cell
          RanM<-runif(1,0,1)
          if (BSecondDC %in% black  && RanM>IntravasationProb){
            otherDenran<-c(1:length(Dens[[2]]))
            if (length(otherDenran)>1){
              otherDenran<-otherDenran[-Denran]
              Denran1<-sample(otherDenran, 1)
              CellTables[[MAX+ options[2]]][2,5]=Dens[[2]][Denran1]
              Lat1[Dens[[2]][Denran1]]=1
              Dens=Density(Dens[[2]][Denran1],Lat1)
              CellTables[[MAX+ options[2]]][2,6]= Dens[[1]]
              if (AA[CellTables[[MAX+ options[2]]][2,5]]==0){CellTables[[MAX+ options[2]]][2,7]=0 }else{ CellTables[[MAX+ options[2]]][2,7]=1}
              if (CellTables[[MAX+ options[2]]][2,7]==0){ CellTables[[MAX+ options[2]]][2,9]=0}else{CellTables[[MAX+ options[2]]][2,9]=DT-DTinR -1}
              CellTables[[MAX+ options[2]]][2,10]= sum(Dens[[1]],CellTables[[MAX+ options[2]]][2,7],CellTables[[MAX+ options[2]]][2,8],CellTables[[MAX+ options[2]]][2,9])
              CellTables[[MAX+ options[2]]][2,11]= NA
              CellTables[[MAX+ options[2]]][2,12]= NA
              CellTables[[MAX+ options[2]]][2,13]= NA
              CellTables[[MAX+ options[2]]][2,14]= round(sqrt((CellTables[[MAX+ options[2]]][2,3]- CellTables[[g]][Ma-1,3])^2 + (CellTables[[MAX+ options[2]]][2,4]- CellTables[[g]][Ma-1,4])^2),digits=2)
              CellTables[[MAX+ options[2]]][2,15]= NA
              CellTables[[MAX+ options[2]]][2,16]= 1
            }

          }
          if (!BSecondDC %in% black) {
            CellTables[[MAX+ options[2]]][2,5]=Dens[[2]][Denran]
            Lat1[Dens[[2]][Denran]]=1
            Dens=Density(Dens[[2]][Denran],Lat1)
            CellTables[[MAX+ options[2]]][2,6]= Dens[[1]]
            if (AA[CellTables[[MAX+ options[2]]][2,5]]==0){CellTables[[MAX+ options[2]]][2,7]=0 }else{ CellTables[[MAX+ options[2]]][2,7]=1}
            if (CellTables[[MAX+ options[2]]][2,7]==0){ CellTables[[MAX+ options[2]]][2,9]=0}else{CellTables[[MAX+ options[2]]][2,9]=DT-DTinR -1}
            CellTables[[MAX+ options[2]]][2,10]= sum(Dens[[1]],CellTables[[MAX+ options[2]]][2,7],CellTables[[MAX+ options[2]]][2,8],CellTables[[MAX+ options[2]]][2,9])
            CellTables[[MAX+ options[2]]][2,11]= NA
            CellTables[[MAX+ options[2]]][2,12]= NA
            CellTables[[MAX+ options[2]]][2,13]= NA
            CellTables[[MAX+ options[2]]][2,14]= round(sqrt((CellTables[[MAX+ options[2]]][2,3]- CellTables[[g]][Ma-1,3])^2 + (CellTables[[MAX+ options[2]]][2,4]- CellTables[[g]][Ma-1,4])^2),digits=2)
            CellTables[[MAX+ options[2]]][2,15]= NA
            CellTables[[MAX+ options[2]]][2,16]= 1
          }

          if (BSecondDC %in% black  && RanM<=IntravasationProb){
            MetastaticCells<-c(MetastaticCells,MAX+ options[2])
            CellTables[[MAX+ options[2]]][2,5]=NA
            CellTables[[MAX+ options[2]]][2,16]=0.9999
            cdDF[step:SimTime,MAX+ options[2]]=0.9999

          }

          cdDF[step,MAX+ options[2]]=CellTables[[MAX+ options[2]]][2,16]
          phylog<-phylogenicDF[g,]
          len.phylogenic<- length(phylog[!is.na(phylog)])
          num.parents<- len.phylogenic +1
          phylogenicDF[MAX+ options[2],1:num.parents]<-c(phylog[!is.na(phylog)],MAX+ options[2])
        }
      }
    }

    ifold<- as.numeric(names(CellState[which(CellState==1)]))                         # to know if we have undivided cells

    ifold<-c(ifold,NewOld)
    if (length(ifold) >0){
      for (o in ifold){
        doneSteps<-c()
        doneSteps<-CellTables[[o]][,2]
        doneSteps<-doneSteps[!is.na(doneSteps)]
        doneSteps<-c(doneSteps,step)
        Ma<-match(step,doneSteps)                                            # to match the simulation steps with the cell steps

        TEST01 <- tryCatch(
          {!is.null(CellTables[[o]][Ma,8]) && length(CellTables[[o]][Ma,8]) > 0  && CellTables[[o]][Ma,8]>= (DT *PeriodForDeath)},
          error = function(e) {FALSE})
        if(TEST01){next}           #### to remove the cell that has spent 431 steps without dividing

        CellTables[[o]][Ma,1]=o
        CellTables[[o]][Ma,2]=step

        x1=CellTables[[o]][Ma-2,3]
        y1=CellTables[[o]][Ma-2,4]
        x2=CellTables[[o]][Ma-1,3]
        y2=CellTables[[o]][Ma-1,4]
        speed1=CellTables[[o]][Ma-1,14]
        if (speed1==0){speed1=1.4}
        speed2=CellTables[[o]][Ma,13] * TimeInterval
        B=CellTables[[o]][Ma-1,5]
        ran=runif(1,0,1)

        if (ran<=(1-PER)){                                               #### if it is smaller or =
          BM<-BackMovement(x1,y1,x2,y2,speed1,speed2,B,Lat1)
          if (length(BM[[1]])>0){
            BMran<-ceiling(runif(1,0,length(BM[[1]])))
            BB<-BM[[1]][BMran]
            NX<-BM[[2]][BMran]
            NY<-BM[[3]][BMran]
          }else{
            BB<-B
            NX<-CellTables[[o]][Ma-1,3]
            NY<-CellTables[[o]][Ma-1,4]
          }
        }

        if (ran>(1-PER)){
          FM<-ForMovement(x1,y1,x2,y2,speed1,speed2,B,Lat1)
          if (length(FM[[1]])>0){
            FMran<-ceiling(runif(1,0,length(FM[[1]])))
            BB<-FM[[1]][FMran]
            NX<-FM[[2]][FMran]
            NY<-FM[[3]][FMran]
          }else{
            BB<-B
            NX<-CellTables[[o]][Ma-1,3]
            NY<-CellTables[[o]][Ma-1,4]
          }
        }
        ############
        RanM<-runif(1,0,1)                       ##########The cell will stay in place
        if (BB %in% black  && RanM>IntravasationProb){
          BB<-B
          NX<-CellTables[[o]][Ma-1,3]
          NY<-CellTables[[o]][Ma-1,4]
          CellTables[[o]][Ma,3]= NX
          CellTables[[o]][Ma,4]= NY
          CellTables[[o]][Ma,5]= BB
          Lat1[B]=0
          Lat1[BB]=1

          DenS=Density(BB,Lat1)
          CellTables[[o]][Ma,6]= DenS[[1]]
          if (AA[CellTables[[o]][Ma,5]]==0){
            CellTables[[o]][Ma,7]=0
          }else{
            CellTables[[o]][Ma,7]=1
          }

          ranRed<-runif(1,0,1)                                                     ################# Reducing the chance that a cell will leave the red zone
          if (CellTables[[o]][Ma,7]==0 && CellTables[[o]][Ma-1,7]==1 && ranRed<=LeavingRzone){
            CellTables[[o]][Ma,3]= CellTables[[o]][Ma-1,3]
            CellTables[[o]][Ma,4]= CellTables[[o]][Ma-1,4]
            CellTables[[o]][Ma,5]= CellTables[[o]][Ma-1,5]
            Lat1[B]=1
            Lat1[BB]=0
            CellTables[[o]][Ma,6]= CellTables[[o]][Ma-1,6]
            CellTables[[o]][Ma,7]= CellTables[[o]][Ma-1,7]
          }

          if (CellTables[[o]][Ma,7]==0){
            CellTables[[o]][Ma,9]=0
          }else{
            CellTables[[o]][Ma,9]= DT-DTinR - 1
          }

          if (sum(CellTables[[o]][2:Ma,7])>=G1 && CellTables[[o]][Ma,7]==0  &&  CellTables[[o]][Ma,8]>=DTinR  &&  Ma <DT){
            WWW=which(CellTables[[o]][,7]==1)
            stepID<-CellTables[[o]][,8]
            stepID=stepID[!is.na(stepID)]
            if (length(WWW)>0){
              ones<-stepID[WWW]
              for (i in ones){
                if(i+((G1)-1)>Ma){
                  break
                }else{
                  if (sum(CellTables[[o]][i:(i+((G1)-1)),7])==G1){
                    CellTables[[o]][Ma,9]=DT-DTinR
                    G1Division<-c(G1Division,o)
                    break
                  }
                }
              }
            }
          }

          CellTables[[o]][Ma,10]=sum(CellTables[[o]][Ma,6],CellTables[[o]][Ma,7],CellTables[[o]][Ma,8],CellTables[[o]][Ma,9])
          CellTables[[o]][Ma,11]=round(ran,digits=2)
          CellTables[[o]][Ma,12]=round(ranRed,digits=2)
          CellTables[[o]][Ma,14]= round(sqrt((CellTables[[o]][Ma,3]- CellTables[[o]][Ma-1,3])^2 + (CellTables[[o]][Ma,4]- CellTables[[o]][Ma-1,4])^2),digits=2)

          if (ran<=(1-PER)){
            CellTables[[o]][Ma,15]=length(BM[[1]])
          }else{
            CellTables[[o]][Ma,15]=length(FM[[1]])
          }

          PRAN<-runif(1,0,1)
          if (CellTables[[o]][Ma,10]>= DT && CellTables[[o]][Ma,6]!=-2 && PRAN<=0.999){
            CellTables[[o]][Ma,16]=2
          }
          else if (CellTables[[o]][Ma,10]< DT && CellTables[[o]][Ma,10]> G1 &&  PRAN<=0.001){
            CellTables[[o]][Ma,16]=2
          }else{
            CellTables[[o]][Ma,16]=1
          }

          CellTables[[o]][Ma,17]=round(PRAN,digits=4)
        }
        ####################### The cell will move

        if (!BB %in% black){
          CellTables[[o]][Ma,3]= NX
          CellTables[[o]][Ma,4]= NY
          CellTables[[o]][Ma,5]= BB
          Lat1[B]=0
          Lat1[BB]=1

          DenS=Density(BB,Lat1)
          CellTables[[o]][Ma,6]= DenS[[1]]
          if (AA[CellTables[[o]][Ma,5]]==0){
            CellTables[[o]][Ma,7]=0
          }else{
            CellTables[[o]][Ma,7]=1
          }

          ranRed<-runif(1,0,1)                                                     ################# Reducing the chance that a cell will leave the red zone
          if (CellTables[[o]][Ma,7]==0 && CellTables[[o]][Ma-1,7]==1 && ranRed<=LeavingRzone){
            CellTables[[o]][Ma,3]= CellTables[[o]][Ma-1,3]
            CellTables[[o]][Ma,4]= CellTables[[o]][Ma-1,4]
            CellTables[[o]][Ma,5]= CellTables[[o]][Ma-1,5]
            Lat1[B]=1
            Lat1[BB]=0
            CellTables[[o]][Ma,6]= CellTables[[o]][Ma-1,6]
            CellTables[[o]][Ma,7]= CellTables[[o]][Ma-1,7]
          }

          if (CellTables[[o]][Ma,7]==0){
            CellTables[[o]][Ma,9]=0
          }else{
            CellTables[[o]][Ma,9]= DT-DTinR - 1
          }

          if (sum(CellTables[[o]][2:Ma,7])>=G1 && CellTables[[o]][Ma,7]==0  &&  CellTables[[o]][Ma,8]>=DTinR  &&  Ma <DT){
            WWW=which(CellTables[[o]][,7]==1)
            stepID<-CellTables[[o]][,8]
            stepID=stepID[!is.na(stepID)]
            if (length(WWW)>0){
              ones<-stepID[WWW]
              for (i in ones){
                if(i+((G1)-1)>Ma){
                  break
                }else{
                  if (sum(CellTables[[o]][i:(i+((G1)-1)),7])==G1){
                    CellTables[[o]][Ma,9]=DT-DTinR
                    G1Division<-c(G1Division,o)
                    break
                  }
                }
              }
            }
          }

          CellTables[[o]][Ma,10]=sum(CellTables[[o]][Ma,6],CellTables[[o]][Ma,7],CellTables[[o]][Ma,8],CellTables[[o]][Ma,9])
          CellTables[[o]][Ma,11]=round(ran,digits=2)
          CellTables[[o]][Ma,12]=round(ranRed,digits=2)
          CellTables[[o]][Ma,14]= round(sqrt((CellTables[[o]][Ma,3]- CellTables[[o]][Ma-1,3])^2 + (CellTables[[o]][Ma,4]- CellTables[[o]][Ma-1,4])^2),digits=2)

          if (ran<=(1-PER)){
            CellTables[[o]][Ma,15]=length(BM[[1]])
          }else{
            CellTables[[o]][Ma,15]=length(FM[[1]])
          }

          PRAN<-runif(1,0,1)
          if (CellTables[[o]][Ma,10]>= DT && CellTables[[o]][Ma,6]!=-2 && PRAN<=0.999){
            CellTables[[o]][Ma,16]=2
          }
          else if (CellTables[[o]][Ma,10]< DT && CellTables[[o]][Ma,10]> G1 &&  PRAN<=0.001){
            CellTables[[o]][Ma,16]=2
          }else{
            CellTables[[o]][Ma,16]=1
          }

          CellTables[[o]][Ma,17]=round(PRAN,digits=4)

        }
        ############################ The cell will go away
        if (BB %in% black  && RanM<=IntravasationProb){
          MetastaticCells<-c(MetastaticCells,o)
          CellTables[[o]][Ma,5]=NA
          CellTables[[o]][Ma,16]=0.9999
          cdDF[step:SimTime,o]=0.9999
          CellTables[[o]][Ma,10]=CellTables[[o]][Ma-1,10]

        }

        ########################
        cdDF[step,o]=CellTables[[o]][Ma,16]
        mainDF[step,o]=CellTables[[o]][Ma,5]
      }

    }
    LatticeTable[step,]<-Lat1
    CelDet<-c()                                        ############ computing cell death per step
    for (i in 1:length(CellTables)){
      NumStep<-CellTables[[i]][,16]
      lenNumStep<-length(NumStep[!is.na(NumStep)])
      if (lenNumStep>=(PeriodForDeath * DT)){CelDet[i]=1}else{CelDet[i]=NA}
    }
    CellDeath[step]<-length(CelDet[!is.na(CelDet)])

    CelSen<-c()                                ############ computing cell  senescence per step
    for (i in 1:length(CellTables)){
      NumStep<-CellTables[[i]][,16]
      lenNumStep<-length(NumStep[!is.na(NumStep)])
      if (lenNumStep>(DT * PeriodForSenescence) && lenNumStep<(PeriodForDeath * DT)){CelSen[i]=1}else{CelSen[i]=NA}
    }

    CellSenescence[step]<-length(CelSen[!is.na(CelSen)])

    main<-mainDF[step,]                                  ############### computing total hours in red zone per step
    m1200<-main[!is.na(main)]
    NRm1200<-c()
    for (n in m1200){
      if (AA[n]==0){
        NRm1200<-c(NRm1200,n)
      }
    }
    RCellNumPerStep[step]<- length(m1200)- length(NRm1200)    # all - not red
    if (length(MetastaticCells)>0){
      MetastatisTable[step,1:length(MetastaticCells)]=MetastaticCells      ### filling the MetastatisTable
    }
  }

  cat("Number of Metastatic Cells = ",length(MetastaticCells), "\n")

  ######################### To get the maximum cellID in the end of the steps
  len1200<-cdDF[SimTime,]
  FinalnumCell<-length(len1200[!is.na(len1200)])
  cat("Number of Cells after ", SimTime, "hours = ",FinalnumCell, "\n")

  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
  totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
  cat("Total Number of Cells generated after ", SimTime, "hours = ",totalnumCell, "\n")
  #################Saving tables
  RanNum<-round(runif(1,100,10000),digits=2)
  ################Saving plots
  TS=SimTime

  RCellNum=RCellNumPerStep

  divisions<-c(1:20)
  CellNumberPerDivision<-c()
  for(i in divisions){
    len17<-phylogenicDF[,i]
    len17<-len17[!is.na(len17)]
    noDub<-len17[!duplicated(len17)]
    CellNumberPerDivision[i]<-length(noDub)
  }
  write.csv(CellNumberPerDivision,file=paste0("CellNumberPerDivision","PR=",PER,"V=",V,"R",RanNum,".csv"))


  ########################## computing the average speed of all cells
  AveSpe<-c()
  for (i in 1:totalnumCell){
    SPEED<-CellTables[[i]][,14]
    SPEED<-SPEED[!is.na(SPEED)]
    if (length(SPEED)>100){
      for(s in SPEED){
        if (s>250){                            ##### to correct the Periodic boundary conditions effect on speed
          mat<-match(s,SPEED)
          SPEED[mat]=4000-s
        }
      }
      AveSpe[i]<-mean(SPEED)
    }else{
      AveSpe[i]<-NA
    }
  }
  FinalAveSpe<-round(mean(AveSpe,na.rm=T),digits=3)


  ########################## computing number of metastatic cells per step
  MCells<-c()
  for (i in 1:SimTime){
    LenMCells<-MetastatisTable[i,]
    LenMCells<-LenMCells[!is.na(LenMCells)]
    MCells[i]<-length(LenMCells)
  }



  ################################################## Computing the number of immotile cells and dividing cells last five steps
  allCells<-c()
  NonDividingCells<-c()
  NotMotile<-c()
  L5S<-c()
  S2155<-mainDF[SimTime-5,]
  S2155<-S2155[!is.na(S2155)]                    # the ID of the first cell in Step=2155
  TotalCells<-totalnumCell                       # total number of cells in the end of the simulation
  for (i in S2155[1]:TotalCells){
    L5S<-mainDF[(SimTime-5):SimTime,i]
    L5S<-L5S[!is.na(L5S)]
    if (length(L5S)>0){
      allCells<-c(allCells,length(L5S))            # to get the length of cell tracks last 5 steps
      if (length(L5S)==6){NonDividingCells<-c(NonDividingCells,1)}
      if (length(L5S)==6 && mainDF[SimTime,i]==mainDF[(SimTime-1),i] && mainDF[(SimTime-1),i]==mainDF[(SimTime-2),i] && mainDF[(SimTime-2),i]==mainDF[(SimTime-3),i]&& mainDF[(SimTime-3),i]==mainDF[(SimTime-4),i]&& mainDF[(SimTime-4),i]==mainDF[(SimTime-5),i]){
        NotMotile<-c(NotMotile,1)
      }
    }
  }

  NotMotileCells<-length(NotMotile)
  DividingCells=length(allCells) - length(NonDividingCells)


  simsumDF<-SimSum(CellTables,cdDF,mainDF,SimTime,AA,CellDeath,RCellNum,LatticeTable,CellSenescence,FinalAveSpe,MCells,V,PER,G1Division,DividingCells,NotMotileCells,RanNum)

  write.csv(phylogenicDF,file=paste0("phylogenicDF","PR=",PER,"Speed=",V,"R",RanNum,".csv"))
  #write.csv(simsumDF,file=paste0("simsumDF","PR=",PER,"Speed=",V,"R",RanNum,".csv"))
  jpeg(paste0("G1B Cells"," PR=",PER," Speed=",V," Number of G1B cells=",length(G1Division),"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
  G1B.Cells=c(1:length(G1Division))
  CellID=G1Division
  LensimsumDF<-length(simsumDF[,1])
  plot(G1B.Cells,CellID,type="p",ylim=c(1,simsumDF[LensimsumDF,5]),main=paste0("G1Bcells=",length(G1Division),", PR=",PER,", Speed=",V,", Tumor cells=",simsumDF[LensimsumDF,5]),cex.main=1)

  for (i in 1:length(simsumDF[,1])){
    abline(h=simsumDF[i,5],col="red")
  }
  dev.off()

  # This seems to raise an error!!!
  # check
  PlottingCellsADLastStep(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows)
  jpeg(paste0("CellNumberPerDivision-",s,"PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
  plot(divisions,CellNumberPerDivision,type="o",ylab="Number of cells undergone cell division",ylim=c(0,30000),
       xlab="Cell division number",pch=16,main=paste0("Cell Division Curve","  (Velecity=",V," ,  PR=",PER,")"))
  dev.off()
  # check
  UIplot(TS,LatticeTable,RanNum,V,PER)

  if (exportSP== TRUE) {
    #check
    PlottingCellsAD(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows)

  }
  if (exportRdata== TRUE) {
    ScdDF<-cdDF[,1:totalnumCell]
    SmainDF<-mainDF[,1:totalnumCell]
    save(ScdDF, file=paste0("cdDF-","PR_",PER,"velecity_",V,"R",RanNum,".RData"))
    save(SmainDF, file=paste0("mainDF-","PR_",PER,"velecity_",V,"R",RanNum,".RData"))

  }
  return(list(simsumDF))
}



# Wrapper function
# to Simulate()
RunCellmigSimulation <- function(
  TimeInterval,            # in hours
  DT,                      # in hours
  DTinR,                   # in hours
  G1,                      # in hours
  SimTime,                 # in hours

  minSpe,                  # um/h
  maxSpe,                  # um/h
  desired_meanSpe,         # um/h
  desired_sdSpe,           # um/h

  VesselsList,             # vessel data

  PER,
  PeriodForSenescence,
  PeriodForDeath,
  LeavingRzone,
  IntravasationProb,
  SetCellN = NULL,
  exportSP=FALSE,
  exportRdata=FALSE)
{

  # Prepare / import
  LatDF <- VesselsList$LatDF
  VesselsXY <- VesselsList$VesselsXY
  BloVesN <- VesselsList$vesselParams$BloVesN
  Rrows <- VesselsList$vesselParams$Rrows

  # define V
  V = desired_meanSpe

  # impute EstCellN
  divN= ceiling(SimTime/DT)+1
  estCellN=1

  for (i in 1:divN){
    estCellN=c(estCellN,2^i)
  }
  EstCellN<-round(sum(estCellN)+(10/PeriodForDeath) *SimTime)

  #bypass
  if (!is.null(SetCellN) && is.numeric(SetCellN)) {
    EstCellN <- SetCellN
  }

  if (EstCellN > 100000){
    EstCellN = 100000
  }

  # main df and other dfs
  mainDF<-data.frame(matrix(NA, nrow = SimTime, ncol = EstCellN))
  rownames(mainDF)<-c(1:SimTime)
  colnames(mainDF)<-c(1:EstCellN)

  speedTable<-data.frame()
  for (i in 1:1000){                  ########## generating 1000 speed profiles
    o <- optim(c(V, desired_sdSpe), eval_function)
    Speed <- rtruncnorm(n = DT, a = minSpe, b = maxSpe, mean = o$par[1], sd = o$par[2])
    SpeedT<-list()
    for (d in 1:PeriodForDeath){
      SpeedT[[d]]<- Speed
    }
    SpeedT <- do.call(c, SpeedT)
    LenSpeedT<-length(SpeedT)
    speedTable[1:LenSpeedT,i]= round(SpeedT,digits=2)
  }

  CellTables <- vector(mode = "list", length = EstCellN)
  for (n in 1:EstCellN){
    CellTables [[n]]<- data.frame(matrix(NA, nrow = (DT *PeriodForDeath), ncol = 17))

    ran=ceiling(runif(1,0,1000))            #########  selecting randomly the speed profile from the speedTable
    CellTables [[n]][2:LenSpeedT,13]= speedTable[1:(LenSpeedT-1),ran]
    CellTables [[n]][2:LenSpeedT,8]=c(1:(LenSpeedT-1))
    colnames(CellTables [[n]])<-c("ID","Steps","X","Y","B","Density","Red","CellStep","NotRed","Sum","RanDR","RanRed","Pre.Distance (um)","Actual Distance (um)","Movement options","CD")
  }
  CellTables[[1]]<-CellTables[[1]][-1,]
  CellTables[[1]][1,13]=NA

  cdDF<-data.frame(matrix(NA, nrow = SimTime, ncol = EstCellN))
  rownames(cdDF)<-c(1:SimTime)
  colnames(cdDF)<-c(1:EstCellN)

  phylogenicDF<-data.frame(matrix(NA, nrow = EstCellN, ncol = 25))
  MetastatisTable<- data.frame(matrix(NA, nrow = SimTime, ncol = (SimTime *3)))

  # Run
  y <- SimulationScript(
    TimeInterval,
    DT,
    DTinR,
    VesselsXY,
    PER,
    LatDF,
    PeriodForSenescence,
    PeriodForDeath,
    LeavingRzone,
    IntravasationProb,
    mainDF,
    cdDF,
    phylogenicDF,
    CellTables,
    speedTable,
    MetastatisTable,
    V,
    SimTime,
    exportSP=FALSE,
    exportRdata=FALSE,
    Rrows = Rrows)

  return(y)
}




#' Initialize Blood vessels
#'
#' Generate a random map of blood vessels to be used in the modeling of cell movement with cellmigRation
#'
#' @param mar integer, Minimum distance between the blood vessels and the sides of the tumor files
#' @param BloVesN integer, Number of blood vessels
#' @param Rrows integer, A number between 1 and 5 determining number of red boxes in the red zone. Specifically:
#'   \item{Rrows= 1 that means that there are 9 boxes in total (8 red + .lumen) }
#'   \item{Rrows= 2 that means that there are 25 boxes in total (24 red + .lumen) }
#'   \item{Rrows= 3 that means that there are 49 boxes in total (48 red + .lumen) }
#'   \item{Rrows= 4 that means that there are 81 boxes in total (80 red + .lumen) }
#'   \item{Rrows= 5 that means that there are 121 boxes in total (120 red + .lumen) }
#' @param minBloVesSep integer, Minimum distance between blood vessels'
#'
#' @param VesselsXY2CSVFile filename used for saving the VesselsXY results (CSV file)
#' @param LatDF2CSVFile filename used for saving the LatDF results (CSV file)
#'
#' @param plot logical, shall a plot be generated to visualize the distribution of vessels on the map
#' @param verbose logical, shall messages be printed to inform user about progress
#' @param seed integer, sets the seed for the random operation
#'
#' @return list including three elements:
#'
#'   \item{VesselsXY, data.frame including all vessel positions on the map }
#'   \item{LatDF, data.frame with 2 columns, i.e. lattice and LatRed }
#'   \item{vesselParams, list including the values of all parameters used for building the vessels map }
#'
#' @export
initializeBloodVessels <- function(mar= 100, BloVesN= 36,  Rrows= 4,
                                   minBloVesSep= 400,
                                   VesselsXY2CSVFile = NULL,
                                   LatDF2CSVFile = NULL,
                                   plot = TRUE, verbose = TRUE,
                                   seed = 12345){

  if (verbose)
    message("Generating vessels")
  # inner function
  inside.circle = function(center, radius, p, plot = FALSE){
    if ((length(center)) !=2){stop ("First argument must be a vector of length two.")}
    if ( radius<=0 ) {stop ("Second argument is not a positive number. Radius must be positive.")}
    if ( (length(radius)) !=1 ) {stop ("Second argument must contain only 1 element.") }

    p.1 = length(p[,1])
    x.coord = numeric(p.1)
    y.coord = numeric(p.1)
    for (i in 1:p.1){
      if ((p[i,1] - center[1])^2 + (p[i,2] - center[2])^2 <= radius^2){
        x.coord[i] = p[i,1]
        y.coord[i] = p[i,2]
      }
    }
    x.coord = x.coord[x.coord != 0]
    y.coord = y.coord[y.coord != 0]
    coordinates = cbind(x.coord,y.coord)
    theta = seq(0,2*pi,length = 2000)
    a = center[1]
    b = center[2]
    x = a + (radius*cos(theta))
    y = b + (radius*sin(theta))

    #Plot cirlce and add points
    if (plot) {
      plot(x,y, type = 'l', xlim = c(0,4000),
           ylim = c(0,4000), ylab = "Y", xlab = "X",col = "red")
      points(p[,1], p[,2], col = "blue", cex=1, pch=19)
      #The points inside the cirlce are red.
      points(x.coord,y.coord, col = "black",cex=1, pch=19)
    }
    return(coordinates)
  }


  # setSeed
  if (is.null(seed)) {
    seed <- sample(x = 1:999999, size = 1)
  } else if (is.na(seed)) {
    seed <- sample(x = 1:999999, size = 1)
  }
  set.seed(seed = seed)

  #initialize DF
  df<-data.frame()
  df[1,1]= runif(1,mar,4000-mar)
  df[1,2]= runif(1,mar,4000-mar)
  for (i in 2:BloVesN){
    if (verbose)
      message(".", appendLF = FALSE)
    df[i,1]= runif(1,mar,4000-mar)
    df[i,2]= runif(1,mar,4000-mar)
    inV<-c()
    for(j in 1:(i-1)){
      p=matrix(c(df[i,1],df[i,2]),1,2)
      cen=c(df[j,1],df[j,2])
      r= minBloVesSep + (Rrows * 20 )
      pp<-inside.circle(cen,r,p)
      inV[j]=pp[1]
    }
    while(length(inV[!is.na(inV)])>0){
      df[i,1]= runif(1,mar,4000-mar)
      df[i,2]= runif(1,mar,4000-mar)
      inV<-c()
      for(j in 1:(i-1)){
        p=matrix(c(df[i,1],df[i,2]),1,2)
        cen=c(df[j,1],df[j,2])
        r= minBloVesSep + (Rrows * 20 )
        pp<-inside.circle(cen,r,p)
        inV[j]=pp[1]
      }
    }
  }

  # convert the X and  Y into boxes
  for (i in 1:BloVesN){
    df[i,3]=(floor(df[i,2]/20)*200)+ceiling(df[i,1]/20)
  }

  colnames(df)<-c("x","y", ".lumen")
  #write.csv(df,file="VesselsXY.csv")
  #save(df,file="VesselsXY.Rdata")



  # Drawing squares
  xleft<-df[,1]- (10 + (Rrows * 20))
  ybottom<-df[,2]-(10 + (Rrows * 20))
  xright<-df[,1]+ (10 + (Rrows * 20))
  ytop<-df[,2]+(10 + (Rrows * 20))

  xleftb<-df[,1]-10
  ybottomb<-df[,2]-10
  xrightb<-df[,1]+10
  ytopb<-df[,2]+10

  # plot
  if (plot) {
    plot(0,0,type = "n", xlim = c(0,4000), ylim = c(0,4000),xlab="X",ylab="Y")
    for (i in 1:BloVesN){
      rect(xleft, ybottom, xright, ytop,col="red",border="red")
    }
    for (i in 1:BloVesN){
      rect(xleftb, ybottomb, xrightb, ytopb,col="black",border="black")
    }
  }


  # making a lattice
  VESdf<-df
  lattice<-c(rep(0,40000))
  length(lattice)
  LatVes<-c()
  for (i in 1:BloVesN){
    LatVes[i]= (floor(VESdf[i,2]/20)*200)+ceiling(VESdf[i,1]/20)
  }
  # print (LatVes)

  for (i in LatVes){    # Now the lattice has the .lumen of the vessels
    lattice[i]=1
  }

  # making a lattice for the red areas
  Red=c()
  TotalRow=c(1:Rrows)
  if (length (TotalRow)<1 | length (TotalRow)> 5)  {
    stop( "Rrows has to be a number between 1 and 5")
  }


  floors1<-c(-200,0,200,NA,NA,NA,NA,NA,NA,NA,NA)
  floors2<-c(-400,-200,0,200,400,NA,NA,NA,NA,NA,NA)
  floors3<-c(-600,-400,-200,0,200,400,600,NA,NA,NA,NA)
  floors4<-c(-800,-600,-400,-200,0,200,400,600,800,NA,NA)
  floors5<-c(-1000,-800,-600,-400,-200,0,200,400,600,800,1000)
  floors<- data.frame(floors1,floors2,floors3,floors4,floors5)

  cDF=floors[,Rrows]
  Floors=cDF[complete.cases(cDF)]

  for (i in LatVes){
    for (n in Floors){
      if (Rrows==1) {
        columns<- c(i+n,i+n+1,i+n-1)
      } else if (Rrows==2) {
        columns<- c(i+n,i+n+1,i+n+2,i+n-1,i+n-2)
      } else if (Rrows==3) {
        columns<- c(i+n,i+n+1,i+n+2,i+n+3,i+n-1,i+n-2,i+n-3)
      } else if (Rrows==4) {
        columns<- c(i+n,i+n+1,i+n+2,i+n+3,i+n+4,i+n-1,i+n-2,i+n-3,i+n-4)
      } else if (Rrows==5) {
        columns<- c(i+n,i+n+1,i+n+2,i+n+3,i+n+4,i+n+5,i+n-1,i+n-2,i+n-3,i+n-4,i+n-5)
      }

      Red=c(Red,columns)
    }
  }

  # removing the .lumena from the red zones
  del<-c()
  for (i in LatVes){
    del<-c(del,which(Red==i))
  }
  RED<-Red[-del]

  #length(Red)
  #length(RED)

  # labeling the boxes of the red zones in the LatRed with "R"
  # generating a lattice with red boxes
  LatRed<-c(rep(0,40000))
  for (i in RED){
    LatRed[i]="R"
  }

  # Making a data frame that has lattice and LatRed
  LatDF<-data.frame(lattice = lattice,LatRed = factor(LatRed))

  if (!is.null(VesselsXY2CSVFile)) {
    write.csv(df,file=VesselsXY2CSVFile, row.names = FALSE)
  }
  if (!is.null(LatDF2CSVFile)) {
    write.csv(LatDF,file=LatDF2CSVFile, row.names = FALSE)
  }


  if (verbose)
    message("Done!")

  return(list(VesselsXY = df, LatDF = LatDF,
              vesselParams = list(mar = mar,
                                  BloVesN = BloVesN,
                                  Rrows = Rrows,
                                  minBloVesSep = minBloVesSep,
                                  seed = seed)))

}


#' Load Blood vessels
#'
#' Load a map of blood vessels to be used in the modeling of cell movement with cellmigRation
#'
#' @param VesselsXY_CSVFile filename including the VesselsXY results (CSV file)
#' @param LatDF_CSVFile filename including the LatDF results (CSV file)
#' @param plot logical, shall a plot be generated to visualize the distribution of vessels on the map
#'
#' @return list including three elements:
#'
#'   \item{VesselsXY, data.frame including all vessel positions on the map }
#'   \item{LatDF, data.frame with 2 columns, i.e. lattice and LatRed }
#'   \item{vesselParams, list including the values of all parameters used for building the vessels map }
#'
#' @export
loadBloodVessels <- function(VesselsXY_CSVFile = NULL,
                             LatDF_CSVFile = NULL,
                             plot = TRUE) {

  # Check
  stopifnot(
    file.exists(VesselsXY_CSVFile),
    file.exists(LatDF_CSVFile))

  # Load files
  Ves <- read.csv(VesselsXY_CSVFile, as.is = TRUE, header = TRUE)
  Lat <- read.csv(LatDF_CSVFile, as.is = TRUE, header = TRUE)

  # COmpute ori params
  VNum <- as.numeric(table(Lat$LatRed))[2]
  NCel <- as.numeric(table(Lat$lattice))[2]
  MCoef <- as.numeric(VNum / NCel)

  # pre-def
  tmp <- data.frame(R = c(1,  2,  3,  4,   5),
                    V = c(8, 24, 48, 80, 120))
  idx <- which.min(abs(tmp$V - MCoef))
  rrow <- tmp[idx, 'R']


  # Plot
  # Drawing squares
  df <- Ves
  Rrows <- rrow
  xleft<-df[,1]- (10 + (Rrows * 20))
  ybottom<-df[,2]-(10 + (Rrows * 20))
  xright<-df[,1]+ (10 + (Rrows * 20))
  ytop<-df[,2]+(10 + (Rrows * 20))

  xleftb<-df[,1]-10
  ybottomb<-df[,2]-10
  xrightb<-df[,1]+10
  ytopb<-df[,2]+10

  # plot
  if (plot) {
    plot(0,0,type = "n", xlim = c(0,4000), ylim = c(0,4000),xlab="X",ylab="Y")
    for (i in 1:VNum){
      rect(xleft, ybottom, xright, ytop,col="red",border="red")
    }
    for (i in 1:VNum){
      rect(xleftb, ybottomb, xrightb, ytopb,col="black",border="black")
    }
  }

  # return
  return(list(VesselsXY = Ves, LatDF = Lat,
              vesselParams = list(mar = NA,
                                  BloVesN = NCel,
                                  Rrows = rrow,
                                  minBloVesSep = NA,
                                  seed = NA)))
}



# --
# Modeling Functions


eval_function <- function(mean_sd){
  mean <- mean_sd[1]
  sd <- mean_sd[2]
  sample <- rtruncnorm(n = DT, a = minSpe, b = maxSpe, mean = mean, sd = sd)
  mean_diff <-abs((desired_meanSpe - mean(sample))/desired_meanSpe)
  sd_diff <- abs((desired_sdSpe - sd(sample))/desired_sdSpe)
  mean_diff + sd_diff
}

Density<-function(B,Lat){
  Den<-c()
  bot<-c(2:199)
  top<-c(39802:39999)
  left<-c()
  right<-c()
  for (r in 1:198){ left<-c(left,1+(r*200))}
  for (r in 1:198){ right<-c(right,200+r*200)}
  loc<-c()

  if (B %in% bot){loc<-c(B+199,B+200,B+201,B-1,B+1,B+39799,B+39800,B+39801)}
  else if  (B %in% top){loc<-c(B-1,B+1,B-201,B-200,B-199,B-39799,B-39800,B-39801)}
  else if  (B %in% left){loc<-c(B+200,B+201,B+1,B-200,B-199,B+199,B+399,B-1)}
  else if  (B %in% right){loc<-c(B+199,B+200,B-1,B-201,B-200,B+1,B-199,B-399)}

  else if  (B ==1){loc<-c(201,202,2,39801,39802,200,400,40000)}
  else if  (B ==200){loc<-c(400,399,199,201,1,40000,39999,39801)}
  else if  (B ==39801){loc<-c(39802,1,2,200,39601,39602,40000,39800)}
  else if  (B ==40000){loc<-c(39999,39799,39800,39801,39601,200,199,1)
  }else{loc=c(B+200,B+201,B+199,B-200,B-201,B-199,B+1,B-1)}
  #for (p in loc){PP<-c(PP,Lat[p])}
  Dens=sum(Lat[loc])
  if (Dens<=2){ Den=2}
  if (Dens<=4 & Dens>2){ Den=1}
  if (Dens<=6 & Dens>4){ Den=0}
  if (Dens==7){ Den=-1}
  if (Dens==8){ Den=-2}

  DenPos<-Lat[loc]     # to know which boxes are empty
  return (list(Den,loc[DenPos==0],DenPos,loc))
}

ForMovement<-function(x1,y1,x2,y2,speed1,speed2,B,Lat){
  alp2=pi/2
  alp1= acos(speed1/sqrt(speed1^2 +  speed2^2))          # initial data
  u=x2-x1
  v=y2-y1
  if (u==0 && v==0){
    u=1
    v=1
  }
  a3=sqrt(u^2+v^2)
  alp3=pi-(alp1+alp2)
  a2=a3*sin(alp2)/sin(alp3)
  RHS1=x1*u+y1*v+a2*a3*cos(alp1)
  RHS2=y2*u-x2*v+a2*a3*sin(alp1)
  x3=(1/a3^2)*(u*RHS1-v*RHS2)
  y3=(1/a3^2)*(v*RHS1+u*RHS2)
  Xr= x2
  Yr= y2
  Xs= x3
  Ys= y3
  An=seq(184,356, by=4)
  AngVx<-c()
  AngVy<-c()
  AvX<-c()
  AvY<-c()
  EmpB<-c()
  AvB<-c()

  for(i in An){
    Nx=Xr + cos(i*pi/180)*(Xs - Xr) - sin(i*pi/180)*(Ys-Yr)
    AngVx[match(i,An)]=Nx
    Ny=Yr + cos(i*pi/180)*(Ys - Yr) + sin(i*pi/180)*(Xs-Xr)
    AngVy[match(i,An)]=Ny
  }
  AngVx<- na.omit(AngVx)
  for (i in 1:length(AngVx)){AngVx[i]<-margin(AngVx[i])}
  AngVy<- na.omit(AngVy)
  for (i in 1:length(AngVy)){AngVy[i]<-margin(AngVy[i])}

  for (i in 1:length(An)){
    EmpB[i]<-(floor(AngVy[i]/20)*200)+ceiling (AngVx[i]/20)
  }
  for (z in 1:length(An)){
    if (EmpB[z]== B){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])

    }
    if (EmpB[z]!= B & Lat[EmpB[z]]== 0){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])
    }
  }
  return(list(v1=AvB,v2=AvX,v3=AvY))
}

BackMovement<-function(x1,y1,x2,y2,speed1,speed2,B,Lat){
  alp2=pi/2
  alp1= acos(speed1/sqrt(speed1^2 +  speed2^2))          # initial data
  u=x2-x1
  v=y2-y1
  if (u==0 && v==0){
    u=1
    v=1
  }
  a3=sqrt(u^2+v^2)
  alp3=pi-(alp1+alp2)
  a2=a3*sin(alp2)/sin(alp3)
  RHS1=x1*u+y1*v+a2*a3*cos(alp1)
  RHS2=y2*u-x2*v+a2*a3*sin(alp1)
  x3=(1/a3^2)*(u*RHS1-v*RHS2)
  y3=(1/a3^2)*(v*RHS1+u*RHS2)
  Xr= x2
  Yr= y2
  Xs= x3
  Ys= y3
  An=-1*seq(180,360, by=4)
  AngVx<-c()
  AngVy<-c()
  AvX<-c()
  AvY<-c()
  EmpB<-c()
  AvB<-c()

  for(i in An){
    Nx=Xr + cos(i*pi/180)*(Xs - Xr) - sin(i*pi/180)*(Ys-Yr)
    AngVx[match(i,An)]=Nx
    Ny=Yr + cos(i*pi/180)*(Ys - Yr) + sin(i*pi/180)*(Xs-Xr)
    AngVy[match(i,An)]=Ny
  }
  AngVx<- na.omit(AngVx)
  for (i in 1:length(AngVx)){AngVx[i]<-margin(AngVx[i])}
  AngVy<- na.omit(AngVy)
  for (i in 1:length(AngVy)){AngVy[i]<-margin(AngVy[i])}

  for (i in 1:length(An)){
    EmpB[i]<-(floor(AngVy[i]/20)*200)+ceiling(AngVx[i]/20)
  }
  for (z in 1:length(An)){
    if (EmpB[z]== B){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])

    }
    if (EmpB[z]!= B & Lat[EmpB[z]]== 0){
      AvB= c(AvB, EmpB[z])
      AvX=c(AvX,AngVx[z])
      AvY=c(AvY,AngVy[z])
    }
  }
  return(list(v1=AvB,v2=AvX,v3=AvY))
}

marginXY<-function(NX,NY){
  if(NX>4000){
    NNX=NX-4000
  }else if(NX<0){
    NNX=4000-abs(NX)
  }else{
    NNX=NX
  }

  if(NY>4000){
    NNY=NY-4000
  }else if(NY<0){
    NNY=4000-abs(NY)
  }else{NNY=NY
  }
  return(c(NNX,NNY))
}

margin<-function(object){
  if(object>4000){
    object=object-4000
  }else if(object<0){
    object=4000-abs(object)
  }else{
    object=object
  }
  return(object)
}


NewDauLoc<-function(object){
  xy<-c()
  y=((ceiling(object/200))*20)-10
  x=((abs(object - (floor(object/200))*200))* 20 )-10
  if (x<=0){x=3990}
  xy<-c(x,y)
  return(xy)
}


PhylogenticTree<- function (phylogenicDF,TargetCell){
  steps<-c(1:TargetCell)
  RevSteps<-rev(steps)
  PHYLOGENIC<-c()
  for (i in RevSteps){
    if (TargetCell %in% phylogenicDF[i,]){
      PHYLOGENIC=phylogenicDF[i,]
      PHYLOGENIC=PHYLOGENIC[!is.na(PHYLOGENIC)]
      break
    }
  }
  return(PHYLOGENIC)
}

SimSum<-function(CellTables,cdDF,mainDF,SimTime,AA,CellDeath,RCellNum,
                 LatticeTable,CellSenescence,FinalAveSpe,MCells,V,PER,
                 G1Division,DividingCells,NotMotileCells,RanNum){
  simsumDF<-data.frame()
  Days<- seq(1,SimTime,by=24)
  totalSteps<-c(Days,SimTime)
  RCellNumPerStep<-c()
  CellNumPerStep<-c()
  totalnumCells<-c()
  CellDeathDay<-c()
  CellSenescenceDay<-c()
  Metastasis<-c()

  TTR=cumsum(RCellNum)
  TTIRZ<-TTR[totalSteps]
  for (i in totalSteps){              ### computing the number of cells in red per step
    main<-mainDF[i,]
    m1200<-main[!is.na(main)]
    NRm1200<-c()
    for (n in m1200){
      if (AA[n]==0){
        NRm1200<-c(NRm1200,n)
      }
    }
    RCellNumPerStep<- c(RCellNumPerStep,(length(m1200)- length(NRm1200)))
    CellNumPerStep<- c(CellNumPerStep,length(m1200))
    CellDeathDay<-c(CellDeathDay,CellDeath[i])
    CellSenescenceDay<-c(CellSenescenceDay,CellSenescence[i])
    Metastasis<-c(Metastasis,MCells[i])

    len1200<-cdDF[i,]         ########## computing the ID of the last cell
    FinalnumCell<-length(len1200[!is.na(len1200)])
    totalnumCell<-len1200[!is.na(len1200)]
    nameDF<-rbind(len1200,names(len1200))
    rownames(nameDF)<-NULL
    colnames(nameDF)<-nameDF[2,]
    w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
    totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
    totalnumCells<-c(totalnumCells,totalnumCell)
  }
  visits<-c()
  MaxCell<-c(1:totalnumCell)    ##### to check every cell
  for(d in MaxCell){
    targetCell<-CellTables[[d]][,7]
    if (length(targetCell[!is.na(targetCell)])>0 && "1" %in% targetCell[!is.na(targetCell)]){ visits[d]=d}else{visits[d]=NA}
  }
  visits=visits[!is.na(visits)]

  Day<-totalSteps        ############## computing uniformity index
  UIvalues<-c()
  vDays<- Day[-c(1:12)]                  # Not looking at the UI for the first 12 days since only very few cels are there
  for (i in vDays){
    TS=i
    UI<-UniformityIndex(TS,LatticeTable)
    UIvalues<-c(UIvalues,UI)
  }
  UIVALUES<-c(rep(0, 12),UIvalues)

  simsumDF[1:length(totalSteps),1]<-totalSteps
  simsumDF[1:length(totalSteps),2]<-CellNumPerStep
  simsumDF[1:length(totalSteps),3]<-RCellNumPerStep
  simsumDF[1:length(totalSteps),4]<-round(RCellNumPerStep/CellNumPerStep,digits=2)
  simsumDF[1:length(totalSteps),5]<- totalnumCells
  simsumDF[1:length(totalSteps),6]<- CellDeathDay
  simsumDF[1:length(totalSteps),7]<- CellSenescenceDay

  simsumDF[1:length(totalSteps),8]<-PER
  simsumDF[1:length(totalSteps),9]<-V
  simsumDF[1:length(totalSteps),10]<-TTIRZ
  simsumDF[1:length(totalSteps),11]<- UIVALUES
  simsumDF[1:length(totalSteps),12]<- c(rep(NA,(length(totalSteps)-1)),length(visits))
  simsumDF[1:length(totalSteps),13]<- Metastasis
  simsumDF[1:length(totalSteps),14]<- c(rep(NA,(length(totalSteps)-1)),FinalAveSpe)
  simsumDF[1:length(totalSteps),15]<- c(rep(NA,(length(totalSteps)-1)),length(G1Division))
  simsumDF[1:length(totalSteps),16]<- c(rep(NA,(length(totalSteps)-1)),NotMotileCells)
  simsumDF[1:length(totalSteps),17]<- c(rep(NA,(length(totalSteps)-1)),DividingCells)

  colnames(simsumDF)<-c("Steps","Number of cells","Number of cells in RedZone","RedZone cells ratio","MaxCellID",
                        "Cumulative Dead Cells","Cumulative Senescent cells","Persistence","Speed(um/h)","Cumulative time in RedZone",
                        "Uniformity Index","Total Number of RedZone Visitors","Cumulative Metastatic cells ","Average speed of all cells","Cumulative G1B cells",
                        "L5S Immotile Cells","L5S Dividing Cells")
  write.csv(simsumDF,file=paste0("simsumDF","R",RanNum,".csv"))
  return(simsumDF)
}


UniformityIndex<-function(TS,LatticeTable){
  lat1<-LatticeTable[TS,]
  smallBoxes<- c(0:49)

  oneOne<-seq(1,9801,by=200)
  oneOneBoxes<-c()
  for (i in oneOne){
    for (n in smallBoxes){
      oneOneBoxes<-c(oneOneBoxes,oneOne[match(i,oneOne)] + n)
    }}
  twoOneBoxes<- oneOneBoxes + 50
  threeOneBoxes<- oneOneBoxes + 100
  fourOneBoxes<- oneOneBoxes + 150
  oneTwo<-seq(10001,19801,by=200)
  oneTwoBoxes<-c()
  for (i in oneTwo){
    for (n in smallBoxes){
      oneTwoBoxes<-c(oneTwoBoxes,oneTwo[match(i,oneTwo)] + n)
    }}
  twoTwoBoxes<- oneTwoBoxes + 50
  threeTwoBoxes<- oneTwoBoxes + 100
  fourTwoBoxes<- oneTwoBoxes+ 150

  oneThree<-seq(20001,29801,by=200)
  oneThreeBoxes<-c()
  for (i in oneThree){
    for (n in smallBoxes){
      oneThreeBoxes<-c(oneThreeBoxes,oneThree[match(i,oneThree)] + n)
    }}

  twoThreeBoxes<- oneThreeBoxes + 50
  threeThreeBoxes<- oneThreeBoxes + 100
  fourThreeBoxes<- oneThreeBoxes + 150

  oneFour<-seq(30001,39801,by=200)
  oneFourBoxes<-c()
  for (i in oneFour){
    for (n in smallBoxes){
      oneFourBoxes<-c(oneFourBoxes,oneFour[match(i,oneFour)] + n)
    }}
  twoFourBoxes<- oneFourBoxes + 50
  threeFourBoxes<- oneFourBoxes + 100
  fourFourBoxes<- oneFourBoxes + 150

  LAT11<-lat1[oneOneBoxes]
  LAT11<-LAT11[LAT11==1]
  LAT11<-length(LAT11)

  LAT21<-lat1[twoOneBoxes]
  LAT21<-LAT21[LAT21==1]
  LAT21<-length(LAT21)

  LAT31<-lat1[threeOneBoxes]
  LAT31<-LAT31[LAT31==1]
  LAT31<-length(LAT31)

  LAT41<-lat1[fourOneBoxes]
  LAT41<-LAT41[LAT41==1]
  LAT41<-length(LAT41)



  LAT12<-lat1[oneTwoBoxes]
  LAT12<-LAT12[LAT12==1]
  LAT12<-length(LAT12)

  LAT22<-lat1[twoTwoBoxes]
  LAT22<-LAT22[LAT22==1]
  LAT22<-length(LAT22)

  LAT32<-lat1[threeTwoBoxes]
  LAT32<-LAT32[LAT32==1]
  LAT32<-length(LAT32)

  LAT42<-lat1[fourTwoBoxes]
  LAT42<-LAT42[LAT42==1]
  LAT42<-length(LAT42)

  LAT13<-lat1[oneThreeBoxes]
  LAT13<-LAT13[LAT13==1]
  LAT13<-length(LAT13)

  LAT23<-lat1[twoThreeBoxes]
  LAT23<-LAT23[LAT23==1]
  LAT23<-length(LAT23)

  LAT33<-lat1[threeThreeBoxes]
  LAT33<-LAT33[LAT33==1]
  LAT33<-length(LAT33)

  LAT43<-lat1[fourThreeBoxes]
  LAT43<-LAT43[LAT43==1]
  LAT43<-length(LAT43)

  LAT14<-lat1[oneFourBoxes]
  LAT14<-LAT14[LAT14==1]
  LAT14<-length(LAT14)

  LAT24<-lat1[twoFourBoxes]
  LAT24<-LAT24[LAT24==1]
  LAT24<-length(LAT24)

  LAT34<-lat1[threeFourBoxes]
  LAT34<-LAT34[LAT34==1]
  LAT34<-length(LAT34)

  LAT44<-lat1[fourFourBoxes]
  LAT44<-LAT44[LAT44==1]
  LAT44<-length(LAT44)
  allLATs<-c(LAT11,LAT21,LAT31,LAT41,
             LAT12,LAT22,LAT32,LAT42,
             LAT13,LAT23,LAT33,LAT43,
             LAT14,LAT24,LAT34,LAT44)

  meanallLATs<-mean(allLATs)
  sdallLATs<-sd(allLATs)
  CV=c()
  CV=round(sdallLATs/meanallLATs,digits=2)
  if (CV<=1){
    UI=(1 -CV )*100
  }else{
    UI=0
  }
  return(UI)
}


UIplot<-function(TS,LatticeTable,RanNum,V,PER){
  lat1<-LatticeTable[TS,]
  smallBoxes<- c(0:49)

  oneOne<-seq(1,9801,by=200)
  oneOneBoxes<-c()
  for (i in oneOne){
    for (n in smallBoxes){
      oneOneBoxes<-c(oneOneBoxes,oneOne[match(i,oneOne)] + n)
    }}
  twoOneBoxes<- oneOneBoxes + 50
  threeOneBoxes<- oneOneBoxes + 100
  fourOneBoxes<- oneOneBoxes + 150
  oneTwo<-seq(10001,19801,by=200)
  oneTwoBoxes<-c()
  for (i in oneTwo){
    for (n in smallBoxes){
      oneTwoBoxes<-c(oneTwoBoxes,oneTwo[match(i,oneTwo)] + n)
    }}
  twoTwoBoxes<- oneTwoBoxes + 50
  threeTwoBoxes<- oneTwoBoxes + 100
  fourTwoBoxes<- oneTwoBoxes+ 150

  oneThree<-seq(20001,29801,by=200)
  oneThreeBoxes<-c()
  for (i in oneThree){
    for (n in smallBoxes){
      oneThreeBoxes<-c(oneThreeBoxes,oneThree[match(i,oneThree)] + n)
    }}

  twoThreeBoxes<- oneThreeBoxes + 50
  threeThreeBoxes<- oneThreeBoxes + 100
  fourThreeBoxes<- oneThreeBoxes + 150

  oneFour<-seq(30001,39801,by=200)
  oneFourBoxes<-c()
  for (i in oneFour){
    for (n in smallBoxes){
      oneFourBoxes<-c(oneFourBoxes,oneFour[match(i,oneFour)] + n)
    }}
  twoFourBoxes<- oneFourBoxes + 50
  threeFourBoxes<- oneFourBoxes + 100
  fourFourBoxes<- oneFourBoxes + 150

  LAT11<-lat1[oneOneBoxes]
  LAT11<-LAT11[LAT11==1]
  LAT11<-length(LAT11)

  LAT21<-lat1[twoOneBoxes]
  LAT21<-LAT21[LAT21==1]
  LAT21<-length(LAT21)

  LAT31<-lat1[threeOneBoxes]
  LAT31<-LAT31[LAT31==1]
  LAT31<-length(LAT31)

  LAT41<-lat1[fourOneBoxes]
  LAT41<-LAT41[LAT41==1]
  LAT41<-length(LAT41)



  LAT12<-lat1[oneTwoBoxes]
  LAT12<-LAT12[LAT12==1]
  LAT12<-length(LAT12)

  LAT22<-lat1[twoTwoBoxes]
  LAT22<-LAT22[LAT22==1]
  LAT22<-length(LAT22)

  LAT32<-lat1[threeTwoBoxes]
  LAT32<-LAT32[LAT32==1]
  LAT32<-length(LAT32)

  LAT42<-lat1[fourTwoBoxes]
  LAT42<-LAT42[LAT42==1]
  LAT42<-length(LAT42)

  LAT13<-lat1[oneThreeBoxes]
  LAT13<-LAT13[LAT13==1]
  LAT13<-length(LAT13)

  LAT23<-lat1[twoThreeBoxes]
  LAT23<-LAT23[LAT23==1]
  LAT23<-length(LAT23)

  LAT33<-lat1[threeThreeBoxes]
  LAT33<-LAT33[LAT33==1]
  LAT33<-length(LAT33)

  LAT43<-lat1[fourThreeBoxes]
  LAT43<-LAT43[LAT43==1]
  LAT43<-length(LAT43)

  LAT14<-lat1[oneFourBoxes]
  LAT14<-LAT14[LAT14==1]
  LAT14<-length(LAT14)

  LAT24<-lat1[twoFourBoxes]
  LAT24<-LAT24[LAT24==1]
  LAT24<-length(LAT24)

  LAT34<-lat1[threeFourBoxes]
  LAT34<-LAT34[LAT34==1]
  LAT34<-length(LAT34)

  LAT44<-lat1[fourFourBoxes]
  LAT44<-LAT44[LAT44==1]
  LAT44<-length(LAT44)
  allLATs<-c(LAT11,LAT21,LAT31,LAT41,
             LAT12,LAT22,LAT32,LAT42,
             LAT13,LAT23,LAT33,LAT43,
             LAT14,LAT24,LAT34,LAT44)

  meanallLATs<-mean(allLATs)
  sdallLATs<-sd(allLATs)
  CV<- round(sdallLATs/meanallLATs,digits=2)
  if (CV<=1){ UI=(1 - CV)*100
  }else{
    UI=0
  }

  jpeg(paste0("UniformityPlot-","PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),
       width = 6, height = 6, units = 'in', res = 500)
  barplot(allLATs, main=paste0("UniformityPlot-","PR=",PER,"velecity=",V,"R",RanNum),
          xlab="Square Quadrat",ylab="Number of cells per Square Quadrat",
          sub=(paste0("Uniformity Index =",UI,"%"))  )
  dev.off()

}


PlottingCellsADLastStep<- function(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows){

  BloVesN <- nrow(VesselsXY)

  len1200<-cdDF[SimTime,]         ########## computing the ID of the last cell
  FinalnumCell<-length(len1200[!is.na(len1200)])
  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
  totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
  color<-c()
  nickColors <- function(n, h = c(120,480), l = c(.60,.70), s = c(.8,1), alpha = 1){
    return (alpha(hex(HLS(seq(h[1],h[2],length.out = n), seq(l[1],l[2],length.out = n), seq(s[1],s[2],length.out=n))), alpha))
  }
  n=totalnumCell
  color = nickColors(n)

  s=SimTime
  len1200<-mainDF[s,]         ########## computing the ID of the cell ID in each step
  FinalnumCell<-length(len1200[!is.na(len1200)])
  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  CellID<-as.numeric(nameDF[2,][!is.na(nameDF[1,])])

  lenID<-mainDF[s,]
  x=c(rep(NA,length(lenID[!is.na(lenID)])))
  y=c(rep(NA,length(lenID[!is.na(lenID)])))

  X=c()
  Y=c()
  boxXY<-data.frame(CellID,lenID[!is.na(lenID)],x,y)
  for (t in 1:length(CellID)){
    Nloc<- NewDauLoc(boxXY[t,2])              # the x and y of each cell
    boxXY[t,3]<-Nloc[1]
    boxXY[t,4]<-Nloc[2]
  }

  X=boxXY[,3]
  Y=boxXY[,4]
  black<-VesselsXY[,3]

  blackTable<-data.frame()
  blackTableX<-c()
  blackTableY<-c()
  for (i in black){
    df<- NewDauLoc(i)              # the x and y of each lumin (black boxes)
    blackTableX<-c(blackTableX,df[1])
    blackTableY<-c(blackTableY,df[2])
  }
  blackTable[1:BloVesN,1]<-blackTableX
  blackTable[1:BloVesN,2]<-blackTableY

  xleft<-blackTableX-(10 + (Rrows * 20))
  ybottom<-blackTableY-(10 + (Rrows * 20))
  xright<-blackTableX+(10 + (Rrows * 20))
  ytop<-blackTableY+(10 + (Rrows * 20))

  xleftb<-blackTableX- 10
  ybottomb<-blackTableY-10
  xrightb<-blackTableX+10
  ytopb<-blackTableY+10


  jpeg(paste0("Track plot-",s,"PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
  plot (0,0, type = "n",xlim=c(0,4000),ylim=c(0,4000),xlab="X",ylab="Y",main=paste0(s," Steps   ",length(X), "cells","  Velecity=",V,"   PR=",PER))

  for (i in 1:BloVesN){
    rect(xleft, ybottom, xright, ytop,col=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3),border=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3))
  }
  for (i in 1:BloVesN){
    rect(xleftb, ybottomb, xrightb, ytopb,col="black",border=NULL)
  }
  points(X,Y, col=color[boxXY[,1]],pch=20,cex=0.3)

  dev.off()

}

PlottingCellsAD<- function(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows){

  BloVesN <- nrow(VesselsXY)

  len1200<-mainDF[SimTime,]         ########## computing the ID of the last cell
  FinalnumCell<-length(len1200[!is.na(len1200)])
  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
  totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
  color<-c()
  nickColors <- function(n, h = c(120,480), l = c(.60,.70), s = c(.8,1), alpha = 1){
    return (alpha(hex(HLS(seq(h[1],h[2],length.out = n), seq(l[1],l[2],length.out = n), seq(s[1],s[2],length.out=n))), alpha))
  }
  n=totalnumCell
  color = nickColors(n)

  for (s in 1:SimTime){
    len1200<-mainDF[s,]         ########## computing the ID of the cell ID in each step
    FinalnumCell<-length(len1200[!is.na(len1200)])
    totalnumCell<-len1200[!is.na(len1200)]
    nameDF<-rbind(len1200,names(len1200))
    rownames(nameDF)<-NULL
    colnames(nameDF)<-nameDF[2,]
    CellID<-as.numeric(nameDF[2,][!is.na(nameDF[1,])])

    lenID<-mainDF[s,]
    x=c(rep(NA,length(lenID[!is.na(lenID)])))
    y=c(rep(NA,length(lenID[!is.na(lenID)])))

    X=c()
    Y=c()
    boxXY<-data.frame(CellID,lenID[!is.na(lenID)],x,y)
    for (t in 1:length(CellID)){
      Nloc<- NewDauLoc(boxXY[t,2])              # the x and y of each cell
      boxXY[t,3]<-Nloc[1]
      boxXY[t,4]<-Nloc[2]
    }

    X=boxXY[,3]
    Y=boxXY[,4]
    black<-VesselsXY[,3]


    blackTable<-data.frame()
    blackTableX<-c()
    blackTableY<-c()
    for (i in black){
      df<- NewDauLoc(i)              # the x and y of each lumin (black boxes)
      blackTableX<-c(blackTableX,df[1])
      blackTableY<-c(blackTableY,df[2])
    }
    blackTable[1:BloVesN,1]<-blackTableX
    blackTable[1:BloVesN,2]<-blackTableY

    xleft<-blackTableX-(10 + (Rrows * 20))
    ybottom<-blackTableY-(10 + (Rrows * 20))
    xright<-blackTableX+(10 + (Rrows * 20))
    ytop<-blackTableY+(10 + (Rrows * 20))


    xleftb<-blackTableX- 10
    ybottomb<-blackTableY-10
    xrightb<-blackTableX+10
    ytopb<-blackTableY+10

    jpeg(paste0("Track plot-",s,"PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
    plot (0,0, type = "n",xlim=c(0,4000),ylim=c(0,4000),xlab="X",ylab="Y",main=paste0(s," Steps   ",length(X), "cells","  Velecity=",V,"   PR=",PER))

    for (i in 1:BloVesN){
      rect(xleft, ybottom, xright, ytop,col=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3),border=rgb(red = 1, green = 0.3, blue = 0.1, alpha = 0.3))
    }
    for (i in 1:BloVesN){
      rect(xleftb, ybottomb, xrightb, ytopb,col="black",border=NULL)
    }
    points(X,Y, col=color[boxXY[,1]],pch=20,cex=0.3)

    dev.off()
  }
}

SimulationScript<- function(TimeInterval,DT,DTinR,VesselsXY,PER,LatDF,
                            PeriodForSenescence,PeriodForDeath,LeavingRzone,IntravasationProb,
                            mainDF,cdDF,phylogenicDF,CellTables,speedTable,MetastatisTable,
                            V,SimTime,exportSP=FALSE,exportRdata=FALSE, Rrows){
  AA<-LatDF[,2]
  LLC1<-as.numeric(LatDF[,1])         #for the first cell
  LatC1=LLC1
  LatticeTableC1<-c()
  for (i in 1:DT){
    LatticeTableC1<-rbind(LatticeTableC1,LatC1)
  }
  dim(LatticeTableC1)

  ####################################
  ###########  First cell  ###########
  ####################################

  #¤¤¤¤¤¤¤¤¤¤¤¤¤ First step
  CellTables[[1]][1,1]=1
  CellTables[[1]][1,2]=1
  Rx= runif(1,100,3900)
  Ry= runif(1,100,3900)
  B=(floor(Ry/20)*200)+ceiling(Rx/20)
  #############a while loop to avoid having the first cell starting from a lumin
  if (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
    while (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
      Rx= floor(runif(1,100,3900))
      Ry= floor(runif(1,100,3900))
      B=(floor(Ry/20)*200)+ceiling(Rx/20)
    }
  }

  CellTables[[1]][1,3]=Rx
  CellTables[[1]][1,4]=Ry
  CellTables[[1]][1,5]= (floor(Ry/20)*200)+ceiling(Rx/20)
  CellTables[[1]][1,6]=2
  if (AA[CellTables[[1]][1,5]]==0){
    CellTables[[1]][1,7]=0
  }else{
    CellTables[[1]][1,7]=1
  }
  CellTables[[1]][1,8]=1
  if (CellTables[[1]][1,7]==0){
    CellTables[[1]][1,9]=0
  }else{
    CellTables[[1]][1,9]=DT-DTinR -1
  }
  CellTables[[1]][1,10]=sum(CellTables[[1]][1,6],CellTables[[1]][1,7],CellTables[[1]][1,8],CellTables[[1]][1,9])
  mainDF[1,1]=CellTables[[1]][1,5]

  mainDF[1:5,1:5]
  CellTables[[1]][1:5,]
  cdDF[1,1]=1



  #¤¤¤¤¤¤¤¤¤¤¤¤¤ Second step ############################
  CellTables[[1]][2,1]=1
  CellTables[[1]][2,2]=2
  x=CellTables[[1]][1,3]
  y=CellTables[[1]][1,4]
  ran=round(runif(1,0,1),digits=2)
  dist=CellTables[[1]][2,13] * TimeInterval
  NX<-x + (cos(ran * 2* pi) *dist)
  NY<-y + (sin(ran * 2* pi) *dist)

  mar=marginXY(NX,NY)
  NX=mar[1]
  NY=mar[2]

  if (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
    while (LLC1[(floor(Ry/20)*200)+ceiling(Rx/20)]==1){
      x=CellTables[[1]][1,3]
      y=CellTables[[1]][1,4]
      ran=round(runif(1,0,1),digits=2)
      dist=CellTables[[1]][2,13] * TimeInterval
      NX<-x + (cos(ran * 2* pi) *dist)
      NY<-y + (sin(ran * 2* pi) *dist)
      mar=marginXY(NX,NY)
      NX=mar[1]
      NY=mar[2]
    }
  }

  CellTables[[1]][2,3]= NX
  CellTables[[1]][2,4]= NY
  CellTables[[1]][2,5]= (floor(Ry/20)*200)+ceiling(Rx/20)
  CellTables[[1]][2,6]=2
  if (AA[CellTables[[1]][2,5]]==0){
    CellTables[[1]][2,7]=0
  }else{
    CellTables[[1]][2,7]=1
  }
  CellTables[[1]][2,8]=2
  if (CellTables[[1]][2,7]==0){
    CellTables[[1]][2,9]=0
  }else{
    CellTables[[1]][2,9]=DT-DTinR -1
  }
  CellTables[[1]][2,10]=sum(CellTables[[1]][2,6],CellTables[[1]][2,7],CellTables[[1]][2,8],CellTables[[1]][2,9])
  CellTables[[1]][2,14]= round(sqrt((CellTables[[1]][2,3]- CellTables[[1]][1,3])^2 + (CellTables[[1]][2,4]- CellTables[[1]][1,4])^2),digits=2)

  mainDF[2,1]=CellTables[[1]][2,5]
  cdDF[2,1]=1
  G1Division<-c()


  #¤¤¤¤¤¤¤¤¤¤¤¤¤remaining steps
  for (s in 3:DT){
    Lat=LatticeTableC1[s,]
    CellTables[[1]][s,1]=1
    CellTables[[1]][s,2]=s

    x1=CellTables[[1]][s-2,3]
    y1=CellTables[[1]][s-2,4]

    x2=CellTables[[1]][s-1,3]
    y2=CellTables[[1]][s-1,4]

    speed1=CellTables[[1]][s-1,14]
    if (speed1==0){speed1=1.4}
    speed2=CellTables[[1]][s,13] * TimeInterval

    B=CellTables[[1]][s-1,5]
    Lat[B]=1
    ran=runif(1,0,1)

    if (ran<=(1-PER)){
      BM<-BackMovement(x1,y1,x2,y2,speed1,speed2,B,Lat)
      if (length(BM[[1]])>0){
        BMran<-ceiling(runif(1,0,length(BM[[1]])))
        BB<-BM[[1]][BMran]
        NX<-BM[[2]][BMran]
        NY<-BM[[3]][BMran]
      }else{
        BB<-B
        NX<-CellTables[[1]][s-1,3]
        NY<-CellTables[[1]][s-1,4]
      }
    }

    if (ran>(1-PER)){
      FM<-ForMovement(x1,y1,x2,y2,speed1,speed2,B,Lat)
      if (length(FM[[1]])>0){
        FMran<-ceiling(runif(1,0,length(FM[[1]])))
        BB<-FM[[1]][FMran]
        NX<-FM[[2]][FMran]
        NY<-FM[[3]][FMran]
      }else{
        BB<-B
        NX<-CellTables[[1]][s-1,3]
        NY<-CellTables[[1]][s-1,4]
      }
    }
    CellTables[[1]][s,3]= NX
    CellTables[[1]][s,4]= NY
    CellTables[[1]][s,5]= BB
    CellTables[[1]][s,6]=2
    if (AA[CellTables[[1]][s,5]]==0){
      CellTables[[1]][s,7]=0
    }else{
      CellTables[[1]][s,7]=1
    }

    ranRed<-runif(1,0,1)                                                     ################# Reducing the chance that a cell will leave the red zone
    if (CellTables[[1]][s,7]==0 && CellTables[[1]][s-1,7]==1 && ranRed<=LeavingRzone){
      CellTables[[1]][s,3]= CellTables[[1]][s-1,3]
      CellTables[[1]][s,4]= CellTables[[1]][s-1,4]
      CellTables[[1]][s,5]= CellTables[[1]][s-1,5]
      CellTables[[1]][s,6]= CellTables[[1]][s-1,6]
      CellTables[[1]][s,7]= CellTables[[1]][s-1,7]
    }

    CellTables[[1]][s,8]=s
    if (CellTables[[1]][s,7]==0){
      CellTables[[1]][s,9]=0
    }else{
      CellTables[[1]][s,9]= DT-DTinR - 1
    }

    if (sum(CellTables[[1]][2:s,7])>=G1 && CellTables[[1]][s,7]==0  &&  CellTables[[1]][s,8]>=DTinR  &&  s <DT){
      WWW=which(CellTables[[1]][,7]==1)
      stepID<-CellTables[[1]][,8]
      stepID=stepID[!is.na(stepID)]
      if (length(WWW)>0){
        ones<-stepID[WWW]
        for (i in ones){
          if(i+((G1)-1)>s){
            break
          }else{
            if (sum(CellTables[[1]][i:(i+((G1)-1)),7])==G1){
              CellTables[[1]][s,9]=DT-DTinR
              G1Division<-c(G1Division,1)
              break
            }
          }
        }
      }
    }

    CellTables[[1]][s,10]=sum(CellTables[[1]][s,6],CellTables[[1]][s,7],CellTables[[1]][s,8],CellTables[[1]][s,9])
    CellTables[[1]][s,11]=round(ran,digits=2)
    CellTables[[1]][s,12]=round(ranRed,digits=2)
    CellTables[[1]][s,14]= round(sqrt((CellTables[[1]][s,3]- CellTables[[1]][s-1,3])^2 + (CellTables[[1]][s,4]- CellTables[[1]][s-1,4])^2),digits=2)

    if (ran<=(1-PER)){
      CellTables[[1]][s,15]=length(BM[[1]])
    }else{
      CellTables[[1]][s,15]=length(FM[[1]])
    }

    PRAN<-runif(1,0,1)
    if (CellTables[[1]][s,10]>= DT & PRAN<=0.999){
      CellTables[[1]][s,16]=2
    }
    else if (CellTables[[1]][s,10]< DT && CellTables[[1]][s,10]> G1 && PRAN<=0.001){
      CellTables[[1]][s,16]=2
    }else{
      CellTables[[1]][s,16]=1
    }
    CellTables[[1]][s,17]=round(PRAN,digits=4)
    cdDF[s,1]=CellTables[[1]][s,16]
    mainDF[s,1]=CellTables[[1]][s,5]
    if (CellTables[[1]][s,16]==2){break}
  }
  phylogenicDF[1,1]=1


  ##################################################################
  ####################################### All other cells ##########
  ##################################################################

  target<-CellTables[[1]][,5]
  LenTar<-length(target[!is.na(target)])    # step of the first cell division
  RCellNumPerStep<-c()
  RCellNumPerStep=CellTables[[1]][1:LenTar,7]
  RCellNumPerStep<-c(RCellNumPerStep,rep(NA,(SimTime-LenTar)))

  CellDeath<-c()
  CellDeath<-c(rep(0,LenTar))
  CellDeath<-c(CellDeath,rep(NA,(SimTime-LenTar)))

  CellSenescence<-c()
  CellSenescence<-c(rep(0,LenTar))

  MetastaticCells<-c()
  cat("The first cell devided after ",LenTar, " steps","\n")

  ###########################################  Generating the LatticeTable ##################
  num=round(SimTime/5)
  black<-VesselsXY[,3]
  LL<-as.numeric(LatDF[,1])
  LL[black]=0
  Lat=LL
  LatticeTable<-c()
  for (i in 1:num){
    LatticeTable<-rbind(LatticeTable,Lat)
  }
  LatticeTable2<-LatticeTable
  for (i in 1:4){
    LatticeTable2<-rbind(LatticeTable2,LatticeTable)
  }
  LatticeTable<-rbind(LatticeTableC1[1:(LenTar-1),],LatticeTable2)

  for (step in (LenTar+1):SimTime){
    print(step)
    Lat=LatticeTable[step-1,]
    Lat1<-Lat
    TargetCDrow<-cdDF[step-1,]
    NetTargetCDrow<-which(!is.na(TargetCDrow))
    CellsPerStep<-as.numeric(names(TargetCDrow[NetTargetCDrow]))                        # to know the ID of the cell
    CellState<-cdDF[step-1,]
    CellNum<-length(CellsPerStep)
    HowManyCells<-ceiling(sum(CellState,na.rm = T))      ################## to convert the 0.9999 of a metastatic cell into 1
    MAX<-CellsPerStep[which.max(CellsPerStep)]
    MIN<-CellsPerStep[which.min(CellsPerStep)]

    NewCells<-MAX + (HowManyCells-CellNum)*2                ##### to know how many new cells are there
    NewCellsID<-c(MAX :NewCells)
    if (length(NewCellsID) >1){
      NewCellsID<-NewCellsID[-1]
      cellIDs<-c(CellState[which(CellState!=2)],NewCellsID)
      gone<-as.numeric(names(CellState[which(CellState==2)]))
      go<-c()
      ifold<-c()
      NewOld<-c()                                      # to save the cells that will divide but they are not able anymore to divide
      for(g in gone){go=c(go,rep(g,2))}

      for(g in gone){
        doneSteps<-c()
        doneSteps<-CellTables[[g]][,2]
        doneSteps<-doneSteps[!is.na(doneSteps)]
        doneSteps<-c(doneSteps,step)
        Ma<-match(step,doneSteps)                    # to match the simulation steps with the cell steps

        Dens<-Density(CellTables[[g]][Ma-1,5],Lat1)        # checking the density
        if (Dens[[1]]==-2  || (Dens[[1]]==-1 && Dens[[2]]  %in% black) ){
          NewOld<-c(NewOld,g)
          go<- go[!go %in% g]
          next
        }else{
          options<-which(go==g)      # for every g there are 2 options    options[[1]] and options[[2]]

          ######################## first daughter cell
          mainDF[step,MAX+ options[1]]= CellTables[[g]][Ma-1,5]
          CellTables[[MAX+ options[1]]][1,1]=MAX+ options[1]
          CellTables[[MAX+ options[1]]][1,2]=step-1
          CellTables[[MAX+ options[1]]][1,3]= CellTables[[g]][Ma-1,3]
          CellTables[[MAX+ options[1]]][1,4]= CellTables[[g]][Ma-1,4]

          CellTables[[MAX+ options[1]]][2,1]=MAX+ options[1]
          CellTables[[MAX+ options[1]]][2,2]=step
          Denran<-ceiling(runif(1,0,length(Dens[[2]])))
          Nloc1<-NewDauLoc(CellTables[[g]][Ma-1,5])                       # the x and y of the first daughter cell is
          CellTables[[MAX+ options[1]]][2,3]= Nloc1[1]
          CellTables[[MAX+ options[1]]][2,4]= Nloc1[2]
          CellTables[[MAX+ options[1]]][2,5]= CellTables[[g]][Ma-1,5]
          CellTables[[MAX+ options[1]]][2,6]= Dens[[1]]
          CellTables[[MAX+ options[1]]][2,7]= CellTables[[g]][Ma-1,7]
          CellTables[[MAX+ options[1]]][2,9]= CellTables[[g]][Ma-1,9]
          CellTables[[MAX+ options[1]]][2,10]= sum(Dens[[1]],CellTables[[g]][Ma-1,7],CellTables[[MAX+ options[1]]][2,8],CellTables[[g]][Ma-1,9])
          CellTables[[MAX+ options[1]]][2,11]= NA
          CellTables[[MAX+ options[1]]][2,12]= NA
          CellTables[[MAX+ options[1]]][2,13]= NA
          CellTables[[MAX+ options[1]]][2,14]= round(sqrt((CellTables[[MAX+ options[1]]][2,3]- CellTables[[g]][Ma-1,3])^2 + (CellTables[[MAX+ options[1]]][2,4]- CellTables[[g]][Ma-1,4])^2),digits=2)
          CellTables[[MAX+ options[1]]][2,15]= NA

          CellTables[[MAX+ options[1]]][2,16]=1
          cdDF[step,MAX+ options[1]]=CellTables[[MAX+ options[1]]][2,16]
          phylog<-phylogenicDF[g,]
          len.phylogenic<- length(phylog[!is.na(phylog)])
          num.parents<- len.phylogenic +1
          phylogenicDF[MAX+ options[1],1:num.parents]<-c(phylog[!is.na(phylog)],MAX+ options[1])

          ####################### second daughter cell
          mainDF[step,MAX+ options[2]]= Dens[[2]][Denran]
          CellTables[[MAX+ options[2]]][1,1]=MAX+ options[2]
          CellTables[[MAX+ options[2]]][1,2]=step-1
          CellTables[[MAX+ options[2]]][1,3]= CellTables[[g]][Ma-1,3]
          CellTables[[MAX+ options[2]]][1,4]= CellTables[[g]][Ma-1,4]

          CellTables[[MAX+ options[2]]][2,1]=MAX+ options[2]
          CellTables[[MAX+ options[2]]][2,2]=step
          Nloc2<- NewDauLoc(Dens[[2]][Denran])              # the x and y of the second daughter cell is
          CellTables[[MAX+ options[2]]][2,3]= Nloc2[1]
          CellTables[[MAX+ options[2]]][2,4]= Nloc2[2]

          BSecondDC=Dens[[2]][Denran]  ############# to get the box of the second daughter cell
          RanM<-runif(1,0,1)
          if (BSecondDC %in% black  && RanM>IntravasationProb){
            otherDenran<-c(1:length(Dens[[2]]))
            if (length(otherDenran)>1){
              otherDenran<-otherDenran[-Denran]
              Denran1<-sample(otherDenran, 1)
              CellTables[[MAX+ options[2]]][2,5]=Dens[[2]][Denran1]
              Lat1[Dens[[2]][Denran1]]=1
              Dens=Density(Dens[[2]][Denran1],Lat1)
              CellTables[[MAX+ options[2]]][2,6]= Dens[[1]]
              if (AA[CellTables[[MAX+ options[2]]][2,5]]==0){CellTables[[MAX+ options[2]]][2,7]=0 }else{ CellTables[[MAX+ options[2]]][2,7]=1}
              if (CellTables[[MAX+ options[2]]][2,7]==0){ CellTables[[MAX+ options[2]]][2,9]=0}else{CellTables[[MAX+ options[2]]][2,9]=DT-DTinR -1}
              CellTables[[MAX+ options[2]]][2,10]= sum(Dens[[1]],CellTables[[MAX+ options[2]]][2,7],CellTables[[MAX+ options[2]]][2,8],CellTables[[MAX+ options[2]]][2,9])
              CellTables[[MAX+ options[2]]][2,11]= NA
              CellTables[[MAX+ options[2]]][2,12]= NA
              CellTables[[MAX+ options[2]]][2,13]= NA
              CellTables[[MAX+ options[2]]][2,14]= round(sqrt((CellTables[[MAX+ options[2]]][2,3]- CellTables[[g]][Ma-1,3])^2 + (CellTables[[MAX+ options[2]]][2,4]- CellTables[[g]][Ma-1,4])^2),digits=2)
              CellTables[[MAX+ options[2]]][2,15]= NA
              CellTables[[MAX+ options[2]]][2,16]= 1
            }

          }
          if (!BSecondDC %in% black) {
            CellTables[[MAX+ options[2]]][2,5]=Dens[[2]][Denran]
            Lat1[Dens[[2]][Denran]]=1
            Dens=Density(Dens[[2]][Denran],Lat1)
            CellTables[[MAX+ options[2]]][2,6]= Dens[[1]]
            if (AA[CellTables[[MAX+ options[2]]][2,5]]==0){CellTables[[MAX+ options[2]]][2,7]=0 }else{ CellTables[[MAX+ options[2]]][2,7]=1}
            if (CellTables[[MAX+ options[2]]][2,7]==0){ CellTables[[MAX+ options[2]]][2,9]=0}else{CellTables[[MAX+ options[2]]][2,9]=DT-DTinR -1}
            CellTables[[MAX+ options[2]]][2,10]= sum(Dens[[1]],CellTables[[MAX+ options[2]]][2,7],CellTables[[MAX+ options[2]]][2,8],CellTables[[MAX+ options[2]]][2,9])
            CellTables[[MAX+ options[2]]][2,11]= NA
            CellTables[[MAX+ options[2]]][2,12]= NA
            CellTables[[MAX+ options[2]]][2,13]= NA
            CellTables[[MAX+ options[2]]][2,14]= round(sqrt((CellTables[[MAX+ options[2]]][2,3]- CellTables[[g]][Ma-1,3])^2 + (CellTables[[MAX+ options[2]]][2,4]- CellTables[[g]][Ma-1,4])^2),digits=2)
            CellTables[[MAX+ options[2]]][2,15]= NA
            CellTables[[MAX+ options[2]]][2,16]= 1
          }

          if (BSecondDC %in% black  && RanM<=IntravasationProb){
            MetastaticCells<-c(MetastaticCells,MAX+ options[2])
            CellTables[[MAX+ options[2]]][2,5]=NA
            CellTables[[MAX+ options[2]]][2,16]=0.9999
            cdDF[step:SimTime,MAX+ options[2]]=0.9999

          }

          cdDF[step,MAX+ options[2]]=CellTables[[MAX+ options[2]]][2,16]
          phylog<-phylogenicDF[g,]
          len.phylogenic<- length(phylog[!is.na(phylog)])
          num.parents<- len.phylogenic +1
          phylogenicDF[MAX+ options[2],1:num.parents]<-c(phylog[!is.na(phylog)],MAX+ options[2])
        }
      }
    }

    ifold<- as.numeric(names(CellState[which(CellState==1)]))                         # to know if we have undivided cells

    ifold<-c(ifold,NewOld)
    if (length(ifold) >0){
      for (o in ifold){
        doneSteps<-c()
        doneSteps<-CellTables[[o]][,2]
        doneSteps<-doneSteps[!is.na(doneSteps)]
        doneSteps<-c(doneSteps,step)
        Ma<-match(step,doneSteps)                                            # to match the simulation steps with the cell steps


        #if(!is.null(CellTables[[o]][Ma,8]) && length(CellTables[[o]][Ma,8]) > 0  && CellTables[[o]][Ma,8]>= (DT *PeriodForDeath)){next}           #### to remove the cell that has spent 431 steps without dividing
        TEST02 <- tryCatch(
          {!is.null(CellTables[[o]][Ma,8]) && length(CellTables[[o]][Ma,8]) > 0  && CellTables[[o]][Ma,8]>= (DT *PeriodForDeath)},
          error = function(e) {FALSE})
        if(TEST02){next}           #### to remove the cell that has spent 431 steps without dividing


        CellTables[[o]][Ma,1]=o
        CellTables[[o]][Ma,2]=step

        x1=CellTables[[o]][Ma-2,3]
        y1=CellTables[[o]][Ma-2,4]
        x2=CellTables[[o]][Ma-1,3]
        y2=CellTables[[o]][Ma-1,4]
        speed1=CellTables[[o]][Ma-1,14]
        if (speed1==0){speed1=1.4}
        speed2=CellTables[[o]][Ma,13] * TimeInterval
        B=CellTables[[o]][Ma-1,5]
        ran=runif(1,0,1)

        if (ran<=(1-PER)){                                               #### if it is smaller or =
          BM<-BackMovement(x1,y1,x2,y2,speed1,speed2,B,Lat1)
          if (length(BM[[1]])>0){
            BMran<-ceiling(runif(1,0,length(BM[[1]])))
            BB<-BM[[1]][BMran]
            NX<-BM[[2]][BMran]
            NY<-BM[[3]][BMran]
          }else{
            BB<-B
            NX<-CellTables[[o]][Ma-1,3]
            NY<-CellTables[[o]][Ma-1,4]
          }
        }

        if (ran>(1-PER)){
          FM<-ForMovement(x1,y1,x2,y2,speed1,speed2,B,Lat1)
          if (length(FM[[1]])>0){
            FMran<-ceiling(runif(1,0,length(FM[[1]])))
            BB<-FM[[1]][FMran]
            NX<-FM[[2]][FMran]
            NY<-FM[[3]][FMran]
          }else{
            BB<-B
            NX<-CellTables[[o]][Ma-1,3]
            NY<-CellTables[[o]][Ma-1,4]
          }
        }
        ############
        RanM<-runif(1,0,1)                       ##########The cell will stay in place
        if (BB %in% black  && RanM>IntravasationProb){
          BB<-B
          NX<-CellTables[[o]][Ma-1,3]
          NY<-CellTables[[o]][Ma-1,4]
          CellTables[[o]][Ma,3]= NX
          CellTables[[o]][Ma,4]= NY
          CellTables[[o]][Ma,5]= BB
          Lat1[B]=0
          Lat1[BB]=1

          DenS=Density(BB,Lat1)
          CellTables[[o]][Ma,6]= DenS[[1]]
          if (AA[CellTables[[o]][Ma,5]]==0){
            CellTables[[o]][Ma,7]=0
          }else{
            CellTables[[o]][Ma,7]=1
          }

          ranRed<-runif(1,0,1)                                                     ################# Reducing the chance that a cell will leave the red zone
          if (CellTables[[o]][Ma,7]==0 && CellTables[[o]][Ma-1,7]==1 && ranRed<=LeavingRzone){
            CellTables[[o]][Ma,3]= CellTables[[o]][Ma-1,3]
            CellTables[[o]][Ma,4]= CellTables[[o]][Ma-1,4]
            CellTables[[o]][Ma,5]= CellTables[[o]][Ma-1,5]
            Lat1[B]=1
            Lat1[BB]=0
            CellTables[[o]][Ma,6]= CellTables[[o]][Ma-1,6]
            CellTables[[o]][Ma,7]= CellTables[[o]][Ma-1,7]
          }

          if (CellTables[[o]][Ma,7]==0){
            CellTables[[o]][Ma,9]=0
          }else{
            CellTables[[o]][Ma,9]= DT-DTinR - 1
          }

          if (sum(CellTables[[o]][2:Ma,7])>=G1 && CellTables[[o]][Ma,7]==0  &&  CellTables[[o]][Ma,8]>=DTinR  &&  Ma <DT){
            WWW=which(CellTables[[o]][,7]==1)
            stepID<-CellTables[[o]][,8]
            stepID=stepID[!is.na(stepID)]
            if (length(WWW)>0){
              ones<-stepID[WWW]
              for (i in ones){
                if(i+((G1)-1)>Ma){
                  break
                }else{
                  if (sum(CellTables[[o]][i:(i+((G1)-1)),7])==G1){
                    CellTables[[o]][Ma,9]=DT-DTinR
                    G1Division<-c(G1Division,o)
                    break
                  }
                }
              }
            }
          }

          CellTables[[o]][Ma,10]=sum(CellTables[[o]][Ma,6],CellTables[[o]][Ma,7],CellTables[[o]][Ma,8],CellTables[[o]][Ma,9])
          CellTables[[o]][Ma,11]=round(ran,digits=2)
          CellTables[[o]][Ma,12]=round(ranRed,digits=2)
          CellTables[[o]][Ma,14]= round(sqrt((CellTables[[o]][Ma,3]- CellTables[[o]][Ma-1,3])^2 + (CellTables[[o]][Ma,4]- CellTables[[o]][Ma-1,4])^2),digits=2)

          if (ran<=(1-PER)){
            CellTables[[o]][Ma,15]=length(BM[[1]])
          }else{
            CellTables[[o]][Ma,15]=length(FM[[1]])
          }

          PRAN<-runif(1,0,1)
          if (CellTables[[o]][Ma,10]>= DT && CellTables[[o]][Ma,6]!=-2 && PRAN<=0.999){
            CellTables[[o]][Ma,16]=2
          }
          else if (CellTables[[o]][Ma,10]< DT && CellTables[[o]][Ma,10]> G1 &&  PRAN<=0.001){
            CellTables[[o]][Ma,16]=2
          }else{
            CellTables[[o]][Ma,16]=1
          }

          CellTables[[o]][Ma,17]=round(PRAN,digits=4)
        }
        ####################### The cell will move

        if (!BB %in% black){
          CellTables[[o]][Ma,3]= NX
          CellTables[[o]][Ma,4]= NY
          CellTables[[o]][Ma,5]= BB
          Lat1[B]=0
          Lat1[BB]=1

          DenS=Density(BB,Lat1)
          CellTables[[o]][Ma,6]= DenS[[1]]
          if (AA[CellTables[[o]][Ma,5]]==0){
            CellTables[[o]][Ma,7]=0
          }else{
            CellTables[[o]][Ma,7]=1
          }

          ranRed<-runif(1,0,1)                                                     ################# Reducing the chance that a cell will leave the red zone
          if (CellTables[[o]][Ma,7]==0 && CellTables[[o]][Ma-1,7]==1 && ranRed<=LeavingRzone){
            CellTables[[o]][Ma,3]= CellTables[[o]][Ma-1,3]
            CellTables[[o]][Ma,4]= CellTables[[o]][Ma-1,4]
            CellTables[[o]][Ma,5]= CellTables[[o]][Ma-1,5]
            Lat1[B]=1
            Lat1[BB]=0
            CellTables[[o]][Ma,6]= CellTables[[o]][Ma-1,6]
            CellTables[[o]][Ma,7]= CellTables[[o]][Ma-1,7]
          }

          if (CellTables[[o]][Ma,7]==0){
            CellTables[[o]][Ma,9]=0
          }else{
            CellTables[[o]][Ma,9]= DT-DTinR - 1
          }

          if (sum(CellTables[[o]][2:Ma,7])>=G1 && CellTables[[o]][Ma,7]==0  &&  CellTables[[o]][Ma,8]>=DTinR  &&  Ma <DT){
            WWW=which(CellTables[[o]][,7]==1)
            stepID<-CellTables[[o]][,8]
            stepID=stepID[!is.na(stepID)]
            if (length(WWW)>0){
              ones<-stepID[WWW]
              for (i in ones){
                if(i+((G1)-1)>Ma){
                  break
                }else{
                  if (sum(CellTables[[o]][i:(i+((G1)-1)),7])==G1){
                    CellTables[[o]][Ma,9]=DT-DTinR
                    G1Division<-c(G1Division,o)
                    break
                  }
                }
              }
            }
          }

          CellTables[[o]][Ma,10]=sum(CellTables[[o]][Ma,6],CellTables[[o]][Ma,7],CellTables[[o]][Ma,8],CellTables[[o]][Ma,9])
          CellTables[[o]][Ma,11]=round(ran,digits=2)
          CellTables[[o]][Ma,12]=round(ranRed,digits=2)
          CellTables[[o]][Ma,14]= round(sqrt((CellTables[[o]][Ma,3]- CellTables[[o]][Ma-1,3])^2 + (CellTables[[o]][Ma,4]- CellTables[[o]][Ma-1,4])^2),digits=2)

          if (ran<=(1-PER)){
            CellTables[[o]][Ma,15]=length(BM[[1]])
          }else{
            CellTables[[o]][Ma,15]=length(FM[[1]])
          }

          PRAN<-runif(1,0,1)
          if (CellTables[[o]][Ma,10]>= DT && CellTables[[o]][Ma,6]!=-2 && PRAN<=0.999){
            CellTables[[o]][Ma,16]=2
          } else if (CellTables[[o]][Ma,10]< DT && CellTables[[o]][Ma,10]> G1 &&  PRAN<=0.001){
            CellTables[[o]][Ma,16]=2
          } else{
            CellTables[[o]][Ma,16]=1
          }

          CellTables[[o]][Ma,17]=round(PRAN,digits=4)

        }
        ############################ The cell will go away
        if (BB %in% black  && RanM<=IntravasationProb){
          MetastaticCells<-c(MetastaticCells,o)
          CellTables[[o]][Ma,5]=NA
          CellTables[[o]][Ma,16]=0.9999
          cdDF[step:SimTime,o]=0.9999
          CellTables[[o]][Ma,10]=CellTables[[o]][Ma-1,10]

        }

        ########################
        cdDF[step,o]=CellTables[[o]][Ma,16]
        mainDF[step,o]=CellTables[[o]][Ma,5]
      }

    }
    LatticeTable[step,]<-Lat1
    CelDet<-c()                                        ############ computing cell death per step
    for (i in 1:length(CellTables)){
      NumStep<-CellTables[[i]][,16]
      lenNumStep<-length(NumStep[!is.na(NumStep)])
      if (lenNumStep>=(PeriodForDeath * DT)){CelDet[i]=1}else{CelDet[i]=NA}
    }
    CellDeath[step]<-length(CelDet[!is.na(CelDet)])

    CelSen<-c()                                ############ computing cell  senescence per step
    for (i in 1:length(CellTables)){
      NumStep<-CellTables[[i]][,16]
      lenNumStep<-length(NumStep[!is.na(NumStep)])
      if (lenNumStep>(DT * PeriodForSenescence) && lenNumStep<(PeriodForDeath * DT)){CelSen[i]=1}else{CelSen[i]=NA}
    }

    CellSenescence[step]<-length(CelSen[!is.na(CelSen)])

    main<-mainDF[step,]                                  ############### computing total hours in red zone per step
    m1200<-main[!is.na(main)]
    NRm1200<-c()
    for (n in m1200){
      if (AA[n]==0){
        NRm1200<-c(NRm1200,n)
      }
    }
    RCellNumPerStep[step]<- length(m1200)- length(NRm1200)    # all - not red
    if (length(MetastaticCells)>0){
      MetastatisTable[step,1:length(MetastaticCells)]=MetastaticCells      ### filling the MetastatisTable
    }
  }

  cat("Number of Metastatic Cells = ",length(MetastaticCells), "\n")

  ######################### To get the maximum cellID in the end of the steps
  len1200<-cdDF[SimTime,]
  FinalnumCell<-length(len1200[!is.na(len1200)])
  cat("Number of Cells after ", SimTime, "hours = ",FinalnumCell, "\n")

  totalnumCell<-len1200[!is.na(len1200)]
  nameDF<-rbind(len1200,names(len1200))
  rownames(nameDF)<-NULL
  colnames(nameDF)<-nameDF[2,]
  w<-which.max(nameDF[2,][!is.na(nameDF[1,])])
  totalnumCell<-as.numeric(nameDF[2,][!is.na(nameDF[1,])][w])
  cat("Total Number of Cells generated after ", SimTime, "hours = ",totalnumCell, "\n")
  #################Saving tables
  RanNum<-round(runif(1,100,10000),digits=2)
  ################Saving plots
  TS=SimTime

  RCellNum=RCellNumPerStep

  divisions<-c(1:20)
  CellNumberPerDivision<-c()
  for(i in divisions){
    len17<-phylogenicDF[,i]
    len17<-len17[!is.na(len17)]
    noDub<-len17[!duplicated(len17)]
    CellNumberPerDivision[i]<-length(noDub)
  }
  write.csv(CellNumberPerDivision,file=paste0("CellNumberPerDivision","PR=",PER,"V=",V,"R",RanNum,".csv"))


  ########################## computing the average speed of all cells
  AveSpe<-c()
  for (i in 1:totalnumCell){
    SPEED<-CellTables[[i]][,14]
    SPEED<-SPEED[!is.na(SPEED)]
    if (length(SPEED)>100){
      for(s in SPEED){
        if (s>250){                            ##### to correct the Periodic boundary conditions effect on speed
          mat<-match(s,SPEED)
          SPEED[mat]=4000-s
        }
      }
      AveSpe[i]<-mean(SPEED)
    }else{
      AveSpe[i]<-NA
    }
  }
  FinalAveSpe<-round(mean(AveSpe,na.rm=T),digits=3)


  ########################## computing number of metastatic cells per step
  MCells<-c()
  for (i in 1:SimTime){
    LenMCells<-MetastatisTable[i,]
    LenMCells<-LenMCells[!is.na(LenMCells)]
    MCells[i]<-length(LenMCells)
  }



  ################################################## Computing the number of immotile cells and dividing cells last five steps
  allCells<-c()
  NonDividingCells<-c()
  NotMotile<-c()
  L5S<-c()
  S2155<-mainDF[SimTime-5,]
  S2155<-S2155[!is.na(S2155)]                    # the ID of the first cell in Step=2155
  TotalCells<-totalnumCell                       # total number of cells in the end of the simulation
  for (i in S2155[1]:TotalCells){
    L5S<-mainDF[(SimTime-5):SimTime,i]
    L5S<-L5S[!is.na(L5S)]
    if (length(L5S)>0){
      allCells<-c(allCells,length(L5S))            # to get the length of cell tracks last 5 steps
      if (length(L5S)==6){NonDividingCells<-c(NonDividingCells,1)}
      if (length(L5S)==6 && mainDF[SimTime,i]==mainDF[(SimTime-1),i] && mainDF[(SimTime-1),i]==mainDF[(SimTime-2),i] && mainDF[(SimTime-2),i]==mainDF[(SimTime-3),i]&& mainDF[(SimTime-3),i]==mainDF[(SimTime-4),i]&& mainDF[(SimTime-4),i]==mainDF[(SimTime-5),i]){
        NotMotile<-c(NotMotile,1)
      }
    }
  }

  NotMotileCells<-length(NotMotile)
  DividingCells=length(allCells) - length(NonDividingCells)


  simsumDF<-SimSum(CellTables,cdDF,mainDF,SimTime,AA,CellDeath,RCellNum,LatticeTable,CellSenescence,FinalAveSpe,MCells,V,PER,G1Division,DividingCells,NotMotileCells,RanNum)

  write.csv(phylogenicDF,file=paste0("phylogenicDF","PR=",PER,"Speed=",V,"R",RanNum,".csv"))
  #write.csv(simsumDF,file=paste0("simsumDF","PR=",PER,"Speed=",V,"R",RanNum,".csv"))
  jpeg(paste0("G1B Cells"," PR=",PER," Speed=",V," Number of G1B cells=",length(G1Division),"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
  G1B.Cells=c(1:length(G1Division))
  CellID=G1Division
  LensimsumDF<-length(simsumDF[,1])
  plot(G1B.Cells,CellID,type="p",ylim=c(1,simsumDF[LensimsumDF,5]),main=paste0("G1Bcells=",length(G1Division),", PR=",PER,", Speed=",V,", Tumor cells=",simsumDF[LensimsumDF,5]),cex.main=1)

  for (i in 1:length(simsumDF[,1])){
    abline(h=simsumDF[i,5],col="red")
  }
  dev.off()

  # This seems to raise an error!!!
  # check
  PlottingCellsADLastStep(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows)
  jpeg(paste0("CellNumberPerDivision-",s,"PR=",PER,"velecity=",V,"R",RanNum,".jpeg"),width = 6, height = 6, units = 'in', res = 500)
  plot(divisions,CellNumberPerDivision,type="o",ylab="Number of cells undergone cell division",ylim=c(0,30000),
       xlab="Cell division number",pch=16,main=paste0("Cell Division Curve","  (Velecity=",V," ,  PR=",PER,")"))
  dev.off()
  # check
  UIplot(TS,LatticeTable,RanNum,V,PER)

  if (exportSP== TRUE) {
    #check
    PlottingCellsAD(cdDF,mainDF,SimTime,RanNum,PER,V,VesselsXY,Rrows)

  }
  if (exportRdata== TRUE) {
    ScdDF<-cdDF[,1:totalnumCell]
    SmainDF<-mainDF[,1:totalnumCell]
    save(ScdDF, file=paste0("cdDF-","PR_",PER,"velecity_",V,"R",RanNum,".RData"))
    save(SmainDF, file=paste0("mainDF-","PR_",PER,"velecity_",V,"R",RanNum,".RData"))

  }
  return(list(simsumDF))
}



# Wrapper function
# to Simulate()
RunCellmigSimulation <- function(
  TimeInterval,            # in hours
  DT,                      # in hours
  DTinR,                   # in hours
  G1,                      # in hours
  SimTime,                 # in hours

  minSpe,                  # um/h
  maxSpe,                  # um/h
  desired_meanSpe,         # um/h
  desired_sdSpe,           # um/h

  VesselsList,             # vessel data

  PER,
  PeriodForSenescence,
  PeriodForDeath,
  LeavingRzone,
  IntravasationProb,

  SetCellN = NULL,
  exportSP=FALSE,
  exportRdata=FALSE)
{

  # Prepare / import
  LatDF <- VesselsList$LatDF
  VesselsXY <- VesselsList$VesselsXY
  BloVesN <- VesselsList$vesselParams$BloVesN
  Rrows <- VesselsList$vesselParams$Rrows

  # define V
  V = desired_meanSpe

  # impute EstCellN
  divN= ceiling(SimTime/DT)+1
  estCellN=1

  for (i in 1:divN){
    estCellN=c(estCellN,2^i)
  }
  EstCellN<-round(sum(estCellN)+(10/PeriodForDeath) *SimTime)

  #bypass
  if (!is.null(SetCellN) && is.numeric(SetCellN)) {
    EstCellN <- SetCellN
  }


  if (EstCellN > 100000){
    EstCellN = 100000
  }

  # main df and other dfs
  mainDF<-data.frame(matrix(NA, nrow = SimTime, ncol = EstCellN))
  rownames(mainDF)<-c(1:SimTime)
  colnames(mainDF)<-c(1:EstCellN)

  speedTable<-data.frame()
  for (i in 1:1000){                  ########## generating 1000 speed profiles
    o <- optim(c(V, desired_sdSpe), eval_function)
    Speed <- rtruncnorm(n = DT, a = minSpe, b = maxSpe, mean = o$par[1], sd = o$par[2])
    SpeedT<-list()
    for (d in 1:PeriodForDeath){
      SpeedT[[d]]<- Speed
    }
    SpeedT <- do.call(c, SpeedT)
    LenSpeedT<-length(SpeedT)
    speedTable[1:LenSpeedT,i]= round(SpeedT,digits=2)
  }

  CellTables <- vector(mode = "list", length = EstCellN)
  for (n in 1:EstCellN){
    CellTables [[n]]<- data.frame(matrix(NA, nrow = (DT *PeriodForDeath), ncol = 17))

    ran=ceiling(runif(1,0,1000))            #########  selecting randomly the speed profile from the speedTable
    CellTables [[n]][2:LenSpeedT,13]= speedTable[1:(LenSpeedT-1),ran]
    CellTables [[n]][2:LenSpeedT,8]=c(1:(LenSpeedT-1))
    colnames(CellTables [[n]])<-c("ID","Steps","X","Y","B","Density","Red","CellStep","NotRed","Sum","RanDR","RanRed","Pre.Distance (um)","Actual Distance (um)","Movement options","CD")
  }
  CellTables[[1]]<-CellTables[[1]][-1,]
  CellTables[[1]][1,13]=NA

  cdDF<-data.frame(matrix(NA, nrow = SimTime, ncol = EstCellN))
  rownames(cdDF)<-c(1:SimTime)
  colnames(cdDF)<-c(1:EstCellN)

  phylogenicDF<-data.frame(matrix(NA, nrow = EstCellN, ncol = 25))
  MetastatisTable<- data.frame(matrix(NA, nrow = SimTime, ncol = (SimTime *3)))

  # Run
  y <- SimulationScript(
    TimeInterval,
    DT,
    DTinR,
    VesselsXY,
    PER,
    LatDF,
    PeriodForSenescence,
    PeriodForDeath,
    LeavingRzone,
    IntravasationProb,
    mainDF,
    cdDF,
    phylogenicDF,
    CellTables,
    speedTable,
    MetastatisTable,
    V,
    SimTime,
    #SetCellN = SetCellN,
    exportSP=FALSE,
    exportRdata=FALSE,
    Rrows=Rrows)

  return(y)
}


