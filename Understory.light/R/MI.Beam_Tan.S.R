
MI.Beam_Tan.S=function(minx, maxx, miny, maxy, b, seq, maxtan,mintan,mi)
{
  library(sp)
  library(gstat)
  library(tcltk)
  xgrid = seq(minx, maxx, length.out = seq + 1)
  ygrid = seq(miny, maxy, length.out = seq + 1)
  basexy = expand.grid(xgrid, ygrid)
  basexy[, 3] = 0
  colnames(basexy) <- c("x", "y", "Size")
  a = dplyr::distinct(basexy)
  

Tan.dependence.Beam_S.mult=function(a,b,maxtan,mintan,mi)
{
  library(tcltk)
  
  Tan.dependence.Beam_S.single=function(a,b,maxtan,mintan,mi)
  {

    angle=pi-pi/mi
    #px=sin(angle)*baix
    px=sin(angle)/mintan
    #py=cos(angle)*baix
    py=cos(angle)/mintan
    r_=(px^2+py^2-(1/maxtan)^2)/(2*py+2/maxtan)
    r=r_+(1/maxtan)
    
    
    
    Neighbourhood.single1 = function(a, b) {
      
      c = b[, 1:2]
      for (i in 1:nrow(b)) {
        c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
        d = (c[, 1] + c[, 2])^(1/2)
      }
      d = cbind(b, d)
      colnames(d) = c("x", "y", "Size", "Distance")
      d
    }
    
    Nei.tree = Neighbourhood.single1(a, b)
    Nei.tree[, 3] = Nei.tree[, 3] - a[, 3]
    Nei.tree = subset(Nei.tree, Size > 0)
    Nei.tree
    
    Nei.tree=unique(Nei.tree)
    row=c(which(Nei.tree[,1]==a[1,1]))
    Nei.tree[row,1]=a[1,1]+0.01
    a = a[rep(1, nrow(Nei.tree)), ]
    Weight = Nei.tree[, 1:2] - a[, 1:2]
    Weight[, 3] = atan(Weight[, 2]/abs(Weight[, 1]))
    Weight[, 3] = Weight[, 3]
    Weight[, 3] = (pi/2) + Weight[, 3]
    Weight[, 4] = ((pi) - mi * Weight[, 3])/(pi)
    Weight[which(Weight[, 4] < 0), 4] = 0
    Nei.tree$Weight=Weight[,4]
    Nei.tree$angle=pi-Weight[, 3]
    Nei.tree=subset(Nei.tree,Weight>0)
    Nei.tree$angle_o=pi-(asin(sin(Nei.tree$angle)/r*r_)+Nei.tree$angle)
    Nei.tree$threshold=(sin(Nei.tree$angle_o)*r/sin(Nei.tree$angle))
    Nei.tree$threshold=1/Nei.tree$threshold
    Nei.tree$tangent=Nei.tree$Size/Nei.tree$Distance
    Nei.tree=subset(Nei.tree,Nei.tree$threshold<Nei.tree$tangent)
    Nei.tree_left=subset(Nei.tree,Nei.tree$x<a[1,1])
    Nei.tree_left$angle=2*pi-Nei.tree_left$angle
    Nei.tree_right=subset(Nei.tree,! Nei.tree$x<a[1,1])
    Nei.tree_right$angle= Nei.tree_right$angle
    Nei.tree=rbind(Nei.tree_left,Nei.tree_right)
    Nei.tree$Wei_tangent=Nei.tree$tangent*Nei.tree$threshold/maxtan
    if(nrow(Nei.tree)==0){
      tangent=0
    }
    
    if(nrow(Nei.tree)>0){
      tangent=sum(Nei.tree$Wei_tangen)
    }
    tangent=tangent
   # tangent=(atan(tangent/(MI *pi))/pi * 2)
    
    Tan.dependence.Beam_S=tangent
    
    if(nrow(Nei.tree)==0)
    {Tan.dependence.Beam_S=0}
    outcome =Tan.dependence.Beam_S
    outcome =list(a=a[1,],Nei.tree=Nei.tree,Tan.dependence.Beam_S=Tan.dependence.Beam_S)
    outcome
  }
  
  d = matrix(NA, nrow(a), 3)
  pb = tkProgressBar("", "Percent complete %", 0, 100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Beam_S.single(a[j,],b,maxtan,mintan,mi)$Tan.dependence.Beam_S))
    info = sprintf("Percent complete %d%%", round(j * 100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)",
                                                  info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  colnames(d) = c("x", "y", "Tan_dependence_Beam_S")
  rownames(d) = 1:nrow(a)
  d = as.data.frame(d)
  d
}

data =Tan.dependence.Beam_S.mult(a,b,maxtan,mintan,mi)
Tan.dependence.Beam_S=data$Tan_dependence_Beam_S
Rank = rank(Tan.dependence.Beam_S)
Percentile_Rank_of_Tan.dependence.Beam_S = Rank/length(Rank)
nls = nls(Percentile_Rank_of_Tan.dependence.Beam_S ~ (atan(Tan.dependence.Beam_S/(MI)-pi)+pi/2)/pi, start = list(MI = 1), algorithm = "port")
res = summary(nls)
MI = res$parameters[1, 1]
Revise_Tan.dependence.Beam_S = (atan(Tan.dependence.Beam_S/(MI)-pi)+pi/2)/pi
plot(Percentile_Rank_of_Tan.dependence.Beam_S,Revise_Tan.dependence.Beam_S,
     main = "Q-Q plot", xlim = c(0, 1), ylim = c(0,
                                                 1)) + abline(0, 1, col = "red4")
R_square = 1 - sum((Percentile_Rank_of_Tan.dependence.Beam_S -
                      Revise_Tan.dependence.Beam_S)^2)/sum((Percentile_Rank_of_Tan.dependence.Beam_S - mean(Percentile_Rank_of_Tan.dependence.Beam_S))^2)
output = list(MI = MI, R_square = R_square)
output
}
