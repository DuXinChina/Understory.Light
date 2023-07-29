

MI.Beam_Tan.S=function (minx, maxx, miny, maxy, b, seq, tan, mi) 
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

Tan.dependence.Beam_S.mult = function(a, b, tan, mi, MI) {
  
  
  density_box=function(entity,bandwidth)
  {
    if (length(entity)== 1){
      bandwidth=(2*pi/mi)/10
    }
    
    if (missing(bandwidth)) {
      bandwidth=density(entity)$bw
    }
    min_x=min(entity)
    max_x=max(entity)
    x=seq(min_x, max_x,length.out=1000)
    y=rep(0, length(x))
    box=function(entity){
      for (i in 1:length(x)) {
        if (abs(x[i] - entity) <= bandwidth/2) {
          y[i] <<- y[i] + 1
        }
      }
    }
    for(i in 1:length(entity))
    {
      box(entity[i])
    }
    y=y/length(entity)
    result=cbind(x,y)
    colnames(result)=c("scale","density")
    result
  }
  
  
  Tan.dependence.Beam_S.single = function(a, b, tan, mi) {
    Neighbourhood.single1 = function(a, b, tan) {
      c = b[, 1:2]
      for (i in 1:nrow(b)) {
        c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
        d = (c[, 1] + c[, 2])^(1/2)
      }
      d = cbind(b, d)
      d$tangent = (d$size - a[, 3])/d$d
      d = subset(d, tangent > tan)
      colnames(d) = c("x", "y", "Size", 
                      "Distance", "tangent")
      d = subset(d, Distance > 0)
      d
    }
    Nei.tree = Neighbourhood.single1(a, b, tan)
    Nei.tree[, 3] = Nei.tree[, 3] - a[, 3]
    Nei.tree = subset(Nei.tree, Size > 0)
    Nei.tree
    
    if (nrow(Nei.tree)>0)
    {
      Nei.tree=unique(Nei.tree)
      a = a[rep(1, nrow(Nei.tree)), ]
      Weight = Nei.tree[, 1:2] - a[, 1:2]
      Weight[, 3] = atan(Weight[, 2]/abs(Weight[, 1]))
      Weight[, 3] = Weight[, 3]
      Weight[, 3] = (pi/2) + Weight[, 3]
      Weight[, 4] = ((pi) - mi * Weight[, 3])/(pi)
      Weight[which(Weight[, 4] < 0), 4] = 0
      Nei.tree$Weight=Weight[,4]
      Nei.tree$angle=Weight[, 3]
      Nei.tree_left=subset(Nei.tree,Nei.tree$x-a[1,1]<0)
      Nei.tree_left$angle=pi-Nei.tree_left$angle
      Nei.tree_right=subset(Nei.tree,Nei.tree$x-a[1,1]>0)
      Nei.tree_right$angle=pi+Nei.tree_right$angle
      Nei.tree=rbind(Nei.tree_left,Nei.tree_right)
      Nei.tree=subset(Nei.tree,Weight>0)
      if(nrow(Nei.tree)==0){
        Dispersion=0
        tangent=0
      }
      if(nrow(Nei.tree)>0){  
        y=density_box(Nei.tree$angle)
        y=subset(y,y[,2]>0)
        
        Dispersion=nrow(y)/1000
        tangent=sum(Nei.tree$tangent)
      }
      tangent=tangent
     # tangent=(atan(tangent/(MI *pi))/pi * 2)
     # tangent=tangent*Dispersion
      Tan.dependence.Beam_S=tangent
    }
    if(nrow(Nei.tree)==0)
    {Tan.dependence.Beam_S=0}
    outcome =Tan.dependence.Beam_S
    outcome
  }
  
  
  d = matrix(NA, nrow(a), 3)
  pb = tkProgressBar("", "Percent complete %", 
                     0, 100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Beam_S.single(a[j, 
    ], b, tan, mi)))
    info = sprintf("Percent complete %d%%", round(j *  100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)", info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  colnames(d) = c("x", "y", "Tan_dependence_Beam_S")
  rownames(d) = 1:nrow(a)
  d = as.data.frame(d)
  d
}

data = Tan.dependence.Beam_S.mult(a, b, tan, mi) 
Tan_dependence_Beam_S=data$Tan_dependence_Beam_S
Rank = rank(Tan_dependence_Beam_S)
Percentile_Rank_of_Tan_dependence_Beam_S = Rank/length(Rank)
nls = nls(Percentile_Rank_of_Tan_dependence_Beam_S ~ (atan(Tan_dependence_Beam_S/(MI * 
                                                                                        pi))/pi * 2), start = list(MI = 1), algorithm = "port")
res = summary(nls)
MI = res$parameters[1, 1]
Revise_Tan_dependence_Beam_S = (atan(Tan_dependence_Beam_S/(MI * 
                                                                  pi))/pi * 2)
plot(Percentile_Rank_of_Tan_dependence_Beam_S,Revise_Tan_dependence_Beam_S,  
     main = "Q-Q plot", xlim = c(0, 1), ylim = c(0, 
                                                 1)) + abline(0, 1, col = "red4")
R_square = 1 - sum((Percentile_Rank_of_Tan_dependence_Beam_S - 
                      Revise_Tan_dependence_Beam_S)^2)/sum((Percentile_Rank_of_Tan_dependence_Beam_S - mean(Percentile_Rank_of_Tan_dependence_Beam_S))^2)
output = list(MI = MI, R_square = R_square)
output
}
