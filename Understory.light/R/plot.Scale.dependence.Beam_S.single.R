


plot.Scale.dependence.Beam_S.single=function(a,b,aaix,baix,mi,MI,wei,bandwidth)
{
Scale.dependence.Beam_S.single=function(a,b,aaix,baix,mi,MI,bandwidth)
{
  density_box=function(entity,bandwidth)
  {
    if (missing(bandwidth)) {
      bandwidth=pi/5
    }
    min_x=min(entity)
    max_x=max(entity)
    x=seq(pi-pi/mi,pi+pi/mi,length.out=1000)
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
    colnames(result)=c("angle","density")
    result
  }
  
  angle=pi-pi/mi
  px=sin(angle)*baix
  py=cos(angle)*baix
  r_=(px^2+py^2-aaix^2)/(2*py+2*aaix)
  r=r_+aaix
  p_o=data.frame(x=a[,1],y=a[,2]+r_)
  
  
  Neighbourhood.single1 = function(a, b) {

    c = b[, 1:2]
    c_= b[,1:2]
    for (i in 1:nrow(b)) {
      c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
      d = (c[, 1] + c[, 2])^(1/2)
    }
    
    for (i in 1:nrow(b)) {
      c_[i, ] = (b[i, 1:2] - p_o[1, 1:2])^2
      d_ = (c_[, 1] + c_[, 2])^(1/2)
    }
        
    d = cbind(b, d, d_)
    colnames(d) = c("x", "y", "Size", "Distance", "Distance_O")
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
    Nei.tree$angle=Weight[, 3]
    Nei.tree=subset(Nei.tree,Weight>0)
    Nei.tree=subset(Nei.tree,Nei.tree$Distance_O<r)
    Nei.tree_left=subset(Nei.tree, x< a[1,1])
    Nei.tree_left$angle=pi+Nei.tree_left$angle
    Nei.tree_right=subset(Nei.tree,!Nei.tree$x< a[1,1])
    Nei.tree_right$angle=pi-Nei.tree_right$angle
    Nei.tree=rbind(Nei.tree_left,Nei.tree_right)
    Nei.tree$tangent=Nei.tree$Size/Nei.tree$Distance
     if(nrow(Nei.tree)==0){
      Dispersion=0
      tangent=0
    }
   
     if(nrow(Nei.tree)>0){
      y=density_box(Nei.tree$angle,bandwidth)
      plot(y,type="l")
      y=subset(y,y[,2]>0)
      Dispersion=nrow(y)/1000
      tangent=sum(Nei.tree$tangent)
    }
    tangent=tangent
    tangent=(atan(tangent/(wei*MI)-pi)+pi/2)/pi
    tangent=tangent*Dispersion
    Scale.dependence.Beam_S=tangent
  
  if(nrow(Nei.tree)==0)
  {Scale.dependence.Beam_S=0}
  outcome =Scale.dependence.Beam_S
  outcome =list(a=a[1,],Nei.tree=Nei.tree,Scale.dependence.Beam_S=Scale.dependence.Beam_S)
  outcome
}
angle=pi-pi/mi
px=sin(angle)*baix
py=cos(angle)*baix
r_=(px^2+py^2-aaix^2)/(2*py+2*aaix)
r=r_+aaix
p_o=data.frame(x=a[,1],y=a[,2]+r_,r=r)

Ne = Scale.dependence.Beam_S.single(a,b,aaix,baix,mi,MI,bandwidth)
big = subset(Ne$Nei.tree, Size > a[1, 3])
small = subset(Ne$Nei.tree, Size <= a[1, 3])
n = nrow(big)
if (n == 0) {
  center = cbind(a[, 1:3], 1)
  colnames(center) = c("x", "y", "size", "group")
  unitg = center
  Neigh = subset(unitg, size < 0)
  colnames(Neigh) = c("x", "y", "size","group")
  colnames(small) = c("x", "y", "size")
  max = 10
}
if (n != 0) {
  center = cbind(cbind(rep(Ne$a[, 1], each = n), rep(Ne$a[,2], each = n), rep(Ne$a[, 3], each = n)), c(1:n))
  center = as.data.frame(center)
  colnames(center) = c("x", "y", "size","group")
  Neigh = cbind(big, c(1:n))
  Neigh = Neigh[, c(1, 2, 3, 9)]
  colnames(Neigh) = c("x", "y", "size","group")
  colnames(small) = c("x", "y", "size")
  unitg = rbind(center, Neigh)
  max = max(abs(Ne$Nei.tree[, 1] - Ne$a[, 1]), abs(Ne$Nei.tree[, 2] - Ne$a[, 2]))
}

ang=(pi-(pi/mi))
pointangx=sin(ang)*(baix)
pointangy=cos(ang)*(baix)
pointangy_=pointangy

max=max(max,abs(pointangx))

pointang=data.frame(x=c(a[1,1],a[1,1],pointangx+a[1,1],-pointangx+a[1,1]),y=c(a[1,2],a[1,2],pointangy+a[1,2],pointangy+a[1,2]),group=c(1,2,1,2))


c = b[, 1:2]
for (i in 1:nrow(b)) {
  c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
  d = (c[, 1] + c[, 2])^(1/2)
}
d = cbind(b, d)
d=subset(d,d<max)

seq=seq(0,2*pi,length.out=5000)
bou_x=sin(seq)*r
bou_y=cos(seq)*r+r_
bou=cbind(bou_x,bou_y)
bou=subset(bou,bou_y<pointangy)
if (nrow(bou)>30)
{
bou=bou[c(round(seq(1,nrow(bou),nrow(bou)/30)),nrow(bou)),]
}
bou1=cbind(bou[-nrow(bou),],1:(nrow(bou)-1))
bou2=cbind(bou[-1,],1:(nrow(bou)-1))
bou=rbind(bou1,bou2)
bou=as.data.frame(bou)
colnames(bou)=c("x","y","group")
bou$x=bou$x+a[1,1]
bou$y=bou$y+a[1,2]

library("ggplot2")
library(ggforce)
p = ggplot()+theme_bw()
p = p + geom_line(data = pointang, aes(x = x, y = y, group = group),linetype = 2, size = 1,color = "red")
p = p + geom_line(data = unitg, aes(x = x, y = y, group = group),linetype = 2, size = 1)
p = p + geom_point(data = Ne$a, aes(x, y), size = 13/max(unitg[,3]) * a[, 3] + 2, color = "red4", alpha = 0.7)
p = p + geom_point(data = small, aes(x, y), size = 13/max(unitg[,3]) * small[, 3] + 2, color = "grey", alpha = 0.7)
p = p + geom_point(data = d, aes(x, y), size = 13/max(unitg[,3]) * d[, 3] + 2, color = "grey", alpha = 0.7)
p = p + geom_point(data = Neigh, aes(x, y), size = 13/max(unitg[,3]) * Neigh[, 3] + 2, color = "green4", alpha = 0.7)
p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = Neigh,aes(x, y), size = 2) + geom_point(data = small, aes(x,y), size = 2) + geom_point(data = d, aes(x, y), size = 2)
p = p + lims(x = c(Ne$a[, 1] - 1.1 * max, Ne$a[, 1] + 1.1 *max), y = c(Ne$a[, 2] - 1.1 * max, Ne$a[, 2] + 1.1 * max))
p = p + geom_line(data = bou, aes(x = x, y = y, group = group),linetype = 4, size = 1)
p = p + annotate("text", x = Ne$a[, 1] - max + 0.3 *(max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("Scale_Beam_S=",round(Ne$Scale.dependence.Beam_S, 3)))
p = p + theme_bw()
p
}
