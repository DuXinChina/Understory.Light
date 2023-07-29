
plot.Tan.dependence.Beam_S.single=function(a, b, tan, mi,MI,bandwidth)
{
Tan.dependence.Beam_S.single = function(a, b, tan, mi,MI) {
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
      y=density_box(Nei.tree$angle,bandwidth)
      print(plot(y,type="l"))
      y=subset(y,y[,2]>0)

      Dispersion=nrow(y)/1000
      tangent=sum(Nei.tree$tangent)
    }
    tangent=tangent
    tangent=(atan(tangent/(MI *pi))/pi * 2)
    tangent=tangent*Dispersion
    Tan.dependence.Beam_S=tangent
  }
  if(nrow(Nei.tree)==0)
  {Tan.dependence.Beam_S=0}
  outcome =list(a=a[1,],Nei.tree=Nei.tree,Tan.dependence.Beam_S=Tan.dependence.Beam_S)
  outcome
}

Ne = Tan.dependence.Beam_S.single(a, b, tan=tan, mi=mi,MI=MI)
big = subset(Ne$Nei.tree, Size > a[, 3])
small = subset(Ne$Nei.tree, Size <= a[, 3])
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
  Neigh = Neigh[, c(1, 2, 3, 8)]
  colnames(Neigh) = c("x", "y", "size","group")
  colnames(small) = c("x", "y", "size")
  unitg = rbind(center, Neigh)
  max = max(abs(Ne$Nei.tree[, 1] - Ne$a[, 1]), abs(Ne$Nei.tree[, 2] - Ne$a[, 2]))
}

ang=(pi-(pi/mi))
pointangx=tan(ang)*(1.1*max)
pointangx=min(1.1*max,abs(pointangx))
pointangy=1/tan(ang)*(1.1*max)
pointangy_=pointangy
pointangy=min(1.1*max,abs(pointangy))
if(pointangy_<0)
{
  pointangy=-1*pointangy
}
pointang=data.frame(x=c(a[,1],a[,1],pointangx+a[,1],-pointangx+a[,1]),y=c(a[,2],a[,2],pointangy+a[,2],pointangy+a[,2]),group=c(1,2,1,2))


c = b[, 1:2]
for (i in 1:nrow(b)) {
  c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
  d = (c[, 1] + c[, 2])^(1/2)
}
d = cbind(b, d)
d=subset(d,d<max)

library("ggplot2")
p = ggplot()+theme_bw()
p = p + geom_line(data = pointang, aes(x = x, y = y, group = group),linetype = 2, size = 1,color = "red")
p = p + geom_line(data = unitg, aes(x = x, y = y, group = group),linetype = 2, size = 1)
p = p + geom_point(data = Ne$a, aes(x, y), size = 13/max(unitg[,3]) * a[, 3] + 2, color = "red4", alpha = 0.7)
p = p + geom_point(data = small, aes(x, y), size = 13/max(unitg[,3]) * small[, 3] + 2, color = "grey", alpha = 0.7)
p = p + geom_point(data = d, aes(x, y), size = 13/max(unitg[,3]) * d[, 3] + 2, color = "grey", alpha = 0.7)
p = p + geom_point(data = Neigh, aes(x, y), size = 13/max(unitg[,3]) * Neigh[, 3] + 2, color = "green4", alpha = 0.7)
p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = Neigh,aes(x, y), size = 2) + geom_point(data = small, aes(x,y), size = 2) + geom_point(data = d, aes(x, y), size = 2)
p = p + lims(x = c(Ne$a[, 1] - 1.1 * max, Ne$a[, 1] + 1.1 *max), y = c(Ne$a[, 2] - 1.1 * max, Ne$a[, 2] + 1.1 * max))
p = p + annotate("text", x = Ne$a[, 1] - max + 0.3 *(max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("Tan_Beam_S=",round(Ne$Tan.dependence.Beam_S, 3)))
p = p + theme_bw()
suppressMessages(print(p))
}
