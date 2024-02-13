#' @title 绘制单一林下位点的加权进界邻体直射荫蔽度的示意图
#'
#' @param a 林下研究位点的横坐标、纵坐标与高度
#' @param b 样地内林木的横坐标与树高
#' @param maxtan 最大临界值
#' @param mintan 最小临界值
#' @param mi 太阳方位角范围参数，pi/mi为与180°太阳方位角的夹角
#' @param MI 校正系数
#' @param wei 线性校正系数
#' @param bandwidth  林木分布太阳方位角核概率密度函数的带宽
#'
#' @return 单一林下位点的加权进界邻体直射荫蔽度的示意图
#' @export
#'
#' @examples
#' a=data.frame(x=75,y=75,size=0)
#' b=Understory.light::b
#' Understory.light::plot.Tan.dependence.Beam_S.single(a=a, b=b, maxtan=2, mintan=0.5, mi=1.5, MI=17.8, wei=0.5, bandwidth=0.5)
#'
plot.Tan.dependence.Beam_S.single=
function(a,b,maxtan,mintan,mi,MI,wei,bandwidth)
{
  Tan.dependence.Beam_S.single=function(a,b,maxtan,mintan,mi,MI,bandwidth)
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
    if (cos(angle)/mintan< (-1/maxtan))
    {
      warning(print("mi too big and mintan too small, wrong"))
      }
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
    Nei.tree$tangent=Nei.tree$tangen*(Nei.tree$threshold/maxtan)
    if(nrow(Nei.tree)==0){
      Dispersion=0
      tangent=0
    }

    if(nrow(Nei.tree)>0){
      y=density_box(Nei.tree$angle,bandwidth)
      plot(y,type="l")
      y=subset(y,y[,2]>0)
      Dispersion=nrow(y)/1000
      tangent=sum(Nei.tree$tangen)
    }
    tangent=tangent
    tangent=(atan(tangent/(wei*MI)-pi)+pi/2)/pi
    tangent=tangent*Dispersion
    Tan.dependence.Beam_S=tangent

    if(nrow(Nei.tree)==0)
    {Tan.dependence.Beam_S=0}
    outcome =Tan.dependence.Beam_S
    outcome =list(a=a[1,],Nei.tree=Nei.tree,Tan.dependence.Beam_S=Tan.dependence.Beam_S)
    outcome
  }


  Ne = Tan.dependence.Beam_S.single(a,b,maxtan,mintan,mi,MI,bandwidth)
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
    Neigh = Neigh[, c(1, 2, 3, 10)]
    colnames(Neigh) = c("x", "y", "size","group")
    colnames(small) = c("x", "y", "size")
    unitg = rbind(center, Neigh)
    max = max(abs(Ne$Nei.tree[, 1] - Ne$a[, 1]), abs(Ne$Nei.tree[, 2] - Ne$a[, 2]))
  }

  ang=(pi-(pi/mi))
  pointangx=sin(ang)*(max(Ne$Nei.tree$Distance))
  pointangy=cos(ang)*(max(Ne$Nei.tree$Distance))
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

  library("ggplot2")
  p = ggplot()+theme_bw()
  p = p + geom_line(data = pointang, aes(x = x, y = y, group = group),linetype = 2, size = 1,color = "red")
  p = p + geom_line(data = unitg, aes(x = x, y = y, group = group),linetype = 2, size = 1)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 13/max(unitg[,3]) * a[1, 3] + 2, color = "red4", alpha = 0.7)
  p = p + geom_point(data = small, aes(x, y), size = 13/max(unitg[,3]) * small[, 3] + 2, color = "grey", alpha = 0.7)
  p = p + geom_point(data = d, aes(x, y), size = 13/max(unitg[,3]) * d[, 3] + 2, color = "grey", alpha = 0.7)
  p = p + geom_point(data = Neigh, aes(x, y), size = 13/max(unitg[,3]) * Neigh[, 3] + 2, color = "green4", alpha = 0.7)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = Neigh,aes(x, y), size = 2) + geom_point(data = small, aes(x,y), size = 2) + geom_point(data = d, aes(x, y), size = 2)
  p = p + lims(x = c(Ne$a[, 1] - 1.1 * max, Ne$a[, 1] + 1.1 *max), y = c(Ne$a[, 2] - 1.1 * max, Ne$a[, 2] + 1.1 * max))

  p = p + annotate("text", x = Ne$a[, 1] - max + 0.3 *(max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("Tan_Beam_S=",round(Ne$Tan.dependence.Beam_S, 3)))
  p = p + theme_bw()
  p
}

