#' @title 计算林下多个位点的加权进界邻体直射荫蔽度
#'
#' @param a 林下研究位点的横坐标、纵坐标与高度
#' @param b 样地内林木的横坐标与树高
#' @param maxtan 最大临界值
#' @param mintan 最小临界值
#' @param mi 太阳方位角范围参数，pi/mi为与180°太阳方位角的夹角
#' @param MI 校正系数
#' @param wei 线性校正系数
#' @param bandwidth 林木分布太阳方位角核概率密度函数的带宽
#'
#' @return 林下多个位点的加权进界邻体直射荫蔽度
#' @export
#'
#' @examples
#' a=data.frame(x=c(25,35,40),y=c(25,21,10),size=c(0,1,2))
#' b=Understory.Light::b
#' Understory.light::Tan.dependence.Beam_S.mult(a=a, b=b, maxtan=2, mintan=0.5, mi=1.5, MI=17.8, wei=0.5, bandwidth=0.5)
Tan.dependence.Beam_S.mult=function(a,b,maxtan,mintan,mi,MI,wei,bandwidth)
{
  library(tcltk)

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
   # Nei.tree$tangent=Nei.tree$tangen*(Nei.tree$threshold/maxtan)
    if(nrow(Nei.tree)==0){
      Dispersion=0
      tangent=0
    }

    if(nrow(Nei.tree)>0){
      y=density_box(Nei.tree$angle,bandwidth)
      # plot(y,type="l")
      y=subset(y,y[,2]>0)
      Dispersion=nrow(y)/1000
      tangent=sum(Nei.tree$tangen)
    }
    tangent=tangent
    tangent=(atan(tangent/(wei*MI)-pi)+pi/2)/pi
    #tangent=(atan(tangent/(MI *pi))/pi * 2)
    tangent=tangent*Dispersion
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
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Beam_S.single(a[j,],b,maxtan,mintan,mi,MI,bandwidth)$Tan.dependence.Beam_S))
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


