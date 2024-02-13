#' @title 绘制尺度化加权直射直射荫蔽度的Kring插值图
#'
#' @param minx 绘图区横坐标最小值
#' @param maxx 绘图区横坐标最大值
#' @param miny 绘图区纵坐标最小值
#' @param maxy 绘图区纵坐标最大值
#' @param b 林分中林木的坐标与树高
#' @param seq 绘制插值图的空间分辨率
#' @param aaix 最小尺度范围
#' @param baix 最大尺度范围
#' @param mi 太阳方位角范围参数，pi/mi为与180°太阳方位角的夹角
#' @param MI 校正系数
#' @param wei 线性校正系数
#' @param bandwidth 林木分布太阳方位角核概率密度函数的带宽
#'
#' @return 尺度化加权直射直射荫蔽度的Kring插值图
#' @export
#'
#' @examples
#' b=Understory.light::b
#' Understory.light::plot.Scale.dependence.Beam_S.Krig(minx=50, maxx=100, miny=50, maxy=100, b=b, seq=20, aaix=6.2, baix=14, mi=1.5, MI=11.2, wei=0.5, bandwidth=0.5)
#'
plot.Scale.dependence.Beam_S.Krig=function(minx,maxx,miny,maxy,b,seq,aaix,baix,mi,MI,wei,bandwidth)
{
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)

  Scale.dependence.Beam_S.mult=function(a,b,aaix,baix,mi,MI,wei,bandwidth)
  {
    library(tcltk)

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
       # plot(y,type="l")
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
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("", "Percent complete %", 0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Scale.dependence.Beam_S.single(a[j,],b,aaix,baix,mi,MI,bandwidth)$Scale.dependence.Beam_S))
      info = sprintf("Percent complete %d%%", round(j * 100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)",
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Scale_dependence_Beam_S")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  xgrid = seq(minx, maxx, length.out = seq + 1)
  ygrid = seq(miny, maxy, length.out = seq + 1)
  basexy = expand.grid(xgrid, ygrid)
  xgrid.cen = seq(minx + 0.5 * (maxx - minx) - 0.5/seq * (maxx -
                                                            minx), minx + 0.5 * (maxx - minx) + 0.5/seq * (maxx -
                                                                                                             minx), length.out = 3)
  ygrid.cen = seq(miny + 0.5 * (maxy - miny) - 0.5/seq * (maxy -
                                                            miny), miny + 0.5 * (maxy - miny) + 0.5/seq * (maxy -
                                                                                                             miny), length.out = 3)
  basexy.cen = expand.grid(xgrid.cen, ygrid.cen)
  xgrid.left.top = seq(minx + 0.25 * (maxx - minx) - 0.5/seq *
                         (maxx - minx), minx + 0.25 * (maxx - minx) + 0.5/seq *
                         (maxx - minx), length.out = 3)
  ygrid.left.top = seq(miny + 0.75 * (maxy - miny) - 0.5/seq *
                         (maxy - miny), miny + 0.75 * (maxy - miny) + 0.5/seq *
                         (maxy - miny), length.out = 3)
  basexy.left.top = expand.grid(xgrid.left.top, ygrid.left.top)
  xgrid.left.bottom = seq(minx + 0.25 * (maxx - minx) - 0.5/seq *
                            (maxx - minx), minx + 0.25 * (maxx - minx) + 0.5/seq *
                            (maxx - minx), length.out = 3)
  ygrid.left.bottom = seq(miny + 0.25 * (maxy - miny) - 0.5/seq *
                            (maxy - miny), miny + 0.25 * (maxy - miny) + 0.5/seq *
                            (maxy - miny), length.out = 3)
  basexy.left.bottom = expand.grid(xgrid.left.bottom, ygrid.left.bottom)
  xgrid.right.top = seq(minx + 0.75 * (maxx - minx) - 0.5/seq *
                          (maxx - minx), minx + 0.75 * (maxx - minx) + 0.5/seq *
                          (maxx - minx), length.out = 3)
  ygrid.right.top = seq(miny + 0.75 * (maxy - miny) - 0.5/seq *
                          (maxy - miny), miny + 0.75 * (maxy - miny) + 0.5/seq *
                          (maxy - miny), length.out = 3)
  basexy.right.top = expand.grid(xgrid.right.top, ygrid.right.top)
  xgrid.right.bottom = seq(minx + 0.75 * (maxx - minx) - 0.5/seq *
                             (maxx - minx), minx + 0.75 * (maxx - minx) + 0.5/seq *
                             (maxx - minx), length.out = 3)
  ygrid.right.bottom = seq(miny + 0.25 * (maxy - miny) - 0.5/seq *
                             (maxy - miny), miny + 0.25 * (maxy - miny) + 0.5/seq *
                             (maxy - miny), length.out = 3)
  basexy.right.bottom = expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  xgrid = seq(minx, maxx, length.out = 200)
  ygrid = seq(miny, maxy, length.out = 200)
  basexy1 = expand.grid(xgrid, ygrid)
  Basexy1 = basexy1
  basexy = rbind(basexy, basexy.cen, basexy.left.top, basexy.left.bottom,
                 basexy.right.top, basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x + y
  gridded(basexy1) <- TRUE
  basexy[, 3] = 0
  colnames(basexy) <- c("x", "y", "Size")
  basexy = dplyr::distinct(basexy)
  data =Scale.dependence.Beam_S.mult(basexy,b,aaix,baix,mi,MI,wei,bandwidth)
  data1 = data
  coordinates(data) <- c("x", "y")
  spplot(data, "Scale_dependence_Beam_S")
  vgm1 <- variogram(Scale_dependence_Beam_S ~ 1, data)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1, vgm("Sph"))
  print(m)
  sd = subset(vgm1$dist, vgm1$dist < m$range[2])
  spre = m$psill[1] + (m$psill[2]) * ((3 * sd)/(2 * m$range[2]) -
                                        sd^3/(2 * (m$range[2])^3))
  bd = subset(vgm1$dist, vgm1$dist > m$range[2])
  bpre = rep((m$psill[1] + m$psill[2]), length(bd))
  pre = rbind(as.matrix(spre), as.matrix(bpre))
  Coefficient_of_Determination = 1 - sum((pre - vgm1$gamma)^2)/sum((vgm1$gamma -
                                                                      mean(vgm1$gamma))^2)
  print(paste("Coefficient_of_Determination=", Coefficient_of_Determination))
  p1 = plot(vgm1, model = m)
  print(p1)
  krige_res <- krige(Scale_dependence_Beam_S ~ 1, data, basexy1,
                     model = m)
  z = krige_res$var1.pred
  Basexyz = cbind(Basexy1, z)
  colnames(Basexyz) = c("x", "y", "Vaule")
  p2 = ggplot() + geom_raster(data = Basexyz, aes(x = x, y = y,
                                                  fill = Vaule)) + theme_bw() + scale_fill_gradientn(limits = c(min(Basexyz$Vaule),
                                                                                                                1), colours = terrain.colors(10))
  p2 = p2 + labs(title = "Scale dependence Beam_S") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,
                                                                         0))
  print(p2)
}
