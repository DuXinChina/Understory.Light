#' @title  绘制尺度化加权散射荫蔽度的Kring插值图
#'
#' @param minx 绘图区域横坐标的最小值
#' @param maxx 绘图区域横坐标的最大值
#' @param miny 绘图区域纵坐标的最小值
#' @param maxy 绘图区域纵坐标的最大值
#' @param b 样地内林木的横坐标与树高
#' @param seq 抽样数据的空间分辨率
#' @param scale 局域尺度圆半径
#' @param MI 校正系数
#' @param wei 线性校正系数
#' @param scale2 背景尺度圆的半径，可以省略
#'
#' @return 绘制尺度化加权散射荫蔽度的Kring插值图
#' @export
#'
#' @examples
#' b=Understory.light::b
#' Understory.light::plot.Scale.dependence.Diffuse_S.Krig(minx=50, maxx=100, miny=50, maxy=100, b=b, seq=20, scale=5.4, MI=10.813,wei=0.5)
#' Understory.light::plot.Scale.dependence.Diffuse_S.Krig(minx=0, maxx=50, miny=0, maxy=50, b=b, seq=20, scale=5.4, MI=10.813, wei=0.5, scale2=15.5)
#'
plot.Scale.dependence.Diffuse_S.Krig=function (minx, maxx, miny, maxy, b, seq, scale, MI,wei,scale2)
{
library(sp)
library(gstat)
library(tcltk)
library(ggplot2)
if(missing(scale2))
{scale2=NA}
Scale.dependence.Diffuse_S.mult=function (a, b, scale, MI,wei,scale2)
{
  if(missing(scale2))
  {scale2=NA}
  Scale.dependence.Diffuse_S.single=function (a, b, scale,scale2)
  {
    if(missing(scale2))
    {scale2=NA}
    Scale.dependence.S.single = function(a, b, scale, MI) {
      Neighbourhood.single1 = function(a, b, scale) {
        c = b[, 1:2]
        for (i in 1:nrow(b)) {
          c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
          d = (c[, 1] + c[, 2])^(1/2)
        }
        d = cbind(b, d)
        if (is.na(scale2)==F)
        {d2=subset(d, d <= scale2 & d > scale)
        d2$tangent=d2$size/d2$d
        Scale2sum<<-sum(d2$tangent)
        }
        d = subset(d, d > 0)
        d = subset(d, d < scale)
        colnames(d) = c("x", "Y", "Size",
                        "Distance")
        d
      }
      Nei.tree = Neighbourhood.single1(a, b, scale)
      n = nrow(Nei.tree)
      Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - a[, 3])
      colnames(Nei.dif.high) = c("x", "Y", "Size",
                                 "Distance", "Size_dif ")
      n.dif = subset(Nei.dif.high, Nei.dif.high$Size_dif >
                       0)
      n.dif = nrow(n.dif)
      if (n.dif == 0) {
        Scale_dependence_S = 0
        Neighbourhood = Nei.dif.high
      }
      if (n.dif != 0) {
        Neighbourhood = Nei.dif.high
        Neighbourhood1 = Neighbourhood
        Neighbourhood[which(Neighbourhood[, 5] <= 0), 5] = NA
        Neighbourhood1[which(Neighbourhood1[, 5] <= 0), 5] = 0
        Neighbourhood1[which(Neighbourhood1[, 4] == 0), 4] = 1e-08
        S1 = Neighbourhood1[, 5]/Neighbourhood1[, 4]
        S2 = as.data.frame(S1)
        S = sum(S2, na.rm = T)
        Scale_dependence_S = S
      }
      Scale_dependence_S = (atan(Scale_dependence_S/(wei*MI)-pi)+pi/2)/pi
      outcome = list(a = a, Neighbourhood = Neighbourhood,
                     Scale_dependence_S = Scale_dependence_S)
      outcome
    }
    Scale.dependence.Simpson = function(a, b, scale) {
      Neighbourhood.single1 = function(a, b, scale) {
        c = b[, 1:2]
        for (i in 1:nrow(b)) {
          c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
          d = (c[, 1] + c[, 2])^(1/2)
        }
        d = cbind(b, d)



        d = subset(d, d <= scale)
        colnames(d) = c("x", "y", "Size",
                        "Distance")
        d
      }
      Nei.tree = Neighbourhood.single1(a, b, scale)
      Nei.tree
      Nei.tree[, 3] = Nei.tree[, 3] - a[, 3]
      Nei.tree = subset(Nei.tree, Size > 0)
      Nei.tree = subset(Nei.tree, Distance > 0)
      Nei.tree
      Simpn = sum(Nei.tree[, 3])
      Nei.tree1 = subset(Nei.tree, x >= a[, 1] & y > a[, 2])
      Nei.tree2 = subset(Nei.tree, x > a[, 1] & y <= a[, 2])
      Nei.tree3 = subset(Nei.tree, x <= a[, 1] & y < a[, 2])
      Nei.tree4 = subset(Nei.tree, x < a[, 1] & y >= a[, 2])
      Simpn1 = sum(Nei.tree1[, 3])
      Simpn2 = sum(Nei.tree2[, 3])
      Simpn3 = sum(Nei.tree3[, 3])
      Simpn4 = sum(Nei.tree4[, 3])
      simp_wei = sum((Simpn1/Simpn)^2, (Simpn2/Simpn)^2, (Simpn3/Simpn)^2,
                     (Simpn4/Simpn)^2)
      if (is.nan(simp_wei) == T)
        (simp_wei = 0.25)
      simp_wei
    }
    S = Scale.dependence.S.single(a, b, scale, MI)
    Levins_Simpson = Scale.dependence.Simpson(a, b, scale)
    a = S$a
    Neighbourhood = S$Neighbourhood
    Scale_dependence_S = S$Scale_dependence_S
    Scale_dependence_Wei_S = Scale_dependence_S * (0.25/Levins_Simpson)
    if (is.na(scale2)==F)
    {
      Ratio=1/(((MI*pi))/((Scale2sum)))
      Scale_dependence_Wei_S=(atan((Ratio+Scale_dependence_Wei_S)/(0.5/pi*(scale2/scale))-pi)+pi/2)/pi
    }

    outcome = list(a = a, Neighbourhood = Neighbourhood, Scale_dependence_S = Scale_dependence_S,
                   Levins_Simpson = Levins_Simpson, Scale_dependence_Wei_S = Scale_dependence_Wei_S)
    outcome
  }
  d = matrix(NA, nrow(a), 3)
  pb = tkProgressBar("", "Percent complete %",
                     0, 100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    res=Scale.dependence.Diffuse_S.single(a[j,
    ], b,scale,scale2)

    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(res$Scale_dependence_Wei_S))
    info = sprintf("Percent complete %d%%", round(j *
                                                    100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)",
                                                  info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  colnames(d) = c("x", "y", "Scale_dependence_Wei_S")
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
data =  Scale.dependence.Diffuse_S.mult(basexy, b, scale, MI,wei,scale2)
data1 = data
coordinates(data) <- c("x", "y")
spplot(data, "Scale_dependence_Wei_S")
vgm1 <- variogram(Scale_dependence_Wei_S ~ 1, data)
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
krige_res <- krige(Scale_dependence_Wei_S ~ 1, data, basexy1,
                   model = m)
z = krige_res$var1.pred
Basexyz = cbind(Basexy1, z)
colnames(Basexyz) = c("x", "y", "Vaule")
p2 = ggplot() + geom_raster(data = Basexyz, aes(x = x, y = y,
                                                fill = Vaule)) + theme_bw() + scale_fill_gradientn(limits = c(min(Basexyz$Vaule),
                                                                                                              1), colours = terrain.colors(10))
p2 = p2 + labs(title = "Scale dependence Wei_S") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,
                                                                       0))
print(p2)
}
