#' @title  绘制林下开敞度的kring插值图
#'
#' @param minx 绘图区域横坐标的最小值
#' @param maxx 绘图区域横坐标的最大值
#' @param miny 绘图区域纵坐标的最小值
#' @param maxy 绘图区域纵坐标的最大值
#' @param b 林分内林木的坐标及树高
#' @param seq 绘制插值图时抽样点的空间分辨率
#'
#' @return 林下开敞度的kring插值图与半变异函数
#' @export
#'
#' @examples
#' b=Understory.light::b
#' Understory.light::plot.K.Krig(minx=50,maxx=100,miny=50,maxy=100,b,seq=20)
plot.K.Krig=function (minx, maxx, miny, maxy, b, seq)
{
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  K.mult = function(a, b) {
    library(tcltk)
    Neighbourhood.K.single = function(a, b) {
      b1 = subset(b, x >= a[, 1] & y > a[, 2] & size >
                    a[, 3])
      b2 = subset(b, x > a[, 1] & y <= a[, 2] & size >
                    a[, 3])
      b3 = subset(b, x <= a[, 1] & y < a[, 2] & size >
                    a[, 3])
      b4 = subset(b, x < a[, 1] & y >= a[, 2] & size >
                    a[, 3])
      Neighbourhood.K.single1 = function(a, b) {
        b1 = subset(b, x >= a[, 1] & y > a[, 2] & size >
                      a[, 3])
        b2 = subset(b, x > a[, 1] & y <= a[, 2] & size >
                      a[, 3])
        b3 = subset(b, x <= a[, 1] & y < a[, 2] & size >
                      a[, 3])
        b4 = subset(b, x < a[, 1] & y >= a[, 2] & size >
                      a[, 3])
        Neighbourhood.single = function(a, b) {
          c = b[, 1:2]
          for (i in 1:nrow(b)) {
            c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
            d = (c[, 1] + c[, 2])^(1/2)
          }
          d = cbind(b, d)
          d = subset(d, d > 0)
          d = d[order(d[, 4])[1], 1:4]
          colnames(d) = c("x", "Y", "Size",
                          "Distance")
          d
        }
        b1.Nei.tree = Neighbourhood.single(a, b1)
        b2.Nei.tree = Neighbourhood.single(a, b2)
        b3.Nei.tree = Neighbourhood.single(a, b3)
        b4.Nei.tree = Neighbourhood.single(a, b4)
        Nei.tree = rbind(b1.Nei.tree, b2.Nei.tree, b3.Nei.tree,
                         b4.Nei.tree)
        Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] -
                               a[, 3])
        colnames(Nei.dif.high) = c("x", "Y",
                                   "Size", "Distance", "Size_dif ")
        Neighbourhood = Nei.dif.high
        k1 = Neighbourhood[, 4]/Neighbourhood[, 5]
        k2 = as.data.frame(k1)
        k2[, 1][is.infinite(k2[, 1])] = NA
        K = sum(k2, na.rm = T)
        key = 2
        outcome = list(a = a, Neighbourhood = Neighbourhood,
                       K = K, key = key)
        outcome
      }
      if (nrow(b1) != 0 & nrow(b2) != 0 & nrow(b3) != 0 &
          nrow(b4) != 0) {
        Neighbourhood.K.single1(a, b)
      }
      else {
        key1 = 1
        key = 1
        out = list(key = key, key1 = key1)
        out
      }
    }
    d = matrix(NA, nrow(a), 4)
    pb <- tkProgressBar("进度", "已完成 %",
                        0, 100)
    star_time <- Sys.time()
    for (j in 1:nrow(a)) {
      KEY = Neighbourhood.K.single(a[j, ], b)
      sw = KEY$key
      if (sw == 1) {
        warning(paste("检查第", c(j), "行"))
      }
      else {
        d[j, ] = cbind(as.matrix(a[j, ]), as.matrix(Neighbourhood.K.single(a[j,
        ], b)$K))
      }
      info <- sprintf("已完成 %d%%", round(j *
                                          100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("进度 (%s)",
                                                    info), info)
    }
    end_time <- Sys.time()
    close(pb)
    run_time <- end_time - star_time
    colnames(d) = c("X", "Y", "Size", "K")
    rownames(d) = 1:nrow(a)
    if (is.na(sum(d)) == TRUE) {
      print("请确保以参考林木为中心的四个象限内均有高于其的立木")
      d
    }
    else {
      d
    }
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
  data = K.mult(basexy, b)
  colnames(data) = c("x", "y", "size", "K")
  data = data[complete.cases(data), ]
  data = as.data.frame(data)
  data1 = data
  coordinates(data) <- c("x", "y")
  spplot(data, "K")
  vgm1 <- variogram(K ~ 1, data)
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
  krige_res <- krige(K ~ 1, data, basexy1, model = m)
  z = krige_res$var1.pred
  Basexyz = cbind(Basexy1, z)
  colnames(Basexyz) = c("x", "y", "Vaule")
  p2 = ggplot() + geom_raster(data = Basexyz, aes(x = x, y = y,
                                                  fill = Vaule)) + theme_bw() + scale_fill_gradientn(limits = c(0,
                                                                                                                max(Basexyz$Vaule)), colours = terrain.colors(10))
  p2 = p2 + labs(title = "K") + scale_x_continuous(expand = c(0,
                                                              0)) + scale_y_continuous(expand = c(0, 0))
  print(p2)
}
