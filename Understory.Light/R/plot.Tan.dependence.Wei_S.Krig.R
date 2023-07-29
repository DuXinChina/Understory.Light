plot.Tan.dependence.Wei_S.Krig=function (minx, maxx, miny, maxy, b, seq, tan, MI)
{
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)
  Tan.dependence.Wei_S.mult = function(a, b, tan) {
    library(tcltk)
    Tan.dependence.Wei_S.single = function(a, b, tan) {
      Tan.dependence.S.single = function(a, b, tan) {
        Neighbourhood.single1 = function(a, b, tan) {
          c = b[, 1:2]
          for (i in 1:nrow(b)) {
            c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
            d = (c[, 1] + c[, 2])^(1/2)
          }
          d = cbind(b, d)
          d$tangent = (d$size - a[, 3])/d$d
          d = subset(d, tangent > tan)
          colnames(d) = c("x", "Y", "Size",
                          "Distance", "tangent")
          d
        }
        Nei.tree = Neighbourhood.single1(a, b, tan)
        n.dif = nrow(Nei.tree)
        if (n.dif == 0) {
          Tan_dependence_S = 0
          Neighbourhood = Nei.tree
        }
        if (n.dif != 0) {
          Neighbourhood = Nei.tree
          Neighbourhood1 = Neighbourhood
          Neighbourhood[which(Neighbourhood[, 5] <= 0),
                        5] = NA
          Neighbourhood1[which(Neighbourhood1[, 5] <=
                                 0), 5] = 0
          S1 = Neighbourhood1[, 5]
          S2 = as.data.frame(S1)
          S = sum(S2, na.rm = T)
          Tan_dependence_S = S
        }
        Tan_dependence_S = (atan(Tan_dependence_S/(MI *
                                                     pi))/pi * 2)
        if (is.nan(Tan_dependence_S) == T)
          (Tan_dependence_S = 1)
        outcome = list(a = a, Neighbourhood = Neighbourhood,
                       Tan_dependence_S = Tan_dependence_S)
        outcome
      }
      Tan.dependence.Simpson = function(a, b, tan) {
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
        Simpn = sum(Nei.tree[, 3])
        Nei.tree = subset(Nei.tree, Distance > 0)
        Nei.tree1 = subset(Nei.tree, x >= a[, 1] & y >
                             a[, 2])
        Nei.tree2 = subset(Nei.tree, x > a[, 1] & y <=
                             a[, 2])
        Nei.tree3 = subset(Nei.tree, x <= a[, 1] & y <
                             a[, 2])
        Nei.tree4 = subset(Nei.tree, x < a[, 1] & y >=
                             a[, 2])
        Simpn1 = sum(Nei.tree1[, 3])
        Simpn2 = sum(Nei.tree2[, 3])
        Simpn3 = sum(Nei.tree3[, 3])
        Simpn4 = sum(Nei.tree4[, 3])
        simp_wei = sum((Simpn1/Simpn)^2, (Simpn2/Simpn)^2,
                       (Simpn3/Simpn)^2, (Simpn4/Simpn)^2)
        if (is.nan(simp_wei) == T)
          (simp_wei = 0.25)
        simp_wei
      }
      S = Tan.dependence.S.single(a, b, tan)
      Levins_Simpson = Tan.dependence.Simpson(a, b, tan)
      a = S$a
      Neighbourhood = S$Neighbourhood
      Tan_dependence_S = S$Tan_dependence_S
      Tan_dependence_Wei_S = Tan_dependence_S * (0.25/Levins_Simpson)
      outcome = list(a = a, Neighbourhood = Neighbourhood,
                     Tan_dependence_S = Tan_dependence_S, Levins_Simpson = Levins_Simpson,
                     Tan_dependence_Wei_S = Tan_dependence_Wei_S)
      outcome
    }
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("", "Percent complete %",
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Wei_S.single(a[j,
      ], b, tan)$Tan_dependence_Wei_S))
      info = sprintf("Percent complete %d%%", round(j *
                                                      100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress  (%s)",
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Tan_dependence_Wei_S")
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
  data = Tan.dependence.Wei_S.mult(basexy, b, tan)
  data1 = data
  coordinates(data) <- c("x", "y")
  spplot(data, "Tan_dependence_Wei_S")
  vgm1 <- variogram(Tan_dependence_Wei_S ~ 1, data)
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
  krige_res <- krige(Tan_dependence_Wei_S ~ 1, data, basexy1,
                     model = m)
  z = krige_res$var1.pred
  Basexyz = cbind(Basexy1, z)
  colnames(Basexyz) = c("x", "y", "Vaule")
  p2 = ggplot() + geom_raster(data = Basexyz, aes(x = x, y = y,
                                                  fill = Vaule)) + theme_bw() + scale_fill_gradientn(limits = c(min(Basexyz$Vaule),
                                                                                                                1), colours = terrain.colors(10))
  p2 = p2 + labs(title = "The Tan dependence Wei_S") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,
                                                                         0))
  print(p2)
}
