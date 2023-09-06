plot.K.single=function (a, b)
{
  library(ggplot2)
  b1 = subset(b, x >= a[, 1] & y > a[, 2] & size > a[, 3])
  b2 = subset(b, x > a[, 1] & y <= a[, 2] & size > a[, 3])
  b3 = subset(b, x <= a[, 1] & y < a[, 2] & size > a[, 3])
  b4 = subset(b, x < a[, 1] & y >= a[, 2] & size > a[, 3])
  options(warn = -1)
  K.single = function(a, b) {
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
      Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] - a[,
                                                       3])
      colnames(Nei.dif.high) = c("x", "y",
                                 "Size", "Distance", "Size_dif ")
      Neighbourhood = Nei.dif.high
      k1 = Neighbourhood[, 4]/Neighbourhood[, 5]
      k2 = as.data.frame(k1)
      k2[, 1][is.infinite(k2[, 1])] = NA
      K = sum(k2, na.rm = T)
      outcome = list(a = a, Neighbourhood = Neighbourhood,
                     K = K)
      outcome
    }
    if (nrow(b1) != 0 & nrow(b2) != 0 & nrow(b3) != 0 & nrow(b4) !=
        0) {
      Neighbourhood.K.single1(a, b)
    }
    else {
    }
  }
  options(warn = 1)
  Ne = K.single(a, b)
  n=nrow(Ne$Neighbourhood)
  if (nrow(b1) != 0 & nrow(b2) != 0 & nrow(b3) != 0 & nrow(b4) !=
      0) {
    center = cbind(cbind(rep(a[, 1], each = n), rep(a[, 2],
                                                    each = n), rep(a[, 3], each = n)), c(1:n))
    center = as.data.frame(center)
    colnames(center) = c("x", "y", "size",
                         "group")
    Neigh = cbind(Ne$Neighbourhood[, 1:3], c(1:n))
    colnames(Neigh) = c("x", "y", "size",
                        "group")
    unitg = rbind(center, Neigh)
    max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[,
                                                                           2] - Ne$a[, 2]))
    p = ggplot()
    p = p + geom_hline(data = Ne$a, aes(yintercept = y),
                       linetype = 5, col = "red") + geom_vline(data = Ne$a,
                                                               aes(xintercept = x), linetype = 5, col = "red")
    p = p + geom_line(data = unitg, aes(x = x, y = y, group = group),
                      linetype = 2, size = 1) + geom_point(data = Ne$a,
                                                           aes(x, y), size = 13/max(unitg[, 3]) * a[, 3] + 2,
                                                           color = "red4", alpha = 0.7) + geom_point(data = Ne$Neighbourhood,
                                                                                                     aes(x, y), size = 13/max(unitg[, 3]) * Ne$Neighbourhood[,
                                                                                                                                                             3] + 2, color = "green4", alpha = 0.7)
    p = p + geom_point(data = Ne$a, aes(x, y), size = 2) +
      geom_point(data = Ne$Neighbourhood, aes(x, y), size = 2)
    p = p + lims(x = c(Ne$a[, 1] - max, Ne$a[, 1] + max),
                 y = c(Ne$a[, 2] - max, Ne$a[, 2] + max))
    p = p + annotate("text", x = Ne$a[, 1] - max +
                       0.125 * (max), y = Ne$a[, 2] + max - 0.125 * (max),
                     label = paste0("K=", round(Ne$K, 3)))
    p + theme_bw()
  }
  else {
    warning("确保以参考林木为中心的四个象限内均有高于其的立木")
  }
}
