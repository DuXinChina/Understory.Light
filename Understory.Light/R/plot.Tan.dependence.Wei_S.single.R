plot.Tan.dependence.Wei_S.single=function (a, b, tan, MI)
{
  library(ggplot2)
  Tan.dependence.Wei_S.single = function(a, b, tan, MI) {
    Tan.dependence.S.single = function(a, b, tan) {
      Neighbourhood.single1 = function(a, b, tan) {
        c = b[, 1:2]
        for (i in 1:nrow(b)) {
          c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
          d = (c[, 1] + c[, 2])^(1/2)
        }
        d = cbind(b, d)
        d = subset(d, d > 0)
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
        Neighbourhood1[which(Neighbourhood1[, 5] <= 0),
                       5] = 0
        S1 = Neighbourhood1[, 5]
        S2 = as.data.frame(S1)
        S = sum(S2, na.rm = T)
        Tan_dependence_S = S
      }
      Tan_dependence_S = (atan(Tan_dependence_S/(MI * pi))/pi *
                            2)
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
        d = subset(d, d > 0)
        d$tangent = (d$size - a[, 3])/d$d
        d = subset(d, tangent > tan)
        colnames(d) = c("x", "y", "Size",
                        "Distance", "tangent")
        d
      }
      Nei.tree = Neighbourhood.single1(a, b, tan)
      Nei.tree[, 3] = Nei.tree[, 3] - a[, 3]
      Nei.tree = subset(Nei.tree, Size > 0)
      Nei.tree = subset(Nei.tree, Distance > 0)
      Nei.tree
      Simpn = sum(Nei.tree[, 3])
      Nei.tree1 = subset(Nei.tree, x >= a[, 1] & y > a[,
                                                       2])
      Nei.tree2 = subset(Nei.tree, x > a[, 1] & y <= a[,
                                                       2])
      Nei.tree3 = subset(Nei.tree, x <= a[, 1] & y < a[,
                                                       2])
      Nei.tree4 = subset(Nei.tree, x < a[, 1] & y >= a[,
                                                       2])
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
  Ne = Tan.dependence.Wei_S.single(a, b, tan, MI)
  big = subset(Ne$Neighbourhood, Size > a[, 3])
  small = subset(Ne$Neighbourhood, Size <= a[, 3])
  n = nrow(big)
  if (n == 0) {
    center = cbind(a[, 1:3], 1)
    colnames(center) = c("x", "y", "size",
                         "group")
    unitg = center
    Neigh = subset(unitg, size < 0)
    colnames(Neigh) = c("x", "y", "size",
                        "group")
    colnames(small) = c("x", "y", "size")
    max = 10
  }
  if (n != 0) {
    center = cbind(cbind(rep(Ne$a[, 1], each = n), rep(Ne$a[,
                                                            2], each = n), rep(Ne$a[, 3], each = n)), c(1:n))
    center = as.data.frame(center)
    colnames(center) = c("x", "y", "size",
                         "group")
    Neigh = cbind(big, c(1:n))
    Neigh = Neigh[, c(1, 2, 3, 6)]
    colnames(Neigh) = c("x", "y", "size",
                        "group")
    colnames(small) = c("x", "y", "size")
    unitg = rbind(center, Neigh)
    max = max(abs(Ne$Neighbourhood[, 1] - Ne$a[, 1]), abs(Ne$Neighbourhood[,
                                                                           2] - Ne$a[, 2]))
  }
  p = ggplot()
  p = p + geom_hline(data = Ne$a, aes(yintercept = y), linetype = 5,
                     col = "red") + geom_vline(data = Ne$a, aes(xintercept = x),
                                               linetype = 5, col = "red")
  p = p + geom_line(data = unitg, aes(x = x, y = y, group = group),
                    linetype = 2, size = 1)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 13/max(unitg[,
                                                                 3]) * a[, 3] + 2, color = "red4", alpha = 0.7) +
    geom_point(data = Neigh, aes(x, y), size = 13/max(unitg[,
                                                            3]) * Neigh[, 3] + 2, color = "green4", alpha = 0.7)
  p = p + geom_point(data = small, aes(x, y), size = 13/max(unitg[,
                                                                  3]) * small[, 3] + 2, color = "grey", alpha = 0.7)
  p = p + geom_point(data = Ne$a, aes(x, y), size = 2) + geom_point(data = Neigh,
                                                                    aes(x, y), size = 2) + geom_point(data = small, aes(x,
                                                                                                                        y), size = 2)
  p = p + lims(x = c(Ne$a[, 1] - 1.1 * max, Ne$a[, 1] + 1.1 *
                       max), y = c(Ne$a[, 2] - 1.1 * max, Ne$a[, 2] + 1.1 *
                                     max))
  p = p + annotate("text", x = Ne$a[, 1] - max + 0.3 *
                     (max), y = Ne$a[, 2] + max - 0.125 * (max), label = paste0("Tan_dependence_Wei_S=",
                                                                                round(Ne$Tan_dependence_Wei_S, 3)))
  p = p + theme_bw()
  suppressMessages(print(p))
}
