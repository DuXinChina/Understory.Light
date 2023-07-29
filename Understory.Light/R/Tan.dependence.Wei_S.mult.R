Tan.dependence.Wei_S.mult=function (a, b, tan, MI)
{
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
  d = matrix(NA, nrow(a), 3)
  pb = tkProgressBar("", "Percent complete %",
                     0, 100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Wei_S.single(a[j,
    ], b, tan)$Tan_dependence_Wei_S))
    info = sprintf("Percent complete %d%%", round(j *
                                                    100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)",
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
