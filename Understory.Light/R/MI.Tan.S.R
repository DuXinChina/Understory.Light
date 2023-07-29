MI.Tan.S=function (minx, maxx, miny, maxy, b, seq, tan)
{
  library(sp)
  library(gstat)
  library(tcltk)
  xgrid = seq(minx, maxx, length.out = seq + 1)
  ygrid = seq(miny, maxy, length.out = seq + 1)
  basexy = expand.grid(xgrid, ygrid)
  basexy[, 3] = 0
  colnames(basexy) <- c("x", "y", "Size")
  a = dplyr::distinct(basexy)
  Tan.dependence.Wei_S.mult = function(a, b, tan, MI) {
    library(tcltk)
    library(sp)
    library(gstat)
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
          Neighbourhood1[which(Neighbourhood1[, 4] ==
                                 0), 4] = 1e-08
          S1 = Neighbourhood1[, 5]
          S2 = as.data.frame(S1)
          S = sum(S2)
          if (S == Inf)
            S = 1e+08
          Tan_dependence_S = S
        }
        if (is.nan(Tan_dependence_S) == T)
          (Tan_dependence_S = 1)
        outcome = list(a = a, Neighbourhood = Neighbourhood,
                       Tan_dependence_S = Tan_dependence_S)
        outcome
      }
      S = Tan.dependence.S.single(a, b, tan)
      a = S$a
      Neighbourhood = S$Neighbourhood
      Tan_dependence_S = S$Tan_dependence_S
      outcome = list(a = a, Neighbourhood = Neighbourhood,
                     Tan_dependence_S = Tan_dependence_S)
      outcome
    }
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("Progress", "Percent complete %",
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Wei_S.single(a[j,
      ], b, tan)$Tan_dependence_S))
      info = sprintf("Percent complete %d%%", round(j *
                                                      100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)",
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Tan_dependence_S")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  data = Tan.dependence.Wei_S.mult(a, b, tan, MI)
  Tan_dependence_S = data$Tan_dependence_S
  Rank = rank(Tan_dependence_S)
  Percentile_Rank_of_Tan_dependence_S = Rank/length(Rank)
  nls = nls(Percentile_Rank_of_Tan_dependence_S ~ (atan(Tan_dependence_S/(MI *
                                                                            pi))/pi * 2), start = list(MI = 1), algorithm = "port")
  res = summary(nls)
  MI = res$parameters[1, 1]
  Revise_Tan_dependence_S = (atan(Tan_dependence_S/(MI * pi))/pi *
                               2)
  plot(Percentile_Rank_of_Tan_dependence_S, Revise_Tan_dependence_S,
       main = "Q-Q plot", xlim = c(0, 1), ylim = c(0,
                                                   1)) + abline(0, 1, col = "red4")
  R_square = 1 - sum((Percentile_Rank_of_Tan_dependence_S -
                        Revise_Tan_dependence_S)^2)/sum((Percentile_Rank_of_Tan_dependence_S -
                                                           mean(Percentile_Rank_of_Tan_dependence_S))^2)
  output = list(MI = MI, R_square = R_square)
  output
}
