#' @title  计算尺度化加权散射荫蔽度在特定局域尺度下最适的校正系数
#'
#' @param minx 研究区域横坐标的最小值
#' @param maxx 研究区域横坐标的最大值
#' @param miny 研究区域纵坐标的最小值
#' @param maxy 研究区域纵坐标的最大值
#' @param b 林分内林木的坐标及树高
#' @param seq  抽样数据的空间分辨率
#' @param scale 局域尺度圆的半径
#'
#' @return 尺度化加权散射荫蔽度在特定局域尺度下最适的校正系数
#' @export
#'
#' @examples
#' b=Understory.light::b
#' Understory.light::MI.Diffuse_Scale.S(minx=50, maxx=100, miny=50, maxy=100, b, seq=5, scale=5.4)
MI.Diffuse_Scale.S=function (minx, maxx, miny, maxy, b, seq, scale)
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
  Scale.dependence.Wei_S.mult = function(a, b, scale) {
    library(tcltk)
    library(sp)
    library(gstat)
    Scale.dependence.Wei_S.single = function(a, b, scale) {
      Scale.dependence.S.single = function(a, b, scale) {
        Neighbourhood.single1 = function(a, b, scale) {
          c = b[, 1:2]
          for (i in 1:nrow(b)) {
            c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
            d = (c[, 1] + c[, 2])^(1/2)
          }
          d = cbind(b, d)
          d = subset(d, d < scale)
          colnames(d) = c("x", "Y", "Size",
                          "Distance")
          d
        }
        Nei.tree = Neighbourhood.single1(a, b, scale)
        n = nrow(Nei.tree)
        Nei.dif.high = cbind(Nei.tree, Nei.tree[, 3] -
                               a[, 3])
        colnames(Nei.dif.high) = c("x", "Y",
                                   "Size", "Distance", "Size_dif ")
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
          Neighbourhood[which(Neighbourhood[, 5] <= 0),
                        5] = NA
          Neighbourhood1[which(Neighbourhood1[, 5] <=
                                 0), 5] = 0
          S1 = Neighbourhood1[, 5]/Neighbourhood1[, 4]
          S2 = as.data.frame(S1)
          S = sum(S2)
          if (S == Inf)
            S = 1e+08
          Scale_dependence_S = S
        }
        outcome = list(a = a, Neighbourhood = Neighbourhood,
                       Scale_dependence_S = Scale_dependence_S)
        outcome
      }
      S = Scale.dependence.S.single(a, b, scale)
      a = S$a
      Neighbourhood = S$Neighbourhood
      Scale_dependence_S = S$Scale_dependence_S
      outcome = list(a = a, Neighbourhood = Neighbourhood,
                     Scale_dependence_S = Scale_dependence_S)
      outcome
    }
    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("Progress", "Percent complete %",
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Scale.dependence.Wei_S.single(a[j,
      ], b, scale)$Scale_dependence_S))
      info = sprintf("Percent complete %d%%", round(j *
                                                      100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)",
                                                    info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Scale_dependence_S")
    rownames(d) = 1:nrow(a)
    d = as.data.frame(d)
    d
  }
  data = Scale.dependence.Wei_S.mult(a, b, scale)
  Scale_dependence_S = data$Scale_dependence_S
  Rank = rank(Scale_dependence_S)
  Percentile_Rank_of_Scale_dependence_S = Rank/length(Rank)
  nls = nls(Percentile_Rank_of_Scale_dependence_S ~ (atan(Scale_dependence_S/(MI)-pi)+pi/2)/pi, start = list(MI = 1), algorithm = "port")
  res = summary(nls)
  MI = res$parameters[1, 1]
  Revise_Scale_dependence_S = (atan(Scale_dependence_S/(MI)-pi)+pi/2)/pi
  plot(Percentile_Rank_of_Scale_dependence_S, Revise_Scale_dependence_S,
       main = "Q-Q plot", xlim = c(0, 1), ylim = c(0,
                                                   1)) + abline(0, 1, col = "red4")
  R_square = 1 - sum((Percentile_Rank_of_Scale_dependence_S -
                        Revise_Scale_dependence_S)^2)/sum((Percentile_Rank_of_Scale_dependence_S -
                                                             mean(Percentile_Rank_of_Scale_dependence_S))^2)
  output = list(MI = MI, R_square = R_square)
  output
}
