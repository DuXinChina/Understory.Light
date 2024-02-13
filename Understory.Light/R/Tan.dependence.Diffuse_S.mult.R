#' @title 计算林下多个位点的加权进界邻体散射荫蔽度
#'
#' @param a 林下研究位点的横坐标、纵坐标与高度
#' @param b 样地内林木的横坐标与树高
#' @param tan 局域临界值
#' @param MI 校正系数
#' @param wei 线性校正系数
#' @param tan2 背景临界值，可以省略
#'
#' @return 林下多个位点的加权进界邻体散射荫蔽度
#' @export
#'
#' @examples
#' a=data.frame(x=c(25,35,40),y=c(25,21,10),size=c(0,1,2))
#' b=Understory.Light::b
#' Understory.light::Tan.dependence.Diffuse_S.mult(a=a,b=b,tan=2.8,MI=10.660,wei=0.5)
Tan.dependence.Diffuse_S.mult=function (a, b, tan, MI,wei,tan2)
{
  if(missing(tan2))
  {tan2=NA}
  library(tcltk)
  Tan.dependence.Wei_S.single = function(a, b, tan,tan2) {
    Tan.dependence.S.single = function(a, b, tan,tan2) {
      Neighbourhood.single1 = function(a, b, tan,tan2) {
        c = b[, 1:2]
        for (i in 1:nrow(b)) {
          c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
          d = (c[, 1] + c[, 2])^(1/2)
        }
        d = cbind(b, d)
        d$tangent = (d$size - a[, 3])/d$d
        if (is.na(tan2)==F)
        {d2=subset(d, tangent >= tan2 & tangent < tan)
        ta2sum<<-sum(d2$tangent)
        }

        d = subset(d, tangent > tan)
        tansum<<-sum(d$tangent)
        colnames(d) = c("x", "Y", "Size",
                        "Distance", "tangent")
        d
      }
      Nei.tree = Neighbourhood.single1(a, b, tan,tan2)
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

      Tan_dependence_S = (atan(Tan_dependence_S/(wei*MI)-pi)+pi/2)/pi
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
    S = Tan.dependence.S.single(a, b, tan,tan2)
    Levins_Simpson = Tan.dependence.Simpson(a, b, tan)
    a = S$a
    Neighbourhood = S$Neighbourhood
    Tan_dependence_S = S$Tan_dependence_S
    Tan_dependence_Wei_S = Tan_dependence_S * (0.25/Levins_Simpson)
    if (is.na(tan2)==F)
    {
      #Ratio=((ta2sum)*sqrt((1/tan2)^2-(1/tan)^2))/((MI*pi)/tan)
      Ratio=((ta2sum))/((MI*pi))
      Tan_dependence_Wei_S=(atan((Ratio+Tan_dependence_Wei_S)/(0.5/pi*(tan/tan2))-pi)+pi/2)/pi

      }
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
    res=Tan.dependence.Wei_S.single(a[j,
    ], b, tan,tan2)

    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(res$Tan_dependence_Wei_S))
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
