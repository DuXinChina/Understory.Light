K.mult=function (a, b)
{
  library(tcltk)
  Neighbourhood.K.single = function(a, b) {
    b1 = subset(b, x >= a[, 1] & y > a[, 2] & size > a[,
                                                       3])
    b2 = subset(b, x > a[, 1] & y <= a[, 2] & size > a[,
                                                       3])
    b3 = subset(b, x <= a[, 1] & y < a[, 2] & size > a[,
                                                       3])
    b4 = subset(b, x < a[, 1] & y >= a[, 2] & size > a[,
                                                       3])
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
    if (nrow(b1) != 0 & nrow(b2) != 0 & nrow(b3) != 0 & nrow(b4) !=
        0) {
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
    info <- sprintf("已完成 %d%%", round(j * 100/nrow(a)))
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
