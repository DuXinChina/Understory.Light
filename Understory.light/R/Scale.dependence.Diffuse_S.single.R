Scale.dependence.Diffuse_S.single=function (a, b, scale, MI,wei,scale2) 
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
