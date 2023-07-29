plot.Tan.dependence.Beam_S.Krig=function (minx, maxx, miny, maxy, b, seq, tan, mi, MI, bandwidth)
{
  library(sp)
  library(gstat)
  library(tcltk)
  library(ggplot2)

  Tan.dependence.Beam_S.mult = function(a, b, tan, mi, MI, bandwidth) {

    library(tcltk)
    Tan.dependence.Beam_S.single = function(a, b, tan, mi,MI,bandwidth) {

      density_box=function(entity,bandwidth)
      {
        if (missing(bandwidth)) {
          bandwidth=pi/5
        }
        min_x=min(entity)
        max_x=max(entity)
        x=seq(pi-pi/mi,pi+pi/mi,length.out=1000)
        y=rep(0, length(x))
        box=function(entity){
          for (i in 1:length(x)) {
            if (abs(x[i] - entity) <= bandwidth/2) {
              y[i] <<- y[i] + 1
            }
          }
        }
        for(i in 1:length(entity))
        {
          box(entity[i])
        }
        y=y/length(entity)
        result=cbind(x,y)
        colnames(result)=c("angle","density")
        result
      }


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

      if (nrow(Nei.tree)>0)
      {
        Nei.tree=unique(Nei.tree)
        a = a[rep(1, nrow(Nei.tree)), ]
        Weight = Nei.tree[, 1:2] - a[, 1:2]
        Weight[, 3] = atan(Weight[, 2]/abs(Weight[, 1]))
        Weight[, 3] = Weight[, 3]
        Weight[, 3] = (pi/2) + Weight[, 3]
        Weight[, 4] = ((pi) - mi * Weight[, 3])/(pi)
        Weight[which(Weight[, 4] < 0), 4] = 0
        Nei.tree$Weight=Weight[,4]
        Nei.tree$angle=Weight[, 3]
        Nei.tree_left=subset(Nei.tree,Nei.tree$x-a[1,1]<0)
        Nei.tree_left$angle=pi-Nei.tree_left$angle
        Nei.tree_right=subset(Nei.tree,Nei.tree$x-a[1,1]>0)
        Nei.tree_right$angle=pi+Nei.tree_right$angle
        Nei.tree=rbind(Nei.tree_left,Nei.tree_right)
        Nei.tree=subset(Nei.tree,Weight>0)
        if(nrow(Nei.tree)==0){
          Dispersion=0
          tangent=0
        }
        if(nrow(Nei.tree)>0){
          y=density_box(Nei.tree$angle,bandwidth)
          y=subset(y,y[,2]>0)

          Dispersion=nrow(y)/1000
          tangent=sum(Nei.tree$tangent)
        }
        tangent=tangent
        tangent=(atan(tangent/(MI *pi))/pi * 2)
        tangent=tangent*Dispersion
        Tan.dependence.Beam_S=tangent
      }
      if(nrow(Nei.tree)==0)
      {Tan.dependence.Beam_S=0}
      outcome =Tan.dependence.Beam_S
      outcome
    }





    d = matrix(NA, nrow(a), 3)
    pb = tkProgressBar("", "Percent complete %",
                       0, 100)
    star_time = Sys.time()
    for (j in 1:nrow(a)) {
      d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Tan.dependence.Beam_S.single(a[j,
      ], b, tan, mi, MI, bandwidth)))
      info = sprintf("Percent complete %d%%", round(j *  100/nrow(a)))
      setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)", info), info)
    }
    end_time = Sys.time()
    close(pb)
    run_time = end_time - star_time
    colnames(d) = c("x", "y", "Tan_dependence_Beam_S")
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
  data = Tan.dependence.Beam_S.mult(basexy, b, tan, mi, MI , bandwidth)
  data1 = data
  coordinates(data) <- c("x", "y")
  spplot(data, "Tan_dependence_Beam_S")
  vgm1 <- variogram(Tan_dependence_Beam_S ~ 1, data)
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
  p1 = plot(vgm1, model = m)
  print(p1)
  krige_res <- krige(Tan_dependence_Beam_S ~ 1, data, basexy1,
                     model = m)
  z = krige_res$var1.pred
  Basexyz = cbind(Basexy1, z)
  colnames(Basexyz) = c("x", "y", "Vaule")
  Basexyz[which(Basexyz$Vaule < 0), 3] = 0
  p2 = ggplot() + geom_raster(data = Basexyz, aes(x = x, y = y,
                                                  fill = Vaule)) + theme_bw() + scale_fill_gradientn(limits = c(min(Basexyz$Vaule),
                                                                                                                max(Basexyz$Vaule)), colours = terrain.colors(10))
  p2 = p2 + labs(title = "Tan_dependence_Beam_S") + scale_x_continuous(expand = c(0,
                                                                                  0)) + scale_y_continuous(expand = c(0, 0))
  print(p2)
}
