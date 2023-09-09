Scale.dependence.Beam_S.mult=function(a,b,aaix,baix,mi,MI,wei,bandwidth)
{
  library(tcltk)
  
  Scale.dependence.Beam_S.single=function(a,b,aaix,baix,mi,MI,bandwidth)
  {
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
    
    angle=pi-pi/mi
    px=sin(angle)*baix
    py=cos(angle)*baix
    r_=(px^2+py^2-aaix^2)/(2*py+2*aaix)
    r=r_+aaix
    p_o=data.frame(x=a[,1],y=a[,2]+r_)
    
    
    Neighbourhood.single1 = function(a, b) {
      
      c = b[, 1:2]
      c_= b[,1:2]
      for (i in 1:nrow(b)) {
        c[i, ] = (b[i, 1:2] - a[1, 1:2])^2
        d = (c[, 1] + c[, 2])^(1/2)
      }
      
      for (i in 1:nrow(b)) {
        c_[i, ] = (b[i, 1:2] - p_o[1, 1:2])^2
        d_ = (c_[, 1] + c_[, 2])^(1/2)
      }
      
      d = cbind(b, d, d_)
      colnames(d) = c("x", "y", "Size", "Distance", "Distance_O")
      d
    }
    
    Nei.tree = Neighbourhood.single1(a, b)
    Nei.tree[, 3] = Nei.tree[, 3] - a[, 3]
    Nei.tree = subset(Nei.tree, Size > 0)
    Nei.tree
    
    Nei.tree=unique(Nei.tree)
    row=c(which(Nei.tree[,1]==a[1,1]))
    Nei.tree[row,1]=a[1,1]+0.01
    a = a[rep(1, nrow(Nei.tree)), ]
    Weight = Nei.tree[, 1:2] - a[, 1:2]
    Weight[, 3] = atan(Weight[, 2]/abs(Weight[, 1]))
    Weight[, 3] = Weight[, 3]
    Weight[, 3] = (pi/2) + Weight[, 3]
    Weight[, 4] = ((pi) - mi * Weight[, 3])/(pi)
    Weight[which(Weight[, 4] < 0), 4] = 0
    Nei.tree$Weight=Weight[,4]
    Nei.tree$angle=Weight[, 3]
    Nei.tree=subset(Nei.tree,Weight>0)
    Nei.tree=subset(Nei.tree,Nei.tree$Distance_O<r)
    Nei.tree_left=subset(Nei.tree, x< a[1,1])
    Nei.tree_left$angle=pi+Nei.tree_left$angle
    Nei.tree_right=subset(Nei.tree,!Nei.tree$x< a[1,1])
    Nei.tree_right$angle=pi-Nei.tree_right$angle
    Nei.tree=rbind(Nei.tree_left,Nei.tree_right)
    Nei.tree$tangent=Nei.tree$Size/Nei.tree$Distance
    if(nrow(Nei.tree)==0){
      Dispersion=0
      tangent=0
    }
    
    if(nrow(Nei.tree)>0){
      y=density_box(Nei.tree$angle,bandwidth)
     # plot(y,type="l")
      y=subset(y,y[,2]>0)
      Dispersion=nrow(y)/1000
      tangent=sum(Nei.tree$tangent)
    }
    tangent=tangent
    tangent=(atan(tangent/(wei*MI)-pi)+pi/2)/pi
    tangent=tangent*Dispersion
    Scale.dependence.Beam_S=tangent
    
    if(nrow(Nei.tree)==0)
    {Scale.dependence.Beam_S=0}
    outcome =Scale.dependence.Beam_S
    outcome =list(a=a[1,],Nei.tree=Nei.tree,Scale.dependence.Beam_S=Scale.dependence.Beam_S)
    outcome
  } 
  d = matrix(NA, nrow(a), 3)
  pb = tkProgressBar("", "Percent complete %", 0, 100)
  star_time = Sys.time()
  for (j in 1:nrow(a)) {
    d[j, ] = cbind(as.matrix(a[j, 1:2]), as.matrix(Scale.dependence.Beam_S.single(a[j,],b,aaix,baix,mi,MI,bandwidth)$Scale.dependence.Beam_S))
    info = sprintf("Percent complete %d%%", round(j * 100/nrow(a)))
    setTkProgressBar(pb, j * 100/nrow(a), sprintf("Progress (%s)",
                                                  info), info)
  }
  end_time = Sys.time()
  close(pb)
  run_time = end_time - star_time
  colnames(d) = c("x", "y", "Scale_dependence_Beam_S")
  rownames(d) = 1:nrow(a)
  d = as.data.frame(d)
  d
}

