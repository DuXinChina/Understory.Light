Scale.dependence.Beam_S.single=function(a,b,scale,mi,MI,bandwidth)
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
  Nei.tree$tangent=Nei.tree$Size/Nei.tree$Distance
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

